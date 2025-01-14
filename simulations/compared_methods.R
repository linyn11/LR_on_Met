

knn_class = function(mfd, X.train, y.train, X.ind, k=4)
{
	pdist = apply(X.train, 3, function(p)geo.dist(mfd, p, q=X.ind))

	voted = NULL
	for(i in k){
		i.idx = base::order(pdist)[1:i]
		vote = major_voting(y.train, idx=i.idx)
		voted = c(voted, vote)
	}

	return(voted)
}


major_voting = function(y, idx)
{
	k = length(idx)
	total.votes = sum(y[idx])
	if( total.votes>k/2 ){
		return(1)
	}else if( total.votes<k/2 ){
		return(0)
	}else{
		return(sample(c(0,1), size=1))
	}
}


knn_cv = function(mfd, X.cur, y.cur, Ks, fold=5, seed=1, cv.crt="ACC")
{
	n = length(y.cur)
	set.seed(seed*fold)
	flds = caret::createFolds(1:n, k=fold, list=TRUE, returnTrain=FALSE)

	cv.mat = NULL
	for(i in 1:fold)
	{
		idx = flds[[i]]
		X.val = X.cur[,,idx]
		y.val = y.cur[idx]
		X.train = X.cur[,,-idx]
		y.train = y.cur[-idx]

		y.cv = t(apply(X.val, 3, function(x)knn_class(mfd, X.train, y.train, x, k=Ks)))
		if( cv.crt=="ACC" ){
			cv.mat = rbind(cv.mat, apply(y.cv, 2, function(x)mean(y.val==x) ))
		}else if( cv.crt=="AUC" ){
			cv.mat = rbind(cv.mat, apply(y.cv, 2, function(x)pROC::roc(y.val, x, quiet=TRUE)$auc ))
		}
	}
	cv.means = colMeans(cv.mat)
	names(cv.means) = Ks

	results = list("K"=Ks, "CV.means"=cv.means, "Opt.K"=Ks[which.max(cv.means)])
	return(results)
}


knn_single = function(mfd, X.train, y.train, X.test, y.test, linear.star, Ks=seq(3, 51, by=2), seed=0, 
		cv.crt="ACC", setting="logit", fold=5)
{
	cv.res = knn_cv(mfd, X.train, y.train, Ks=Ks, fold=fold, seed=seed*100, cv.crt=cv.crt)
	K.cv = cv.res$Opt.K
	y.pred = apply(X.test, 3, function(x)knn_class(mfd, X.train, y.train, x, k=K.cv))	# make prediction
	
	# evaluations
	evals = eval_class_method(y.test, linear.star, y.pred, setting=setting)
	return(evals)
}


ksvm_cv_inner = function(mfd, gamma.r, df.train, X.val.mat, y.val)
{
	# define kernel
	m = mfd$dim[1]
	Gk = function(X1, X2)
	{
		X1.mat = matrix(X1, byrow=TRUE, nrow=m, ncol=m)
		X2.mat = matrix(X2, byrow=TRUE, nrow=m, ncol=m)
		d2 = geo.dist(mfd, X1.mat, X2.mat)^2
		kv = exp(-gamma.r*d2)
		return(kv)
	}
	class(Gk)="kernel"

	fit = ksvm(y~., data=df.train, kernel=Gk, scaled=FALSE)
	y.pred = predict(fit, X.val.mat)
	acc = mean(as.numeric(y.val)==as.numeric(y.pred))
	return(acc)
}

#X.cur=X.train;y.cur=y.train
ksvm_cv = function(mfd, X.cur, y.cur, Gammas, fold=5, seed=1)
{
	m = mfd$dim[1]
	n = length(y.cur)
	set.seed(seed*fold)
	flds = caret::createFolds(1:n, k=fold, list=TRUE, returnTrain=FALSE)

	cv.mat = NULL
	for(i in 1:fold)
	{
		idx = flds[[i]]
		X.val = X.cur[,,idx]
		y.val = y.cur[idx]
		X.train = X.cur[,,-idx]
		y.train = y.cur[-idx]

		X.train.mat = matrix(X.train, byrow=TRUE, nrow=dim(X.train)[3], ncol=m*m)
		y.train = as.factor(y.train)
		X.val.mat = matrix(X.val, byrow=TRUE, nrow=dim(X.val)[3], ncol=m*m)
		y.val = as.factor(y.val)
		df.train = data.frame(X.train.mat, y=y.train)

		# accuracy as the criterion
		ACCs = sapply(Gammas, function(gamma.r)ksvm_cv_inner(mfd, gamma.r, df.train, X.val.mat, y.val))
		cv.mat = rbind(cv.mat, ACCs)
	}
	cv.means = colMeans(cv.mat)
	names(cv.means) = Gammas

	results = list("gamma.r"=Gammas, "CV.means"=cv.means, "Opt.gamma.r"=Gammas[which.max(cv.means)])
	return(results)
}


ksvm_single = function(mfd, X.train, y.train, X.test, y.test, linear.star, Gammas=10^seq(-3, 3, by=1), seed=0, 
		setting="logit", fold=5)
{
	cv.res = ksvm_cv(mfd, X.train, y.train, Gammas=Gammas, fold=fold, seed=seed*100)
	gamma.r = cv.res$Opt.gamma.r

	# define kernel
	m = mfd$dim[1]
	Gk = function(X1, X2)
	{
		X1.mat = matrix(X1, byrow=TRUE, nrow=m, ncol=m)
		X2.mat = matrix(X2, byrow=TRUE, nrow=m, ncol=m)
		d2 = geo.dist(mfd, X1.mat, X2.mat)^2
		kv = exp(-gamma.r*d2)
		return(kv)
	}
	class(Gk)="kernel"

	X.train.mat = matrix(X.train, byrow=TRUE, nrow=dim(X.train)[3], ncol=m*m)
	y.train = as.factor(y.train)
	X.val.mat = matrix(X.test, byrow=TRUE, nrow=dim(X.test)[3], ncol=m*m)
	y.val = as.factor(y.test)
	df.train = data.frame(X.train.mat, y=y.train)

	fit = ksvm(y~., data=df.train, kernel=Gk, scaled=FALSE)
	y.pred = as.numeric(predict(fit, X.val.mat))-1

	# evaluations
	evals = eval_class_method(y.test, linear.star, y.pred, setting=setting)
	return(evals)
}


kdc_cv_inner = function(mfd, h.cur, X.train, y.train, X.val, y.val)
{
	# define kernel
	m = mfd$dim[1]
	D = m*(m+1)/2
	gk = function(y, x, h)
	{
		y.mat = matrix(y, byrow=TRUE, nrow=m, ncol=m)
		x.mat = matrix(x, byrow=TRUE, nrow=m, ncol=m)
		d2 = geo.dist(mfd, y.mat, x.mat)^2
		k = 1/sqrt(2*pi) * exp(-d2/(2*h^2))
		return(k)
	}

	X.train.mat = matrix(X.train, byrow=TRUE, nrow=dim(X.train)[3], ncol=m*m)
	X.val.mat = matrix(X.val, byrow=TRUE, nrow=dim(X.val)[3], ncol=m*m)
	y.pred = NULL
	for( i in 1:nrow(X.val.mat))
	{
		y = X.val.mat[i,]
		c.0.kde = sum(apply(X.train.mat[which(y.train==0),], 1, function(x)gk(y, x, h.cur)))
		c.1.kde = sum(apply(X.train.mat[which(y.train==1),], 1, function(x)gk(y, x, h.cur)))
		y.pred.cur = ifelse(c.0.kde>c.1.kde, 0, 1)
		y.pred = c(y.pred, y.pred.cur)
	}

	res = mean(y.pred==y.val)
	return(res)
}


kdc_cv = function(mfd, X.cur, y.cur, Hs, fold=5, seed=1)
{
	n = length(y.cur)
	set.seed(seed*fold)
	flds = caret::createFolds(1:n, k=fold, list=TRUE, returnTrain=FALSE)

	cv.mat = NULL
	for(i in 1:fold)
	{
		idx = flds[[i]]
		X.val = X.cur[,,idx]
		y.val = y.cur[idx]
		X.train = X.cur[,,-idx]
		y.train = y.cur[-idx]

		cv.cur = sapply(Hs, function(h)kdc_cv_inner(mfd, h, X.train, y.train, X.val, y.val))
		cv.mat = rbind(cv.mat, cv.cur)
	}
	cv.means = colMeans(cv.mat)
	names(cv.means) = Hs

	results = list("Hs"=Hs, "CV.means"=cv.means, "Opt.h"=Hs[which.max(cv.means)])
	return(results)
}


#linear.star=data.list$linear.star[test.idx];fold=5;Hs=10^seq(-3, 3, by=1)
kdc_single = function(mfd, X.train, y.train, X.test, y.test, linear.star, Hs=10^seq(-3, 3, by=1), seed=0, 
		setting="logit", fold=5)
{
	cv.res = kdc_cv(mfd, X.train, y.train, Hs=Hs, fold=fold, seed=seed*100)
	h.cv = cv.res$Opt.h

	# define kernel
	m = mfd$dim[1]
	D = m*(m+1)/2
	gk = function(y, x, h)
	{
		y.mat = matrix(y, byrow=TRUE, nrow=m, ncol=m)
		x.mat = matrix(x, byrow=TRUE, nrow=m, ncol=m)
		d2 = geo.dist(mfd, y.mat, x.mat)^2
		k = 1/sqrt(2*pi) * exp(-d2/(2*h^2))
		return(k)
	}

	# kde, an estimator from Pelletier (2005), 
	# with the volume density function equals to 1 everywhere for flat space, 
	# which is the case for the SPD endowed with the Log-Cholesky metric
	X.train.mat = matrix(X.train, byrow=TRUE, nrow=dim(X.train)[3], ncol=m*m)
	X.test.mat = matrix(X.test, byrow=TRUE, nrow=dim(X.test)[3], ncol=m*m)
	y.pred = NULL
	for( i in 1:nrow(X.test.mat))
	{
		y = X.test.mat[i,]
		c.0.kde = sum(apply(X.train.mat[which(y.train==0),], 1, function(x)gk(y, x, h.cv)))
		c.1.kde = sum(apply(X.train.mat[which(y.train==1),], 1, function(x)gk(y, x, h.cv)))
		y.pred.cur = ifelse(c.0.kde>c.1.kde, 0, 1)
		y.pred = c(y.pred, y.pred.cur)
	}

	# evaluations
	evals = eval_class_method(y.test, linear.star, y.pred, setting=setting)
	return(evals)
}


eval_class_method = function(y, linear.star, y.pred, th=0.5, setting="logit")
{
	if( setting=="IP" ){
		prob.star = sigmoid(2*sin(linear.star*pi)+linear.star)
	}else{
		prob.star = sigmoid(linear.star)
	}
	oracle.class = as.numeric(prob.star>th)

	# classification performance
	AUC = pROC::roc(as.vector(y), y.pred, quiet=TRUE)$auc
	AUC.oracle = pROC::roc(as.vector(y), oracle.class, quiet=TRUE)$auc
	AUC.diff = AUC.oracle - AUC

	accuracy = mean(y.pred == y)
	accuracy.oracle = mean(oracle.class == y)
	acc.diff = accuracy.oracle - accuracy

	sen = sensitivity(y, y.pred)
	spe = specificity(y, y.pred)
	sen.oracle = sensitivity(y, oracle.class)
	spe.oracle = specificity(y, oracle.class)
	sen.diff = sen.oracle - sen
	spe.diff = spe.oracle - spe

	res = c(accuracy, acc.diff, AUC, AUC.diff, sen, sen.diff, spe, spe.diff)
	names(res) = c("Accuracy", "Accuracy-diff", "AUC", "AUC-diff", "Sensitivity", "Sensitivity-diff", 
			"Specificity", "Specificity-diff")
	return(res)
}


#data.list=NULL;cv.crt="ACC"
method_single = function(mfd, data.list=NULL, seed=0, test.ratio=0.2, cv.crt="ACC", setting="logit")
{
	file.compared = paste0("tmp_compared_", seed, ".RData")
	if( file.exists(file.compared) ){
		load(file.compared)
		return(res)
	}

	if( is.null(data.list) ){
		file = paste0("tmp_", seed, ".RData")
		load(file)
	}

	convergence = res.cur[length(res.cur)]
	X = data.list$X
	y = data.list$y

	m = mfd$dim[1]
	n = length(y)
	set.seed(seed)
	test.idx = sample(n, size=round(test.ratio*n))
	X.train = X[,,-test.idx]
	y.train = y[-test.idx]
	X.test = X[,,test.idx]
	y.test = y[test.idx]

	flag = 10000
	while( length(unique(y.test))==1 | length(unique(y.train))==1 ){
		set.seed(seed+flag)
		test.idx = sample(n, size=round(test.ratio*n))
		X.train = X[,,-test.idx]
		y.train = y[-test.idx]
		X.test = X[,,test.idx]
		y.test = y[test.idx]
		flag = flag + 10000
	}

	#data.list$linear.star[test.idx]
	# knn
	res.knn = knn_single(mfd, X.train, y.train, X.test, y.test, linear.star=data.list$linear.star[test.idx], seed=seed*100, cv.crt=cv.crt, setting=setting)
	res.knn = c(res.knn, convergence)
	# kernel svm
	res.ksvm = ksvm_single(mfd, X.train, y.train, X.test, y.test, linear.star=data.list$linear.star[test.idx], seed=seed*100, setting=setting)
	res.ksvm = c(res.ksvm, convergence)
	# kernel density classifier
	res.kdc = kdc_single(mfd, X.train, y.train, X.test, y.test, linear.star=data.list$linear.star[test.idx], seed=seed*100, setting=setting)
	res.kdc = c(res.kdc, convergence)

	res = rbind(res.knn, res.ksvm, res.kdc)
	rownames(res) = c("KNN", "KSVM", "KDC")

	save.list = c("res", "seed")
	save(list=save.list, file=file.compared)
	return(res)
}


output_matlab = function(path, tmp.res.path, seed=0, test.ratio=0.2, data.list=NULL, setting="logit")
{
	setwd(tmp.res.path)
	if( is.null(data.list) ){
		file = paste0("tmp_", seed, ".RData")
		load(file)
	}

	cur.path = paste0(path, "tmp-", seed, "/")
	if(!dir.exists(cur.path)) dir.create(cur.path, recursive=TRUE)
	setwd(cur.path)
	csv.file = paste0("oracle-results-tmp-", seed, ".csv")
	if( file.exists(csv.file) & dir.exists(paste0(cur.path, "train/")) & dir.exists(paste0(cur.path, "test/")) ) return(0)
	setwd(path)

	# split data
	X = data.list$X
	y = data.list$y

	y = y + 1
	n = length(y)
	set.seed(seed)
	test.idx = sample(n, size=round(test.ratio*n))
	X.train = X[,,-test.idx]
	y.train = y[-test.idx]
	X.test = X[,,test.idx]
	y.test = y[test.idx]

	flag = 10000
	while( length(unique(y.test))==1 | length(unique(y.train))==1 ){
		set.seed(seed+flag)
		test.idx = sample(n, size=round(test.ratio*n))
		X.train = X[,,-test.idx]
		y.train = y[-test.idx]
		X.test = X[,,test.idx]
		y.test = y[test.idx]
		flag = flag + 10000
	}

	pos.path.train = paste0(cur.path, "train/1/")
	if(!dir.exists(pos.path.train)) dir.create(pos.path.train, recursive=TRUE)
	neg.path.train = paste0(cur.path, "train/2/")
	if(!dir.exists(neg.path.train)) dir.create(neg.path.train, recursive=TRUE)
	pos.path.val = paste0(cur.path, "val/1/")
	if(!dir.exists(pos.path.val)) dir.create(pos.path.val, recursive=TRUE)
	neg.path.val = paste0(cur.path, "val/2/")
	if(!dir.exists(neg.path.val)) dir.create(neg.path.val, recursive=TRUE)

	setwd(pos.path.train)
	for(i in which(y.train==1)) writeMat(paste0(i, ".mat"), Y1=X.train[,,i])
	setwd(neg.path.train)
	for(i in which(y.train==2)) writeMat(paste0(i, ".mat"), Y1=X.train[,,i])
	setwd(pos.path.val)
	for(i in which(y.test==1)) writeMat(paste0(i, ".mat"), Y1=X.test[,,i])
	setwd(neg.path.val)
	for(i in which(y.test==2)) writeMat(paste0(i, ".mat"), Y1=X.test[,,i])

	# save the oracle classifier
	setwd(cur.path)
	if( setting=="IP" ){
		prob.star = sigmoid(2*sin(data.list$linear.star[test.idx]*pi)+data.list$linear.star[test.idx])
	}else{
		prob.star = sigmoid(data.list$linear.star[test.idx])
	}
	oracle.class = as.numeric(prob.star>0.5)
	# classification performance
	AUC.oracle = pROC::roc(as.vector(y.test-1), oracle.class, quiet=TRUE)$auc
	accuracy.oracle = mean(oracle.class == (y.test-1))
	sen.oracle = sensitivity((y.test-1), oracle.class)
	spe.oracle = specificity((y.test-1), oracle.class)
	res = c(accuracy.oracle, AUC.oracle, sen.oracle, spe.oracle)
	names(res) = c("Accuracy.oracle", "AUC.oracle", "Sensitivity.oracle", "Specificity.oracle")
	convergence = res.cur[length(res.cur)]
	res = c(res, convergence)
	write.csv(t(res), file=csv.file)

	setwd(path)
}





#############################################################################

args = commandArgs(T)
#args = c("null", "null", "500", "5", "1")

# set directory, should be consistent with 'simulation_manifold.R', 'stat_table.R' and './spdnet/task.py'
main.path = "/home/linyn/LR_on_Met/simulations/"
setwd(main.path)


library(CovTools)
library(caret)
library(numDeriv)
library(MatrixManifold)
library(pROC)
library(doSNOW)
library(parallel)
library(expm)
library(kernlab)
source("symmetric.matrix.R")
source("matrix.utility.R")
source("affine.invariant.R")
source("utilities_manifold.R")

res.path = paste0(main.path, "results/")
if( !dir.exists(res.path) ) dir.create(res.path, recursive=TRUE)
setwd(res.path)



mu.types = args[1]
if( mu.types=="null" ) mu.types = c("diag", "AR1")
beta.types = args[2]
if( beta.types=="null" ) beta.types=c("diag", "AR1")
R = as.numeric(args[3])
ncore = as.numeric(args[4])
metrics = "LogCholesky"
if( metrics=="null" ) metrics = c("LogCholesky", "LogEuclidean", "AffineInvariant")
D = 3
case = args[5]
if( case=="null" ){
	settings = c("logit", "IP", "additive")
}else{
	case = as.numeric(case)
	if( case==1 ){
		settings = "logit"
	}else if( case==2 ){
		settings = "IP"
	}else if( case==3 ){
		settings = "additive"
	}
}
mat.save = TRUE
matdata.only = FALSE
if( matdata.only ) mat.save = TRUE

beta.normalize = FALSE
test.ratio = 0.2
beta.0 = 0 # the intercept (scalar) in LR model
div = 50

Ns.total = c(100, 500)
vars = c(1, 4)
Ms = c(3)						# dim of matrix
grad.type = "numerical"
grad.method = "simple"	#"simple", "Richardson"

save.list = c("m", "var", "results", "stats.knn", "stats.ksvm", "stats.kdc", "mu.type", "R", "beta.normalize", 
		"test.ratio", "mfd", "beta.0", "beta.type", "setting", "Ns", "M")

#setting=settings[1];beta.type=beta.types[1];mu.type=mu.types[1];metric=metrics[1];m=Ms[1];var=vars[1];n=100;seed=1
for( setting in settings )
{
for( beta.type in beta.types )
{
for( mu.type in mu.types )
{
	for( metric in metrics )
	{
		if( metric=="Frobenius" ){
			manifold = 'sym'
		}else if( metric=="LogEuclidean" | metric=="LogCholesky" | metric=="AffineInvariant" ){
			manifold = 'spd'
		}else{
			stop("Manifold or metric not matched.")
		}

		for(m in Ms)
		{
			mfd = matrix.manifold(manifold=manifold, metric=metric, dim=m)
			if( mu.type=="diag" ){
				mu = diag(m) # the true frechet mean of Xs
			}else if( mu.type=="random" ){
				set.seed(11)
				mu = rmatrix(mfd, sig=0.25)
			}else if( mu.type=="AR1" ){
				mu = AR_cor(m, rho=0.5)
			}
			
			for(var in vars)
			{
				cur.name = paste0("LR-", setting, "-mfd_", mfd$manifold, "-metric_", mfd$metric, "-R", R, "-m", m, "-var", var, "-mu_", mu.type, "-beta_", beta.type, "-normalize_", beta.normalize, "-GradType_", ifelse(grad.type=="numerical", grad.method, grad.type))
				cur.name.noR = paste0("LR-", setting, "-mfd_", mfd$manifold, "-metric_", mfd$metric, "-m", m, "-var", var, "-mu_", mu.type, "-beta_", beta.type, "-normalize_", beta.normalize, "-GradType_", ifelse(grad.type=="numerical", grad.method, grad.type))
				Rfile = paste0("Compared_methods-", cur.name, ".RData")

				if( !matdata.only ){
					knn.csv = paste0("knn-", cur.name, ".csv")
					ksvm.csv = paste0("ksvm-", cur.name, ".csv")
					kdc.csv = paste0("kdc-", cur.name, ".csv")
					if( file.exists(knn.csv) & file.exists(ksvm.csv) & file.exists(kdc.csv) ) next
				}

				stats.knn = NULL
				stats.ksvm = NULL
				stats.kdc = NULL
				Ns = Ns.total

				#start = proc.time()
				for(n in Ns)
				{
					tmp.res.path = paste0(res.path, "tmps/", cur.name.noR, "-n", n, "/")
					if( !dir.exists(tmp.res.path) ) dir.create(tmp.res.path, recursive=TRUE)

					setwd(tmp.res.path)
					tmp.files = dir(pattern="tmp_[0-9]")
					M = min(R+1000, length(tmp.files))
					print(M)
					
					cat(sprintf("setting: %s, mu.type: %s, beta.type: %s, metric: %s, m: %d, n: %d, var: %.2f\n", setting, mu.type, beta.type, metric, m, n, var))

					if( mat.save ){
					print("Create Matlab Data")

					require(R.matlab)
					matlab.path = paste0(res.path, "tmps_matlab/", cur.name.noR, "-n", n, "/")
					if( !dir.exists(matlab.path) ) dir.create(matlab.path, recursive=TRUE)

					ncore = min(ncore,detectCores())
					cl = snow::makeCluster(ncore, type="SOCK")  
					snow::clusterExport(cl, list=c("output_matlab", "sigmoid", "sensitivity", "specificity"))
					registerDoSNOW(cl)   
					
					progress = function(r) if( r %% div == 0 ) cat(sprintf("task %d is completed.\n", r))
					opts = list(progress=progress)
					start = proc.time()
					result = foreach(i=1:M,
							.packages=c('pROC', 'base', 'MatrixManifold', 'R.matlab'),
							.options.snow=opts) %dopar% {
								output_matlab(matlab.path, tmp.res.path, seed=i, test.ratio=test.ratio, setting=setting)
							}
					snow::stopCluster(cl)
					end = proc.time()
					print(end-start)
					}

					if( !matdata.only ){
					print("Compared Methods")
					setwd(tmp.res.path)
					fun.list = as.vector(lsf.str())
					ncore = min(ncore,detectCores())
					cl = snow::makeCluster(ncore, type="SOCK")  
					snow::clusterExport(cl, list=fun.list)
					registerDoSNOW(cl)   
					
					progress = function(r) if( r %% div == 0 ) cat(sprintf("task %d is completed.\n", r))
					opts = list(progress=progress)

					start = proc.time()
					results = foreach(i=1:M,
							.packages=c('pROC', 'base', 'MatrixManifold', 'numDeriv', 'expm', 'caret', 'kernlab'),
							.options.snow=opts) %dopar% {
								method_single(mfd=mfd, seed=i, test.ratio=test.ratio, setting=setting)
							}
					snow::stopCluster(cl)

					res.knn = t(sapply(results, function(res)res[1,]))
					res.ksvm = t(sapply(results, function(res)res[2,]))
					res.kdc = t(sapply(results, function(res)res[3,]))
					means.knn = colMeans(res.knn[which(res.knn[,ncol(res.knn)]==1),])
					means.ksvm = colMeans(res.ksvm[which(res.ksvm[,ncol(res.ksvm)]==1),])
					means.kdc = colMeans(res.kdc[which(res.kdc[,ncol(res.kdc)]==1),])
					sds.knn = apply(res.knn[which(res.knn[,ncol(res.knn)]==1),], 2, sd)
					sds.ksvm = apply(res.ksvm[which(res.ksvm[,ncol(res.ksvm)]==1),], 2, sd)
					sds.kdc = apply(res.kdc[which(res.kdc[,ncol(res.kdc)]==1),], 2, sd)
					stats.knn = rbind(stats.knn, sapply(1:length(means.knn), function(i)paste0(means.knn[i], "(", sds.knn[i], ")")))
					stats.ksvm = rbind(stats.ksvm, sapply(1:length(means.ksvm), function(i)paste0(means.ksvm[i], "(", sds.ksvm[i], ")")))
					stats.kdc = rbind(stats.kdc, sapply(1:length(means.kdc), function(i)paste0(means.kdc[i], "(", sds.kdc[i], ")")))

					end = proc.time()
					print(end-start)
					}
				}# end of n

				if( !matdata.only ){
					setwd(res.path)

					save(list=save.list, file=Rfile)
					rownames(stats.knn) = Ns.total
					write.csv(stats.knn, file=knn.csv)
					rownames(stats.ksvm) = Ns.total
					write.csv(stats.ksvm, file=ksvm.csv)
					rownames(stats.kdc) = Ns.total
					write.csv(stats.kdc, file=kdc.csv)
				}
			}# end of var
		}# end of m
	}# end of metrics
}# end of mu.types
}# end of beta.types
}# end of settings



