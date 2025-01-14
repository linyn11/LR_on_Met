

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


knn_single = function(mfd, X.train, y.train, X.test, y.test, Ks=seq(3, 51, by=2), seed=0, 
		cv.crt="ACC", fold=5)
{
	cv.res = knn_cv(mfd, X.train, y.train, Ks=Ks, fold=fold, seed=seed, cv.crt=cv.crt)
	K.cv = cv.res$Opt.K
	y.pred = apply(X.test, 3, function(x)knn_class(mfd, X.train, y.train, x, k=K.cv))	# make prediction
	
	# evaluations
	evals = eval_class_method(y.test, y.pred)
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
	acc = mean(y.val==y.pred)
	return(acc)
}


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


ksvm_single = function(mfd, X.train, y.train, X.test, y.test, Gammas=10^seq(-3, 3, by=1), seed=0, fold=5)
{
	cv.res = ksvm_cv(mfd, X.train, y.train, Gammas=Gammas, fold=fold, seed=seed)
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
	evals = eval_class_method(y.test, y.pred)
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
kdc_single = function(mfd, X.train, y.train, X.test, y.test, Hs=10^seq(-3, 3, by=1), seed=0, fold=5)
{
	cv.res = kdc_cv(mfd, X.train, y.train, Hs=Hs, fold=fold, seed=seed)
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
	evals = eval_class_method(y.test, y.pred)
	return(evals)
}


eval_class_method = function(y, y.pred, th=0.5)
{
	# classification performance
	AUC = pROC::roc(as.vector(y), y.pred, quiet=TRUE)$auc
	accuracy = mean(y.pred == y)
	sen = sensitivity(y, y.pred)
	spe = specificity(y, y.pred)

	res = c(accuracy, AUC, sen, spe)
	names(res) = c("Accuracy", "AUC", "Sensitivity", "Specificity")
	
	return(res)
}


method_single = function(mfd, X, y, seed=0, test.ratio=0.2, cv.crt="ACC")
{
	m = mfd$dim[1]
	n = length(y)
	set.seed(seed)
	test.idx = sample(n, size=round(test.ratio*n))

	X.train = X[,,-test.idx]
	y.train = y[-test.idx]
	X.test = X[,,test.idx]
	y.test = y[test.idx]

	# knn
	res.knn = knn_single(mfd, X.train, y.train, X.test, y.test, seed=seed*100, cv.crt=cv.crt)
	# kernel svm
	res.ksvm = ksvm_single(mfd, X.train, y.train, X.test, y.test, seed=seed*100)
	# kernel density classifier
	res.kdc = kdc_single(mfd, X.train, y.train, X.test, y.test, seed=seed*100)

	res = rbind(res.knn, res.ksvm, res.kdc)
	rownames(res) = c("KNN", "KSVM", "KDC")
	return(res)
}


#############################################################################

# set directory, should be consistent with 'real_analysis.R', 'nullity_test.R'
main.path = "/home/linyn/LR_on_Met/real_data/fMRI/"
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


load("fMRI_data.RData")
M = 10	# fold of cross-validation
ncore = 10	# number of CPU cores to use for 'foreach' function

manifold = "spd"
metric = "LogCholesky"	#"LogEuclidean" | "LogCholesky" | "AffineInvariant"
mfd = matrix.manifold(manifold=manifold, metric=metric, dim=dim(RRCM)[1])
fun.list = as.vector(lsf.str())
ncore = min(ncore,detectCores())
cl = snow::makeCluster(ncore, type="SOCK")  
snow::clusterExport(cl, list=fun.list)
registerDoSNOW(cl)  

progress = function(r) if( r %% 1 == 0 ) cat(sprintf("task %d is completed.\n", r))
opts = list(progress=progress)
result = foreach(i=1:M,
		.packages=c('pROC', 'base', 'MatrixManifold', 'numDeriv', 'expm', 'caret', 'kernlab'),
		.options.snow=opts) %dopar% {
			method_single(mfd=mfd, X=RRCM, y=y, seed=i*111, test.ratio=1/M)
		}
snow::stopCluster(cl)

res.knn = t(sapply(result, function(res)res[1,]))
res.ksvm = t(sapply(result, function(res)res[2,]))
res.kdc = t(sapply(result, function(res)res[3,]))
means.knn = colMeans(res.knn)
means.ksvm = colMeans(res.ksvm)
means.kdc = colMeans(res.kdc)
res.others = rbind(means.knn, means.ksvm, means.kdc)

res.LRs = read.csv("realdata_res_LR.csv", header=TRUE)[,-1]
res.spdnet = read.table(paste0(main.path, "spdnet/CV_avg.txt"))
res.spdnet = as.numeric(sapply(unlist(res.spdnet[c(1:4*2)]), function(x)unlist(strsplit(x, split=","))[1]))
res = rbind(res.spdnet, res.others, data.matrix(res.LRs))
rownames(res) = c("spdnet", "knn", "ksvm", "kdc", "LR-BOLD", "proposed")
write.csv(res, file="realdata_methods_comparisons.csv")



