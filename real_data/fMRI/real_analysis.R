

fit_logistic = function(X.ori, y, test.ratio=0.2, seed=0)
{
	# extract diag as input
	X = t(apply(X.ori, MARGIN=3, function(x){diag(x)}))
	data = data.frame(y, X)

	n = length(y)
	set.seed(seed)
	test.idx = sample(n, size=round(test.ratio*n))
	data.train = data[-test.idx,]
	data.test = data[test.idx,]

	fit = glm(y~., family="binomial", data=data.train)
	pred.linear = predict(fit, newdata=data.test, type="link")
	#y.pred = as.numeric(prob.pred>=0.5)

	res.diag = eval_class(data.test$y, pred.linear)
	names(res.diag) = c("ACC", "AUC", "SEN", "SPE")

	# lower-triangular entries
	X = t(apply(X.ori, MARGIN=3, function(x){x[lower.tri(x, diag=T)]}))
	data = data.frame(y, X)

	n = length(y)
	set.seed(seed)
	test.idx = sample(n, size=round(test.ratio*n))
	data.train = data[-test.idx,]
	data.test = data[test.idx,]

	fit = glm(y~., family="binomial", data=data.train)
	pred.linear = predict(fit, newdata=data.test, type="link")
	#y.pred = as.numeric(prob.pred>=0.5)

	res.lower = eval_class(data.test$y, pred.linear)
	names(res.lower) = c("ACC", "AUC", "SEN", "SPE")

	res = rbind(res.diag, res.lower)
	return(res)
}


output_matlab = function(X, y, path, seed=0, test.ratio=0.2)
{
	require(R.matlab)

	y = y + 1
	n = length(y)
	set.seed(seed)
	test.idx = sample(n, size=round(test.ratio*n))
	X.train = X[,,-test.idx]
	y.train = y[-test.idx]
	X.test = X[,,test.idx]
	y.test = y[test.idx]

	cur.path = paste0(path, "matlab-data/CV-", seed/111, "/")
	if(!dir.exists(cur.path)) dir.create(cur.path, recursive=TRUE)

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
	setwd(path)
}

#############################################################################

# set directory, should be consistent with 'compared_methods.R', 'nullity_test.R'
main.path = "C:/Users/Administrator/Desktop/LR_on_Met/real_data/fMRI/"
setwd(main.path)


M = 10	# fold of cross-validation
ncore = 10	# number of CPU cores to use for 'foreach' function

# load dependency
#devtools::install_github("linulysses/matrix-manifold")
library(MatrixManifold)
library(doSNOW)
library(parallel)
source("symmetric.matrix.R")
source("matrix.utility.R")
source("affine.invariant.R")
source("utilities_manifold.R")

# load data
load("fMRI_data.RData")


#####################################################################################
## for classification performance comparisons
## proposed classifier and Logistic-BOLD
#####################################################################################

# prepare matlab-type data for the SPDNet
# this procedure is time-consuming
for(i in 1:M) output_matlab(RRCM, y, main.path, seed=i*111, test.ratio=1/M)

# manifold optimization for the proposed classifier
# this procedure would cost a lot of time
num.iterations = 300
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
		.packages=c('pROC', 'base', 'MatrixManifold', 'numDeriv', 'expm'),
		.options.snow=opts) %dopar% {
			fit_eval(mfd=mfd, X=RRCM, y=y, seed=i*111, num.iterations=num.iterations, trace=TRUE, test.ratio=1/M)
		}
snow::stopCluster(cl)
result.eval = lapply(result, function(x)x$res.class[1:4])
res.eval = do.call(rbind, result.eval)
res.mean.proposed = colMeans(res.eval)


# logistic-BOLD
res.LR.BOLD = NULL
for(i in 1:M) 
{
	res.cur = fit_logistic(X.ori=RRCM, y, seed=i*111, test.ratio=1/M)
	res.LR.BOLD = rbind(res.LR.BOLD, res.cur[1,])
}	
res.mean.LR.BOLD = colMeans(res.LR.BOLD)

res = rbind(res.mean.LR.BOLD, res.mean.proposed)
write.csv(res, file="realdata_res_LR.csv")



#####################################################################################
# for interpolations
#####################################################################################
library(corrplot)

# fit an overall model to get beta.hat
manifold = "spd"
metric = "LogCholesky"	#"LogEuclidean" | "LogCholesky" | "AffineInvariant"
mfd = matrix.manifold(manifold=manifold, metric=metric, dim=dim(RRCM)[1])
#X=RRCM;y=y;alpha=0.8;threshold.abs=1e-2;threshold.rel=1e-4;num.iterations=300;trace=TRUE;type="numerical";grad.method="simple"
opt.GD = gradient_descent(mfd=mfd, X=RRCM, y=y, alpha=0.7, 
				threshold.abs=1e-2, threshold.rel=1e-4, 
				num.iterations=100, trace=TRUE, type="numerical",
				grad.method="simple")
mu.hat = frechet.mean(mfd, RRCM)
beta.hat = opt.GD$beta.hat


#load("realdata_overall_8regions_mot-lang_LogCholesky-LogCholesky-P16-M10-It300-StartBlock0-sliding_TRUE.RData")

# ROC curve 
library(pROC)
fitted.linear = apply(RRCM, MARGIN=3, h, mfd=mfd, mu=mu.hat, beta=beta.hat)
roc.lr = roc(y, sigmoid(fitted.linear))
pdf("ROC-lr.pdf")
plot(roc.lr, legacy.axes=TRUE, cex.lab=2, cex.axis=2, lwd=3, mgp=c(3, 1, 0), mar=c(4, 4.5, 2, 2)+.3)
dev.off()

# geodesic connecting mu and beta
T.th = 0.03	# need to choose a suitable gap T.th to get a better virtualization
Ts = seq(from=-T.th, to=T.th, length.out=9)
geod.mu.beta.hat = lapply(Ts, function(t)geodesic(mfd, p=mu.hat, u=beta.hat, t=t))

# compute odds, log-odds and odds ratio
odds = round(sapply(geod.mu.beta.hat, function(X)exp(h(mfd, mu.hat, beta.hat, X))), 4)
log.odds = round(log(sapply(geod.mu.beta.hat, function(X)exp(h(mfd, mu.hat, beta.hat, X)))), 4)
odds.ratio = round(sapply(geod.mu.beta.hat, function(X)exp(h(mfd, mu.hat, beta.hat, X) / geo.dist(mfd, mu.hat, X))), 3)
print(odds)
print(log.odds)

# prepare to plot
C = mean(unlist(geod.mu.beta.hat))
for(i in 1:length(Ts)) geod.mu.beta.hat[[i]] = geod.mu.beta.hat[[i]]-C
col.neg <- colorRampPalette(c("#0033FF", "#FFFFCC", "#FF3333"))(50)
col.pos <- colorRampPalette(c("#FFFFFF", "#BB4444"))(50)
max.cut = quantile(unlist(geod.mu.beta.hat), probs=c(0.98))
for(i in 1:length(geod.mu.beta.hat))
{
	tmp = geod.mu.beta.hat[[i]]
	tmp[tmp>max.cut]=max.cut
	geod.mu.beta.hat[[i]] = tmp
}
ranges = range(geod.mu.beta.hat)

s = 1
pdf(file="geodesic_mu_beta_hat.pdf", width=length(Ts)*1.7*s, height=2*s)
S = 1
layout.mat = NULL
for(i in 1:length(Ts)) layout.mat = cbind(layout.mat, matrix(i, nrow=S, ncol=S))
layout(mat=layout.mat, heights=c(4, 1), widths=rep(1, ncol(layout.mat)))
for(i in 1:length(Ts))
{
	par(mar=c(1.1, 1, 0, 0))
	if( min(geod.mu.beta.hat[[i]])<0 ){
		corrplot(geod.mu.beta.hat[[i]], method="color", is.corr=FALSE, tl.col="black", cl.pos="n", col=col.neg, col.lim=ranges)
	}else{
		corrplot(geod.mu.beta.hat[[i]], method="color", is.corr=FALSE, tl.col="black", cl.pos="n", col=col.pos, col.lim=ranges)
	}
}
dev.off()


pdf(file="legend_geodesic_mu_beta_hat.pdf", width=3*1.7*s, height=1)
par(mar=c(0, 0, 0, 0))
plot.new()
at = round(seq(from=ranges[1]+C, to=ranges[2]+C, length.out=10), 1)
colorlegend(col.neg, at, align = 'l', cex = 1, xlim = c(0, 1),
              ylim = c(0-0.1, 1-0.1), vertical = FALSE)
dev.off()


