
permut_test_ind = function(mfd, X, Y, mu.hat=NULL, seed=1, alpha.GD=0.8, thr.GD.abs=1e-3, 
			thr.GD.rel=1e-6, num.iterations=600, trace=FALSE, type="numerical", 
			grad.method="simple", beta=diag(1, mfd$dim[1]))
{
	n = length(Y)
	if( is.null(seed) ){
		idx.b = 1:n
	}else{
		set.seed(seed*111)
		idx.b = sample(n, size=n)
	}
	X.b = X[,,idx.b]
	Y.b = Y

	opt.GD = gradient_descent(mfd=mfd, X=X.b, y=Y.b, alpha=alpha.GD, 
					threshold.abs=thr.GD.abs, threshold.rel=thr.GD.rel, 
					num.iterations=num.iterations, trace=trace, type=type,
					grad.method=grad.method, beta=beta)
	beta.hat.b = opt.GD$beta.hat
	hs.b = apply(X.b, MARGIN=3, h, mfd=mfd, mu=mu.hat, beta=beta.hat.b)
	L.b = sum(Y.b*hs.b - log(1+exp(hs.b)))	# likelihood
	G.b = -2*L.b	# deviance

	return(G.b)
}

#############################################################################

# set directory, should be consistent with 'real_analysis.R', 'compared_methods.R'
main.path = "/home/linyn/LR_on_Met/real_data/fMRI/"
setwd(main.path)


library(MatrixManifold)
library(doSNOW)
library(parallel)
source("symmetric.matrix.R")
source("matrix.utility.R")
source("affine.invariant.R")
source("utilities_manifold.R")

# load data
load("fMRI_data.RData")
B = 1000
num.iterations=2
ncore = 20
iters = 1:B
div = 50

# prepare data
X = RRCM
Y = y
n = length(Y)
mfd = matrix.manifold('spd', metric="LogCholesky", dim=c(8,8))
mu.hat = frechet.mean(mfd, RRCM)


# overall measure
print("Compute Overall Deviance")
X.0 = X
Y.0 = Y
G.0 = permut_test_ind(mfd, X.0, Y.0, mu.hat, seed=NULL, num.iterations=num.iterations)

# permuated measures
print("Start Permutations")
fun.list = as.vector(lsf.str())
ncore = min(ncore,detectCores())
cl = snow::makeCluster(ncore, type="SOCK")  
snow::clusterExport(cl, list=fun.list)
registerDoSNOW(cl)  

progress = function(r) if( r %% div == 0 ) cat(sprintf("task %d is completed.\n", r))
opts = list(progress=progress)
result = foreach(i=iters,
		.packages=c('pROC', 'base', 'MatrixManifold', 'numDeriv', 'expm'),
		.options.snow=opts) %dopar% {
			permut_test_ind(mfd, X.0, Y.0, mu.hat, seed=i, num.iterations=num.iterations)
		}
snow::stopCluster(cl)
print("Permutation finished.")
Gs = do.call(c, result)
p.val = mean(Gs<G.0)
print(p.val)

# plot Figure S1 in Supplementary Material
pdf("nullity_test.pdf")
C = 2
par(mar=par()$mar+c(0, 1, 0, 0))
plot(Gs, pch=16, xlab="b", ylab="Deviance", ylim=c(min(Gs, G.0)*0.9, max(Gs)), cex.lab=C, xaxt='n', yaxt='n')
abline(h=G.0, lty=2, lwd=C+1)
box(lwd=C)
axis(2, at=seq(300, 550, by=50), cex.axis=C, lwd.ticks=C)
axis(1, at=c(1, seq(0, 1000, by=200)[-1]), labels=FALSE, lwd.ticks=C, cex.axis=C)
axis(1, at=c(1, seq(0, 1000, by=200)[-1]), cex.axis=C, line=0.2, tick=FALSE)
dev.off()



