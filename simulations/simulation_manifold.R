

args = commandArgs(T)
#args = c("null", "null", "10", "10", "1")

# set directory, should be consistent with 'compared_methods.R', 'stat_table.R' and './spdnet/task.py'
main.path = "/home/linyn/LR_on_Met/simulations/"
setwd(main.path)

library(CovTools)
library(numDeriv)
library(MatrixManifold)
library(pROC)
library(parallel)
library(expm)
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


test.ratio = 0.2
beta.0 = 0 # the intercept (scalar) in LR model
beta.normalize = FALSE
start.R = NULL
end.R = NULL

Ns.total = c(100, 500)
vars = c(1, 4)
Ms = c(3)						# dim of matrix
grad.type = "numerical"
grad.method = "simple"	#"simple", "Richardson"
fun.list = as.vector(lsf.str())

save_list = c("m", "var", "n", "res.list", "stats", "mu.type", "R", "beta.normalize", 
		"test.ratio", "mfd", "beta.0", "beta.type", "setting")

#setting=settings[1];beta.type=beta.types[1];mu.type=mu.types[1];metric=metrics[1];m=Ms[1];var=vars[1];n=Ns[1]
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

				csv.name = paste0(cur.name, ".csv")
				if( file.exists(csv.name) ){
					#save(list=save_list, file=Rfile)
					next
				}

				stats = NULL
				res.list = list()
				Rfile = paste0(cur.name, ".RData")
				if( file.exists(Rfile) ){
					load(Rfile)
					Ns = Ns.total[-c(1:length(res.list))]
				}else{
					Ns = Ns.total
				}

				start = proc.time()
				for(n in Ns)
				{
					tmp.res.path = paste0(res.path, "tmps/", cur.name.noR, "-n", n, "/")
					if( !dir.exists(tmp.res.path) ) dir.create(tmp.res.path, recursive=TRUE)

					R.cur = R
					while( R.cur<=2*R )
					{
						cat(sprintf("setting: %s, mu.type: %s, beta.type: %s, metric: %s, m: %d, n: %d, var: %.2f, R.cur: %d\n", setting, mu.type, beta.type, metric, m, n, var, R.cur))
						res = simu(n, mfd, test.ratio=test.ratio, mu=mu, var=var, R=R.cur, ncore=ncore, grad.method=grad.method, 
							beta.normalize=beta.normalize, grad.type=grad.type, tmp.path=tmp.res.path, fun.list=fun.list,
							start.R=start.R, end.R=end.R, beta.type=beta.type, D=D, setting=setting)
						converged.num = sum(res[,ncol(res)])
						cat(sprintf("Converged: %d.\n", converged.num))
						if( converged.num>=R | !is.null(end.R) ) break

						R.cur = R.cur + min((R-converged.num)*2, 50)
					}

					if( is.null(end.R) ){
						res.list = c(res.list, list(res))
						stats = rbind(stats, out_stat(res, R, outlier.drop=FALSE, only.converged=TRUE))
						save(list=save_list, file=Rfile)
					}
					end = proc.time()
					print(end-start)
				}# end of n
				if( is.null(end.R) ){
					rownames(stats) = Ns.total
					names(res.list) = Ns.total
					save(list=save_list, file=Rfile)
		
					write.csv(stats, file=csv.name)
				}
			}# end of var
		}# end of m
	}# end of metrics
}# end of mu.types
}# end of beta.types
}# end of settings


