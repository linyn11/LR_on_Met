
args = commandArgs(T)

# set directory, should be consistent with 'compared_methods.R', 'simulation_manifold.R' and './spdnet/task.py'
main.path = "/home/linyn/LR_on_Met/simulations/"
setwd(main.path)

res.path = paste0(main.path, "results/")
spdnet.path = paste0(main.path, "spdnet/results/")
setwd(res.path)
## simulation task, estimation or classification. 
## For estimation, the case number could be 1, for generating Table 1 in the main paper and Table S1 in Supplementary Material.
## For classification, the case number could be 1, 2, or 3, for generating Table 2 in the main paper and Tables S2-S4 in Supplementary Material.

task = args[1]
case = as.numeric(args[2])	# case number
R = as.numeric(args[3])		# number of replicates
dig = 3				# the digit for output

if( case==1 ){
	setting = "logit"
}else if( case==2 ){
	setting = "IP"
}else if( case==3 ){
	setting = "additive"
}

if( setting=="logit" ){
	beta.types = c("diag", "AR1")
	mu.types = c("diag", "AR1")
}else{
	beta.types = c("AR1")
	mu.types = c("diag")
}
metric = "LogCholesky"
beta.normalize = "FALSE"
vars = c(1, 4)
m = 3
grad.type = "numerical"
grad.method = "simple"
Ns = c(100, 500)
if( task=="classification" ){
	COLs = c(1, 3, 5, 7, 2)
	COL.names = c("Accuracy", "AUC", "Sensitivity", "Specificity", "EER")
}else if( task=="estimation" ){
	COLs = c(2, 8, 11)
	COL.names = c("d(mu.star, mu.hat)", "d(beta.star, beta.hat)", "RMSE")
}

Res = NULL
for( mu.type in mu.types )
{
	for( beta.type in beta.types )
	{
		for( var in vars )
		{
			for( n in Ns )
			{
				cur.name = paste0("LR-", setting, "-mfd_spd-metric_", metric, "-R", R, "-m", m, "-var", var, "-mu_", mu.type, "-beta_", beta.type, "-normalize_", beta.normalize, "-GradType_", ifelse(grad.type=="numerical", grad.method, grad.type))
				setwd(res.path)

				LR.file = paste0(cur.name, ".csv")
				LR.res = read.csv(LR.file)
				row = which(LR.res[,1]==n)

				if( task=="classification" ){
					LR.res = as.character(LR.res[row, 13:20])
					LR.res.mean = as.numeric(sapply(LR.res, function(x)unlist(strsplit(x, split="\\(|\\)"))[1]))	# take means
					LR.res.sd = as.numeric(sapply(LR.res, function(x)unlist(strsplit(x, split="\\(|\\)"))[2]))	# take sds
					LR.res = sapply(COLs, function(i)paste0(sprintf(paste0("%.", dig, "f"), LR.res.mean[i]), "(", sprintf(paste0("%.", dig, "f"), LR.res.sd[i]), ")"))
				}else if( task=="estimation" ){
					LR.res = as.character(LR.res[row, COLs])
					LR.res.mean = as.numeric(sapply(LR.res, function(x)unlist(strsplit(x, split="\\(|\\)"))[1]))	# take means
					LR.res.sd = as.numeric(sapply(LR.res, function(x)unlist(strsplit(x, split="\\(|\\)"))[2]))	# take sds
					LR.res = sapply(1:length(COLs), function(i)paste0(sprintf(paste0("%.", dig, "f"), LR.res.mean[i]), "(", sprintf(paste0("%.", dig, "f"), LR.res.sd[i]), ")"))

					res.cur = matrix(LR.res, nrow=1)
					colnames(res.cur) = COL.names
					rownames = c("proposed")
					rownames = paste0("mu-", mu.type, "_beta-", beta.type, "_r-", var, "_n-", n, "-", rownames)
					rownames(res.cur) = rownames
				}
				
				if( task=="classification" ){
				knn.file = paste0("knn-", cur.name, ".csv")
				ksvm.file = paste0("ksvm-", cur.name, ".csv")
				kdc.file = paste0("kdc-", cur.name, ".csv")
				spdnet.file = paste0("LR-", setting, "-mfd_spd-metric_", metric, "-m", m, "-var", var, "-mu_", mu.type, "-beta_", beta.type, "-normalize_", beta.normalize, "-GradType_", ifelse(grad.type=="numerical", grad.method, grad.type), "-n", n, ".txt")

				knn.res = read.csv(knn.file)
				row = which(knn.res[,1]==n)
				knn.res = knn.res[row,-c(1,10)]
				knn.res.mean = as.numeric(sapply(knn.res, function(x)unlist(strsplit(x, split="\\(|\\)"))[1]))	# take means
				knn.res.sd = as.numeric(sapply(knn.res, function(x)unlist(strsplit(x, split="\\(|\\)"))[2]))	# take sds
				knn.res = sapply(COLs, function(i)paste0(sprintf(paste0("%.", dig, "f"), knn.res.mean[i]), "(", sprintf(paste0("%.", dig, "f"), knn.res.sd[i]), ")"))


				ksvm.res = read.csv(ksvm.file)
				row = which(ksvm.res[,1]==n)
				ksvm.res = ksvm.res[row,-c(1,10)]
				ksvm.res.mean = as.numeric(sapply(ksvm.res, function(x)unlist(strsplit(x, split="\\(|\\)"))[1]))	# take means
				ksvm.res.sd = as.numeric(sapply(ksvm.res, function(x)unlist(strsplit(x, split="\\(|\\)"))[2]))	# take sds
				ksvm.res = sapply(COLs, function(i)paste0(sprintf(paste0("%.", dig, "f"), ksvm.res.mean[i]), "(", sprintf(paste0("%.", dig, "f"), ksvm.res.sd[i]), ")"))
				
				kdc.res = read.csv(kdc.file)
				row = which(kdc.res[,1]==n)
				kdc.res = kdc.res[row,-c(1,10)]
				kdc.res.mean = as.numeric(sapply(kdc.res, function(x)unlist(strsplit(x, split="\\(|\\)"))[1]))	# take means
				kdc.res.sd = as.numeric(sapply(kdc.res, function(x)unlist(strsplit(x, split="\\(|\\)"))[2]))	# take sds
				kdc.res = sapply(COLs, function(i)paste0(sprintf(paste0("%.", dig, "f"), kdc.res.mean[i]), "(", sprintf(paste0("%.", dig, "f"), kdc.res.sd[i]), ")"))

				setwd(spdnet.path)
				tb = read.table(spdnet.file)[,2*(1:8)]
				spdnet.res.mean = as.numeric(sapply(unlist(tb[1,]), function(x)unlist(strsplit(x, split=","))[1]))
				spdnet.res.sd = as.numeric(sapply(unlist(tb[2,]), function(x)unlist(strsplit(x, split=","))[1]))
				spdnet.res = sapply(COLs, function(i)paste0(sprintf(paste0("%.", dig, "f"), spdnet.res.mean[i]), "(", sprintf(paste0("%.", dig, "f"), spdnet.res.sd[i]), ")"))

				res.cur = rbind(spdnet.res, knn.res, ksvm.res, kdc.res, LR.res)
				colnames(res.cur) = COL.names
				rownames = c("spdnet", "knn", "ksvm", "kdc", "proposed")
				rownames = paste0("mu-", mu.type, "_beta-", beta.type, "_r-", var, "_n-", n, "-", rownames)
				rownames(res.cur) = rownames
				}

				Res = rbind(Res, res.cur)
			}
		}
	}
}

setwd(main.path)
csv.file = paste0("simulations-", task, "-Case", case, ".csv")
write.csv(Res, file=csv.file)

