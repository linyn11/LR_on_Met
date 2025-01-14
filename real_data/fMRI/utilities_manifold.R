
matrix.manifold <- function(manifold,metric,dim)
{
    # ltpd: lower triangular with positive diagonals
    stopifnot(manifold %in% c('ltpd', 'spd','sym')) #,'stiefel'))
    if(manifold == 'spd')
        stopifnot(metric %in% c('LogCholesky','LogEuclidean','AffineInvariant'))
    if(manifold == 'sym')
        stopifnot(metric %in% c('Frobenius'))
#    if(manifold == 'lt') # lower triangular matrix
#        stopifnot(metric %in% c('Frobenius'))
    if(manifold == 'ltpd') # lower triangular matrix with positive diagonals
        stopifnot(metric %in% c('LogCholesky'))
    
    class.name <- paste0(manifold,'_',metric)
    if(length(dim) == 1) dim <- rep(dim,2)
    mfd <- structure(list(dim=dim,manifold=manifold,metric=metric),class=class.name)
    return(mfd)
}


inner_F = function(X, Y)
{
	X = as.matrix(X)
	Y = as.matrix(Y)
	return(sum( as.vector(X) * as.vector(Y) ))
}


norm_mfd = function(mfd, p, q)
{
	sqrt(rie.metric(mfd, p, q, q))
}


cos_mfd = function(mfd, p, X, Y)
{
	if( norm_mfd(mfd, p, X)*norm_mfd(mfd, p, Y)==0 ){
		# zero matrix
		similarity = 0
	}else{
		similarity = rie.metric(mfd, p, X, Y) / ( norm_mfd(mfd, p, X)*norm_mfd(mfd, p, Y) )
	}

	return( similarity )
}


sigmoid = function(z)
{
	#Sigmoid function for logistic regression
	g = 1/(1+exp(-z))
	return(as.vector(g))
}


l2_norm = function(x)
{
	return( sqrt(sum(x^2)) )
}


sensitivity = function(y, y.hat)
{
	TP = sum(y==1 & y.hat==1)
	FN = sum(y==1 & y.hat==0)
	TP / (TP + FN)
}


specificity = function(y, y.hat)
{
	TN = sum(y==0 & y.hat==0)
	FP = sum(y==0 & y.hat==1)
	TN / (TN+FP)
}


eval_class = function(y, fitted.linear, th=0.5)
{
	predicted.prob = sigmoid(fitted.linear)

	# classification performance
	AUC.prob = pROC::roc(as.vector(y), predicted.prob, quiet=TRUE)$auc

	predicted.class = as.numeric(predicted.prob>th)
	accuracy = mean(predicted.class == y)
	sen = sensitivity(y, predicted.class)
	spe = specificity(y, predicted.class)

	res = c(accuracy, AUC.prob, sen, spe)
	return(res)
}

#dig=4;alpha.GD=0.8;thr.GD.abs=1e-3;thr.GD.rel=1e-6;type="numerical";grad.method="simple";beta=diag(1, mfd$dim[1])
fit_eval = function(mfd, X, y, test.ratio=0.2, seed=0, dig=4, alpha.GD=0.8, thr.GD.abs=1e-3, thr.GD.rel=1e-6, num.iterations=600,
			trace=FALSE, type="numerical", grad.method="simple", beta=diag(1, mfd$dim[1]))
{
	n = length(y)
	set.seed(seed)
	test.idx = sample(n, size=round(test.ratio*n))
	X.train = X[,,-test.idx]
	y.train = y[-test.idx]
	mu.hat = frechet.mean(mfd, X.train)
	X.test = X[,,test.idx]
	y.test = y[test.idx]

	# evaluation
	opt.GD = gradient_descent(mfd=mfd, X=X.train, y=y.train, alpha=alpha.GD, 
					threshold.abs=thr.GD.abs, threshold.rel=thr.GD.rel, 
					num.iterations=num.iterations, trace=trace, type=type,
					grad.method=grad.method, beta=beta)
	beta.hat = opt.GD$beta.hat
	fitted.linear = apply(X.test, MARGIN=3, h, mfd=mfd, mu=mu.hat, beta=beta.hat)
	converged = opt.GD$converged

	res.class = eval_class(y.test, fitted.linear)
	res = c(res.class, converged)
	names(res) = c("ACC", "AUC", "SEN", "SPE", "Converged")
	res.list = list(res.class=res, opt.GD=opt.GD)
	#save.list = c("res.list")
	#save(list=save.list, file=paste0(mfd$metric, "_", seed, ".RData"))
	return(res.list)
}

#error.num=-9999;step.type="constant";beta=diag(1, mfd$dim[1])
gradient_descent = function(mfd, X, y, beta=diag(1, mfd$dim[1]), alpha=0.8, num.iterations=600, 
				threshold.abs=1e-3, threshold.rel=1e-6, step.type="constant", trace=FALSE,
				error.num=-9999, type="numerical", grad.method="simple") 
{
	# Look at the values over each iteration
	beta.path = list(beta)
	mu.hat = frechet.mean(mfd, X)

	grad.norms <- beta.diffs <- grad.rel.changes <- NULL
	converged = FALSE
	for(t in 1:num.iterations)
	{
		if(t %% (num.iterations/10) == 0 & trace) print(sprintf('t: %d',t))

		# negative gradient of ell
#		if( mfd$metric=="AffineInvariant" ){
		if( FALSE ){
			E.beta = NULL
			type = "analytic"
		}else{
			E.beta = tryCatch(
				{
					generate_E(mfd$dim[1], mfd, beta)
				},
				error=function(cond){
					return(error.num)
				})
			if( identical(E.beta, error.num) ) break
		}
		grad  = -grad_ell_beta(mfd, beta, mu.hat, X, y, E.beta, type=type, grad.method=grad.method)
		grad.norms = c(grad.norms, sqrt(rie.metric(mfd, p=beta, grad, grad)))
		beta = rie.exp(mfd, p=beta, v=alpha*grad)
		check.beta = tryCatch(
			{
				chol(beta)
			},
			error=function(cond){
				return(error.num)
			})
		if( identical(check.beta, error.num) ){
			beta = beta.path[[t]]
			break
		}
		if( step.type=="shrink" ) alpha = alpha/t

		beta.path = c(beta.path, list(beta))
		beta.diffs = c( beta.diffs, geo.dist(mfd, beta.path[[t]], beta) )

		if( t>1 ){
			grad.rel.change = (grad.norms[t] - grad.norms[t-1]) / grad.norms[t-1]
			grad.rel.changes = c(grad.rel.changes, grad.rel.change)

#			if( grad.norms[t]<threshold.abs ){
			if( grad.norms[t]<threshold.abs | abs(grad.rel.change)<threshold.rel ){
				if( grad.rel.change<0 )	converged = TRUE
				break 
			}
		}
	}

	res = list(converged=converged, beta.hat=beta, beta.path=beta.path, grad.norms=grad.norms, 
			beta.diffs=beta.diffs, grad.rel.changes=grad.rel.changes)
	return(res)
}


normalize = function(mfd, mu, beta, dis=1, eps=1e-10, left=0, right=1)
{
	# move positions to make sure right position is larger than dis, and left position is smaller than dis
	while( TRUE )
	{
		beta.n = geodesic(mfd, p=mu, u=beta, right)
		dis.beta = geo.dist(mfd, mu, beta.n)
	
		if( dis.beta < dis ){
			left = right
			right = 2*right
		}else{
			break
		}
	}
	
	# use binary search to find the normalzed beta within the precision
	while( abs(dis.beta - dis) >= eps )
	{
		mid = (left + right) / 2
	
		beta.n = geodesic(mfd, p=mu, u=beta, t=mid)
		dis.beta = geo.dist(mfd, mu, beta.n)
		#print(sprintf("%.10f", dis.beta))
	
		if( dis.beta > dis ){
			right = mid
		}else{
			left = mid
		}
	}

	return(beta.n)
}


generate_eij = function(m, i, j)
{
	e.ij = matrix(0, m, m)
	e.ij[i,j] = 1
	return(e.ij)
}


generate_xi_tilde = function(m, mfd)
{
	# generate linearly independent matrices in \mathcal{S}_m
	# return as a list, each element is a matrix in \mathcal{S}_m

	if( mfd$manifold != "spd" )  stop("The generation process of Es can only be valid for the SPD manifold.")

	xis.tilde = NULL
	for(i in 1:m)
	{
		for(j in 1:i)
		{
			xi.tilde.ij = generate_eij(m, i, j) + t(generate_eij(m, i, j))
			xis.tilde = c(xis.tilde, as.vector(xi.tilde.ij))
		}
	}

	n = m*(m+1) / 2
	xis.tilde = array(xis.tilde, dim=c(m, m, n))

	return(xis.tilde)
}


generate_E = function(m, mfd, beta)
{
	# generate linearly independent matrices \xis
	E.tilde = generate_xi_tilde(m, mfd)
	n = dim(E.tilde)[3]

	# Gram.Schmidt process for E.tilde
	 # first vector
	xi = E.tilde[,,1]
	xi = xi / sqrt(MatrixManifold::rie.metric(mfd, p=beta, u=xi, v=xi))
	data = NULL
	data = c(data, as.vector(xi))
	E = array(data, dim=c(m, m, 1))
	 # other vectors 
	for( i in 2:n )
	{
		xi = E.tilde[,,i]	
		for( j in 1:(i-1) )
		{
			g.ij = MatrixManifold::rie.metric(mfd, p=beta, u=xi, v=E[,,j])
			g.jj = MatrixManifold::rie.metric(mfd, p=beta, u=E[,,j], v=E[,,j])
			xi = xi - g.ij / g.jj * E[,,j] 
		}
		xi = xi / sqrt(MatrixManifold::rie.metric(mfd, p=beta, u=xi, v=xi))
		data = c(data, as.vector(xi))
		E = array(data, dim=c(m, m, i))
	}

	return(E)
}


grad_ell_beta = function(mfd, beta, mu, X, y, E.beta, type="numerical", grad.method="simple")
{
#	if( mfd$metric=="AffineInvariant" & type=="analytic"  ){
	if( FALSE  ){
		inds = sapply(1:n, function(i)ind_term_AIM(mfd, mu, beta, X[,,i], y[i]))
		grad = mean(inds)
	}else{
		# generalized Fourier coefficients with respect to the basis E
		M = dim(E.beta)[3]
		C = sapply(1:M, function(j)D_ell_beta_xi(mfd, mu, beta, X, y, E.beta[,,j], type=type, grad.method=grad.method))
	
		# gradient of \ell at beta
		grad = matrix(0, nrow=mfd$dim[1], ncol=mfd$dim[2])
		for(j in 1:M) grad = grad + C[j]*E.beta[,,j]
	}

	return(grad)
}



D_ell_beta_xi = function(mfd, mu, beta, X, y, xi, L.mu=NULL, L.beta=NULL, D.beta.xi=NULL, type="numerical", grad.method="simple")
{
	# (X, y) are the data pairs consists of all the single pair (Xi, yi)
	n = length(y)

	if( type=="numerical" ){
		# the likelihood term
		inds = sapply(1:n, function(i)ind_term_numerical(mfd, X[,,i], y[i], xi, mu, beta, method=grad.method))
	}
	dell.0 = mean(inds)

	return(dell.0)
}


ind_term_AIM = function(mfd, mu, beta, X, y)
{
#	grad.inner = ***
	grad.inner = 0
	grad = beta %*% grad.inner %*% beta
	return(grad)
}

#library(numDeriv)
#numDeriv::grad

ind_term_numerical = function(mfd, X, y, xi, mu, beta, method="simple", method.args=list(eps=1e-6))
{
	logit.py = h(mfd, mu, beta, X)
	py = sigmoid(logit.py)

	fun = function(t)
	{
		gamma.t = geodesic(mfd, beta, xi, t)
		h.t = h(mfd, mu, gamma.t, X)	
		h.t
	}

	dh.0 = numDeriv::grad(fun, x=0, method=method, method.args=method.args)
	res = (py - y) * dh.0

	return(res)
}



h = function(mfd, mu, beta, X, beta.0=0)
{
	# intermediates
	Log.mu.X = rie.log(mfd, p=mu, q=X)
	Log.mu.beta = rie.log(mfd, p=mu, q=beta)
	d.X = geo.dist(mfd, mu, X)
	d.beta = geo.dist(mfd, mu, beta)
	cos.beta.X = cos_mfd(mfd, mu, Log.mu.X, Log.mu.beta)

	# relation part
	res = beta.0 + d.X * d.beta * cos.beta.X
	return(res)
}

simu_worker = function(n, mfd, beta.0=0, test.ratio=0.2, mu=diag(mfd$dim[1]), var=1, seed=0, 
			dig=4, beta.normalize=TRUE, trace=FALSE, grad.type="numerical", tmp.path=NULL,
			grad.method="simple")
{
	if( !is.null(tmp.path) ){
		tmp.name = paste0(tmp.path, "tmp_", seed, ".RData")
		if( file.exists(tmp.name) ){
			load(tmp.name)
			return(res.cur)
		}
	}

	if( trace ) print("Data Generating.")
	data.list = data_generation(n, mfd, beta.0=beta.0, mu=mu, var=var, seed=seed, beta.normalize=beta.normalize)
	beta.mu.dist = geo.dist(mfd, data.list$beta.star, mu)
	mu.mu.hat.dist = geo.dist(mfd, frechet.mean(mfd, data.list$X), mu)
	if( trace ) print("Optimizing.")
	res.cur = fit_eval(mfd=mfd, X=data.list$X, y=data.list$y, linear.star=data.list$linear.star, beta.star=data.list$beta.star, 
				mu=mu, beta.0=beta.0, test.ratio=test.ratio, seed=seed, dig=dig, trace=trace, type=grad.type,
				grad.method=grad.method)
	res.cur = c(mu.mu.hat.dist, beta.mu.dist, res.cur)
	
	if( !is.null(tmp.path) ){
		save.list = c("res.cur", "seed")
		save(list=save.list, file=tmp.name)
	}
	return(res.cur)
}


simu = function(n, mfd, beta.0=0, R=100, mu=diag(mfd$dim[1]), var=1, test.ratio=0.2, ncore=1, ndiv=R/5, 
		dig=4, beta.normalize=TRUE, trace=FALSE, grad.type="numerical", tmp.path=NULL, grad.method="simple",
		fun.list=NULL, start.R=NULL, end.R=NULL)
{
	colNames = c("Mu-Mu.hat-Dist", "Beta-Mu-Dist", "Beta-Mu-hat-Dist",
			"Beta-Dist", "Beta-Abs-Diff-Mu", "Beta-Abs-Diff-Muhat",
			"Linear Part-RMSE", "Probability-RMSE", 
			"Accuracy", "Accuracy-diff", 
			"AUC", "AUC-diff", 
			"Sensitivity", "Sensitivity-diff", 
			"Specificity", "Specificity-diff", 
			"Converged")
	if( is.null(end.R) ){
		loop.R = 1:R
	}else{
		loop.R = (start.R+1):end.R
		ndiv = length(loop.R) / 5
	}
	print("loop.R: ")
	print(loop.R)

	if(ncore > 1){
		require(doSNOW)
		ncore = min(ncore,detectCores())
        
		cat(sprintf('ncore=%d is used\n', ncore))
        
		cl = snow::makeCluster(ncore, type="SOCK")  
		snow::clusterExport(cl, list=fun.list)
		registerDoSNOW(cl)  
        
		progress = function(r) if( r %% ndiv == 0 ) cat(sprintf("task %d is completed.\n", r))
		opts = list(progress=progress)
		result = foreach(i=loop.R,
				.packages=c('pROC', 'base', 'MatrixManifold', 'numDeriv', 'expm'),
				.options.snow=opts) %dopar% {
					simu_worker(n, mfd, beta.0=beta.0, test.ratio=test.ratio, mu=mu, var=var, seed=i, 
							dig=dig, beta.normalize=beta.normalize, trace=trace, grad.type=grad.type,
							tmp.path=tmp.path, grad.method=grad.method)
				}
		snow::stopCluster(cl)
		res = do.call(rbind, result)
	}else{
		res = matrix(0, length(loop.R), length(colNames))
		for(i in loop.R)
		{
			print(i)
			tmp = simu_worker(n, mfd, beta.0=beta.0, test.ratio=test.ratio, mu=mu, var=var, seed=i, 
						dig=dig, beta.normalize=beta.normalize, trace=trace, grad.type=grad.type,
						tmp.path=tmp.path, grad.method=grad.method)
			res[i,] = tmp 
		}
	}
	colnames(res) = colNames 

	return(res)
}


rm_outlier = function(res)
{
	# only consider outliers in beta-rmse and linear-part-rmse
	rm.idx = unique(unlist(apply(res[,3:6], 2, function(x){which(x==boxplot.stats(x)$out)})))
	if( length(rm.idx)>0 )  res = res[-rm.idx,]

	return(res)
}


out_stat = function(res, R=NULL, outlier.drop=TRUE, dig=3, only.converged=TRUE)
{
	if( outlier.drop )  res = rm_outlier(res)
	if( only.converged ){
		converged.mean = mean(res[,ncol(res)])
		converged.sd = sd(res[,ncol(res)])
		res = res[which(res[,ncol(res)]==1),]
		converged.num = nrow(res)
	}

	if( is.null(R) ) R = nrow(res)
	res = res[1:R,]

	means = colMeans(res)
	sds = apply(res, 2, sd)
	if( only.converged ){
		means[length(means)] = converged.mean
		sds[length(sds)] = converged.sd
	}
	stats = sapply(1:length(means), function(i){paste0(round(means[i], dig=dig), "(", round(sds[i], dig=dig), ")")})
	names(stats) = names(means)

	if( only.converged ){
		stats = c(stats, paste0(converged.num, "(", R, ")"))
		names(stats)[length(stats)] = "converged.num"
	}

	return(stats)
}


