

AR_cor = function(p, rho)
{
	Sigma = diag(0.5, p)
	for(i in 1:(p-1))
	{
		for(j in (i+1):p) Sigma[i,j] = rho^(abs(i-j))
	}
	Sigma = (t(Sigma) + Sigma)

	return(Sigma)
}



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


eval_numeric = function(mfd, beta.hat, beta.star, fitted.linear, linear.star, mu.hat, mu, beta.0=0)
{
	n = length(fitted.linear)

	predicted.prob = sigmoid(fitted.linear)
	prob.star = sigmoid(linear.star)

	beta.mu.hat.dist = geo.dist(mfd, beta.hat, mu.hat)
	beta.dist = geo.dist(mfd, beta.hat, beta.star)
	beta.abs_diff.root_mu = abs(geo.dist(mfd, mu, beta.hat) - geo.dist(mfd, mu, beta.star))
	beta.abs_diff.root_muhat = abs(geo.dist(mfd, mu.hat, beta.hat) - geo.dist(mfd, mu.hat, beta.star))

	linear.rmse = l2_norm(as.vector(fitted.linear - linear.star)/n)
	prob.rmse = l2_norm(as.vector(predicted.prob - prob.star)/n)

	res = c(beta.mu.hat.dist, beta.dist, beta.abs_diff.root_mu, beta.abs_diff.root_muhat, linear.rmse, prob.rmse)
	return(res)
}


eval_class = function(y, fitted.linear, linear.star, th=0.5, setting="logit")
{
	predicted.prob = sigmoid(fitted.linear)
	if( setting=="IP" ){
		prob.star = sigmoid(2*sin(linear.star*pi)+linear.star)
	}else{
		prob.star = sigmoid(linear.star)
	}

	# classification performance
	AUC.prob = pROC::roc(as.vector(y), predicted.prob, quiet=TRUE)$auc
	AUC.oracle = pROC::roc(as.vector(y), prob.star, quiet=TRUE)$auc
	AUC.diff = AUC.oracle - AUC.prob

	predicted.class = as.numeric(predicted.prob>th)
	oracle.class = as.numeric(prob.star>th)	

	accuracy = mean(predicted.class == y)
	accuracy.oracle = mean(oracle.class == y)
	acc.diff = accuracy.oracle - accuracy

	sen = sensitivity(y, predicted.class)
	spe = specificity(y, predicted.class)
	sen.oracle = sensitivity(y, oracle.class)
	spe.oracle = specificity(y, oracle.class)
	sen.diff = sen.oracle - sen
	spe.diff = spe.oracle - spe

	res = c(accuracy, acc.diff, AUC.prob, AUC.diff, sen, sen.diff, spe, spe.diff)
	return(res)
}


eval = function(mfd, y, fitted.linear, linear.star, beta.hat, beta.star, mu.hat, mu, dig=4, beta.0=0, setting="logit")
{
	res.num = eval_numeric(mfd, beta.hat, beta.star, fitted.linear, linear.star, mu.hat=mu.hat, mu=mu, beta.0=beta.0)
	res.class = eval_class(y, fitted.linear, linear.star, setting=setting)
	res = round(c(res.num, res.class), dig=dig)

	return(res)
}


fit_eval = function(mfd, X, y, linear.star, beta.star, mu, test.ratio=0.2, 
			seed=0, dig=4, beta.0=0, alpha.GD=0.8, thr.GD.abs=1e-3, thr.GD.rel=1e-6, num.iterations=500,
			trace=FALSE, type="numerical", grad.method="simple", setting="logit")
{
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

	# evaluation
	opt.GD = gradient_descent(mfd=mfd, X=X.train, y=y.train, alpha=alpha.GD, 
					threshold.abs=thr.GD.abs, threshold.rel=thr.GD.rel, 
					num.iterations=num.iterations, trace=trace, type=type,
					grad.method=grad.method)
	beta.hat = opt.GD$beta.hat
	mu.hat = frechet.mean(mfd, X.train)
	fitted.linear = apply(X.test, MARGIN=3, h, mfd=mfd, mu=mu.hat, beta=beta.hat)
	converged = opt.GD$converged

	res = eval(mfd, y.test, fitted.linear, linear.star[test.idx], beta.hat, beta.star, 
			dig=dig, mu.hat=mu.hat, mu=mu, beta.0=beta.0, setting=setting)
	res = c(res, converged)
	return(res)
}


gradient_descent = function(mfd, X, y, beta=diag(1, mfd$dim[1]), alpha=0.8, num.iterations=500, 
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


data_generation = function(n, mfd, beta.0=0, mu=diag(mfd$dim[1]), var=1, seed=0, beta.normalize=TRUE, 
			error.num=-9999, beta.type="diag", D=2, setting="logit")
{
	# generate Xs and compute sample frechet mean
	set.seed(1+1000*seed)
	X = rmatrix(mfd, n=n, mu=mu, sig=sqrt(var))

	# generate beta
	if( beta.type=="random" ){
		set.seed(23+2000*seed)
		beta.star = D*rmatrix(mfd) # \beta \in a given mfd
	}else if( beta.type=="diag" ){
		beta.star = D*diag(mfd$dim[1])
	}else if( beta.type=="AR1" ){
		beta.star = D*AR_cor(mfd$dim[1], rho=0.5)
	}
	if( beta.normalize ) beta.star = normalize(mfd, mu, beta.star)# normalize beta s.t. d.beta=1

	# generate y
	linear.star = NULL
	for(i in 1:n)
	{
		if( setting=="additive" ){
			#mu, mfd
			#h.cur = compute(X[,,i], beta.star, beta.0=beta.0)
			log.X.cur = rie.log(mfd, p=mu, q=X[,,i])
			log.beta.star = rie.log(mfd, p=mu, q=beta.star)
			X.cur.vec = log.X.cur[lower.tri(log.X.cur, diag=TRUE)]
			beta.star.vec = log.beta.star[lower.tri(log.beta.star, diag=TRUE)]
			h.cur = beta.star.vec[1]*sin(pi*X.cur.vec[1]) + beta.star.vec[2]*(X.cur.vec[2])^2 + beta.star.vec[3]*exp(X.cur.vec[3]) + sum(beta.star.vec[4:6]*X.cur.vec[4:6])
		}else{
			h.cur = error.num
			round = 1
			while( h.cur==error.num )
			{
				h.cur = tryCatch(
					{
						h(mfd, mu=mu, beta=beta.star, X=X[,,i], beta.0=beta.0)
					},
					error=function(cond){
						return(error.num)
					})
				if( h.cur != error.num ) break
				set.seed(i+11*seed+10000*round)
				X[,,i] = rmatrix(mfd, n=1, mu=mu, sig=sqrt(var))
				round = round + 1
			}
		}
		linear.star = c(linear.star, h.cur)
	}
#	linear.star2 = sapply(1:n, function(i){h(mfd, mu=mu, beta=beta.star, X=X[,,i], beta.0=beta.0)})

	if( setting=="IP" ){
		prob.star = sigmoid(2*sin(linear.star*pi)+linear.star)
	}else{
		prob.star = sigmoid(linear.star)
	}
	set.seed(45+3000*seed)
	y = rbinom(n, size=1, prob=prob.star)

	res = list(X=X, y=y, beta.star=beta.star, linear.star=linear.star, prob.star=prob.star)
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


simu_worker = function(n, mfd, beta.0=0, test.ratio=0.2, mu=diag(mfd$dim[1]), var=1, seed=0, 
			dig=4, beta.normalize=TRUE, trace=FALSE, grad.type="numerical", tmp.path=NULL,
			grad.method="simple", beta.type="diag", D=2, setting="logit")
{
	if( !is.null(tmp.path) ){
		tmp.name = paste0(tmp.path, "tmp_", seed, ".RData")
		if( file.exists(tmp.name) ){
			load(tmp.name)

			if( !("prob.star" %in% names(data.list)) ){
				if( setting=="IP" ){
					prob.star = sigmoid(2*sin(data.list$linear.star*pi)+data.list$linear.star)
				}else{
					prob.star = sigmoid(data.list$linear.star)
				}
				data.list = c(data.list, prob.star=list(prob.star))

				save.list = c("res.cur", "seed", "data.list")
				save(list=save.list, file=tmp.name)
			}

			return(res.cur)
		}
	}

	if( trace ) print("Data Generating.")
	data.list = data_generation(n, mfd, beta.0=beta.0, mu=mu, var=var, seed=seed, beta.normalize=beta.normalize, beta.type=beta.type, D=D, setting=setting)
	beta.mu.dist = geo.dist(mfd, data.list$beta.star, mu)
	mu.mu.hat.dist = geo.dist(mfd, frechet.mean(mfd, data.list$X), mu)
	var.X.oracle = mean(apply(data.list$X, 3, function(x)geo.dist(mfd, x, mu)^2))
	var.X.hat = mean(apply(data.list$X, 3, function(x)geo.dist(mfd, x, frechet.mean(mfd, data.list$X))^2))
	var.X.diff = abs(var.X.hat - var.X.oracle)	

	if( trace ) print("Optimizing.")
	res.cur = fit_eval(mfd=mfd, X=data.list$X, y=data.list$y, linear.star=data.list$linear.star, beta.star=data.list$beta.star, 
				mu=mu, beta.0=beta.0, test.ratio=test.ratio, seed=seed, dig=dig, trace=trace, type=grad.type,
				grad.method=grad.method, setting=setting)
	res.cur = c(mu.mu.hat.dist, var.X.oracle, var.X.hat, var.X.diff, beta.mu.dist, res.cur)
	
	if( !is.null(tmp.path) ){
		save.list = c("res.cur", "seed", "data.list")
		save(list=save.list, file=tmp.name)
	}
	return(res.cur)
}


simu = function(n, mfd, beta.0=0, R=100, mu=diag(mfd$dim[1]), var=1, test.ratio=0.2, ncore=1, ndiv=R/5, 
		dig=4, beta.normalize=TRUE, trace=FALSE, grad.type="numerical", tmp.path=NULL, grad.method="simple",
		fun.list=NULL, start.R=NULL, end.R=NULL, beta.type="diag", D=2, setting="logit")
{
	colNames = c("Mu-Mu.hat-Dist", 
			"Var-X-Oracle", "Var-X-Hat", "Var-Diff",
			"Beta-Mu-Dist", "Beta-Mu-hat-Dist",
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
#		if( mfd$metric=="AffineInvariant" ){
#			loop.R = setdiff(loop.R, c(400, 409, 412, 425, 426, 427, 434, 435))
#		}
		ndiv = length(loop.R) / 5
	}
	cat(sprintf("loop.R: from-%d, to-%d\n", loop.R[1], loop.R[length(loop.R)]))

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
							tmp.path=tmp.path, grad.method=grad.method, beta.type=beta.type, D=D, setting=setting)
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
						tmp.path=tmp.path, grad.method=grad.method, beta.type=beta.type, D=D, setting=setting)
			res[i,] = tmp 
		}
	}
	colnames(res) = colNames 

	return(res)
}






#####################################################################################################

rm_outlier = function(res)
{
	# only consider outliers in beta-rmse and linear-part-rmse
	rm.idx = unique(unlist(apply(res[,3:6], 2, function(x){which(x==boxplot.stats(x)$out)})))
	if( length(rm.idx)>0 )  res = res[-rm.idx,]

	return(res)
}


out_stat = function(res, R=NULL, outlier.drop=TRUE, dig=5, only.converged=TRUE)
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

