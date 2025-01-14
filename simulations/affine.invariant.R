
rie.exp.spd_AffineInvariant <- function(mfd,p,v,...)
{
    sqrp = expm::sqrtm(p)
    make.sym(sqrp %*% expm::expm(solve(p)%*%v) %*% sqrp)
}


rie.log.spd_AffineInvariant <- function(mfd,p,q,...)
{
    sqrp = expm::sqrtm(p)
    make.sym(sqrp %*% expm::logm(solve(p)%*%q) %*% sqrp)
}


geo.dist.spd_AffineInvariant <- function(mfd,p,q,...)
{
    sqrp.inv = solve(expm::sqrtm(p))
    sqrt( sum( expm::logm(sqrp.inv%*%q%*%sqrp.inv)^2 ) )
}


rie.metric.spd_AffineInvariant <- function(mfd,p,u,v,...)
{
    p.inv = solve(p)
    sum( (p.inv%*%u) * (p.inv%*%v) )
}


geodesic.spd_AffineInvariant <- function(mfd,p,u,t,...)
{
    d = mfd$dim[1]
    if(length(t) > 1)
    {
        R <- array(0,c(d,d,length(t)))
        for(i in 1:length(t))
            R[,,i] <- make.sym(rie.exp(mfd,p,t[i]*u))
    }
    else
    {
        R <- make.sym(rie.exp(mfd,p,t*u))
    }
    
    return(R)
}


rtvecor.spd_AffineInvariant <- function(mfd,n=1,sig=1,drop=T)
{
    return(rsym(d=mfd$dim[1],n=n,sig=sig,drop=drop))
}


rmatrix.spd_AffineInvariant <- function(mfd,n=1,mu=NULL,sig=1,drop=T)
{
    if(is.null(mu)) mu <- diag(rep(1,mfd$dim[1]))
    stopifnot(is.spd(mu))
    
    S <- rtvecor(mfd,n=n,sig,drop=F)
    
    R <- array(0,c(mfd$dim[1],mfd$dim[1],n))
    for(i in 1:n)
    {
        R[,,i] <- rie.exp(mfd,mu,S[,,i])
    }
    
    if(n==1 && drop) return(as.matrix(R[,,1]))
    else return(R)
}


frechet.mean.spd_AffineInvariant <- function(mfd,S)
{
    if( !is.array(S) ) stop('S must be an array for the Affine Invariant metric.')

    M = CovTools::CovMean(S, method="AIRM")
    return(M)
}
