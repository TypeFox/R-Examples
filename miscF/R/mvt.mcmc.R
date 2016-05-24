mvtmcmc <- function(niter, lambda, Sigma.inv, Sigma0.inv, X, prior.Mu0,
                      prior.p, V.inv, v, prior.lower.v, prior.upper.v){
    .Call("mvtmcmc", niter, lambda, Sigma.inv, Sigma0.inv, X, prior.Mu0,
                      prior.p, V.inv, v, prior.lower.v, prior.upper.v)
}

mvt.mcmc <- function(X, niter, prior.lower.v, prior.upper.v,
                     prior.Mu0=rep(0, ncol(X)),
                     prior.Sigma0=diag(10000, ncol(X)),
                     prior.p=ncol(X), prior.V=diag(1, ncol(X)),
                     initial.v=NULL, initial.Sigma=NULL){
        
    if(!is.matrix(X) || ncol(X) < 2 || nrow(X) < 3){
        stop("The input observations of the 'mvt.mcmc' has to be a matrix
             of more than three rows and one column.")
    }
    if(prior.lower.v<0 || prior.upper.v<0 || prior.lower.v>prior.upper.v){
        stop("The bounds of the degrees of freedom of mvt have to be
             positive and the lower bound has to be smaller than the upper.")
    }
    if(!is.vector(prior.Mu0) || length(prior.Mu0)!=ncol(X)){
        stop("The mean vector of multivariate normal prior of the
             location of mvt has to be a vector of the same
             length as the number of columns of observations.")
    }
    if(!is.matrix(prior.Sigma0) || any(eigen(prior.Sigma0)$values<=0) || any(dim(prior.Sigma0)!=ncol(X))){
        stop("The variance matrix of multivariate normal
             prior of the location of mvt has to be a positive definite
             matrix with the number of rows and columns equal to
             the number of columns of observations.")
    }
    if(prior.p < ncol(X)){
        stop("The degrees of freedom of wishart prior of inverse of the
             scale matrix of mvt has be equal to or larger than the
             number of columns of observations.")
    }
    if(!is.matrix(prior.V) || any(eigen(prior.V)$values<=0) ||  any(dim(prior.V)!=ncol(X))){
        stop("The scale matrix of wishart prior of inverse of the scale
             matrix of mvt has to be a positive definite
             matrix with the number of rows and columns equal to
             the number of columns of observations.")
    }
    if(!is.null(initial.v) && (initial.v < prior.lower.v || initial.v > prior.upper.v)){
        stop("The initial value of the df of the mvt has to be within the bounds
             (prior.lower.v, prior.upper.v).")
    }
    if(!is.null(initial.Sigma) &&
           (!is.matrix(initial.Sigma) || any(eigen(initial.Sigma)$values<=0) || any(dim(initial.Sigma)!=ncol(X)))){
        stop("The initial value of the scale matrix of mvt has to be
             a positive definite matrix with the number of rows and columns
             equal to the number of columns of observations.")
    }
    
    
    n <- nrow(X)            
    d <- ncol(X)             
    
    Sigma0.inv <- solve(prior.Sigma0)
    V.inv <- solve(prior.V)
    
    #assign initial values obtained from ecme
    if(is.null(initial.v) || is.null(initial.Sigma)){
        ecme <- mvt.ecme(X, prior.lower.v, prior.upper.v)
        if(is.null(initial.v)){
            initial.v <- ecme$v
        }
        if(is.null(initial.Sigma)){
            initial.Sigma <- ecme$Sigma
        }
    }
    v <- initial.v
    Sigma <- initial.Sigma
    lambda <- rgamma(n, v/2, v/2)
    Sigma.inv <- solve(Sigma)
    
    #assign matrice for saving results
    Mu.save <- matrix(0, niter, d)
    Sigma.save <- array(0, dim=c(d,d,niter))
    v.save <- rep(0, niter)

    result <- mvtmcmc(niter, lambda, Sigma.inv, Sigma0.inv, X, prior.Mu0,
                      prior.p, V.inv, v, prior.lower.v, prior.upper.v)
    Mu.save <- matrix(result$Mu.save, nrow=niter, byrow=TRUE) 
    Sigma.save <- array(result$Sigma.save, dim=c(d, d, niter))
    v.save <- result$v

    list(Mu.save=Mu.save, Sigma.save=Sigma.save, v.save=v.save)

}
