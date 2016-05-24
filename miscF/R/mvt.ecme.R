mvt.ecme <- function(X, lower.v, upper.v, err=1e-4){

    if(!is.matrix(X) || ncol(X) < 2 || nrow(X) < 3){
        stop("The input observations of the 'mvt.mcmc' has to be a matrix
              of more than three rows and one column.")
    }
    if(lower.v<0 || upper.v<0 || lower.v>upper.v){
        stop("The bounds of the degrees of freedom of mvt have to be
              positive and the lower bound has to be smaller than the upper.")
    }

	
	n <- nrow(X)               
	p <- ncol(X)              

    v <- 7
	Mu <- colMeans(X)
	S <- matrix(0, p, p)
	for(i in 1:n){
		S <- S + (X[i,]-Mu) %*% t(X[i,]-Mu)
	}
	S <- S/n
	Sigma <- S * (v-2)/v

    flag <- 0
    while(flag==0){
        #E step
        wi <- apply(X, 1, function(x) t(x - Mu) %*% solve(Sigma) %*% (x - Mu))
        wi <- (p + v) / (wi + v)
        #M step
        Mu <- colSums(X * wi) / sum(wi)
        Sigma <- rowSums(sapply(1:n, function(i)
                                (X[i,]-Mu) %*% t(X[i,]-Mu) * wi[i]))
        Sigma <- Sigma / n
        Sigma <- matrix(Sigma, p, p)

        f.lower <-  -digamma(lower.v/2) + log(lower.v/2) +
                    mean(log(wi) - wi) + 1 +
                    digamma((p+lower.v)/2) - log((p+lower.v)/2)
        f.upper <-  -digamma(upper.v/2) + log(upper.v/2) +
                    mean(log(wi) - wi) + 1 +
                    digamma((p+upper.v)/2) - log((p+upper.v)/2)
        if(f.lower * f.upper < 0){
            v.old <- v
            v <- uniroot(function(v,...)
                     -digamma(v/2) + log(v/2) + mean(log(wi) - wi) + 1 +
                     digamma((p+v)/2) - log((p+v)/2),
                     interval=c(lower.v, upper.v), wi=wi, p=p,
                     f.lower=f.lower, f.upper=f.upper, tol = 1e-4)$root
        
            if(abs(v.old-v)/v.old < err){
                flag <- 1
            }
        }
        else{
            flag <- 1
        }
    }
	list(Mu=Mu, Sigma=Sigma, v=v)
}

