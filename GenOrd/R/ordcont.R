ordcont<-function (marginal, Sigma, support = list(), Spearman = FALSE, 
    epsilon = 1e-06, maxit = 100) 
{
        if (!all(unlist(lapply(marginal, function(x) (sort(x)==x & min(x)>0 & max(x)<1))))) stop("Error in assigning marginal distributions!")
    if(!isSymmetric(Sigma) | min(eigen(Sigma)$values)<0 | !all(diag(Sigma)==1)) stop("Correlation matrix not valid!")
    len <- length(support)
    k <- length(marginal)
    niter<-matrix(0,k,k)
    kj <- numeric(k)
    for (i in 1:k) {
        kj[i] <- length(marginal[[i]]) + 1
		if (len == 0) {
            support[[i]] <- 1:kj[i]
        }

    }
    for (g in 2:k) {
        if (det(Sigma[1:g, 1:g]) <= 0) {
            stop("Main minor number ", g, " is not positive!", 
                "\n")
        }
    }
    mcmin <- corrcheck(marginal, support, Spearman)[[1]]
    mcmax <- corrcheck(marginal, support, Spearman)[[2]]
    if (sum(mcmin <= Sigma & Sigma <= mcmax) != k^2) {
        stop("Some correlation coefficients are not feasible!", 
            "\n", "Please use function corrcheck to get lower and upper bounds!", 
            "\n")
    }
    Sigma0 <- Sigma
    Sigmaord <- Sigma
    Sigmaord <- contord(marginal, Sigma, support, Spearman)
    Sigmaold <- Sigma
    Sigmaordold <- Sigmaord

    
        for (q in 1:(k - 1)) {
            for (r in (q + 1):k) {
                if (Sigma0[q, r] == 0) {
                  Sigma[q, r] <- 0
                }
                else {
		      it <- 0
			while (max(abs(Sigmaord[q, r] - Sigma0[q, r]) > epsilon) & it < maxit) {
                  if (Sigma0[q, r] * (Sigma0[q, r]/Sigmaordold[q, 
                    r]) >= 1) {
                    Sigma[q, r] <- Sigmaold[q, r] * (1 + 0.1 * 
                      (1 - Sigmaold[q, r]) * sign(Sigma0[q, r] - 
                      Sigmaord[q, r]))
                  }
                  else {
                    Sigma[q, r] <- Sigmaold[q, r] * (Sigma0[q, 
                      r]/Sigmaord[q, r])
                  }
                  
                Sigma[r, q] <- Sigma[q, r]
      	    Sigmaord[r, q] <- contord(list(marginal[[q]],marginal[[r]]), matrix(c(1,Sigma[q, r],Sigma[q, r],1),2,2), list(support[[q]],support[[r]]), Spearman)[2]
		    Sigmaord[q, r]<-Sigmaord[r, q]
        	    Sigmaold[q, r]<- Sigma[q, r]
		    Sigmaold[r, q]<-Sigmaold[q, r]
                it <- it + 1
                }
           niter[q,r]<-it
           niter[r,q]<-it
            }
               
        }
	            
}
		

if(eigen(Sigma)$values[k]<=0)
{
stop("The Sigma matrix is not coherent with the given margins or
it is impossible to find a feasible correlation matrix for MVN
ensuring Sigma for the given margins")
}
    emax <- max(abs(Sigmaord - Sigma0))
    list(SigmaC = Sigma, SigmaO = Sigmaord, Sigma = Sigma0, niter = niter, 
        maxerr = emax)
}
