# Main ssd function
ssd <- function(p1, p2, k, ratio=1, alpha=0.05, beta=0.2, cc=0.02, d=0.2, r=0.3, m, scheme='M2', Niter=500) {
       
       # Original method function
       original_method <- function(p1, p2, lambda0) {
                          ceiling( 2*lambda0/sum(((p1-p2)^2)/(0.5*(p1+p2))) )
       }
       
       # Minimum difference method function
       minimum_method <- function(p1, p2, lambda0,cc) {
                         absolute_difference <- abs(p1-p2)
                         absolute_difference[ absolute_difference<cc ] <- cc
                         ceiling( 2*lambda0/sum( (absolute_difference^2)/(0.5*(p1+p2)) ) )
       }
       
       # Correction method function
       correction_method_D <- function(p1, p2, tj) {
                              k <- length(p1)
                              zj <- rnorm(2*k, mean=0, sd=1)
                              vB <- 2*(p1-p2)*sqrt(tj)*zj[1:k]/(0.5*(p1+p2))
                              vC <- ((p1-p2)^2)*sqrt(tj)*zj[(k+1):(2*k)]/(2*(0.5*(p1+p2)))
                              sum(vB)-sum(vC)
       }
         
       correction_method <- function(p1, p2, n0, lambda0, Niter){
                            k <- length(p1)
                            a <- sum( ((p1-p2)^2)/(0.5*(p1+p2)) )
                            A <- 2*lambda0*(1/sqrt(n0))/(a^2)
                            tj <- p1*(1-p1)+p2*(1-p2)
                            ntemp <- replicate(Niter, correction_method_D(p1, p2, tj))
                            ceiling( original_method(p1, p2, lambda0) + A*mean(ntemp) )
       }
                   
       # Parametric bootstrap method function
       generate_bootstrap <- function(p1, p2, n0, lambda0) {
                             estimate_p1 <- rmultinom(1,n0,p1)/n0
                             estimate_p2 <- rmultinom(1,n0,p2)/n0
                             estimate_p <- 0.5*(estimate_p1+estimate_p2) 
                             while ( isTRUE( all.equal(0, min(estimate_p)) ) ) {
                                   estimate_p1 <- rmultinom(1,n0,p1)/n0
                                   estimate_p2 <- rmultinom(1,n0,p2)/n0
                                   estimate_p <- 0.5*(estimate_p1+estimate_p2)
                             }
                             original_method(estimate_p1, estimate_p2, lambda0) 
       }
                             
       bootstrap_method <- function(p1, p2, n0, Niter, lambda0) {
                           bootstrap_result <- replicate(Niter, generate_bootstrap(p1, p2, n0, lambda0))
                           c( ceiling(mean(bootstrap_result)), ceiling(median(bootstrap_result)), ceiling(quantile(bootstrap_result,0.75)), ceiling(quantile(bootstrap_result,0.80)) )
       }
       
       # Unbalanced Parametric bootstrap method function
       unbalanced_generate_bootstrap <- function(p1, p2, n1, n2, lambda0) {
                                        estimate_p1 <- rmultinom(1,n1,p1)/n1
                                        estimate_p2 <- rmultinom(1,n2,p2)/n2
                                        estimate_p <- 0.5*(estimate_p1+estimate_p2)
                                        while ( isTRUE( all.equal(0, min(estimate_p)) ) ) {
                                              estimate_p1 <- rmultinom(1,n1,p1)/n1
                                              estimate_p2 <- rmultinom(1,n2,p2)/n2
                                              estimate_p <- 0.5*(estimate_p1+estimate_p2)
                                        }
                                        unbalanced_original_method(estimate_p1, estimate_p2, n1, n2, lambda0)
       }
                                        
       unbalanced_bootstrap_method <- function(p1, p2, n1, n2, Niter, lambda0) {
                                      bootstrap_result <- replicate(Niter, unbalanced_generate_bootstrap(p1, p2, n1, n2, lambda0))
                                      c( ceiling(mean(bootstrap_result)), ceiling(median(bootstrap_result)), ceiling(quantile(bootstrap_result,0.75)), ceiling(quantile(bootstrap_result,0.80)) )
       }
       
       # Scheme M1 function
       M1_method <- function(r, d, k, lambda) {
                    ceiling( 2*lambda/(r*d*k) )
       }
       
       # unbalanced original method function
       unbalanced_original_method <- function(p1, p2, n1, n2, lambda0) {
                                     r <- n1/n2
                                     nn1 <- ceiling( (1+r)*lambda0/sum(((p1-p2)^2)/(0.5*(p1+p2))) )
                                     nn2 <- ceiling(nn1/r) 
                                     c(nn1,nn2)
       }
       
       # Main function
       if ( !missing(p1) ) {
          if ( (min(p1)<0) | (max(p1)>1) ) {
             stop( "the specified 'p1' is not valid" )
          }
       }
       
       if ( !missing(p2) ) {
          if ( (min(p2)<0) | (max(p2)>1) ) {
             stop( "the specified 'p2' is not valid" )
          }
       }
       
       if ( !missing(k) ) {
          if ( k<0 ) {
             stop( "the specified 'k' is not valid" )
          } else {
             k <- ceiling(k)
          }
       }
       
       if ( !missing(ratio) ) {
          if ( ratio<0 ) {
             stop( "the specified 'ratio' is not valid" ) 
          }
       }
       
       if ( !missing(alpha) ) {
          if ( (alpha<0) | (alpha>1) ) {
             stop( "the specified 'alpha' is not valid" )
          }
       }
       
       if ( !missing(beta) ) {
          if ( (beta<0) | (beta>1) ) {
             stop( "the specified 'beta' is not valid" )
          }
       }
       
       if ( !missing(d) ) {
          if ( (d<0) | (d>1) ) {
             stop( "the specified 'd' is not valid" )
          }
       }
       
       if ( !missing(r) ) {
          if ( r<0 ) {
             stop( "the specified 'r' is not valid" )
          }
       }
       
       if ( !missing(scheme) ) {
          if ( !(scheme %in% c('M1','M2')) ) {
             stop( "the specified 'scheme' is not valid" )
          }
       }
       
       if ( !missing(Niter) ) {
          if ( Niter<0 ) {
             stop( "the specified 'Niter' is not valid" )
          } else {
             Niter <- ceiling(Niter)
          }
       }
       
       if ( !missing(cc) ) {
          if ( (cc<0) | (cc>1) ) {
             stop( "the specified 'cc' is not valid" )
          }
       }
       
       if ( !missing(m) ) {
          if ( m<0 ) {
             stop( "the specified 'm' is not valid" )
          } else {
             m <- ceiling(m)
          }
       }
       
       if ( !missing(m) & !missing(ratio) ) {
          n2 <- ceiling(m/ratio)
       }
       
       kk <- 0
       if ( !missing(p1) & !missing(p2) ) {
          if ( !isTRUE( all.equal(length(p1), length(p2)) ) ) {
             stop( "the specified 'p1' and 'p2' are not valid" )
          } else {
             kk <- length(p1)
          }
       }
       if ( !missing(k) ) {
          kk <- k
       }
       
       if ( !isTRUE( all.equal(0, kk) ) ) {
          lambda <- 1
          ncchi <- qchisq(beta, df=kk-1, ncp=lambda)
          cchi <- qchisq(1-alpha, df=kk-1)
          while ( ncchi<cchi ) {
                lambda <- lambda+0.01
                ncchi <- qchisq(beta, df=kk-1, ncp=lambda)
                cchi <- qchisq(1-alpha, df=kk-1)
          }
       } else {
          stop( "either 'p1' and 'p2' or 'k' needs to be specified" )
       }
       
       if ( scheme=='M2' ) {
          if ( !missing(p1) & !missing(p2) ) {
             if ( isTRUE( all.equal(1, ratio) ) ) {
                cat( "The calculated sample sizes under scheme M2 are following:\n" )
                cat( paste( "Original Method:",original_method(p1, p2, lambda),"\n",sep=' ' ) )
                cat( paste( "Minimum Difference Method:",minimum_method(p1, p2, lambda, cc),paste("(cc=",cc,")",sep=''),"\n",sep=' ' ) )
                if ( !missing(m) ) {      
                   cat( paste( "Correction Method:",correction_method(p1, p2, m, lambda, Niter),"\n",sep=' ' ) )
                   bootstrap_result <- bootstrap_method(p1, p2, m, Niter, lambda)
                   cat( paste( "Bootstrap Mean Method:",bootstrap_result[1],"\n",sep=' ' ) )
                   cat( paste( "Bootstrap Median Method:",bootstrap_result[2],"\n",sep=' ' ) )
                   cat( paste( "Bootstrap 75% Percentile Method:",bootstrap_result[3],"\n",sep=' ' ) )
                   cat( paste( "Bootstrap 80% Percentile Method:",bootstrap_result[4],"\n",sep=' ' ) )
                }
             }
             if ( !isTRUE( all.equal(1, ratio) ) ) {
                if ( !missing(m) ) {
                   cat( "The calculated sample sizes under scheme M2 are following:\n" )
                   unbalanced_original_result <- unbalanced_original_method(p1, p2, m, n2, lambda)
                   cat( paste( "Original Method:",paste("n1=",unbalanced_original_result[1],";",sep=''),paste("n2=",unbalanced_original_result[2],sep=''),"\n",sep=' ' ) )
                   unbalanced_bootstrap_result <- unbalanced_bootstrap_method(p1, p2, m, n2, Niter, lambda)
                   cat( paste( "Bootstrap Mean Method:",paste("n1=",unbalanced_bootstrap_result[1],";",sep=''),paste("n2=",ceiling(unbalanced_bootstrap_result[1]/ratio),sep=''),"\n",sep=' ' ) )
                   cat( paste( "Bootstrap Median Method:",paste("n1=",unbalanced_bootstrap_result[2],";",sep=''),paste("n2=",ceiling(unbalanced_bootstrap_result[2]/ratio),sep=''),"\n",sep=' ' ) )
                   cat( paste( "Bootstrap 75% Quantile Method:",paste("n1=",unbalanced_bootstrap_result[3],";",sep=''),paste("n2=",ceiling(unbalanced_bootstrap_result[3]/ratio),sep=''),"\n",sep=' ' ) )
                   cat( paste( "Bootstrap 80% Quantile Method:",paste("n1=",unbalanced_bootstrap_result[4],";",sep=''),paste("n2=",ceiling(unbalanced_bootstrap_result[4]/ratio),sep=''),"\n",sep=' ' ) )
                } else {
                   stop( "'m' needs to be specified" )
                }
             }
          } else {
             stop( "'p1' and 'p2' need to be specified" )
          }
       }
       
       if ( scheme=='M1' ) {
          if ( !missing(k) ) {
             cat( "The calculated sample sizes under scheme M1 are following:\n" )
             cat( paste( "Scheme M1:",M1_method(r, d, k, lambda),"\n",sep=' ' ) )
          } else {
             stop( "the parameter setting is not valid" )
          }
       }
}