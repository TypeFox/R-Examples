nextB <- function(a, b, d, p, f, env)
{
     
     constraintsNorm <- env$constraintsNorm
     
     # note the multiplication by rho. this is the "jumping distribution"
     bnew <- b + matrix(rnorm(env$gNP*env$gNIV),nrow=env$gNP,ncol=env$gNIV)%*%chol(d)*sqrt(env$rho)
     
     # apply constraints here
     # cycle through constraints till no constraints are violated
     
     if(!is.null(constraintsNorm))
     {
          totalConstraintsViolated <- 999
          
          while(totalConstraintsViolated > 0)
          {
               totalConstraintsViolated <- 0
               for(i in 1:length(constraintsNorm))
               {
                    firstpar  <- constraintsNorm[[i]][1]
                    condition <- constraintsNorm[[i]][2]  # 1 = <, 2 = >
                    secondpar <- constraintsNorm[[i]][3]
                    # less than constraint
                    if(condition==1)
                    {
                         if(secondpar==0)
                         {          
                              constraintsViolated <- bnew[,firstpar]>0
                              bnew[constraintsViolated,firstpar] <- 0
                         }
                         if(secondpar!=0)
                         {
                              constraintsViolated <- bnew[,firstpar]>bnew[,secondpar]
                              bnew[constraintsViolated,c(firstpar,secondpar)] <- bnew[constraintsViolated,secondpar]
                         }
                    } # if 1
                    if(condition==2)
                    {
                         if(secondpar==0)
                         {          
                              constraintsViolated <- bnew[,firstpar]<0
                              bnew[constraintsViolated,firstpar] <- 0
                         }
                         if(secondpar!=0)
                         {
                              constraintsViolated <- bnew[,firstpar]<bnew[,secondpar]
                              bnew[constraintsViolated,c(firstpar,secondpar)] <- bnew[constraintsViolated,secondpar]
                         }
                    } # if 2
               
                    totalConstraintsViolated <- totalConstraintsViolated + sum(constraintsViolated)
                    
               } # for
               
          } # while
     } # if
     
     bn   <- bnew - matrix(t(a),nrow=env$gNP,ncol=env$gNIV,byrow=T)
     bo   <- b - matrix(t(a),nrow=env$gNP,ncol=env$gNIV,byrow=T)
     
     pnew <- env$likelihood(f,bnew,env)

     # gives the relative density (or the probability of seeing the betas) given current estimates
     # of D and A	
     # the relative densities are important because without any reference to how the betas fit in with the current
     # estimates of A and D, you'd end up with 
     # r.new <- pnew / p * exp(-0.5 * (colSums(t(bn)*(solve(d)%*%t(bn))) -	colSums(t(bo)*(solve(d)%*%t(bo)))))
     r.new <- (log(pnew) + -0.5 * (colSums(t(bn)*(solve(d)%*%t(bn)))) - log(p) - -0.5 * (colSums(t(bo)*(solve(d)%*%t(bo)))))
     
     # if r.new > 1 then we accept the new estimate of beta. if r.new < 1 then we accept the new estimate
     # with probability = r.new
     
     #ind  <- (r.new >= 1) + (r.new < 0)*(matrix(runif(env$gNP),nrow=env$gNP) <= r.new)
     ind  <- (r.new >= 0) + (r.new < 0)*(log(matrix(runif(env$gNP),nrow=env$gNP)) <= r.new)
     nind <- 1*(ind==0)
     
     # this is the acceptance rate. the target for this 0.3 (though Sawtooth allows for the user to specify this).
     i <- colSums(ind)/env$gNP
     
     if(i < env$targetAcceptanceNormal)
     {
          env$rho <- env$rho - env$rho/10
     }
     if(i > env$targetAcceptanceNormal)
     {
          env$rho <- env$rho + env$rho/10
     }

     #if(env$rho<0)
     #{
     #     env$rho <- 0.0001          
     #}
     
     # i've just converted it to matrix form to make the multiplication simpler.
     mind  <- matrix(0,nrow = env$gNP,ncol=env$gNIV)
     mnind <- matrix(0,nrow = env$gNP,ncol=env$gNIV)
     
     mind[,1:env$gNIV]  <- ind
     mnind[,1:env$gNIV] <- 1*(ind == 0)
     
     return(list(mind * bnew + mnind * b, i, ind * pnew + nind * p))
     
}
