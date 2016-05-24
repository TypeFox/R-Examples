GBprob <-
function(parameter,nt) {
      m <- nt + 1
      a <- parameter[1] 
      b <- parameter[2]
      if (b>=0) {
         if (round(b,digits=14)==1) {
            p <- 1 - exp(-a)
            probability <- dbinom(c(0:nt),nt,p,log=FALSE) 
                                          } else {
            vlambda <- c(rep(a,m))
            if (b>0) {            
               vlambda <- vlambda*(c(rep(nt,m)) - c(0:nt))^c(rep(b,m)) } # end if (b>0)
# limiting value for lambda
            lambda.limit <- 745
            vlambda <- sapply(1:m, function(j) 
                 if ((is.finite(vlambda[j])==FALSE) | (vlambda[j]>lambda.limit)) { 
                                        vlambda[j] <- lambda.limit  
                                      } else { vlambda[j] <- vlambda[j] } )
            probability <- EPPMprob(vlambda) 
                                    } # end of if (b==1)
          } else { probability <- rep(NA,m) } # end if (b>=0)
      return(probability)           }
