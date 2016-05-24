Faddyprob.limiting <-
function(parameter,nmax) { 
   nmax1 <- nmax + 1
   vlambda <- rep(0,nmax1) 
# if parameter[1] is Inf (infinity) na's are produced in all elements of vlambda 
# except the first
   if (is.na(parameter[1])==TRUE) { parameter[1] <- Inf }
   if (is.infinite(parameter[1])==TRUE) { probability <- rep(0,nmax1) 
       } else { if (parameter[2]==0) { vlambda <- rep(parameter[1],nmax1)  
                              } else { vbeta <- rep(parameter[2],nmax1)
                                       vnum  <- c(0:nmax)
                                       vlambda <- parameter[1]*exp(parameter[2]*vnum)
                                    } # end if (parameter(2)
         probability <- EPPMprob(vlambda) }
    return(probability)   }
