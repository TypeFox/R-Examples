MatFit <-
function(p, par, matdat, method, control=ls(), itnmax)
                  {
                   results1 <- do.call(optimx, list(par=par,
                                                    fn=eval(as.name(paste('.matlik.',as.character(p),'p',sep=""))),
                                                    matdat=matdat,
                                                    gr=NULL,
                                                    method=method,
                                                    lower = -Inf,
                                                    upper = Inf,
                                                    control=list(),
                                                    hessian=TRUE,
                                                    itnmax=itnmax));
                   results2        <- vector("list",length(method));
                   names(results2) <- method;
                   temp            <- attr(results1, "details");
                      for(i in 1:length(method))
                      {
                       results2[[i]]$Model     <- paste('matlik.',as.character(p),'p',sep="");
                       results2[[i]]$converg   <- NA;
                       results2[[i]]$kkt       <- NA;
                       results2[[i]]$AIC       <- NA;
                       results2[[i]]$par.mle   <- NA;
                       results2[[i]]$num.grads <- NA;
                       results2[[i]]$stdev     <- NA;
                       results2[[i]]$Cor       <- matrix(NA,length(par),length(par));
                       if(length(temp[i,]$ngatend) == length(par) & !any(is.na(temp[i,]$nhatend)) & 1/kappa(temp[i,]$nhatend)>1e-15)
                         {
                          results2[[i]]$converg    <- results1[i,length(par)+5];
                          results2[[i]]$kkt        <- results1[i,(length(par)+6):(length(par)+7)];
                          results2[[i]]$AIC        <- 2*length(par)-2*results1[i,length(par)+1];
                          results2[[i]]$par.mle   <- results1[i,1:length(par)];
                          results2[[i]]$num.grads  <- temp[i,]$ngatend;
                                                v  <- matrix(0,length(par),length(par));
                                                v  <- try(solve(temp[i,]$nhatend));
                          results2[[i]]$stdev      <- sqrt(diag(v));
                          results2[[i]]$Cor        <- cor(v);
                         }
                      }
                  rm(v);
                  rm(temp);
                  return(results2);
                  }
