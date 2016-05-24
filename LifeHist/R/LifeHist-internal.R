.matlik.2p <-
function(pars, matdat)
                  {
                   prop.est  <- 1/(1+exp(log(1/19)*(matdat$pred-pars[1])/(pars[2]-pars[1])));
                   iprop.est <- 1-prop.est;                                
                   negloglik <- -sum(matdat$count*(matdat$mat*log(prop.est)+matdat$juv*log(iprop.est)));  
                  }
.matlik.3p <-
function(pars, matdat)
                  {
                   prop.est  <- pars[1]/(1+exp(log(1/19)*(matdat$pred-pars[2])/(pars[3]-pars[2])));
                   iprop.est <- 1-prop.est;
                   negloglik <- -sum(matdat$count*(matdat$mat*log(prop.est)+matdat$juv*log(iprop.est)));
                  }
.vonBer <-
function(par, age, obslen, output, distr)
               {
                loglinf <- par[1]
                logl0   <- par[2]
                logk    <- par[3]
                predlen   <- exp(loglinf)*(1-(1-exp(logl0)/exp(loglinf))*exp(-exp(logk)*age))
                Likel     <- .GroAgeLik(predlen,obslen,distr,par)
                if(output=="predict")
                  {
                   resid             <- unlist(Likel[["DevianceResidual"]])
                   groageexp         <- data.frame(age=age,obslen=obslen,predlen=predlen,resid=resid)
                   return(groageexp)
                  }
                else
                  {
                   if(distr == "apnormal" | distr == "aplnormal")
                     {
                      negsup <- ((length(age)-2)/2)*log(sum(Likel[["Likelihood"]]));
                     }
                   else
                     {
                      negsup <- -sum(Likel[["Likelihood"]]);
                     }
                   return(negsup);
                  }
               }
.schnute1 <-
function(par, age, obslen, output, distr)
               {
                logl1 <- par[1]
                loga1 <- par[2]
                logmu <- par[3]
                logg1 <- par[4]
                logg2 <- par[5]
                predlen <- (exp(logl1)^exp(logg2)+(exp(logmu)^exp(logg2)/(1-exp(logg2))-exp(logl1)^exp(logg2))*(1-exp(-exp(logg1)*(age-exp(loga1)))))^(1/exp(logg2))
                Likel     <- .GroAgeLik(predlen,obslen,distr,par)
                if(output == "predict")
                  {
                   resid             <- unlist(Likel[["DevianceResidual"]])
                   groageexp         <- data.frame(age=age,obslen=obslen,predlen=predlen,resid=resid)
                   return(groageexp)
                  }
                else
                  {
                   if(distr == "apnormal" | distr == "aplnormal")
                     {
                      negsup <- ((length(age)-2)/2)*log(sum(Likel[["Likelihood"]]));
                     }
                   else
                     {
                      negsup <- -sum(Likel[["Likelihood"]]);
                     }
                   return(negsup);
                  }
               }
.schnute2 <-
function(par, age, obslen, output, distr)
               {
                loglinf <- par[1]
                loggam  <- par[2]
                logg1   <- par[3]
                predlen <- exp(loglinf)*exp(-exp(loggam)*exp(-exp(logg1)*age))
                Likel     <- .GroAgeLik(predlen,obslen,distr,par)
                if(output == "predict")
                  {
                   resid             <- unlist(Likel[["DevianceResidual"]])
                   groageexp         <- data.frame(age=age,obslen=obslen,predlen=predlen,resid=resid)
                   return(groageexp)
                  }
                else
                  {
                   if(distr == "apnormal" | distr == "aplnormal")
                     {
                      negsup <- ((length(age)-2)/2)*log(sum(Likel[["Likelihood"]]))
                     }
                   else
                     {
                      negsup <- -sum(Likel[["Likelihood"]])
                     }
                   return(negsup)
                  }
               }
.schnute3 <-
function(par, l0, age, obslen, output, distr)
               {
                logg1 <- par[1]
                logg2 <- par[2]
                predlen <- l0*(1+(exp(logg2)/(1-exp(logg2)))(1-exp(-exp(logg1)*age)))^(1/exp(logg2))
                Likel     <- .GroAgeLik(predlen,obslen,distr,par)
                if(output == "predict")
                  {
                   resid             <- unlist(Likel[["DevianceResidual"]])
                   groageexp         <- data.frame(age=age,obslen=obslen,predlen=predlen,resid=resid)
                   return(groageexp)
                  }
                else
                  {
                   if(distr == "apnormal" | distr == "aplnormal")
                     {
                      negsup <- ((length(age)-2)/2)*log(sum(Likel[["Likelihood"]]))
                     }
                   else
                     {
                      negsup <- -sum(Likel[["Likelihood"]])
                     }
                   return(negsup)
                  }
               }
.GroAgeLik <-
function(obslen,predlen,distr,par)
  {
                  Likel        <- vector("list",4);
                  names(Likel) <- c("Dispersion","Deviance","Likelihood","DevianceResidual");
                  sealen       <- length(obslen);
                  if(distr=='normal')
                    {
                     psi1        <- exp(tail(par,1));
                     dev1        <- (obslen-predlen)^2;
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obslen-predlen)^2;
                     res1        <- sign(obslen-predlen)*sqrt(dev1);
                    }
                  else if(distr=='apnormal')
                    {
                     psi1        <- NA;
                     dev1        <- NA;
                     likcontr1   <- (obslen-predlen)^2;
                     res1        <- obslen-predlen;
                    }
                  else if(distr=='lognormal')
                    {
                     psi1        <- exp(tail(par,1));
                     dev1        <- (log(obslen)-log(predlen))^2;
                     likcontr1   <- -(1/2)*log(2*obslen^2*pi*psi1)-(1/(2*psi1))*(log(obslen)-log(predlen))^2;
                     res1        <- sign(obslen-predlen)*sqrt(dev1);
                    }
                  else if(distr=='aplnormal')
                    {
                     psi1        <- NA;
                     dev1        <- NA;
                     likcontr1   <- (log(obslen)-log(predlen))^2;
                     res1        <- log(obslen)-log(predlen);
                    }
                  else if(distr=='gamma')
                    {
                     psi1        <- exp(tail(par,1));
                     dev1        <- -(2/psi1)*(log(obslen/predlen)-(obslen-predlen)/predlen);
                     likcontr1   <- (1/psi1)*(log(obslen)-log(psi1*predlen))-log(obslen)-sealen*lgamma(1/psi1)-obslen/(psi1*predlen);
                     res1        <- sign(obslen-predlen)*sqrt(dev1);
                    }
                  Likel[["Dispersion"]]       <- psi1
                  Likel[["Deviance"]]         <- dev1
                  Likel[["Likelihood"]]       <- likcontr1
                  Likel[["DevianceResidual"]] <- res1
                  return(Likel);
  }
