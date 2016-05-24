.onAttach <-
function (lib, pkg) {
    pkg.info <- drop(read.dcf(file=system.file("DESCRIPTION", package="CatDyn"),
                              fields=c("Title","Version","Date")))

   packageStartupMessage(paste("--------------------------------------------------------------\n",
			       pkg.info["Title"]),"\n",
    paste(" CatDyn version ", pkg.info["Version"],
              " (built on ", pkg.info["Date"], ") is now loaded\n", sep=""),
			 "--------------------------------------------------------------\n"
    )
}
.psi.nb2 <-
function(psi,predcat1,obscat1) {sum(2*psi*log((psi+predcat1)/(psi+obscat1)) +
                                                     psi/(psi+obscat1) -
                                                     obscat1*(2*psi+obscat1)/(6*psi*(psi+obscat1)^2)) - psi}
.CatDynLik1F <-
function(obscat1,predcat1,distr,par) # <- put par in the arguments
  {
                  Likel        <- vector("list",4);
                  names(Likel) <- c("Dispersion","Deviance","Likelihood","DevianceResidual");
                  sealen       <- length(obscat1);
                  if(distr=='poisson')
                    {
                     psi1        <- 1;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                    }
                  else if(distr=='negbin')
                    {
                     psi1        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                    }
                  else if(distr=='normal')
                    {
                     psi1        <- exp(tail(par,1));
                     dev1        <- (obscat1-predcat1)^2;
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                    }
                  else if(distr=='apnormal')
                    {
                     psi1        <- NA;
                     dev1        <- NA;
                     likcontr1   <- (obscat1-predcat1)^2;
                     res1        <- obscat1-predcat1;
                    }
                  else if(distr=='lognormal')
                    {
                     psi1        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(1/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                    }
                  else if(distr=='aplnormal')
                    {
                     psi1        <- NA;
                     dev1        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));;
                    }
                  else if(distr=='gamma')
                    {
                     psi1        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-(2/psi1)*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                    }
                  Likel[["Dispersion"]]       <- psi1
                  Likel[["Deviance"]]         <- dev1
                  Likel[["Likelihood"]]       <- likcontr1
                  Likel[["DevianceResidual"]] <- res1
                  return(Likel);
  }
.CatDynLik2F <-
function(obscat1,predcat1,obscat2,predcat2,distr,par)
  {
                  Likel        <- vector("list",4);
                  names(Likel) <- c("Dispersion","Deviance","Likelihood","DevianceResidual");
                  sealen       <- length(obscat1);
                  if(sum(distr==c('poisson','poisson')) == 2)
                    {
                     psi1        <- 1;
                     psi2        <- 1;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,2*(obscat2*log(obscat2)-obscat2*log(predcat2)-(obscat2-predcat2)));
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     likcontr2   <- ifelse(predcat2==0,0,obscat2*log(predcat2) - predcat2 - lfactorial(obscat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #1/49
                  else if(sum(distr==c('poisson','negbin')) == 2)
                    {
                     psi1        <- 1;
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     dev2        <- ifelse(obscat2==0,2*log(1+predcat2),
                                                      2*(obscat2*log(obscat2/predcat2)-(obscat2+1)*log((1+obscat2)/(1+predcat2))));
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     likcontr2   <- ifelse(predcat2==0,0,
                                           obscat2*log(1/psi2) +
                                           obscat2*log(predcat2) -
                                           (obscat2 + psi2)*log(1 + predcat2/psi2) +
                                           lgamma(obscat2 + psi2) -
                                           lgamma(obscat2 + 1) -
                                           lgamma(psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #2/49
                  else if(sum(distr==c('poisson','normal')) == 2)
                    {
                     psi1        <- 1;
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     dev2        <- (obscat2-predcat2)^2;
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     likcontr2   <- -(1/2)*log(2*pi*psi2)-(1/(2*psi2))*(obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #3/49
                  else if(sum(distr==c('poisson','apnormal')) == 2)
                    {
                     psi1        <- 1;
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     dev2        <- NA;
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     likcontr2   <- (obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- obscat2-predcat2;
                    }
                  #4/49
                  else if(sum(distr==c('poisson','lognormal')) == 2)
                    {
                     psi1        <- 1;
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,-(1/2)*log(2*obscat2^2*pi*psi2)-(1/(2*psi2))*(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #5/49
                  else if(sum(distr==c('poisson','aplnormal')) == 2)
                    {
                     psi1        <- 1;
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     dev2        <- NA;
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- ifelse(obscat2==0 | predcat2==0,0,log(obscat2)-log(predcat2));
                    }
                  #6/49
                  else if(sum(distr==c('poisson','gamma')) == 2)
                    {
                     psi1        <- 1;
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,2*(obscat1*log(obscat1)-obscat1*log(predcat1)-(obscat1-predcat1)));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(log(obscat2/predcat2)-(obscat2-predcat2)/predcat2));
                     likcontr1   <- ifelse(predcat1==0,0,obscat1*log(predcat1) - predcat1 - lfactorial(obscat1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(1/psi2)*(log(obscat2)-log(psi2*predcat2))-log(obscat2)-sealen*lgamma(1/psi2)-obscat2/(psi2*predcat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #7/49
                  else if(sum(distr==c('negbin','poisson')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- 1;
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,2*(obscat2*log(obscat2)-obscat2*log(predcat2)-(obscat2-predcat2)));
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     likcontr2   <- ifelse(predcat2==0,0,obscat2*log(predcat2) - predcat2 - lfactorial(obscat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #8/49
                  else if(sum(distr==c('negbin','negbin')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     dev2        <- ifelse(obscat2==0,2*log(1+predcat2),
                                                      2*(obscat2*log(obscat2/predcat2)-(obscat2+1)*log((1+obscat2)/(1+predcat2))));
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     likcontr2   <- ifelse(predcat2==0,0,
                                           obscat2*log(1/psi2) +
                                           obscat2*log(predcat2) -
                                           (obscat2 + psi2)*log(1 + predcat2/psi2) +
                                           lgamma(obscat2 + psi2) -
                                           lgamma(obscat2 + 1) -
                                           lgamma(psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #9/49
                  else if(sum(distr==c('negbin', 'normal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     dev2        <- ((obscat2-predcat2)^2);
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     likcontr2   <- -(1/2)*log(2*pi*psi2)-(1/(2*psi2))*(obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #10/49
                  else if(sum(distr==c('negbin', 'apnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     dev2        <- obscat2-predcat2;
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     likcontr2   <- dev2^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- dev2;
                    }
                  #11/49
                  else if(sum(distr==c('negbin', 'lognormal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,-(1/2)*log(2*obscat2^2*pi*psi2)-(1/(2*psi2))*(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #12/49
                  else if(sum(distr==c('negbin', 'aplnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     dev2        <- NA;
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- log(obscat2)-log(predcat2);
                    }
                  #13/49
                  else if(sum(distr==c('negbin', 'gamma')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0,2*log(1+predcat1),
                                                      2*(obscat1*log(obscat1/predcat1)-(obscat1+1)*log((1+obscat1)/(1+predcat1))));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(log(obscat2/predcat2)-(obscat2-predcat2)/predcat2));
                     likcontr1   <- ifelse(predcat1==0,0,
                                           obscat1*log(1/psi1) +
                                           obscat1*log(predcat1) -
                                           (obscat1 + psi1)*log(1 + predcat1/psi1) +
                                           lgamma(obscat1 + psi1) -
                                           lgamma(obscat1 + 1) -
                                           lgamma(psi1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(1/psi2)*(log(obscat2)-log(psi2*predcat2))-log(obscat2)-sealen*lgamma(1/psi2)-obscat2/(psi2*predcat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #14/49
                  else if(sum(distr==c('normal','poisson')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- 1;
                     dev1        <- ((obscat1-predcat1)^2)/sealen;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,2*(obscat2*log(obscat2)-obscat2*log(predcat2)-(obscat2-predcat2)));
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     likcontr2   <- ifelse(predcat2==0,0,obscat2*log(predcat2) - predcat2 - lfactorial(obscat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #15/49
                  else if(sum(distr==c('normal', 'negbin')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ((obscat1-predcat1)^2);
                     dev2        <- ifelse(obscat2==0,2*log(1+predcat2),
                                                      2*(obscat2*log(obscat2/predcat2)-(obscat2+1)*log((1+obscat2)/(1+predcat2))));
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     likcontr2   <- ifelse(predcat2==0,0,
                                           obscat2*log(1/psi2) +
                                           obscat2*log(predcat2) -
                                           (obscat2 + psi2)*log(1 + predcat2/psi2) +
                                           lgamma(obscat2 + psi2) -
                                           lgamma(obscat2 + 1) -
                                           lgamma(psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #16/49
                  else if(sum(distr==c('normal','normal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ((obscat1-predcat1)^2);
                     dev2        <- ((obscat2-predcat2)^2);
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     likcontr2   <- -(1/2)*log(2*pi*psi2)-(1/(2*psi2))*(obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #17/49
                  else if(sum(distr==c('normal','apnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ((obscat1-predcat1)^2);
                     dev2        <- NA;
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     likcontr2   <- (obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- obscat2-predcat2;
                    }
                  #18/49
                  else if(sum(distr==c('normal','lognormal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ((obscat1-predcat1)^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,-(1/2)*log(2*obscat2^2*pi*psi2)-(1/(2*psi2))*(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #19/49
                  else if(sum(distr==c('normal','aplnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ((obscat1-predcat1)^2);
                     dev2        <- NA;
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- ifelse(obscat2==0 | predcat2==0,0,log(obscat2)-log(predcat2));
                    }
                  #20/49
                  else if(sum(distr==c('normal','gamma')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ((obscat1-predcat1)^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(log(obscat2/predcat2)-(obscat2-predcat2)/predcat2));
                     likcontr1   <- -(1/2)*log(2*pi*psi1)-(1/(2*psi1))*(obscat1-predcat1)^2;
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(1/psi2)*(log(obscat2)-log(psi2*predcat2))-log(obscat2)-sealen*lgamma(1/psi2)-obscat2/(psi2*predcat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #21/49
                  else if(sum(distr==c('apnormal','poisson')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- 1;
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,2*(obscat2*log(obscat2)-obscat2*log(predcat2)-(obscat2-predcat2)));
                     likcontr1   <- (obscat1-predcat1)^2;
                     likcontr2   <- ifelse(predcat2==0,0,obscat2*log(predcat2) - predcat2 - lfactorial(obscat2));
                     res1        <- obscat1-predcat1;
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #22/49
                  else if(sum(distr==c('apnormal','negbin')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0,2*log(1+predcat2),
                                                      2*(obscat2*log(obscat2/predcat2)-(obscat2+1)*log((1+obscat2)/(1+predcat2))));
                     likcontr1   <- (obscat1-predcat1)^2;
                     likcontr2   <- ifelse(predcat2==0,0,
                                           obscat2*log(1/psi2) +
                                           obscat2*log(predcat2) -
                                           (obscat2 + psi2)*log(1 + predcat1/psi2) +
                                           lgamma(obscat2 + psi2) -
                                           lgamma(obscat2 + 1) -
                                           lgamma(psi2));
                     res1        <- obscat1-predcat1;
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #23/49
                  else if(sum(distr==c('apnormal','normal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ((obscat2-predcat2)^2)/sealen;
                     likcontr1   <- (obscat1-predcat1)^2;
                     likcontr2   <- -(sealen/2)*log(2*pi*psi2)-(1/(2*psi2))*(obscat2-predcat2)^2;
                     res1        <- obscat1-predcat1;
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #24/49
                  else if(sum(distr==c('apnormal','apnormal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- NA;
                     dev1        <- NA;
                     dev2        <- NA;
                     likcontr1   <- (obscat1-predcat1)^2;
                     likcontr2   <- (obscat2-predcat2)^2;
                     res1        <- obscat1-predcat1;
                     res2        <- obscat2-predcat2;
                    }
                  #25/49
                  else if(sum(distr==c('apnormal','lognormal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- (obscat1-predcat1)^2;
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,-(1/2)*log(2*obscat2^2*pi*psi2)-(1/(2*psi2))*(log(obscat2)-log(predcat2))^2);
                     res1        <- obscat1-predcat1;
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #26/49
                  else if(sum(distr==c('apnormal','aplnormal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- NA;
                     dev1        <- NA;
                     dev2        <- NA;
                     likcontr1   <- (obscat1-predcat1)^2;
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     res1        <- obscat1-predcat1;
                     res2        <- ifelse(obscat2==0 | predcat2==0,0,log(obscat2)-log(predcat2));
                    }
                  #27/49
                  else if(sum(distr==c('apnormal','gamma')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(log(obscat2/predcat2)-(obscat2-predcat2)/predcat2));
                     likcontr1   <- (obscat1-predcat1)^2;
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(1/psi2)*(log(obscat2)-log(psi2*predcat2))-log(obscat2)-sealen*lgamma(1/psi2)-obscat2/(psi2*predcat2));
                     res1        <- obscat1-predcat1;
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #28/49
                  else if(sum(distr==c('lognormal', 'poisson')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- 1;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2)/sealen;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,2*(obscat2*log(obscat2)-obscat2*log(predcat2)-(obscat2-predcat2)));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(1/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(predcat2==0,0,obscat2*log(predcat2) - predcat2 - lfactorial(obscat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #29/49
                  else if(sum(distr==c('lognormal','negbin')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ifelse(obscat2==0,2*log(1+predcat2),
                                                      2*(obscat2*log(obscat2/predcat2)-(obscat2+1)*log((1+obscat2)/(1+predcat2))));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(1/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(predcat2==0,0,
                                           obscat2*log(1/psi2) +
                                           obscat2*log(predcat2) -
                                           (obscat2 + psi2)*log(1 + predcat2/psi2) +
                                           lgamma(obscat2 + psi2) -
                                           lgamma(obscat2 + 1) -
                                           lgamma(psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #30/49
                  else if(sum(distr==c('lognormal','normal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ((obscat2-predcat2)^2)/sealen;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(1/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- -(1/2)*log(2*pi*psi2)-(1/(2*psi2))*(obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #31/49
                  else if(sum(distr==c('lognormal','apnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(sealen/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- (obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- obscat2-predcat2;
                    }
                  #32/49
                  else if(sum(distr==c('lognormal','lognormal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(1/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,-(1/2)*log(2*obscat2^2*pi*psi2)-(1/(2*psi2))*(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #33/49
                  else if(sum(distr==c('lognormal','aplnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(1/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- ifelse(obscat2==0 | predcat2==0,0,log(obscat2)-log(predcat2));
                    }
                  #34/49
                  else if(sum(distr==c('lognormal','gamma')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(log(obscat2/predcat2)-(obscat2-predcat2)/predcat2));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,-(1/2)*log(2*obscat1^2*pi*psi1)-(1/(2*psi1))*(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(1/psi2)*(log(obscat2)-log(psi2*predcat2))-log(obscat2)-sealen*lgamma(1/psi2)-obscat2/(psi2*predcat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #35/49
                  else if(sum(distr==c('aplnormal','poisson')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- 1;
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,2*(obscat2*log(obscat2)-obscat2*log(predcat2)-(obscat2-predcat2)));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(predcat2==0,0,obscat2*log(predcat2) - predcat2 - lfactorial(obscat2));
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #36/49
                  else if(sum(distr==c('aplnormal','negbin')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0,2*log(1+predcat2),
                                                      2*(obscat2*log(obscat2/predcat2)-(obscat2+1)*log((1+obscat2)/(1+predcat2))));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(predcat2==0,0,
                                           obscat2*log(1/psi2) +
                                           obscat2*log(predcat2) -
                                           (obscat2 + psi2)*log(1 + predcat2/psi2) +
                                           lgamma(obscat2 + psi2) -
                                           lgamma(obscat2 + 1) -
                                           lgamma(psi2));
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #37/49
                  else if(sum(distr==c('aplnormal','normal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ((obscat2-predcat2)^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- -(sealen/2)*log(2*pi*psi2)-(1/(2*psi2))*(obscat2-predcat2)^2;
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #38/49
                  else if(sum(distr==c('aplnormal','apnormal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- NA;
                     dev1        <- NA;
                     dev2        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- (obscat2-predcat2)^2;
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                     res2        <- obscat2-predcat2;
                    }
                  #39/49
                  else if(sum(distr==c('aplnormal','lognormal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,-(1/2)*log(2*obscat2^2*pi*psi2)-(1/(2*psi2))*(log(obscat2)-log(predcat2))^2);
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #40/49
                  else if(sum(distr==c('aplnormal','aplnormal')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- NA;
                     dev1        <- NA;
                     dev2        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                     res2        <- ifelse(obscat2==0 | predcat2==0,0,log(obscat2)-log(predcat2));
                    }
                  #41/49
                  else if(sum(distr==c('aplnormal','gamma')) == 2)
                    {
                     psi1        <- NA;
                     psi2        <- exp(tail(par,1));
                     dev1        <- NA;
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(log(obscat2/predcat2)-(obscat2-predcat2)/predcat2));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(log(obscat1)-log(predcat1))^2);
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(1/psi2)*(log(obscat2)-log(psi2*predcat2))-log(obscat2)-sealen*lgamma(1/psi2)-obscat2/(psi2*predcat2));
                     res1        <- ifelse(obscat1==0 | predcat1==0,0,log(obscat1)-log(predcat1));
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #42/49
                  else if(sum(distr==c('gamma','poisson')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- 1;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,2*(obscat2*log(obscat2)-obscat2*log(predcat2)-(obscat2-predcat2)));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     likcontr2   <- ifelse(predcat2==0,0,obscat2*log(predcat2) - predcat2 - lfactorial(obscat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #43/49
                  else if(sum(distr==c('gamma','negbin')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     dev2        <- ifelse(obscat2==0,2*log(1+predcat2),
                                                      2*(obscat2*log(obscat2/predcat2)-(obscat2+1)*log((1+obscat2)/(1+predcat2))));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     likcontr2   <- ifelse(predcat2==0,0,
                                           obscat2*log(1/psi2) +
                                           obscat2*log(predcat2) -
                                           (obscat2 + psi2)*log(1 + predcat2/psi2) +
                                           lgamma(obscat2 + psi2) -
                                           lgamma(obscat2 + 1) -
                                           lgamma(psi2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #44/49
                  else if(sum(distr==c('gamma','normal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     dev2        <- ((obscat2-predcat2)^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     likcontr2   <- -(sealen/2)*log(2*pi*psi2)-(1/(2*psi2))*(obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #45/49
                  else if(sum(distr==c('gamma','apnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     dev2        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     likcontr2   <- (obscat2-predcat2)^2;
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- obscat2-predcat2;
                    }
                  #46/49
                  else if(sum( distr==c('gamma','lognormal')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,-(1/2)*log(2*obscat2^2*pi*psi2)-(1/(2*psi2))*(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #47/49
                  else if(sum(distr==c('gamma','aplnormal')) == 2)
                    {
                     psi1        <- exp(tail(par,1));
                     psi2        <- NA;
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     dev2        <- NA;
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(log(obscat2)-log(predcat2))^2);
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- ifelse(obscat2==0 | predcat2==0,0,log(obscat2)-log(predcat2));
                    }
                  #48/49
                  else if(sum(distr==c('gamma','gamma')) == 2)
                    {
                     psi1        <- exp(tail(par,2)[1]);
                     psi2        <- exp(tail(par,1));
                     dev1        <- ifelse(obscat1==0 | predcat1==0,0,-2*(log(obscat1/predcat1)-(obscat1-predcat1)/predcat1));
                     dev2        <- ifelse(obscat2==0 | predcat2==0,0,-2*(log(obscat2/predcat2)-(obscat2-predcat2)/predcat2));
                     likcontr1   <- ifelse(obscat1==0 | predcat1==0,0,(1/psi1)*(log(obscat1)-log(psi1*predcat1))-log(obscat1)-sealen*lgamma(1/psi1)-obscat1/(psi1*predcat1));
                     likcontr2   <- ifelse(obscat2==0 | predcat2==0,0,(1/psi2)*(log(obscat2)-log(psi2*predcat2))-log(obscat2)-sealen*lgamma(1/psi2)-obscat2/(psi2*predcat2));
                     res1        <- sign(obscat1-predcat1)*sqrt(dev1);
                     res2        <- sign(obscat2-predcat2)*sqrt(dev2);
                    }
                  #49/49
                  Likel[[1]] <- c(psi1,psi2)
                  Likel[[2]] <- data.frame(fleet1=dev1,fleet2=dev2)
                  Likel[[3]] <- data.frame(fleet1=likcontr1,fleet2=likcontr2)
                  Likel[[4]] <- data.frame(fleet1=res1,fleet2=res2)
                  return(Likel);
  }
.CatDynExp1F.Res <-
function(properties,nstep,obsmbm1,ppp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
  {
                  if(properties$Units[2] == "ind")
                    {
                     biom <- NA
                    }
                  else
                    {
                     biom <- if(properties$Units[4] == "bill" && properties$Units[3] == "g")
                               {
                                obsmbm1*1e-6*nstep*1e9;
                               }
                             else if(properties$Units[4] == "bill" && properties$Units[3] == "kg")
                               {
                                obsmbm1*1e-3*nstep*1e9;
                               }
                             else if(properties$Units[4] == "mill" && properties$Units[3] == "g")
                               {
                                obsmbm1*1e-6*nstep*1e6;
                               }
                             else if(properties$Units[4] == "mill" && properties$Units[3] == "kg")
                               {
                               obsmbm1*1e-3*nstep*1e6;
                               }
                             else if(properties$Units[4] == "thou" && properties$Units[3] == "g")
                               {
                               obsmbm1*1e-6*nstep*1e3;
                               }
                             else if(properties$Units[4] == "thou" && properties$Units[3] == "kg")
                               {
                               obsmbm1*1e-3*nstep*1e3;
                               }
                             else if(properties$Units[4] == "hund" && properties$Units[3] == "g")
                               {
                               obsmbm1*1e-6*nstep*1e2;
                               }
                             else
                               {
                               obsmbm1*1e-3*nstep*1e2;
                               }
                    }
                  catdynexp                         <- vector("list",2);
                  names(catdynexp)                  <- c("Properties","Model");
                  catdynexp$Properties              <- properties;
                  catdynexp$Model                   <- vector("list",5);
                  names(catdynexp$Model)            <- c("Type","Dates","Distribution","Parameters","Results");
                  catdynexp$Model$Type              <- ppp;
                  catdynexp$Model$Dates             <- dates;
                  catdynexp$Model$Distribution      <- distr;
                  catdynexp$Model$Parameters        <- exp(par);
                  if(ppp==0)
                    {
                     if(distr == "poisson" | distr == "apnormal" | distr == "aplnormal")
                       {
                        names(catdynexp$Model$Parameters) <- c("M","N0","k","alpha","beta")
                       }
                     else
                       {
                        names(catdynexp$Model$Parameters) <- c("M","N0","k","alpha","beta","psi")
                       }
                    }
                  else if(ppp>0)
                    {
                     if(distr == "poisson" | distr == "apnormal" | distr == "aplnormal")
                       {
                        names(catdynexp$Model$Parameters) <- c("M","N0",paste("P",1:ppp,sep=""),"k","alpha","beta")
                       }
                     else
                       {
                        names(catdynexp$Model$Parameters) <- c("M","N0",paste("P",1:ppp,sep=""),"k","alpha","beta","psi")
                       }
                    }
                  else
                    {
                     ppp                               <- abs(ppp);
                     if(distr == "poisson" | distr == "apnormal" | distr == "aplnormal")
                       {
                        names(catdynexp$Model$Parameters) <- c("M","N0",paste("P",1:ppp,sep=""),"k","alpha","beta")
                       }
                     else
                       {
                        names(catdynexp$Model$Parameters) <- c("M","N0",paste("P",1:ppp,sep=""),"k","alpha","beta","psi")
                       }
                    }
                  catdynexp$Model$Results           <- data.frame(period=ts.start:ts.end,
                                                                  obseff1=obseff1,
                                                                  obscat1=obscat1,
                                                                  modcat=predcat1,
                                                                  deviance=Likel[[2]],
                                                                  likcontr=Likel[[3]],
                                                                  devresid=Likel[[4]],
                                                                  npred=nstep,
                                                                  biompred=biom);
                  names(catdynexp$Model$Results)    <- c(paste("Period (",properties$Units[1],")",sep=""),
                                                         paste("Effort ",properties$Fleets[1,1]," Fleet (",properties$Fleets[1,2],")",sep=""),
                                                         paste("Obs. Catch ",properties$Fleets[1,1]," Fleet (",properties$Units[4],")",sep=""),
                                                         paste("Pred. Catch ",properties$Fleets[1,1]," Fleet (",properties$Units[4],")",sep=""),
                                                         "Deviance",
                                                         "Likelihood",
                                                         "Deviance Residuals",
                                                         paste("Pred. Abundance (",properties$Units[4],")",sep=""),
                                                         "Pred. Biomass (tonnes)")
                  return(catdynexp);
  }
.CatDynExp2F.Res <-
function(properties,nstep,obsmbm1,obsmbm2,ppp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
  {
                  sdistr.set <- c("poisson","apnormal","aplnormal");
                  fleet.name <- properties$Fleets$Fleet
                  if(properties$Units[2] == "ind")
                    {
                     biom <- NA
                    }
                  else
                    {
                     biom <- if(properties$Units[4] == "bill" && properties$Units[3] == "g")
                               {
                                (obsmbm1+obsmbm2)*0.5*1e-6*nstep*1e9;
                               }
                             else if(properties$Units[4] == "bill" && properties$Units[3] == "kg")
                               {
                                (obsmbm1+obsmbm2)*0.5*1e-3*nstep*1e9;
                               }
                             else if(properties$Units[4] == "mill" && properties$Units[3] == "g")
                               {
                                (obsmbm1+obsmbm2)*0.5*1e-6*nstep*1e6;
                               }
                             else if(properties$Units[4] == "mill" && properties$Units[3] == "kg")
                               {
                               (obsmbm1+obsmbm2)*0.5*1e-3*nstep*1e6;
                               }
                             else if(properties$Units[4] == "thou" && properties$Units[3] == "g")
                               {
                               (obsmbm1+obsmbm2)*0.5*1e-6*nstep*1e3;
                               }
                             else if(properties$Units[4] == "thou" && properties$Units[3] == "kg")
                               {
                               (obsmbm1+obsmbm2)*0.5*1e-3*nstep*1e3;
                               }
                             else if(properties$Units[4] == "hund" && properties$Units[3] == "g")
                               {
                               (obsmbm1+obsmbm2)*0.5*1e-6*nstep*1e2;
                               }
                             else
                               {
                               (obsmbm1+obsmbm2)*0.5*1e-3*nstep*1e2;
                               }
                    }
                  catdynexp                         <- vector("list",2);
                  names(catdynexp)                  <- c("Properties","Model");
                  catdynexp$Properties              <- properties;
                  catdynexp$Model                   <- vector("list",5);
                  names(catdynexp$Model)            <- c("Type","Dates","Distribution","Parameters","Results");
                  catdynexp$Model$Type              <- ppp;
                  catdynexp$Model$Dates             <- dates;
                  catdynexp$Model$Distribution      <- distr;
                  catdynexp$Model$Parameters        <- exp(par);
                  if(ppp == c(0,0))
                    {
                     if(sum(distr%in%sdistr.set) == 2)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""));
                       }
                     if(sum(distr%in%sdistr.set) == 1)
                       {
                        if(distr[1]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi2.",fleet.name[2],sep=""));
                          }
                        if(distr[2]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi1.",fleet.name[1],sep=""));
                          }
                       }
                     if(sum(distr%in%sdistr.set) == 0)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                paste("psi1.",fleet.name[1],sep=""),paste("psi2.",fleet.name[2],sep=""));
                       }
                    }
                  else if(ppp == c(0,1))
                    {
                     if(sum(distr%in%sdistr.set) == 2)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""));
                       }
                     if(sum(distr%in%sdistr.set) == 1)
                       {
                        if(distr[1]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi2.",fleet.name[2],sep=""));
                          }
                        if(distr[2]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi1.",fleet.name[1],sep=""));
                          }
                       }
                     if(sum(distr%in%sdistr.set) == 0)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                paste("psi1.",fleet.name[1],sep=""),paste("psi2.",fleet.name[2],sep=""));
                       }
                    }
                  else if(ppp == c(0,2))
                    {
                     if(sum(distr%in%sdistr.set) == 2)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""));
                       }
                     if(sum(distr%in%sdistr.set) == 1)
                       {
                        if(distr[1]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi2.",fleet.name[2],sep=""));
                          }
                        if(distr[2]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi1.",fleet.name[1],sep=""));
                          }
                       }
                     if(sum(distr%in%sdistr.set) == 0)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                paste("psi1.",fleet.name[1],sep=""),paste("psi2.",fleet.name[2],sep=""));
                       }
                    }
                  else if(ppp == c(0,3))
                    {
                     if(sum(distr%in%sdistr.set) == 2)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""));
                       }
                     if(sum(distr%in%sdistr.set) == 1)
                       {
                        if(distr[1]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi2.",fleet.name[2],sep=""));
                          }
                        if(distr[2]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi1.",fleet.name[1],sep=""));
                          }
                       }
                     if(sum(distr%in%sdistr.set) == 0)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                paste("psi1.",fleet.name[1],sep=""),paste("psi2.",fleet.name[2],sep=""));
                       }
                    }
                  else if(ppp == c(0,4))
                    {
                     if(sum(distr%in%sdistr.set) == 2)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""));
                       }
                     if(sum(distr%in%sdistr.set) == 1)
                       {
                        if(distr[1]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi2.",fleet.name[2],sep=""));
                          }
                        if(distr[2]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi1.",fleet.name[1],sep=""));
                          }
                       }
                     if(sum(distr%in%sdistr.set) == 0)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                paste("psi1.",fleet.name[1],sep=""),paste("psi2.",fleet.name[2],sep=""));
                       }
                    }
                  else if(ppp == c(0,5))
                    {
                     if(sum(distr%in%sdistr.set) == 2)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","P5.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""));
                       }
                     if(sum(distr%in%sdistr.set) == 1)
                       {
                        if(distr[1]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","P5.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi2.",fleet.name[2],sep=""));
                          }
                        if(distr[2]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","P5.","k.","alpha.","beta."),
                                                         c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                   paste("psi1.",fleet.name[1],sep=""));
                          }
                       }
                     if(sum(distr%in%sdistr.set) == 0)
                       {
                        par.names <- c("M","N0",paste(c("k.","alpha.","beta.","P1.","P2.","P3.","P4.","P5.","k.","alpha.","beta."),
                                                      c(rep(fleet.name[1],3+ppp[1]),rep(fleet.name[2],3+ppp[2])),sep=""),
                                                paste("psi1.",fleet.name[1],sep=""),paste("psi2.",fleet.name[2],sep=""));
                       }
                    }
                  else
                    {
                     if(sum(distr%in%sdistr.set) == 2)
                       {
                        par.names <- c("M","N0",paste("P",1:ppp[1],".",fleet.name[1],sep=""),
                                                paste(c("k.","alpha.","beta."),fleet.name[1],sep=""),
                                                paste("P",1:ppp[2],".",fleet.name[2],sep=""),
                                                paste(c("k.","alpha.","beta."),fleet.name[2],sep=""));
                       }
                     if(sum(distr%in%sdistr.set) == 1)
                       {
                        if(distr[1]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste("P",1:ppp[1],".",fleet.name[1],sep=""),
                                                   paste(c("k.","alpha.","beta."),fleet.name[1],sep=""),
                                                   paste("P",1:ppp[2],".",fleet.name[2],sep=""),
                                                   paste(c("k.","alpha.","beta."),fleet.name[2],sep=""),
                                                   paste("psi2.",fleet.name[2],sep=""));
                          }
                        if(distr[2]%in%sdistr.set)
                          {
                           par.names <- c("M","N0",paste("P",1:ppp[1],".",fleet.name[1],sep=""),
                                                   paste(c("k.","alpha.","beta."),fleet.name[1],sep=""),
                                                   paste("P",1:ppp[2],".",fleet.name[2],sep=""),
                                                   paste(c("k.","alpha.","beta."),fleet.name[2],sep=""),
                                                   paste("psi1.",fleet.name[1],sep=""));
                          }
                       }
                     if(sum(distr%in%sdistr.set) == 0)
                       {
                        par.names <- c("M","N0",paste("P",1:ppp[1],".",fleet.name[1],sep=""),
                                                paste(c("k.","alpha.","beta."),fleet.name[1],sep=""),
                                                paste("P",1:ppp[2],".",fleet.name[2],sep=""),
                                                paste(c("k.","alpha.","beta."),fleet.name[2],sep=""),
                                                paste("psi.",fleet.name[1],sep=""),paste("psi.",fleet.name[2],sep=""));
                       }
                    }
                  names(catdynexp$Model$Parameters) <- par.names;
                  catdynexp$Model$Results           <- data.frame(period=ts.start:ts.end,
                                                                  obseff1=obseff1,
                                                                  obscat1=obscat1,
                                                                  modcat1=predcat1,
                                                                  deviance1=Likel[[2]][1],
                                                                  likcontr1=Likel[[3]][1],
                                                                  devresid1=Likel[[4]][1],
                                                                  obseff2=obseff2,
                                                                  obscat2=obscat2,
                                                                  modcat2=predcat2,
                                                                  deviance2=Likel[[2]][2],
                                                                  likcontr2=Likel[[3]][2],
                                                                  devresid2=Likel[[4]][2],
                                                                  npred=nstep,
                                                                  biompred=biom);
                  names(catdynexp$Model$Results) <- c(paste("Period (",properties$Units[1],")",sep=""),
                                                      paste("Effort ",properties$Fleets[1,1]," Fleet (",properties$Fleets[1,2],")",sep=""),
                                                      paste("Obs. Catch ",properties$Fleets[1,1]," Fleet (",properties$Units[4],")",sep=""),
                                                      paste("Pred. Catch ",properties$Fleets[1,1]," Fleet (",properties$Units[4],")",sep=""),
                                                      paste("Deviance ",properties$Fleets[1,1]," Fleet (",properties$Units[4],")",sep=""),
                                                      paste("Likelihood ",properties$Fleets[1,1]," Fleet (",properties$Units[4],")",sep=""),
                                                      paste("Deviance Residuals ",properties$Fleets[1,1]," Fleet (",properties$Units[4],")",sep=""),
                                                      paste("Effort ",properties$Fleets[2,1]," Fleet (",properties$Fleets[2,2],")",sep=""),
                                                      paste("Obs. Catch ",properties$Fleets[2,1]," Fleet (",properties$Units[4],")",sep=""),
                                                      paste("Pred. Catch ",properties$Fleets[2,1]," Fleet (",properties$Units[4],")",sep=""),
                                                      paste("Deviance ",properties$Fleets[2,1]," Fleet (",properties$Units[4],")",sep=""),
                                                      paste("Likelihood ",properties$Fleets[2,1]," Fleet (",properties$Units[4],")",sep=""),
                                                      paste("Deviance Residuals ",properties$Fleets[2,1]," Fleet (",properties$Units[4],")",sep=""),
                                                      paste("Pred. Abundance (",properties$Units[4],")",sep=""),
                                                      "Pred. Biomass (tonnes)")
                  return(catdynexp);
  }
.CDMN0P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=0)
  {
                  ts.start  <- head(dates,1);
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logscale  <- par[3];
                  logalpha  <- par[4];
                  logbeta   <- par[5];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                    mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                    nstep[i] <- exp(logN0)*exp(-exp(logM)*i) - mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN1P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=1)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logscale  <- par[4];
                  logalpha  <- par[5];
                  logbeta   <- par[6];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                    mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                    nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) -
                                mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN2P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=2)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logscale  <- par[5];
                  logalpha  <- par[6];
                  logbeta   <- par[7];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                    mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                    nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) -
                                mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN3P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=3)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                   ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logscale  <- par[6];
                  logalpha  <- par[7];
                  logbeta   <- par[8];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen){
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN4P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=4)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logscale  <- par[7];
                  logalpha  <- par[8];
                  logbeta   <- par[9];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN5P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=5)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logscale  <- par[8];
                  logalpha  <- par[9];
                  logbeta   <- par[10];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN6P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=6)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logscale  <- par[9];
                  logalpha  <- par[10];
                  logbeta   <- par[11];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN7P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=7)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logscale  <- par[10];
                  logalpha  <- par[11];
                  logbeta   <- par[12];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN8P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=8)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logscale  <- par[11];
                  logalpha  <- par[12];
                  logbeta   <- par[13];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN9P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=9)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logscale  <- par[12];
                  logalpha  <- par[13];
                  logbeta   <- par[14];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN10P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=10)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logscale  <- par[13];
                  logalpha  <- par[14];
                  logbeta   <- par[15];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN11P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=11)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logscale  <- par[14];
                  logalpha  <- par[15];
                  logbeta   <- par[16];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN12P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=12)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logscale  <- par[15];
                  logalpha  <- par[16];
                  logbeta   <- par[17];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN13P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=13)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logscale  <- par[16];
                  logalpha  <- par[17];
                  logbeta   <- par[18];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN14P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=14)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logscale  <- par[17];
                  logalpha  <- par[18];
                  logbeta   <- par[19];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN15P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=15)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logscale  <- par[18];
                  logalpha  <- par[19];
                  logbeta   <- par[20];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN16P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=16)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.P16    <- dates[17];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logscale  <- par[19];
                  logalpha  <- par[20];
                  logbeta   <- par[21];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN17P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=17)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.P16    <- dates[17];
                  ts.P17    <- dates[18];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logscale  <- par[20];
                  logalpha  <- par[21];
                  logbeta   <- par[22];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start)+1)) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN18P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=18)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.P16    <- dates[17];
                  ts.P17    <- dates[18];
                  ts.P18    <- dates[19];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logP18    <- par[20];
                  logscale  <- par[21];
                  logalpha  <- par[22];
                  logbeta   <- par[23];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start)+1)) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start)+1)) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN19P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=19)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.P16    <- dates[17];
                  ts.P17    <- dates[18];
                  ts.P18    <- dates[19];
                  ts.P19    <- dates[20];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logP18    <- par[20];
                  logP19    <- par[21];
                  logscale  <- par[22];
                  logalpha  <- par[23];
                  logbeta   <- par[24];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start)+1)) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start)+1)) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start)+1)) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN20P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=20)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.P2     <- dates[3];
                  ts.P3     <- dates[4];
                  ts.P4     <- dates[5];
                  ts.P5     <- dates[6];
                  ts.P6     <- dates[7];
                  ts.P7     <- dates[8];
                  ts.P8     <- dates[9];
                  ts.P9     <- dates[10];
                  ts.P10    <- dates[11];
                  ts.P11    <- dates[12];
                  ts.P12    <- dates[13];
                  ts.P13    <- dates[14];
                  ts.P14    <- dates[15];
                  ts.P15    <- dates[16];
                  ts.P16    <- dates[17];
                  ts.P17    <- dates[18];
                  ts.P18    <- dates[19];
                  ts.P19    <- dates[20];
                  ts.P20    <- dates[21];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start), 0, 1);
                  ind.P20   <- ifelse(1:sealen < (ts.P20-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logP18    <- par[20];
                  logP19    <- par[21];
                  logP20    <- par[22];
                  logscale  <- par[23];
                  logalpha  <- par[24];
                  logbeta   <- par[25];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start)+1)) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start)+1)) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start)+1)) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start)+1)) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start)+1)) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start)+1)) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start)+1)) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start)+1)) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start)+1)) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start)+1)) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start)+1)) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start)+1)) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start)+1)) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start)+1)) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start)+1)) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start)+1)) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start)+1)) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start)+1)) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start)+1)) +
                                 ind.P20[i]*exp(logP20)*exp(-exp(logM)*(i-(ts.P20-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha));
                  effn1     <- nstep^(exp(logbeta));
                  predcat1  <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT1P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-1)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logscale  <- par[4];
                  logalpha  <- par[5];
                  logbeta   <- par[6];
                  logtetha  <- par[7];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i]    <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT2P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-2)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logscale  <- par[5];
                  logalpha  <- par[6];
                  logbeta   <- par[7];
                  logtetha  <- par[8];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT3P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-3)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logscale  <- par[6];
                  logalpha  <- par[7];
                  logbeta   <- par[8];
                  logtetha  <- par[9];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT4P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-4)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logscale  <- par[7];
                  logalpha  <- par[8];
                  logbeta   <- par[9];
                  logtetha  <- par[10];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT5P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-5)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logscale  <- par[8];
                  logalpha  <- par[9];
                  logbeta   <- par[10];
                  logtetha  <- par[11];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)))
                                          - ind.N5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT6P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-6)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logscale  <- par[9];
                  logalpha  <- par[10];
                  logbeta   <- par[11];
                  logtetha  <- par[12];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)))
                                          - ind.N5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1)))
                                          - ind.N6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT7P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-7)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logscale  <- par[10];
                  logalpha  <- par[11];
                  logbeta   <- par[12];
                  logtetha  <- par[13];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)))
                                          - ind.N5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1)))
                                          - ind.N6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1)))
                                          - ind.N7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT8P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-8)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logscale  <- par[11];
                  logalpha  <- par[12];
                  logbeta   <- par[13];
                  logtetha  <- par[14];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)))
                                          - ind.N5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1)))
                                          - ind.N6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1)))
                                          - ind.N7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1)))
                                          - ind.N8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT9P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-9)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logscale  <- par[12];
                  logalpha  <- par[13];
                  logbeta   <- par[14];
                  logtetha  <- par[15];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)))
                                          - ind.N5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1)))
                                          - ind.N6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1)))
                                          - ind.N7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1)))
                                          - ind.N8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1)))
                                          - ind.N9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT10P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-10)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logscale  <- par[13];
                  logalpha  <- par[14];
                  logbeta   <- par[15];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)))
                                          - ind.N5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1)))
                                          - ind.N6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1)))
                                          - ind.N7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1)))
                                          - ind.N8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1)))
                                          - ind.N9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1)))
                                          - ind.N10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT11P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-11)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logscale  <- par[14];
                  logalpha  <- par[15];
                  logbeta   <- par[16];
                  logtetha  <- par[17];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)))
                                          - ind.N5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1)))
                                          - ind.N6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1)))
                                          - ind.N7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1)))
                                          - ind.N8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1)))
                                          - ind.N9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1)))
                                          - ind.N10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1)))
                                          - ind.N11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT12P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-12)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logscale  <- par[15];
                  logalpha  <- par[16];
                  logbeta   <- par[17];
                  logtetha  <- par[18];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)))
                                          - ind.N5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1)))
                                          - ind.N6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1)))
                                          - ind.N7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1)))
                                          - ind.N8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1)))
                                          - ind.N9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1)))
                                          - ind.N10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1)))
                                          - ind.N11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1)))
                                          - ind.N12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT13P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-13)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logscale  <- par[16];
                  logalpha  <- par[17];
                  logbeta   <- par[18];
                  logtetha  <- par[19];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)))
                                          - ind.N5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1)))
                                          - ind.N6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1)))
                                          - ind.N7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1)))
                                          - ind.N8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1)))
                                          - ind.N9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1)))
                                          - ind.N10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1)))
                                          - ind.N11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1)))
                                          - ind.N12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1)))
                                          - ind.N13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT14P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-14)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logscale  <- par[17];
                  logalpha  <- par[18];
                  logbeta   <- par[19];
                  logtetha  <- par[20];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)))
                                          - ind.N5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1)))
                                          - ind.N6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1)))
                                          - ind.N7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1)))
                                          - ind.N8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1)))
                                          - ind.N9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1)))
                                          - ind.N10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1)))
                                          - ind.N11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1)))
                                          - ind.N12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1)))
                                          - ind.N13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1)))
                                          - ind.N14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT15P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-15)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logscale  <- par[18];
                  logalpha  <- par[19];
                  logbeta   <- par[20];
                  logtetha  <- par[21];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)))
                                          - ind.N5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1)))
                                          - ind.N6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1)))
                                          - ind.N7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1)))
                                          - ind.N8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1)))
                                          - ind.N9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1)))
                                          - ind.N10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1)))
                                          - ind.N11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1)))
                                          - ind.N12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1)))
                                          - ind.N13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1)))
                                          - ind.N14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1)))
                                          - ind.N15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT16P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-16)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.P16    <- dates[32];
                  ts.N16    <- dates[33];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start+1), 0, 1);
                  ind.N16   <- ifelse(1:sealen < (ts.N16-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logscale  <- par[19];
                  logalpha  <- par[20];
                  logbeta   <- par[21];
                  logtetha  <- par[22];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1))) -
                                mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)))
                                          - ind.N5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1)))
                                          - ind.N6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1)))
                                          - ind.N7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1)))
                                          - ind.N8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1)))
                                          - ind.N9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1)))
                                          - ind.N10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1)))
                                          - ind.N11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1)))
                                          - ind.N12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1)))
                                          - ind.N13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1)))
                                          - ind.N14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1)))
                                          - ind.N15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1)))
                                          - ind.N16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT17P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-17)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.P16    <- dates[32];
                  ts.N16    <- dates[33];
                  ts.P17    <- dates[34];
                  ts.N17    <- dates[35];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start+1), 0, 1);
                  ind.N16   <- ifelse(1:sealen < (ts.N16-ts.start+1), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start+1), 0, 1);
                  ind.N17   <- ifelse(1:sealen < (ts.N17-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logscale  <- par[20];
                  logalpha  <- par[21];
                  logbeta   <- par[22];
                  logtetha  <- par[23];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1))) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)))
                                          - ind.N5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1)))
                                          - ind.N6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1)))
                                          - ind.N7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1)))
                                          - ind.N8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1)))
                                          - ind.N9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1)))
                                          - ind.N10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1)))
                                          - ind.N11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1)))
                                          - ind.N12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1)))
                                          - ind.N13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1)))
                                          - ind.N14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1)))
                                          - ind.N15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1)))
                                          - ind.N16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1)))
                                          - ind.N17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT18P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-18)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.P16    <- dates[32];
                  ts.N16    <- dates[33];
                  ts.P17    <- dates[34];
                  ts.N17    <- dates[35];
                  ts.P18    <- dates[36];
                  ts.N18    <- dates[37];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start+1), 0, 1);
                  ind.N16   <- ifelse(1:sealen < (ts.N16-ts.start+1), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start+1), 0, 1);
                  ind.N17   <- ifelse(1:sealen < (ts.N17-ts.start+1), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start+1), 0, 1);
                  ind.N18   <- ifelse(1:sealen < (ts.N18-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logP18    <- par[20];
                  logscale  <- par[21];
                  logalpha  <- par[22];
                  logbeta   <- par[23];
                  logtetha  <- par[24];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1))) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1))) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)))
                                          - ind.N5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1)))
                                          - ind.N6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1)))
                                          - ind.N7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1)))
                                          - ind.N8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1)))
                                          - ind.N9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1)))
                                          - ind.N10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1)))
                                          - ind.N11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1)))
                                          - ind.N12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1)))
                                          - ind.N13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1)))
                                          - ind.N14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1)))
                                          - ind.N15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1)))
                                          - ind.N16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1)))
                                          - ind.N17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1)))
                                          - ind.N18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT19P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-19)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.P16    <- dates[32];
                  ts.N16    <- dates[33];
                  ts.P17    <- dates[34];
                  ts.N17    <- dates[35];
                  ts.P18    <- dates[36];
                  ts.N18    <- dates[37];
                  ts.P19    <- dates[38];
                  ts.N19    <- dates[39];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start+1), 0, 1);
                  ind.N16   <- ifelse(1:sealen < (ts.N16-ts.start+1), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start+1), 0, 1);
                  ind.N17   <- ifelse(1:sealen < (ts.N17-ts.start+1), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start+1), 0, 1);
                  ind.N18   <- ifelse(1:sealen < (ts.N18-ts.start+1), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start+1), 0, 1);
                  ind.N19   <- ifelse(1:sealen < (ts.N19-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logP18    <- par[20];
                  logP19    <- par[21];
                  logscale  <- par[22];
                  logalpha  <- par[23];
                  logbeta   <- par[24];
                  logtetha  <- par[25];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1))) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1))) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start+1))) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)))
                                          - ind.N5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1)))
                                          - ind.N6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1)))
                                          - ind.N7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1)))
                                          - ind.N8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1)))
                                          - ind.N9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1)))
                                          - ind.N10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1)))
                                          - ind.N11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1)))
                                          - ind.N12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1)))
                                          - ind.N13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1)))
                                          - ind.N14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1)))
                                          - ind.N15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1)))
                                          - ind.N16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1)))
                                          - ind.N17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1)))
                                          - ind.N18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start+1)))
                                          - ind.N19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMNT20P <-
function(par,dates,obscat1,obseff1,obsmbm1,distr,properties,output,pp=-20)
  {
                  ts.start  <- head(dates,1);
                  ts.P1     <- dates[2];
                  ts.N1     <- dates[3];
                  ts.P2     <- dates[4];
                  ts.N2     <- dates[5];
                  ts.P3     <- dates[6];
                  ts.N3     <- dates[7];
                  ts.P4     <- dates[8];
                  ts.N4     <- dates[9];
                  ts.P5     <- dates[10];
                  ts.N5     <- dates[11];
                  ts.P6     <- dates[12];
                  ts.N6     <- dates[13];
                  ts.P7     <- dates[14];
                  ts.N7     <- dates[15];
                  ts.P8     <- dates[16];
                  ts.N8     <- dates[17];
                  ts.P9     <- dates[18];
                  ts.N9     <- dates[19];
                  ts.P10    <- dates[20];
                  ts.N10    <- dates[21];
                  ts.P11    <- dates[22];
                  ts.N11    <- dates[23];
                  ts.P12    <- dates[24];
                  ts.N12    <- dates[25];
                  ts.P13    <- dates[26];
                  ts.N13    <- dates[27];
                  ts.P14    <- dates[28];
                  ts.N14    <- dates[29];
                  ts.P15    <- dates[30];
                  ts.N15    <- dates[31];
                  ts.P16    <- dates[32];
                  ts.N16    <- dates[33];
                  ts.P17    <- dates[34];
                  ts.N17    <- dates[35];
                  ts.P18    <- dates[36];
                  ts.N18    <- dates[37];
                  ts.P19    <- dates[38];
                  ts.N19    <- dates[39];
                  ts.P20    <- dates[40];
                  ts.N20    <- dates[41];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  resn1     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind.P1    <- ifelse(1:sealen < (ts.P1-ts.start+1), 0, 1);
                  ind.N1    <- ifelse(1:sealen < (ts.N1-ts.start+1), 0, 1);
                  ind.P2    <- ifelse(1:sealen < (ts.P2-ts.start+1), 0, 1);
                  ind.N2    <- ifelse(1:sealen < (ts.N2-ts.start+1), 0, 1);
                  ind.P3    <- ifelse(1:sealen < (ts.P3-ts.start+1), 0, 1);
                  ind.N3    <- ifelse(1:sealen < (ts.N3-ts.start+1), 0, 1);
                  ind.P4    <- ifelse(1:sealen < (ts.P4-ts.start+1), 0, 1);
                  ind.N4    <- ifelse(1:sealen < (ts.N4-ts.start+1), 0, 1);
                  ind.P5    <- ifelse(1:sealen < (ts.P5-ts.start+1), 0, 1);
                  ind.N5    <- ifelse(1:sealen < (ts.N5-ts.start+1), 0, 1);
                  ind.P6    <- ifelse(1:sealen < (ts.P6-ts.start+1), 0, 1);
                  ind.N6    <- ifelse(1:sealen < (ts.N6-ts.start+1), 0, 1);
                  ind.P7    <- ifelse(1:sealen < (ts.P7-ts.start+1), 0, 1);
                  ind.N7    <- ifelse(1:sealen < (ts.N7-ts.start+1), 0, 1);
                  ind.P8    <- ifelse(1:sealen < (ts.P8-ts.start+1), 0, 1);
                  ind.N8    <- ifelse(1:sealen < (ts.N8-ts.start+1), 0, 1);
                  ind.P9    <- ifelse(1:sealen < (ts.P9-ts.start+1), 0, 1);
                  ind.N9    <- ifelse(1:sealen < (ts.N9-ts.start+1), 0, 1);
                  ind.P10   <- ifelse(1:sealen < (ts.P10-ts.start+1), 0, 1);
                  ind.N10   <- ifelse(1:sealen < (ts.N10-ts.start+1), 0, 1);
                  ind.P11   <- ifelse(1:sealen < (ts.P11-ts.start+1), 0, 1);
                  ind.N11   <- ifelse(1:sealen < (ts.N11-ts.start+1), 0, 1);
                  ind.P12   <- ifelse(1:sealen < (ts.P12-ts.start+1), 0, 1);
                  ind.N12   <- ifelse(1:sealen < (ts.N12-ts.start+1), 0, 1);
                  ind.P13   <- ifelse(1:sealen < (ts.P13-ts.start+1), 0, 1);
                  ind.N13   <- ifelse(1:sealen < (ts.N13-ts.start+1), 0, 1);
                  ind.P14   <- ifelse(1:sealen < (ts.P14-ts.start+1), 0, 1);
                  ind.N14   <- ifelse(1:sealen < (ts.N14-ts.start+1), 0, 1);
                  ind.P15   <- ifelse(1:sealen < (ts.P15-ts.start+1), 0, 1);
                  ind.N15   <- ifelse(1:sealen < (ts.N15-ts.start+1), 0, 1);
                  ind.P16   <- ifelse(1:sealen < (ts.P16-ts.start+1), 0, 1);
                  ind.N16   <- ifelse(1:sealen < (ts.N16-ts.start+1), 0, 1);
                  ind.P17   <- ifelse(1:sealen < (ts.P17-ts.start+1), 0, 1);
                  ind.N17   <- ifelse(1:sealen < (ts.N17-ts.start+1), 0, 1);
                  ind.P18   <- ifelse(1:sealen < (ts.P18-ts.start+1), 0, 1);
                  ind.N18   <- ifelse(1:sealen < (ts.N18-ts.start+1), 0, 1);
                  ind.P19   <- ifelse(1:sealen < (ts.P19-ts.start+1), 0, 1);
                  ind.N19   <- ifelse(1:sealen < (ts.N19-ts.start+1), 0, 1);
                  ind.P20   <- ifelse(1:sealen < (ts.P20-ts.start+1), 0, 1);
                  ind.N20   <- ifelse(1:sealen < (ts.N20-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1     <- par[3];
                  logP2     <- par[4];
                  logP3     <- par[5];
                  logP4     <- par[6];
                  logP5     <- par[7];
                  logP6     <- par[8];
                  logP7     <- par[9];
                  logP8     <- par[10];
                  logP9     <- par[11];
                  logP10    <- par[12];
                  logP11    <- par[13];
                  logP12    <- par[14];
                  logP13    <- par[15];
                  logP14    <- par[16];
                  logP15    <- par[17];
                  logP16    <- par[18];
                  logP17    <- par[19];
                  logP18    <- par[20];
                  logP19    <- par[21];
                  logP20    <- par[22];
                  logscale  <- par[23];
                  logalpha  <- par[24];
                  logbeta   <- par[25];
                  logtetha  <- par[26];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  resn1[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind.P1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1))) +
                                 ind.P2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1))) +
                                 ind.P3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1))) +
                                 ind.P4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1))) +
                                 ind.P5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1))) +
                                 ind.P6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1))) +
                                 ind.P7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1))) +
                                 ind.P8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1))) +
                                 ind.P9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1))) +
                                 ind.P10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1))) +
                                 ind.P11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1))) +
                                 ind.P12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1))) +
                                 ind.P13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1))) +
                                 ind.P14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1))) +
                                 ind.P15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1))) +
                                 ind.P16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1))) +
                                 ind.P17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1))) +
                                 ind.P18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start+1))) +
                                 ind.P19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start+1))) -
                                 ind.P20[i]*exp(logP20)*exp(-exp(logM)*(i-(ts.P20-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                     resn1[i] <- nstep[i] - ind.N1[i]*exp(logP1)*exp(-exp(logM)*(i-(ts.P1-ts.start+1)))
                                          - ind.N2[i]*exp(logP2)*exp(-exp(logM)*(i-(ts.P2-ts.start+1)))
                                          - ind.N3[i]*exp(logP3)*exp(-exp(logM)*(i-(ts.P3-ts.start+1)))
                                          - ind.N4[i]*exp(logP4)*exp(-exp(logM)*(i-(ts.P4-ts.start+1)))
                                          - ind.N5[i]*exp(logP5)*exp(-exp(logM)*(i-(ts.P5-ts.start+1)))
                                          - ind.N6[i]*exp(logP6)*exp(-exp(logM)*(i-(ts.P6-ts.start+1)))
                                          - ind.N7[i]*exp(logP7)*exp(-exp(logM)*(i-(ts.P7-ts.start+1)))
                                          - ind.N8[i]*exp(logP8)*exp(-exp(logM)*(i-(ts.P8-ts.start+1)))
                                          - ind.N9[i]*exp(logP9)*exp(-exp(logM)*(i-(ts.P9-ts.start+1)))
                                          - ind.N10[i]*exp(logP10)*exp(-exp(logM)*(i-(ts.P10-ts.start+1)))
                                          - ind.N11[i]*exp(logP11)*exp(-exp(logM)*(i-(ts.P11-ts.start+1)))
                                          - ind.N12[i]*exp(logP12)*exp(-exp(logM)*(i-(ts.P12-ts.start+1)))
                                          - ind.N13[i]*exp(logP13)*exp(-exp(logM)*(i-(ts.P13-ts.start+1)))
                                          - ind.N14[i]*exp(logP14)*exp(-exp(logM)*(i-(ts.P14-ts.start+1)))
                                          - ind.N15[i]*exp(logP15)*exp(-exp(logM)*(i-(ts.P15-ts.start+1)))
                                          - ind.N16[i]*exp(logP16)*exp(-exp(logM)*(i-(ts.P16-ts.start+1)))
                                          - ind.N17[i]*exp(logP17)*exp(-exp(logM)*(i-(ts.P17-ts.start+1)))
                                          - ind.N18[i]*exp(logP18)*exp(-exp(logM)*(i-(ts.P18-ts.start+1)))
                                          - ind.N19[i]*exp(logP19)*exp(-exp(logM)*(i-(ts.P19-ts.start+1)))
                                          - ind.N20[i]*exp(logP20)*exp(-exp(logM)*(i-(ts.P20-ts.start+1)));
                    }
                  effeff1  <- obseff1^(exp(logalpha));
                  effn1    <- resn1^(exp(logbeta));
                  predcat1 <- exp(logscale)*(effeff1*effn1)*exp(-exp(logM)/2);
                  Likel    <- .CatDynLik1F(obscat1,predcat1,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp1F.Res(properties,nstep,obsmbm1,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr == "apnormal" | distr == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]]));
                       }
                     else
                       {
                        negsup <- -sum(Likel[["Likelihood"]]);
                       }
                     return(negsup);
                    }
 }
.CDMN0P0P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(0,0))
  {

                  ts.start  <- head(dates,1);
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effeff2   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logscale1 <- par[3];
                  logalpha1 <- par[4];
                  logbeta1  <- par[5];
                  logscale2 <- par[6];
                  logalpha2 <- par[7];
                  logbeta2  <- par[8];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] +mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) - mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN0P1P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(0,1))
  {
                  ts.start  <- head(dates,1);
                  ts2.P1    <- dates[2];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logscale1 <- par[3];
                  logalpha1 <- par[4];
                  logbeta1  <- par[5];
                  logP1F2   <- par[6];
                  logscale2 <- par[7];
                  logalpha2 <- par[8];
                  logbeta2  <- par[9];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN0P2P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(0,2))
  {
                  ts.start  <- head(dates,1);
                  ts2.P1    <- dates[2];
                  ts2.P2    <- dates[3];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logscale1 <- par[3];
                  logalpha1 <- par[4];
                  logbeta1  <- par[5];
                  logP1F2   <- par[6];
                  logP2F2   <- par[7];
                  logscale2 <- par[8];
                  logalpha2 <- par[9];
                  logbeta2  <- par[10];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN0P3P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(0,3))
  {
                  ts.start  <- head(dates,1);
                  ts2.P1    <- dates[2];
                  ts2.P2    <- dates[3];
                  ts2.P3    <- dates[4];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logscale1 <- par[3];
                  logalpha1 <- par[4];
                  logbeta1  <- par[5];
                  logP1F2   <- par[6];
                  logP2F2   <- par[7];
                  logP3F2   <- par[8];
                  logscale2 <- par[9];
                  logalpha2 <- par[10];
                  logbeta2  <- par[11];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN0P4P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(0,4))
  {
                  ts.start  <- head(dates,1);
                  ts2.P1    <- dates[2];
                  ts2.P2    <- dates[3];
                  ts2.P3    <- dates[4];
                  ts2.P4    <- dates[5];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logscale1 <- par[3];
                  logalpha1 <- par[4];
                  logbeta1  <- par[5];
                  logP1F2   <- par[6];
                  logP2F2   <- par[7];
                  logP3F2   <- par[8];
                  logP4F2   <- par[9];
                  logscale2 <- par[10];
                  logalpha2 <- par[11];
                  logbeta2  <- par[12];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN0P5P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(0,5))
  {
                  ts.start  <- head(dates,1);
                  ts2.P1    <- dates[2];
                  ts2.P2    <- dates[3];
                  ts2.P3    <- dates[4];
                  ts2.P4    <- dates[5];
                  ts2.P5    <- dates[6];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logscale1 <- par[3];
                  logalpha1 <- par[4];
                  logbeta1  <- par[5];
                  logP1F2   <- par[6];
                  logP2F2   <- par[7];
                  logP3F2   <- par[8];
                  logP4F2   <- par[9];
                  logP5F2   <- par[10];
                  logscale2 <- par[11];
                  logalpha2 <- par[12];
                  logbeta2  <- par[13];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                  ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN1P1P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(1,1))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts2.P1    <- dates[3];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start+1), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start+1), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logscale1 <- par[4];
                  logalpha1 <- par[5];
                  logbeta1  <- par[6];
                  logP1F2   <- par[7];
                  logscale2 <- par[8];
                  logalpha2 <- par[9];
                  logbeta2  <- par[10];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start+1))) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start+1))) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN1P2P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(1,2))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts2.P1    <- dates[3];
                  ts2.P2    <- dates[4];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logscale1 <- par[4];
                  logalpha1 <- par[5];
                  logbeta1  <- par[6];
                  logP1F2   <- par[7];
                  logP2F2   <- par[8];
                  logscale2 <- par[9];
                  logalpha2 <- par[10];
                  logbeta2  <- par[11];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN1P3P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(1,3))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts2.P1    <- dates[3];
                  ts2.P2    <- dates[4];
                  ts2.P3    <- dates[5];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logscale1 <- par[4];
                  logalpha1 <- par[5];
                  logbeta1  <- par[6];
                  logP1F2   <- par[7];
                  logP2F2   <- par[8];
                  logP3F2   <- par[9];
                  logscale2 <- par[10];
                  logalpha2 <- par[11];
                  logbeta2  <- par[12];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN1P4P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(1,4))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts2.P1    <- dates[3];
                  ts2.P2    <- dates[4];
                  ts2.P3    <- dates[5];
                  ts2.P4    <- dates[6];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logscale1 <- par[4];
                  logalpha1 <- par[5];
                  logbeta1  <- par[6];
                  logP1F2   <- par[7];
                  logP2F2   <- par[8];
                  logP3F2   <- par[9];
                  logP4F2   <- par[10];
                  logscale2 <- par[11];
                  logalpha2 <- par[12];
                  logbeta2  <- par[13];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN1P5P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(1,5))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts2.P1    <- dates[3];
                  ts2.P2    <- dates[4];
                  ts2.P3    <- dates[5];
                  ts2.P4    <- dates[6];
                  ts2.P5    <- dates[7];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logscale1 <- par[4];
                  logalpha1 <- par[5];
                  logbeta1  <- par[6];
                  logP1F2   <- par[7];
                  logP2F2   <- par[8];
                  logP3F2   <- par[9];
                  logP4F2   <- par[10];
                  logP5F2   <- par[11];
                  logscale2 <- par[12];
                  logalpha2 <- par[13];
                  logbeta2  <- par[14];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                  ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN2P2P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(2,2))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts2.P1    <- dates[4];
                  ts2.P2    <- dates[5];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logscale1 <- par[5];
                  logalpha1 <- par[6];
                  logbeta1  <- par[7];
                  logP1F2   <- par[8];
                  logP2F2   <- par[9];
                  logscale2 <- par[10];
                  logalpha2 <- par[11];
                  logbeta2  <- par[12];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN2P3P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(2,3))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts2.P1    <- dates[4];
                  ts2.P2    <- dates[5];
                  ts2.P3    <- dates[6];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logscale1 <- par[5];
                  logalpha1 <- par[6];
                  logbeta1  <- par[7];
                  logP1F2   <- par[8];
                  logP2F2   <- par[9];
                  logP3F2   <- par[10];
                  logscale2 <- par[11];
                  logalpha2 <- par[12];
                  logbeta2  <- par[13];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN2P4P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(2,4))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts2.P1    <- dates[4];
                  ts2.P2    <- dates[5];
                  ts2.P3    <- dates[6];
                  ts2.P4    <- dates[7];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logscale1 <- par[5];
                  logalpha1 <- par[6];
                  logbeta1  <- par[7];
                  logP1F2   <- par[8];
                  logP2F2   <- par[9];
                  logP3F2   <- par[10];
                  logP4F2   <- par[11];
                  logscale2 <- par[12];
                  logalpha2 <- par[13];
                  logbeta2  <- par[14];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN2P5P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(2,5))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts2.P1    <- dates[4];
                  ts2.P2    <- dates[5];
                  ts2.P3    <- dates[6];
                  ts2.P4    <- dates[7];
                  ts2.P5    <- dates[8];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logscale1 <- par[5];
                  logalpha1 <- par[6];
                  logbeta1  <- par[7];
                  logP1F2   <- par[8];
                  logP2F2   <- par[9];
                  logP3F2   <- par[10];
                  logP4F2   <- par[11];
                  logP5F2   <- par[12];
                  logscale2 <- par[13];
                  logalpha2 <- par[14];
                  logbeta2  <- par[15];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                  ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN3P3P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(3,3))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts2.P1    <- dates[5];
                  ts2.P2    <- dates[6];
                  ts2.P3    <- dates[7];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logscale1 <- par[6];
                  logalpha1 <- par[7];
                  logbeta1  <- par[8];
                  logP1F2   <- par[9];
                  logP2F2   <- par[10];
                  logP3F2   <- par[11];
                  logscale2 <- par[12];
                  logalpha2 <- par[13];
                  logbeta2  <- par[14];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                    {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                    }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN3P4P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(3,4))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts2.P1    <- dates[5];
                  ts2.P2    <- dates[6];
                  ts2.P3    <- dates[7];
                  ts2.P4    <- dates[8];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logscale1 <- par[6];
                  logalpha1 <- par[7];
                  logbeta1  <- par[8];
                  logP1F2   <- par[9];
                  logP2F2   <- par[10];
                  logP3F2   <- par[11];
                  logP4F2   <- par[12];
                  logscale2 <- par[13];
                  logalpha2 <- par[14];
                  logbeta2  <- par[15];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                  ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN3P5P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(3,5))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts2.P1    <- dates[5];
                  ts2.P2    <- dates[6];
                  ts2.P3    <- dates[7];
                  ts2.P4    <- dates[8];
                  ts2.P5    <- dates[9];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logscale1 <- par[6];
                  logalpha1 <- par[7];
                  logbeta1  <- par[8];
                  logP1F2   <- par[9];
                  logP2F2   <- par[10];
                  logP3F2   <- par[11];
                  logP4F2   <- par[12];
                  logP5F2   <- par[13];
                  logscale2 <- par[14];
                  logalpha2 <- par[15];
                  logbeta2  <- par[16];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                  ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                  ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN4P4P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(4,4))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts2.P1    <- dates[6];
                  ts2.P2    <- dates[7];
                  ts2.P3    <- dates[8];
                  ts2.P4    <- dates[9];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logscale1 <- par[7];
                  logalpha1 <- par[8];
                  logbeta1  <- par[9];
                  logP1F2   <- par[10];
                  logP2F2   <- par[11];
                  logP3F2   <- par[12];
                  logP4F2   <- par[13];
                  logscale2 <- par[14];
                  logalpha2 <- par[15];
                  logbeta2  <- par[16];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                  ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                  ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN4P5P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(4,5))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts2.P1    <- dates[6];
                  ts2.P2    <- dates[7];
                  ts2.P3    <- dates[8];
                  ts2.P4    <- dates[9];
                  ts2.P5    <- dates[10];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logscale1 <- par[7];
                  logalpha1 <- par[8];
                  logbeta1  <- par[9];
                  logP1F2   <- par[10];
                  logP2F2   <- par[11];
                  logP3F2   <- par[12];
                  logP4F2   <- par[13];
                  logP5F2   <- par[14];
                  logscale2 <- par[15];
                  logalpha2 <- par[16];
                  logbeta2  <- par[17];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                  ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                  ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                  ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN5P5P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(5,5))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts2.P1    <- dates[7];
                  ts2.P2    <- dates[8];
                  ts2.P3    <- dates[9];
                  ts2.P4    <- dates[10];
                  ts2.P5    <- dates[11];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logscale1 <- par[8];
                  logalpha1 <- par[9];
                  logbeta1  <- par[10];
                  logP1F2   <- par[11];
                  logP2F2   <- par[12];
                  logP3F2   <- par[13];
                  logP4F2   <- par[14];
                  logP5F2   <- par[15];
                  logscale2 <- par[16];
                  logalpha2 <- par[17];
                  logbeta2  <- par[18];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                  ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                  ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                  ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                  ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN6P6P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(6,6))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts2.P1    <- dates[8];
                  ts2.P2    <- dates[9];
                  ts2.P3    <- dates[10];
                  ts2.P4    <- dates[11];
                  ts2.P5    <- dates[12];
                  ts2.P6    <- dates[13];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logscale1 <- par[9];
                  logalpha1 <- par[10];
                  logbeta1  <- par[11];
                  logP1F2   <- par[12];
                  logP2F2   <- par[13];
                  logP3F2   <- par[14];
                  logP4F2   <- par[15];
                  logP5F2   <- par[16];
                  logP6F2   <- par[17];
                  logscale2 <- par[18];
                  logalpha2 <- par[19];
                  logbeta2  <- par[20];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                      mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                      nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                  ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                  ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                  ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                  ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                  ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                  ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                  ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                  ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                  ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                  ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                  ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                  ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) -
                                  mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN7P7P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(7,7))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts2.P1    <- dates[9];
                  ts2.P2    <- dates[10];
                  ts2.P3    <- dates[11];
                  ts2.P4    <- dates[12];
                  ts2.P5    <- dates[13];
                  ts2.P6    <- dates[14];
                  ts2.P7    <- dates[15];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logscale1 <- par[10];
                  logalpha1 <- par[11];
                  logbeta1  <- par[12];
                  logP1F2   <- par[13];
                  logP2F2   <- par[14];
                  logP3F2   <- par[15];
                  logP4F2   <- par[16];
                  logP5F2   <- par[17];
                  logP6F2   <- par[18];
                  logP7F2   <- par[19];
                  logscale2 <- par[20];
                  logalpha2 <- par[21];
                  logbeta2  <- par[22];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN8P8P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(8,8))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts2.P1    <- dates[10];
                  ts2.P2    <- dates[11];
                  ts2.P3    <- dates[12];
                  ts2.P4    <- dates[13];
                  ts2.P5    <- dates[14];
                  ts2.P6    <- dates[15];
                  ts2.P7    <- dates[16];
                  ts2.P8    <- dates[17];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logscale1 <- par[11];
                  logalpha1 <- par[12];
                  logbeta1  <- par[13];
                  logP1F2   <- par[14];
                  logP2F2   <- par[15];
                  logP3F2   <- par[16];
                  logP4F2   <- par[17];
                  logP5F2   <- par[18];
                  logP6F2   <- par[19];
                  logP7F2   <- par[20];
                  logP8F2   <- par[21];
                  logscale2 <- par[22];
                  logalpha2 <- par[23];
                  logbeta2  <- par[24];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN9P9P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(9,9))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts2.P1    <- dates[11];
                  ts2.P2    <- dates[12];
                  ts2.P3    <- dates[13];
                  ts2.P4    <- dates[14];
                  ts2.P5    <- dates[15];
                  ts2.P6    <- dates[16];
                  ts2.P7    <- dates[17];
                  ts2.P8    <- dates[18];
                  ts2.P9    <- dates[19];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logscale1 <- par[12];
                  logalpha1 <- par[13];
                  logbeta1  <- par[14];
                  logP1F2   <- par[15];
                  logP2F2   <- par[16];
                  logP3F2   <- par[17];
                  logP4F2   <- par[18];
                  logP5F2   <- par[19];
                  logP6F2   <- par[20];
                  logP7F2   <- par[21];
                  logP8F2   <- par[22];
                  logP9F2   <- par[23];
                  logscale2 <- par[24];
                  logalpha2 <- par[25];
                  logbeta2  <- par[26];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN10P10P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(10,10))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts2.P1    <- dates[12];
                  ts2.P2    <- dates[13];
                  ts2.P3    <- dates[14];
                  ts2.P4    <- dates[15];
                  ts2.P5    <- dates[16];
                  ts2.P6    <- dates[17];
                  ts2.P7    <- dates[18];
                  ts2.P8    <- dates[19];
                  ts2.P9    <- dates[20];
                  ts2.P10   <- dates[21];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                   ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logscale1 <- par[13];
                  logalpha1 <- par[14];
                  logbeta1  <- par[15];
                  logP1F2   <- par[16];
                  logP2F2   <- par[17];
                  logP3F2   <- par[18];
                  logP4F2   <- par[19];
                  logP5F2   <- par[20];
                  logP6F2   <- par[21];
                  logP7F2   <- par[22];
                  logP8F2   <- par[23];
                  logP9F2   <- par[24];
                  logP10F2  <- par[25];
                  logscale2 <- par[26];
                  logalpha2 <- par[27];
                  logbeta2  <- par[28];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN11P11P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(11,11))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts2.P1    <- dates[13];
                  ts2.P2    <- dates[14];
                  ts2.P3    <- dates[15];
                  ts2.P4    <- dates[16];
                  ts2.P5    <- dates[17];
                  ts2.P6    <- dates[18];
                  ts2.P7    <- dates[19];
                  ts2.P8    <- dates[20];
                  ts2.P9    <- dates[21];
                  ts2.P10   <- dates[22];
                  ts2.P11   <- dates[23];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logscale1 <- par[14];
                  logalpha1 <- par[15];
                  logbeta1  <- par[16];
                  logP1F2   <- par[17];
                  logP2F2   <- par[18];
                  logP3F2   <- par[19];
                  logP4F2   <- par[20];
                  logP5F2   <- par[21];
                  logP6F2   <- par[22];
                  logP7F2   <- par[23];
                  logP8F2   <- par[24];
                  logP9F2   <- par[25];
                  logP10F2  <- par[26];
                  logP11F2  <- par[27];
                  logscale2 <- par[28];
                  logalpha2 <- par[29];
                  logbeta2  <- par[30];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN12P12P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(12,12))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts2.P1    <- dates[14];
                  ts2.P2    <- dates[15];
                  ts2.P3    <- dates[16];
                  ts2.P4    <- dates[17];
                  ts2.P5    <- dates[18];
                  ts2.P6    <- dates[19];
                  ts2.P7    <- dates[20];
                  ts2.P8    <- dates[21];
                  ts2.P9    <- dates[22];
                  ts2.P10   <- dates[23];
                  ts2.P11   <- dates[24];
                  ts2.P12   <- dates[25];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logscale1 <- par[15];
                  logalpha1 <- par[16];
                  logbeta1  <- par[17];
                  logP1F2   <- par[18];
                  logP2F2   <- par[19];
                  logP3F2   <- par[20];
                  logP4F2   <- par[21];
                  logP5F2   <- par[22];
                  logP6F2   <- par[23];
                  logP7F2   <- par[24];
                  logP8F2   <- par[25];
                  logP9F2   <- par[26];
                  logP10F2  <- par[27];
                  logP11F2  <- par[28];
                  logP12F2  <- par[29];
                  logscale2 <- par[30];
                  logalpha2 <- par[31];
                  logbeta2  <- par[32];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN13P13P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(13,13))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts2.P1    <- dates[15];
                  ts2.P2    <- dates[16];
                  ts2.P3    <- dates[17];
                  ts2.P4    <- dates[18];
                  ts2.P5    <- dates[19];
                  ts2.P6    <- dates[20];
                  ts2.P7    <- dates[21];
                  ts2.P8    <- dates[22];
                  ts2.P9    <- dates[23];
                  ts2.P10   <- dates[24];
                  ts2.P11   <- dates[25];
                  ts2.P12   <- dates[26];
                  ts2.P13   <- dates[27];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logscale1 <- par[16];
                  logalpha1 <- par[17];
                  logbeta1  <- par[18];
                  logP1F2   <- par[19];
                  logP2F2   <- par[20];
                  logP3F2   <- par[21];
                  logP4F2   <- par[22];
                  logP5F2   <- par[23];
                  logP6F2   <- par[24];
                  logP7F2   <- par[25];
                  logP8F2   <- par[26];
                  logP9F2   <- par[27];
                  logP10F2  <- par[28];
                  logP11F2  <- par[29];
                  logP12F2  <- par[30];
                  logP13F2  <- par[31];
                  logscale2 <- par[32];
                  logalpha2 <- par[33];
                  logbeta2  <- par[34];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN14P14P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(14,14))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts2.P1    <- dates[16];
                  ts2.P2    <- dates[17];
                  ts2.P3    <- dates[18];
                  ts2.P4    <- dates[19];
                  ts2.P5    <- dates[20];
                  ts2.P6    <- dates[21];
                  ts2.P7    <- dates[22];
                  ts2.P8    <- dates[23];
                  ts2.P9    <- dates[24];
                  ts2.P10   <- dates[25];
                  ts2.P11   <- dates[26];
                  ts2.P12   <- dates[27];
                  ts2.P13   <- dates[28];
                  ts2.P14   <- dates[29];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logscale1 <- par[17];
                  logalpha1 <- par[18];
                  logbeta1  <- par[19];
                  logP1F2   <- par[20];
                  logP2F2   <- par[21];
                  logP3F2   <- par[22];
                  logP4F2   <- par[23];
                  logP5F2   <- par[24];
                  logP6F2   <- par[25];
                  logP7F2   <- par[26];
                  logP8F2   <- par[27];
                  logP9F2   <- par[28];
                  logP10F2  <- par[29];
                  logP11F2  <- par[30];
                  logP12F2  <- par[31];
                  logP13F2  <- par[32];
                  logP14F2  <- par[33];
                  logscale2 <- par[34];
                  logalpha2 <- par[35];
                  logbeta2  <- par[36];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN15P15P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(15,15))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts2.P1    <- dates[17];
                  ts2.P2    <- dates[18];
                  ts2.P3    <- dates[19];
                  ts2.P4    <- dates[20];
                  ts2.P5    <- dates[21];
                  ts2.P6    <- dates[22];
                  ts2.P7    <- dates[23];
                  ts2.P8    <- dates[24];
                  ts2.P9    <- dates[25];
                  ts2.P10   <- dates[26];
                  ts2.P11   <- dates[27];
                  ts2.P12   <- dates[28];
                  ts2.P13   <- dates[29];
                  ts2.P14   <- dates[30];
                  ts2.P15   <- dates[31];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logscale1 <- par[18];
                  logalpha1 <- par[19];
                  logbeta1  <- par[20];
                  logP1F2   <- par[21];
                  logP2F2   <- par[22];
                  logP3F2   <- par[23];
                  logP4F2   <- par[24];
                  logP5F2   <- par[25];
                  logP6F2   <- par[26];
                  logP7F2   <- par[27];
                  logP8F2   <- par[28];
                  logP9F2   <- par[29];
                  logP10F2  <- par[30];
                  logP11F2  <- par[31];
                  logP12F2  <- par[32];
                  logP13F2  <- par[33];
                  logP14F2  <- par[34];
                  logP15F2  <- par[35];
                  logscale2 <- par[36];
                  logalpha2 <- par[37];
                  logbeta2  <- par[38];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN16P16P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(16,16))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts1.P16   <- dates[17];
                  ts2.P1    <- dates[18];
                  ts2.P2    <- dates[19];
                  ts2.P3    <- dates[20];
                  ts2.P4    <- dates[21];
                  ts2.P5    <- dates[22];
                  ts2.P6    <- dates[23];
                  ts2.P7    <- dates[24];
                  ts2.P8    <- dates[25];
                  ts2.P9    <- dates[26];
                  ts2.P10   <- dates[27];
                  ts2.P11   <- dates[28];
                  ts2.P12   <- dates[29];
                  ts2.P13   <- dates[30];
                  ts2.P14   <- dates[31];
                  ts2.P15   <- dates[32];
                  ts2.P16   <- dates[33];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  ind1.P16  <- ifelse(1:sealen < (ts1.P16-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  ind2.P16  <- ifelse(1:sealen < (ts2.P16-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logP16F1  <- par[18];
                  logscale1 <- par[19];
                  logalpha1 <- par[20];
                  logbeta1  <- par[21];
                  logP1F2   <- par[22];
                  logP2F2   <- par[23];
                  logP3F2   <- par[24];
                  logP4F2   <- par[25];
                  logP5F2   <- par[26];
                  logP6F2   <- par[27];
                  logP7F2   <- par[28];
                  logP8F2   <- par[29];
                  logP9F2   <- par[30];
                  logP10F2  <- par[31];
                  logP11F2  <- par[32];
                  logP12F2  <- par[33];
                  logP13F2  <- par[34];
                  logP14F2  <- par[35];
                  logP15F2  <- par[36];
                  logP16F2  <- par[37];
                  logscale2 <- par[38];
                  logalpha2 <- par[39];
                  logbeta2  <- par[40];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind1.P16[i]*exp(logP16F1)*exp(-exp(logM)*(i-(ts1.P16-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) +
                                 ind2.P16[i]*exp(logP16F2)*exp(-exp(logM)*(i-(ts2.P16-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN17P17P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(17,17))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts1.P16   <- dates[17];
                  ts1.P17   <- dates[18];
                  ts2.P1    <- dates[19];
                  ts2.P2    <- dates[20];
                  ts2.P3    <- dates[21];
                  ts2.P4    <- dates[22];
                  ts2.P5    <- dates[23];
                  ts2.P6    <- dates[24];
                  ts2.P7    <- dates[25];
                  ts2.P8    <- dates[26];
                  ts2.P9    <- dates[27];
                  ts2.P10   <- dates[28];
                  ts2.P11   <- dates[29];
                  ts2.P12   <- dates[30];
                  ts2.P13   <- dates[31];
                  ts2.P14   <- dates[32];
                  ts2.P15   <- dates[33];
                  ts2.P16   <- dates[34];
                  ts2.P17   <- dates[35];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  ind1.P16  <- ifelse(1:sealen < (ts1.P16-ts.start), 0, 1);
                  ind1.P17  <- ifelse(1:sealen < (ts1.P17-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  ind2.P16  <- ifelse(1:sealen < (ts2.P16-ts.start), 0, 1);
                  ind2.P17  <- ifelse(1:sealen < (ts2.P17-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logP16F1  <- par[18];
                  logP17F1  <- par[19];
                  logscale1 <- par[20];
                  logalpha1 <- par[21];
                  logbeta1  <- par[22];
                  logP1F2   <- par[23];
                  logP2F2   <- par[24];
                  logP3F2   <- par[25];
                  logP4F2   <- par[26];
                  logP5F2   <- par[27];
                  logP6F2   <- par[28];
                  logP7F2   <- par[29];
                  logP8F2   <- par[30];
                  logP9F2   <- par[31];
                  logP10F2  <- par[32];
                  logP11F2  <- par[33];
                  logP12F2  <- par[34];
                  logP13F2  <- par[35];
                  logP14F2  <- par[36];
                  logP15F2  <- par[37];
                  logP16F2  <- par[38];
                  logP17F2  <- par[39];
                  logscale2 <- par[40];
                  logalpha2 <- par[41];
                  logbeta2  <- par[42];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind1.P16[i]*exp(logP16F1)*exp(-exp(logM)*(i-(ts1.P16-ts.start)+1)) +
                                 ind1.P17[i]*exp(logP17F1)*exp(-exp(logM)*(i-(ts1.P17-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) +
                                 ind2.P16[i]*exp(logP16F2)*exp(-exp(logM)*(i-(ts2.P16-ts.start)+1)) +
                                 ind2.P17[i]*exp(logP17F2)*exp(-exp(logM)*(i-(ts2.P17-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN18P18P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(18,18))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts1.P16   <- dates[17];
                  ts1.P17   <- dates[18];
                  ts1.P18   <- dates[19];
                  ts2.P1    <- dates[20];
                  ts2.P2    <- dates[21];
                  ts2.P3    <- dates[22];
                  ts2.P4    <- dates[23];
                  ts2.P5    <- dates[24];
                  ts2.P6    <- dates[25];
                  ts2.P7    <- dates[26];
                  ts2.P8    <- dates[27];
                  ts2.P9    <- dates[28];
                  ts2.P10   <- dates[29];
                  ts2.P11   <- dates[30];
                  ts2.P12   <- dates[31];
                  ts2.P13   <- dates[32];
                  ts2.P14   <- dates[33];
                  ts2.P15   <- dates[34];
                  ts2.P16   <- dates[35];
                  ts2.P17   <- dates[36];
                  ts2.P18   <- dates[37];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  ind1.P16  <- ifelse(1:sealen < (ts1.P16-ts.start), 0, 1);
                  ind1.P17  <- ifelse(1:sealen < (ts1.P17-ts.start), 0, 1);
                  ind1.P18  <- ifelse(1:sealen < (ts1.P18-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  ind2.P16  <- ifelse(1:sealen < (ts2.P16-ts.start), 0, 1);
                  ind2.P17  <- ifelse(1:sealen < (ts2.P17-ts.start), 0, 1);
                  ind2.P18  <- ifelse(1:sealen < (ts2.P18-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logP16F1  <- par[18];
                  logP17F1  <- par[19];
                  logP18F1  <- par[20];
                  logscale1 <- par[21];
                  logalpha1 <- par[22];
                  logbeta1  <- par[23];
                  logP1F2   <- par[24];
                  logP2F2   <- par[25];
                  logP3F2   <- par[26];
                  logP4F2   <- par[27];
                  logP5F2   <- par[28];
                  logP6F2   <- par[29];
                  logP7F2   <- par[30];
                  logP8F2   <- par[31];
                  logP9F2   <- par[32];
                  logP10F2  <- par[33];
                  logP11F2  <- par[34];
                  logP12F2  <- par[35];
                  logP13F2  <- par[36];
                  logP14F2  <- par[37];
                  logP15F2  <- par[38];
                  logP16F2  <- par[39];
                  logP17F2  <- par[40];
                  logP18F2  <- par[41];
                  logscale2 <- par[42];
                  logalpha2 <- par[43];
                  logbeta2  <- par[44];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind1.P16[i]*exp(logP16F1)*exp(-exp(logM)*(i-(ts1.P16-ts.start)+1)) +
                                 ind1.P17[i]*exp(logP17F1)*exp(-exp(logM)*(i-(ts1.P17-ts.start)+1)) +
                                 ind1.P18[i]*exp(logP18F1)*exp(-exp(logM)*(i-(ts1.P18-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) +
                                 ind2.P16[i]*exp(logP16F2)*exp(-exp(logM)*(i-(ts2.P16-ts.start)+1)) +
                                 ind2.P17[i]*exp(logP17F2)*exp(-exp(logM)*(i-(ts2.P17-ts.start)+1)) +
                                 ind2.P18[i]*exp(logP18F2)*exp(-exp(logM)*(i-(ts2.P18-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN19P19P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(19,19))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts1.P16   <- dates[17];
                  ts1.P17   <- dates[18];
                  ts1.P18   <- dates[19];
                  ts1.P19   <- dates[20];
                  ts2.P1    <- dates[21];
                  ts2.P2    <- dates[22];
                  ts2.P3    <- dates[23];
                  ts2.P4    <- dates[24];
                  ts2.P5    <- dates[25];
                  ts2.P6    <- dates[26];
                  ts2.P7    <- dates[27];
                  ts2.P8    <- dates[28];
                  ts2.P9    <- dates[29];
                  ts2.P10   <- dates[30];
                  ts2.P11   <- dates[31];
                  ts2.P12   <- dates[32];
                  ts2.P13   <- dates[33];
                  ts2.P14   <- dates[34];
                  ts2.P15   <- dates[35];
                  ts2.P16   <- dates[36];
                  ts2.P17   <- dates[37];
                  ts2.P18   <- dates[38];
                  ts2.P19   <- dates[39];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  ind1.P16  <- ifelse(1:sealen < (ts1.P16-ts.start), 0, 1);
                  ind1.P17  <- ifelse(1:sealen < (ts1.P17-ts.start), 0, 1);
                  ind1.P18  <- ifelse(1:sealen < (ts1.P18-ts.start), 0, 1);
                  ind1.P19  <- ifelse(1:sealen < (ts1.P19-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  ind2.P16  <- ifelse(1:sealen < (ts2.P16-ts.start), 0, 1);
                  ind2.P17  <- ifelse(1:sealen < (ts2.P17-ts.start), 0, 1);
                  ind2.P18  <- ifelse(1:sealen < (ts2.P18-ts.start), 0, 1);
                  ind2.P19  <- ifelse(1:sealen < (ts2.P19-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logP16F1  <- par[18];
                  logP17F1  <- par[19];
                  logP18F1  <- par[20];
                  logP19F1  <- par[21];
                  logscale1 <- par[22];
                  logalpha1 <- par[23];
                  logbeta1  <- par[24];
                  logP1F2   <- par[25];
                  logP2F2   <- par[26];
                  logP3F2   <- par[27];
                  logP4F2   <- par[28];
                  logP5F2   <- par[29];
                  logP6F2   <- par[30];
                  logP7F2   <- par[31];
                  logP8F2   <- par[32];
                  logP9F2   <- par[33];
                  logP10F2  <- par[34];
                  logP11F2  <- par[35];
                  logP12F2  <- par[36];
                  logP13F2  <- par[37];
                  logP14F2  <- par[38];
                  logP15F2  <- par[39];
                  logP16F2  <- par[40];
                  logP17F2  <- par[41];
                  logP18F2  <- par[42];
                  logP19F2  <- par[43];
                  logscale2 <- par[44];
                  logalpha2 <- par[45];
                  logbeta2  <- par[46];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind1.P16[i]*exp(logP16F1)*exp(-exp(logM)*(i-(ts1.P16-ts.start)+1)) +
                                 ind1.P17[i]*exp(logP17F1)*exp(-exp(logM)*(i-(ts1.P17-ts.start)+1)) +
                                 ind1.P18[i]*exp(logP18F1)*exp(-exp(logM)*(i-(ts1.P18-ts.start)+1)) +
                                 ind1.P19[i]*exp(logP19F1)*exp(-exp(logM)*(i-(ts1.P19-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) +
                                 ind2.P16[i]*exp(logP16F2)*exp(-exp(logM)*(i-(ts2.P16-ts.start)+1)) +
                                 ind2.P17[i]*exp(logP17F2)*exp(-exp(logM)*(i-(ts2.P17-ts.start)+1)) +
                                 ind2.P18[i]*exp(logP18F2)*exp(-exp(logM)*(i-(ts2.P18-ts.start)+1)) +
                                 ind2.P19[i]*exp(logP19F2)*exp(-exp(logM)*(i-(ts2.P19-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
.CDMN20P20P <-
function(par,dates,obscat1,obseff1,obsmbm1,obscat2,obseff2,obsmbm2,distr,properties,output,pp=c(20,20))
  {
                  ts.start  <- head(dates,1);
                  ts1.P1    <- dates[2];
                  ts1.P2    <- dates[3];
                  ts1.P3    <- dates[4];
                  ts1.P4    <- dates[5];
                  ts1.P5    <- dates[6];
                  ts1.P6    <- dates[7];
                  ts1.P7    <- dates[8];
                  ts1.P8    <- dates[9];
                  ts1.P9    <- dates[10];
                  ts1.P10   <- dates[11];
                  ts1.P11   <- dates[12];
                  ts1.P12   <- dates[13];
                  ts1.P13   <- dates[14];
                  ts1.P14   <- dates[15];
                  ts1.P15   <- dates[16];
                  ts1.P16   <- dates[17];
                  ts1.P17   <- dates[18];
                  ts1.P18   <- dates[19];
                  ts1.P19   <- dates[20];
                  ts1.P20   <- dates[21];
                  ts2.P1    <- dates[22];
                  ts2.P2    <- dates[23];
                  ts2.P3    <- dates[24];
                  ts2.P4    <- dates[25];
                  ts2.P5    <- dates[26];
                  ts2.P6    <- dates[27];
                  ts2.P7    <- dates[28];
                  ts2.P8    <- dates[29];
                  ts2.P9    <- dates[30];
                  ts2.P10   <- dates[31];
                  ts2.P11   <- dates[32];
                  ts2.P12   <- dates[33];
                  ts2.P13   <- dates[34];
                  ts2.P14   <- dates[35];
                  ts2.P15   <- dates[36];
                  ts2.P16   <- dates[37];
                  ts2.P17   <- dates[38];
                  ts2.P18   <- dates[39];
                  ts2.P19   <- dates[40];
                  ts2.P20   <- dates[41];
                  ts.end    <- tail(dates,1);
                  sealen    <- ts.end-ts.start+1;
                  nstep     <- vector("numeric",sealen);
                  mccum     <- vector("numeric",sealen);
                  effeff1   <- vector("numeric",sealen);
                  effn1     <- vector("numeric",sealen);
                  predcat1  <- vector("numeric",sealen);
                  ind1.P1   <- ifelse(1:sealen < (ts1.P1-ts.start), 0, 1);
                  ind1.P2   <- ifelse(1:sealen < (ts1.P2-ts.start), 0, 1);
                  ind1.P3   <- ifelse(1:sealen < (ts1.P3-ts.start), 0, 1);
                  ind1.P4   <- ifelse(1:sealen < (ts1.P4-ts.start), 0, 1);
                  ind1.P5   <- ifelse(1:sealen < (ts1.P5-ts.start), 0, 1);
                  ind1.P6   <- ifelse(1:sealen < (ts1.P6-ts.start), 0, 1);
                  ind1.P7   <- ifelse(1:sealen < (ts1.P7-ts.start), 0, 1);
                  ind1.P8   <- ifelse(1:sealen < (ts1.P8-ts.start), 0, 1);
                  ind1.P9   <- ifelse(1:sealen < (ts1.P9-ts.start), 0, 1);
                  ind1.P10  <- ifelse(1:sealen < (ts1.P10-ts.start), 0, 1);
                  ind1.P11  <- ifelse(1:sealen < (ts1.P11-ts.start), 0, 1);
                  ind1.P12  <- ifelse(1:sealen < (ts1.P12-ts.start), 0, 1);
                  ind1.P13  <- ifelse(1:sealen < (ts1.P13-ts.start), 0, 1);
                  ind1.P14  <- ifelse(1:sealen < (ts1.P14-ts.start), 0, 1);
                  ind1.P15  <- ifelse(1:sealen < (ts1.P15-ts.start), 0, 1);
                  ind1.P16  <- ifelse(1:sealen < (ts1.P16-ts.start), 0, 1);
                  ind1.P17  <- ifelse(1:sealen < (ts1.P17-ts.start), 0, 1);
                  ind1.P18  <- ifelse(1:sealen < (ts1.P18-ts.start), 0, 1);
                  ind1.P19  <- ifelse(1:sealen < (ts1.P19-ts.start), 0, 1);
                  ind1.P20  <- ifelse(1:sealen < (ts1.P20-ts.start), 0, 1);
                  effeff2   <- vector("numeric",sealen);
                  effn2     <- vector("numeric",sealen);
                  predcat2  <- vector("numeric",sealen);
                  ind2.P1   <- ifelse(1:sealen < (ts2.P1-ts.start), 0, 1);
                  ind2.P2   <- ifelse(1:sealen < (ts2.P2-ts.start), 0, 1);
                  ind2.P3   <- ifelse(1:sealen < (ts2.P3-ts.start), 0, 1);
                  ind2.P4   <- ifelse(1:sealen < (ts2.P4-ts.start), 0, 1);
                  ind2.P5   <- ifelse(1:sealen < (ts2.P5-ts.start), 0, 1);
                  ind2.P6   <- ifelse(1:sealen < (ts2.P6-ts.start), 0, 1);
                  ind2.P7   <- ifelse(1:sealen < (ts2.P7-ts.start), 0, 1);
                  ind2.P8   <- ifelse(1:sealen < (ts2.P8-ts.start), 0, 1);
                  ind2.P9   <- ifelse(1:sealen < (ts2.P9-ts.start), 0, 1);
                  ind2.P10  <- ifelse(1:sealen < (ts2.P10-ts.start), 0, 1);
                  ind2.P11  <- ifelse(1:sealen < (ts2.P11-ts.start), 0, 1);
                  ind2.P12  <- ifelse(1:sealen < (ts2.P12-ts.start), 0, 1);
                  ind2.P13  <- ifelse(1:sealen < (ts2.P13-ts.start), 0, 1);
                  ind2.P14  <- ifelse(1:sealen < (ts2.P14-ts.start), 0, 1);
                  ind2.P15  <- ifelse(1:sealen < (ts2.P15-ts.start), 0, 1);
                  ind2.P16  <- ifelse(1:sealen < (ts2.P16-ts.start), 0, 1);
                  ind2.P17  <- ifelse(1:sealen < (ts2.P17-ts.start), 0, 1);
                  ind2.P18  <- ifelse(1:sealen < (ts2.P18-ts.start), 0, 1);
                  ind2.P19  <- ifelse(1:sealen < (ts2.P19-ts.start), 0, 1);
                  ind2.P20  <- ifelse(1:sealen < (ts2.P20-ts.start), 0, 1);
                  logM      <- par[1];
                  logN0     <- par[2];
                  logP1F1   <- par[3];
                  logP2F1   <- par[4];
                  logP3F1   <- par[5];
                  logP4F1   <- par[6];
                  logP5F1   <- par[7];
                  logP6F1   <- par[8];
                  logP7F1   <- par[9];
                  logP8F1   <- par[10];
                  logP9F1   <- par[11];
                  logP10F1  <- par[12];
                  logP11F1  <- par[13];
                  logP12F1  <- par[14];
                  logP13F1  <- par[15];
                  logP14F1  <- par[16];
                  logP15F1  <- par[17];
                  logP16F1  <- par[18];
                  logP17F1  <- par[19];
                  logP18F1  <- par[20];
                  logP19F1  <- par[21];
                  logP20F1  <- par[22];
                  logscale1 <- par[23];
                  logalpha1 <- par[24];
                  logbeta1  <- par[25];
                  logP1F2   <- par[26];
                  logP2F2   <- par[27];
                  logP3F2   <- par[28];
                  logP4F2   <- par[29];
                  logP5F2   <- par[30];
                  logP6F2   <- par[31];
                  logP7F2   <- par[32];
                  logP8F2   <- par[33];
                  logP9F2   <- par[34];
                  logP10F2  <- par[35];
                  logP11F2  <- par[36];
                  logP12F2  <- par[37];
                  logP13F2  <- par[38];
                  logP14F2  <- par[39];
                  logP15F2  <- par[40];
                  logP16F2  <- par[41];
                  logP17F2  <- par[42];
                  logP18F2  <- par[43];
                  logP19F2  <- par[44];
                  logP20F2  <- par[45];
                  logscale2 <- par[46];
                  logalpha2 <- par[47];
                  logbeta2  <- par[48];
                  mccum[1]  <- 0;
                  nstep[1]  <- exp(logN0)*exp(-exp(logM));
                  for(i in 2:sealen)
                     {
                     mccum[i] <- obscat1[i-1] + obscat2[i-1] + mccum[i-1]*exp(-exp(logM));
                     nstep[i] <- exp(logN0)*exp(-exp(logM)*i) +
                                 ind1.P1[i]*exp(logP1F1)*exp(-exp(logM)*(i-(ts1.P1-ts.start)+1)) +
                                 ind1.P2[i]*exp(logP2F1)*exp(-exp(logM)*(i-(ts1.P2-ts.start)+1)) +
                                 ind1.P3[i]*exp(logP3F1)*exp(-exp(logM)*(i-(ts1.P3-ts.start)+1)) +
                                 ind1.P4[i]*exp(logP4F1)*exp(-exp(logM)*(i-(ts1.P4-ts.start)+1)) +
                                 ind1.P5[i]*exp(logP5F1)*exp(-exp(logM)*(i-(ts1.P5-ts.start)+1)) +
                                 ind1.P6[i]*exp(logP6F1)*exp(-exp(logM)*(i-(ts1.P6-ts.start)+1)) +
                                 ind1.P7[i]*exp(logP7F1)*exp(-exp(logM)*(i-(ts1.P7-ts.start)+1)) +
                                 ind1.P8[i]*exp(logP8F1)*exp(-exp(logM)*(i-(ts1.P8-ts.start)+1)) +
                                 ind1.P9[i]*exp(logP9F1)*exp(-exp(logM)*(i-(ts1.P9-ts.start)+1)) +
                                 ind1.P10[i]*exp(logP10F1)*exp(-exp(logM)*(i-(ts1.P10-ts.start)+1)) +
                                 ind1.P11[i]*exp(logP11F1)*exp(-exp(logM)*(i-(ts1.P11-ts.start)+1)) +
                                 ind1.P12[i]*exp(logP12F1)*exp(-exp(logM)*(i-(ts1.P12-ts.start)+1)) +
                                 ind1.P13[i]*exp(logP13F1)*exp(-exp(logM)*(i-(ts1.P13-ts.start)+1)) +
                                 ind1.P14[i]*exp(logP14F1)*exp(-exp(logM)*(i-(ts1.P14-ts.start)+1)) +
                                 ind1.P15[i]*exp(logP15F1)*exp(-exp(logM)*(i-(ts1.P15-ts.start)+1)) +
                                 ind1.P16[i]*exp(logP16F1)*exp(-exp(logM)*(i-(ts1.P16-ts.start)+1)) +
                                 ind1.P17[i]*exp(logP17F1)*exp(-exp(logM)*(i-(ts1.P17-ts.start)+1)) +
                                 ind1.P18[i]*exp(logP18F1)*exp(-exp(logM)*(i-(ts1.P18-ts.start)+1)) +
                                 ind1.P19[i]*exp(logP19F1)*exp(-exp(logM)*(i-(ts1.P19-ts.start)+1)) +
                                 ind1.P20[i]*exp(logP20F1)*exp(-exp(logM)*(i-(ts1.P20-ts.start)+1)) +
                                 ind2.P1[i]*exp(logP1F2)*exp(-exp(logM)*(i-(ts2.P1-ts.start)+1)) +
                                 ind2.P2[i]*exp(logP2F2)*exp(-exp(logM)*(i-(ts2.P2-ts.start)+1)) +
                                 ind2.P3[i]*exp(logP3F2)*exp(-exp(logM)*(i-(ts2.P3-ts.start)+1)) +
                                 ind2.P4[i]*exp(logP4F2)*exp(-exp(logM)*(i-(ts2.P4-ts.start)+1)) +
                                 ind2.P5[i]*exp(logP5F2)*exp(-exp(logM)*(i-(ts2.P5-ts.start)+1)) +
                                 ind2.P6[i]*exp(logP6F2)*exp(-exp(logM)*(i-(ts2.P6-ts.start)+1)) +
                                 ind2.P7[i]*exp(logP7F2)*exp(-exp(logM)*(i-(ts2.P7-ts.start)+1)) +
                                 ind2.P8[i]*exp(logP8F2)*exp(-exp(logM)*(i-(ts2.P8-ts.start)+1)) +
                                 ind2.P9[i]*exp(logP9F2)*exp(-exp(logM)*(i-(ts2.P9-ts.start)+1)) +
                                 ind2.P10[i]*exp(logP10F2)*exp(-exp(logM)*(i-(ts2.P10-ts.start)+1)) +
                                 ind2.P11[i]*exp(logP11F2)*exp(-exp(logM)*(i-(ts2.P11-ts.start)+1)) +
                                 ind2.P12[i]*exp(logP12F2)*exp(-exp(logM)*(i-(ts2.P12-ts.start)+1)) +
                                 ind2.P13[i]*exp(logP13F2)*exp(-exp(logM)*(i-(ts2.P13-ts.start)+1)) +
                                 ind2.P14[i]*exp(logP14F2)*exp(-exp(logM)*(i-(ts2.P14-ts.start)+1)) +
                                 ind2.P15[i]*exp(logP15F2)*exp(-exp(logM)*(i-(ts2.P15-ts.start)+1)) +
                                 ind2.P16[i]*exp(logP16F2)*exp(-exp(logM)*(i-(ts2.P16-ts.start)+1)) +
                                 ind2.P17[i]*exp(logP17F2)*exp(-exp(logM)*(i-(ts2.P17-ts.start)+1)) +
                                 ind2.P18[i]*exp(logP18F2)*exp(-exp(logM)*(i-(ts2.P18-ts.start)+1)) +
                                 ind2.P19[i]*exp(logP19F2)*exp(-exp(logM)*(i-(ts2.P19-ts.start)+1)) +
                                 ind2.P20[i]*exp(logP20F2)*exp(-exp(logM)*(i-(ts2.P20-ts.start)+1)) -
                                 mccum[i]*exp(-exp(logM)/2);
                     }
                  effeff1   <- obseff1^(exp(logalpha1));
                  effn1     <- nstep^(exp(logbeta1));
                  predcat1  <- exp(logscale1)*(effeff1*effn1)*exp(-exp(logM)/2);
                  effeff2   <- obseff2^(exp(logalpha2));
                  effn2     <- nstep^(exp(logbeta2));
                  predcat2  <- exp(logscale2)*(effeff2*effn2)*exp(-exp(logM)/2);
                  Likel     <- .CatDynLik2F(obscat1,predcat1,obscat2,predcat2,distr,par)
                  if(output == "predict")
                    {
                    catdynexp <- .CatDynExp2F.Res(properties,nstep,obsmbm1,obsmbm2,pp,dates,distr,par,Likel,ts.start,ts.end,obseff1,obscat1,predcat1,obseff2,obscat2,predcat2)
                    class(catdynexp) <- "CatDynExp";
                    return(catdynexp);
                    }
                  else
                    {
                     if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] != "apnormal" & distr[2] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,1])) - sum(Likel[["Likelihood"]][,2]);
                       }

                     else if(distr[2] == "apnormal" | distr[2] == "aplnormal" & distr[1] != "apnormal" & distr[1] != "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*log(sum(Likel[["Likelihood"]][,2])) - sum(Likel[["Likelihood"]][,1]);
                       }
                     else if(distr[1] == "apnormal" | distr[1] == "aplnormal" & distr[2] == "apnormal" | distr[2] == "aplnormal")
                       {
                        negsup <- ((sealen-2)/2)*(log(sum(Likel[["Likelihood"]][,1])) + log(sum(Likel[["Likelihood"]][,2])));
                       }
                     else
                       {
                        negsup <- - sum(Likel[["Likelihood"]][,1]) - sum(Likel[["Likelihood"]][,2])
                       }
                     return(negsup);
                    }
 }
