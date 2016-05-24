GroAgeFit <-
function(x, m, unsex.action=NULL, par, distr, method, control=list(), hessian=TRUE, itnmax)
  {
             distr.set  <- c("normal","apnormal","lognormal","aplnormal","gamma")
             model.set  <- c("vonber","schnute1","schnute2","schnute3")
             unsex.set  <- c("split","males","females","ignore")
             if(x$Properties$Sex == "Total" & is.null(unsex.action))
               {
                stop("When the data includes males, females, and unsexed individuals \n
                      and option for what to do with the unsexed individuals must be provided")
               }
             if(!is.null(unsex.action) & !unsex.action%in%unsex.set)
               {
                stop("'unsex.action' must be 'split', 'males', 'females' or 'ignore'")
               }
             if(class(par) != "numeric")
               {
                stop("'par' must be a numeric vector")
               }
             if(!distr%in%distr.set)
               {
                stop("'distr' must be 'normal','apnormal','lognormal','aplnormal', or 'gamma'")
               }
             if(!m%in%model.set)
               {
                stop("'m' must be 'vonber','schnute1','schnute2', or 'schnute3'")
               }
             if(x$Properties$Sex == "Pooled" | x$Properties$Sex == "Males" | x$Properties$Sex == "Females")
               {
                if(m == "vonber")
                   {
                    if(distr == "apnormal" | distr == "aplnormal")
                      {
                       if(length(par) != 3)
                         {
                          stop("For a Von Bertalanffy model of type 1 with adjusted profile likelihood approximation to the normal \n
                          or lognormal, par must be vector of length 3")
                         }
                      }
                    if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                      {
                       if(length(par) != 4)
                         {
                          stop("For a Von Bertalanffy model of type 1 with exact normal, lognormal or gamma likelihood, \n
                          par must be vector of length 4")
                         }
                      }
                    results1        <- do.call(optimx, list(par=par,
                                                            fn=.vonBer,
                                                            gr=NULL,
                                                            age=x$Data[which(!is.na(x$Data[,2]) & !is.na(x$Data[,3])),2],
                                                            obslen=x$Data[which(!is.na(x$Data[,2]) & !is.na(x$Data[,3])),3],
                                                            distr=distr,
                                                            method=method,
                                                            lower = -Inf,
                                                            upper = Inf,
                                                            control=list(),
                                                            hessian=hessian,
                                                            itnmax=itnmax,
                                                            output="estimate"));
                    results2        <- vector("list",length(method));
                    names(results2) <- method;
                    temp            <- attr(results1, "details")
                    for(i in 1:length(method))
                      {
                       results2[[i]]$Type         <- m;
                       results2[[i]]$Distr        <- distr;
                       par.names                  <- c("Linf","L0","K");
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          par.names <- c(par.names,"Dispersion")
                         }
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev         <- NA;
                       results2[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                       if(length(temp[i,]$ngatend) == length(par.names) & !any(is.na(temp[i,]$nhatend)) & 1/kappa(temp[i,]$nhatend) > 1e-15)
                         {
                          results2[[i]]$converg          <- results1[i,length(par)+5];
                          results2[[i]]$kkt              <- results1[i,(length(par)+6):(length(par)+7)];
                          results2[[i]]$AIC              <- 2*length(par)-2*results1[i,length(par)+1];
                          results2[[i]]$bt.par           <- exp(results1[i,1:length(par)])
                          results2[[i]]$num.grads        <- temp[i,]$ngatend;
                                                v        <- matrix(0,length(par.names),length(par.names));
                                                if(distr == "apnormal" | distr == "aplnormal")
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                           mean=as.numeric(results1[i,1:length(par)]),
                                                                           cov=try(solve(temp[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                                                else
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                           mean=as.numeric(results1[i,1:length(par)]),
                                                                           cov=try(solve(temp[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                          results2[[i]]$bt.stdev         <- sqrt(diag(v));
                          results2[[i]]$Cor              <- cor(v);
                          names(results2[[i]]$num.grads) <- par.names;
                          names(results2[[i]]$bt.par)    <- par.names;
                          names(results2[[i]]$bt.stdev)  <- par.names;
                          colnames(results2[[i]]$Cor)    <- par.names;
                          rownames(results2[[i]]$Cor)    <- par.names;
                         }
                      }
                    rm(v);
                    rm(temp);
                   }
                else if(m == "schnute1")
                   {
                    if(distr == "apnormal" | distr == "aplnormal")
                      {
                       if(length(par) != 5)
                         {
                          stop("For a Schnute model of type 1 with adjusted profile likelihood approximation to the normal or lognormal, par must be vector of length 5")
                         }
                      }
                    if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                      {
                       if(length(par) != 6)
                         {
                          stop("For a Schnute model of type 1 with exact normal, lognormal or gamma likelihood, par must be vector of length 6")
                         }
                      }
                    results1        <- do.call(optimx, list(par=par,
                                                            fn=.schnute1,
                                                            gr=NULL,
                                                            age=x$Data[which(!is.na(x$Data[,2]) & !is.na(x$Data[,3])),2],
                                                            obslen=x$Data[which(!is.na(x$Data[,2]) & !is.na(x$Data[,3])),3],
                                                            distr=distr,
                                                            method=method,
                                                            lower = -Inf,
                                                            upper = Inf,
                                                            control=list(),
                                                            hessian=hessian,
                                                            itnmax=itnmax,
                                                            output="estimate"));
                    results2        <- vector("list",length(method));
                    names(results2) <- method;
                    temp            <- attr(results1, "details")
                    for(i in 1:length(method))
                      {
                       results2[[i]]$Type         <- m;
                       results2[[i]]$Distr        <- distr;
                       par.names                  <- c("L1","a1","mu","g1","g2");
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          par.names <- c(par.names,"Dispersion")
                         }
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev         <- NA;
                       results2[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                       if(length(temp[i,]$ngatend) == length(par.names) & !any(is.na(temp[i,]$nhatend)) & 1/kappa(temp[i,]$nhatend) > 1e-15)
                         {
                          results2[[i]]$converg          <- results1[i,length(par)+5];
                          results2[[i]]$kkt              <- results1[i,(length(par)+6):(length(par)+7)];
                          results2[[i]]$AIC              <- 2*length(par)-2*results1[i,length(par)+1];
                          results2[[i]]$bt.par           <- exp(results1[i,1:length(par)])
                          results2[[i]]$num.grads        <- temp[i,]$ngatend;
                                                v        <- matrix(0,length(par.names),length(par.names));
                                                if(distr == "apnormal" | distr == "aplnormal")
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                           mean=as.numeric(results1[i,1:length(par)]),
                                                                           cov=try(solve(temp[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                                                else
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                           mean=as.numeric(results1[i,1:length(par)]),
                                                                           cov=try(solve(temp[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                          results2[[i]]$bt.stdev         <- sqrt(diag(v));
                          results2[[i]]$Cor              <- cor(v);
                          names(results2[[i]]$num.grads) <- par.names;
                          names(results2[[i]]$bt.par)    <- par.names;
                          names(results2[[i]]$bt.stdev)  <- par.names;
                          colnames(results2[[i]]$Cor)    <- par.names;
                          rownames(results2[[i]]$Cor)    <- par.names;
                         }
                      }
                    rm(v);
                    rm(temp);
                   }
                else if(m == "schnute2")
                   {
                    if(distr == "apnormal" | distr == "aplnormal")
                      {
                       if(length(par) != 3)
                         {
                          stop("For a Schnute model of type 2 with adjusted profile likelihood approximation to the normal or lognormal, par must be vector of length 3")
                         }
                      }
                    if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                      {
                       if(length(par) != 4)
                         {
                          stop("For a Schnute model of type 2 with exact normal, lognormal or gamma likelihood, par must be vector of length 4")
                         }
                      }
                    results1        <- do.call(optimx, list(par=par,
                                                            fn=.schnute2,
                                                            gr=NULL,
                                                            age=x$Data[which(!is.na(x$Data[,2]) & !is.na(x$Data[,3])),2],
                                                            obslen=x$Data[which(!is.na(x$Data[,2]) & !is.na(x$Data[,3])),3],
                                                            distr=distr,
                                                            method=method,
                                                            lower = -Inf,
                                                            upper = Inf,
                                                            control=list(),
                                                            hessian=hessian,
                                                            itnmax=itnmax,
                                                            output="estimate"));
                    results2        <- vector("list",length(method));
                    names(results2) <- method;
                    temp            <- attr(results1, "details")
                    for(i in 1:length(method))
                      {
                       results2[[i]]$Type         <- m;
                       results2[[i]]$Distr        <- distr;
                       par.names                  <- c("Ling","gamma","g1");
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          par.names <- c(par.names,"Dispersion")
                         }
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev         <- NA;
                       results2[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                       if(length(temp[i,]$ngatend) == length(par.names) & !any(is.na(temp[i,]$nhatend)) & 1/kappa(temp[i,]$nhatend) > 1e-15)
                         {
                          results2[[i]]$converg          <- results1[i,length(par)+5];
                          results2[[i]]$kkt              <- results1[i,(length(par)+6):(length(par)+7)];
                          results2[[i]]$AIC              <- 2*length(par)-2*results1[i,length(par)+1];
                          results2[[i]]$bt.par           <- exp(results1[i,1:length(par)])
                          results2[[i]]$num.grads        <- temp[i,]$ngatend;
                                                v        <- matrix(0,length(par.names),length(par.names));
                                                if(distr == "apnormal" | distr == "aplnormal")
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                           mean=as.numeric(results1[i,1:length(par)]),
                                                                           cov=try(solve(temp[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                                                else
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                           mean=as.numeric(results1[i,1:length(par)]),
                                                                           cov=try(solve(temp[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                          results2[[i]]$bt.stdev         <- sqrt(diag(v));
                          results2[[i]]$Cor              <- cor(v);
                          names(results2[[i]]$num.grads) <- par.names;
                          names(results2[[i]]$bt.par)    <- par.names;
                          names(results2[[i]]$bt.stdev)  <- par.names;
                          colnames(results2[[i]]$Cor)    <- par.names;
                          rownames(results2[[i]]$Cor)    <- par.names;
                         }
                      }
                    rm(v);
                    rm(temp);
                   }
                else if(m == "schnute3")
                   {
                    if(distr == "apnormal" | distr == "aplnormal")
                      {
                       if(length(par) != 2)
                         {
                          stop("For a Schnute model of type 2 with adjusted profile likelihood approximation to the normal or lognormal, par must be vector of length 2")
                         }
                      }
                    if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                      {
                       if(length(par) != 3)
                         {
                          stop("For a Schnute model of type 2 with exact normal, lognormal or gamma likelihood, par must be vector of length 3")
                         }
                      }
                    results1        <- do.call(optimx, list(par=par,
                                                            fn=.schnute2,
                                                            gr=NULL,
                                                            age=x$Data[which(!is.na(x$Data[,2]) & !is.na(x$Data[,3])),2],
                                                            obslen=x$Data[which(!is.na(x$Data[,2]) & !is.na(x$Data[,3])),3],
                                                            distr=distr,
                                                            method=method,
                                                            lower = -Inf,
                                                            upper = Inf,
                                                            control=list(),
                                                            hessian=hessian,
                                                            itnmax=itnmax,
                                                            output="estimate"));
                    results2        <- vector("list",length(method));
                    names(results2) <- method;
                    temp            <- attr(results1, "details")
                    for(i in 1:length(method))
                      {
                       results2[[i]]$Type         <- m;
                       results2[[i]]$Distr        <- distr;
                       par.names                  <- c("g1","g2");
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          par.names <- c(par.names,"Dispersion")
                         }
                       results2[[i]]$converg          <- "FAIL";
                       results2[[i]]$kkt              <- NA;
                       results2[[i]]$AIC              <- NA;
                       results2[[i]]$bt.par           <- NA;
                       results2[[i]]$num.grads        <- NA;
                       results2[[i]]$bt.stdev         <- NA;
                       results2[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                       if(length(temp[i,]$ngatend) == length(par.names) & !any(is.na(temp[i,]$nhatend)) & 1/kappa(temp[i,]$nhatend) > 1e-15)
                         {
                          results2[[i]]$converg          <- results1[i,length(par)+5];
                          results2[[i]]$kkt              <- results1[i,(length(par)+6):(length(par)+7)];
                          results2[[i]]$AIC              <- 2*length(par)-2*results1[i,length(par)+1];
                          results2[[i]]$bt.par           <- exp(results1[i,1:length(par)])
                          results2[[i]]$num.grads        <- temp[i,]$ngatend;
                                                v        <- matrix(0,length(par.names),length(par.names));
                                                if(distr == "apnormal" | distr == "aplnormal")
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2)),
                                                                           mean=as.numeric(results1[i,1:length(par)]),
                                                                           cov=try(solve(temp[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                                                else
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                           mean=as.numeric(results1[i,1:length(par)]),
                                                                           cov=try(solve(temp[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                          results2[[i]]$bt.stdev         <- sqrt(diag(v));
                          results2[[i]]$Cor              <- cor(v);
                          names(results2[[i]]$num.grads) <- par.names;
                          names(results2[[i]]$bt.par)    <- par.names;
                          names(results2[[i]]$bt.stdev)  <- par.names;
                          colnames(results2[[i]]$Cor)    <- par.names;
                          rownames(results2[[i]]$Cor)    <- par.names;
                         }
                      }
                    rm(v);
                    rm(temp);
                  }
               }
             else if(x$Properties$Sex == "Both")
               {
                if(m == "vonber")
                   {
                    if(distr == "apnormal" | distr == "aplnormal")
                      {
                       if(length(par) != 6)
                         {
                          stop("For a two-sexes Von Bertalanffy model with adjusted profile likelihood approximation to the normal \n
                                or lognormal, par must be vector of length 6")
                         }
                       par.fem     <- par[1:3];
                       par.mal     <- par[4:6];
                       par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]));
                      }
                    if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                      {
                       if(length(par) != 8)
                         {
                          stop("For a two-sexes Von Bertalanffy model with exact normal, lognormal or gamma likelihood, \n
                                par must be vector of length 8")
                         }
                       par.fem     <- par[c(1:3,7)];
                       par.mal     <- par[c(4:6,8)];
                       par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]), 0.5*(par[7]+par[8]));
                      }
                    age.fem     <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2];
                    obslen.fem  <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3];
                    results.fem <- do.call(optimx, list(par=par.fem,
                                                        fn=.vonBer,
                                                        gr=NULL,
                                                        age=age.fem,
                                                        obslen=obslen.fem,
                                                        distr=distr,
                                                        method=method,
                                                        lower = -Inf,
                                                        upper = Inf,
                                                        control=list(),
                                                        hessian=hessian,
                                                        itnmax=itnmax,
                                                        output="estimate"));
                    age.mal     <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2];
                    obslen.mal  <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3];
                    results.mal <- do.call(optimx, list(par=par.mal,
                                                        fn=.vonBer,
                                                        gr=NULL,
                                                        age=age.mal,
                                                        obslen=obslen.mal,
                                                        distr=distr,
                                                        method=method,
                                                        lower = -Inf,
                                                        upper = Inf,
                                                        control=list(),
                                                        hessian=hessian,
                                                        itnmax=itnmax,
                                                        output="estimate"));
                    age.pld     <- c(age.fem,age.mal);
                    obslen.pld  <- c(obslen.fem,obslen.mal);
                    results.pld <- do.call(optimx, list(par=par.pld,
                                                        fn=.vonBer,
                                                        gr=NULL,
                                                        age=age.pld,
                                                        obslen=obslen.pld,
                                                        distr=distr,
                                                        method=method,
                                                        lower = -Inf,
                                                        upper = Inf,
                                                        control=list(),
                                                        hessian=hessian,
                                                        itnmax=itnmax,
                                                        output="estimate"));
                    rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                    results2                <- vector("list",3);
                    names(results2)         <- c("Females","Males","Pooled");
                    #Females
                    results2$Females        <- vector("list",length(method));
                    names(results2$Females) <- method;
                    temp.fem                <- attr(results.fem, "details")
                    for(i in 1:length(method))
                      {
                       results2$Females[[i]]$Type         <- m;
                       results2$Females[[i]]$Distr        <- distr;
                       par.names                  <- c("Linf","L0","K");
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          par.names <- c(par.names,"Dispersion")
                         }
                       results2$Females[[i]]$converg          <- "FAIL";
                       results2$Females[[i]]$kkt              <- NA;
                       results2$Females[[i]]$AIC              <- NA;
                       results2$Females[[i]]$bt.par           <- NA;
                       results2$Females[[i]]$num.grads        <- NA;
                       results2$Females[[i]]$bt.stdev         <- NA;
                       results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                       if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                         {
                          results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                          results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                          results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                          results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                          results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                v        <- matrix(0,length(par.names),length(par.names));
                                                if(distr == "apnormal" | distr == "aplnormal")
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                           mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                           cov=try(solve(temp.fem[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                                                else
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                           mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                           cov=try(solve(temp.fem[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                          results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                          results2$Females[[i]]$Cor              <- cor(v);
                          names(results2$Females[[i]]$num.grads) <- par.names;
                          names(results2$Females[[i]]$bt.par)    <- par.names;
                          names(results2$Females[[i]]$bt.stdev)  <- par.names;
                          colnames(results2$Females[[i]]$Cor)    <- par.names;
                          rownames(results2$Females[[i]]$Cor)    <- par.names;
                         }
                      }
                    rm(v);
                    rm(temp.fem,par.fem);
                    #Males
                    results2$Males        <- vector("list",length(method));
                    names(results2$Males)   <- method;
                    temp.mal                <- attr(results.mal, "details")
                    for(i in 1:length(method))
                      {
                       results2$Males[[i]]$Type         <- m;
                       results2$Males[[i]]$Distr        <- distr;
                       par.names                  <- c("Linf","L0","K");
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          par.names <- c(par.names,"Dispersion")
                         }
                       results2$Males[[i]]$converg          <- "FAIL";
                       results2$Males[[i]]$kkt              <- NA;
                       results2$Males[[i]]$AIC              <- NA;
                       results2$Males[[i]]$bt.par           <- NA;
                       results2$Males[[i]]$num.grads        <- NA;
                       results2$Males[[i]]$bt.stdev         <- NA;
                       results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                       if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                         {
                          results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                          results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                          results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                          results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                          results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                v        <- matrix(0,length(par.names),length(par.names));
                                                if(distr == "apnormal" | distr == "aplnormal")
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                           mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                           cov=try(solve(temp.mal[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                                                else
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                           mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                           cov=try(solve(temp.mal[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                          results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                          results2$Males[[i]]$Cor              <- cor(v);
                          names(results2$Males[[i]]$num.grads) <- par.names;
                          names(results2$Males[[i]]$bt.par)    <- par.names;
                          names(results2$Males[[i]]$bt.stdev)  <- par.names;
                          colnames(results2$Males[[i]]$Cor)    <- par.names;
                          rownames(results2$Males[[i]]$Cor)    <- par.names;
                         }
                      }
                    rm(v);
                    rm(temp.mal,par.mal);
                    #Pooled
                    results2$Pooled        <- vector("list",length(method));
                    names(results2$Pooled)   <- method;
                    temp.pld                <- attr(results.pld, "details")
                    for(i in 1:length(method))
                      {
                       results2$Pooled[[i]]$Type         <- m;
                       results2$Pooled[[i]]$Distr        <- distr;
                       par.names                  <- c("Linf","L0","K");
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          par.names <- c(par.names,"Dispersion")
                         }
                       results2$Pooled[[i]]$converg          <- "FAIL";
                       results2$Pooled[[i]]$kkt              <- NA;
                       results2$Pooled[[i]]$AIC              <- NA;
                       results2$Pooled[[i]]$bt.par           <- NA;
                       results2$Pooled[[i]]$num.grads        <- NA;
                       results2$Pooled[[i]]$bt.stdev         <- NA;
                       results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                       if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                         {
                          results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                          results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                          results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                          results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                          results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                v        <- matrix(0,length(par.names),length(par.names));
                                                if(distr == "apnormal" | distr == "aplnormal")
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                           mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                           cov=try(solve(temp.pld[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                                                else
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                           mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                           cov=try(solve(temp.pld[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                          results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                          results2$Pooled[[i]]$Cor              <- cor(v);
                          names(results2$Pooled[[i]]$num.grads) <- par.names;
                          names(results2$Pooled[[i]]$bt.par)    <- par.names;
                          names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                          colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                          rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                         }
                      }
                    rm(v);
                    rm(temp.pld,par.pld);
                   }
                else if(m == "schnute1")
                   {
                    if(distr == "apnormal" | distr == "aplnormal")
                      {
                       if(length(par) != 10)
                         {
                          stop("For a two-sexes Schnute model of type 1 with adjusted profile likelihood approximation to the normal \n
                          or lognormal, par must be vector of length 10")
                         }
                       par.fem     <- par[1:5];
                       par.mal     <- par[6:10];
                       par.pld     <- c(0.5*(par[1]+par[6]), 0.5*(par[2]+par[7]), 0.5*(par[3]+par[8]), 0.5*(par[4]+par[9]), 0.5*(par[5]+par[10]));
                      }
                    if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                      {
                       if(length(par) != 12)
                         {
                          stop("For a two-sexes Schnute model of type 1 with exact normal, lognormal or gamma likelihood, \n
                          par must be vector of length 12")
                         }
                       par.fem     <- par[c(1:5,11)];
                       par.mal     <- par[c(6:10,12)];
                       par.pld     <- c(0.5*(par[1]+par[6]), 0.5*(par[2]+par[7]), 0.5*(par[3]+par[8]), 0.5*(par[4]+par[9]), 0.5*(par[5]+par[10]), 0.5*(par[11]+par[12]));
                      }
                    age.fem     <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2];
                    obslen.fem  <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3];
                    results.fem <- do.call(optimx, list(par=par.fem,
                                                        fn=.schnute1,
                                                        gr=NULL,
                                                        age=age.fem,
                                                        obslen=obslen.fem,
                                                        distr=distr,
                                                        method=method,
                                                        lower = -Inf,
                                                        upper = Inf,
                                                        control=list(),
                                                        hessian=hessian,
                                                        itnmax=itnmax,
                                                        output="estimate"));
                    age.mal     <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2];
                    obslen.mal  <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3];
                    results.mal <- do.call(optimx, list(par=par.mal,
                                                        fn=.schnute1,
                                                        gr=NULL,
                                                        age=age.mal,
                                                        obslen=obslen.mal,
                                                        distr=distr,
                                                        method=method,
                                                        lower = -Inf,
                                                        upper = Inf,
                                                        control=list(),
                                                        hessian=hessian,
                                                        itnmax=itnmax,
                                                        output="estimate"));
                    age.pld     <- c(age.fem,age.mal);
                    obslen.pld  <- c(obslen.fem,obslen.mal);
                    results.pld <- do.call(optimx, list(par=par.pld,
                                                        fn=.schnute1,
                                                        gr=NULL,
                                                        age=age.pld,
                                                        obslen=obslen.pld,
                                                        distr=distr,
                                                        method=method,
                                                        lower = -Inf,
                                                        upper = Inf,
                                                        control=list(),
                                                        hessian=hessian,
                                                        itnmax=itnmax,
                                                        output="estimate"));
                    rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                    results2                <- vector("list",3);
                    names(results2)         <- c("Females","Males","Pooled");
                    #Females
                    results2$Females        <- vector("list",length(method));
                    names(results2$Females) <- method;
                    temp.fem                <- attr(results.fem, "details")
                    for(i in 1:length(method))
                      {
                       results2$Females[[i]]$Type         <- m;
                       results2$Females[[i]]$Distr        <- distr;
                       par.names                  <- c("l1","a1","mu","g1","g2");
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          par.names <- c(par.names,"Dispersion")
                         }
                       results2$Females[[i]]$converg          <- "FAIL";
                       results2$Females[[i]]$kkt              <- NA;
                       results2$Females[[i]]$AIC              <- NA;
                       results2$Females[[i]]$bt.par           <- NA;
                       results2$Females[[i]]$num.grads        <- NA;
                       results2$Females[[i]]$bt.stdev         <- NA;
                       results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                       if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                         {
                          results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                          results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                          results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                          results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                          results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                v        <- matrix(0,length(par.names),length(par.names));
                                                if(distr == "apnormal" | distr == "aplnormal")
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                           mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                           cov=try(solve(temp.fem[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                                                else
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                           mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                           cov=try(solve(temp.fem[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                          results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                          results2$Females[[i]]$Cor              <- cor(v);
                          names(results2$Females[[i]]$num.grads) <- par.names;
                          names(results2$Females[[i]]$bt.par)    <- par.names;
                          names(results2$Females[[i]]$bt.stdev)  <- par.names;
                          colnames(results2$Females[[i]]$Cor)    <- par.names;
                          rownames(results2$Females[[i]]$Cor)    <- par.names;
                         }
                      }
                    rm(v);
                    rm(temp.fem,par.fem);
                    #Males
                    results2$Males        <- vector("list",length(method));
                    names(results2$Males)   <- method;
                    temp.mal                <- attr(results.mal, "details")
                    for(i in 1:length(method))
                      {
                       results2$Males[[i]]$Type         <- m;
                       results2$Males[[i]]$Distr        <- distr;
                       par.names                  <- c("l1","a1","mu","g1","g2");
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          par.names <- c(par.names,"Dispersion")
                         }
                       results2$Males[[i]]$converg          <- "FAIL";
                       results2$Males[[i]]$kkt              <- NA;
                       results2$Males[[i]]$AIC              <- NA;
                       results2$Males[[i]]$bt.par           <- NA;
                       results2$Males[[i]]$num.grads        <- NA;
                       results2$Males[[i]]$bt.stdev         <- NA;
                       results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                       if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                         {
                          results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                          results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                          results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                          results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                          results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                v        <- matrix(0,length(par.names),length(par.names));
                                                if(distr == "apnormal" | distr == "aplnormal")
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                           mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                           cov=try(solve(temp.mal[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                                                else
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                           mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                           cov=try(solve(temp.mal[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                          results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                          results2$Males[[i]]$Cor              <- cor(v);
                          names(results2$Males[[i]]$num.grads) <- par.names;
                          names(results2$Males[[i]]$bt.par)    <- par.names;
                          names(results2$Males[[i]]$bt.stdev)  <- par.names;
                          colnames(results2$Males[[i]]$Cor)    <- par.names;
                          rownames(results2$Males[[i]]$Cor)    <- par.names;
                         }
                      }
                    rm(v);
                    rm(temp.mal,par.mal);
                    #Pooled
                    results2$Pooled        <- vector("list",length(method));
                    names(results2$Pooled)   <- method;
                    temp.pld                <- attr(results.pld, "details")
                    for(i in 1:length(method))
                      {
                       results2$Pooled[[i]]$Type         <- m;
                       results2$Pooled[[i]]$Distr        <- distr;
                       par.names                  <- c("l1","a1","mu","g1","g2");
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          par.names <- c(par.names,"Dispersion")
                         }
                       results2$Pooled[[i]]$converg          <- "FAIL";
                       results2$Pooled[[i]]$kkt              <- NA;
                       results2$Pooled[[i]]$AIC              <- NA;
                       results2$Pooled[[i]]$bt.par           <- NA;
                       results2$Pooled[[i]]$num.grads        <- NA;
                       results2$Pooled[[i]]$bt.stdev         <- NA;
                       results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                       if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                         {
                          results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                          results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                          results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                          results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                          results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                v        <- matrix(0,length(par.names),length(par.names));
                                                if(distr == "apnormal" | distr == "aplnormal")
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                           mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                           cov=try(solve(temp.pld[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                                                else
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                           mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                           cov=try(solve(temp.pld[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                          results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                          results2$Pooled[[i]]$Cor              <- cor(v);
                          names(results2$Pooled[[i]]$num.grads) <- par.names;
                          names(results2$Pooled[[i]]$bt.par)    <- par.names;
                          names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                          colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                          rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                         }
                      }
                    rm(v);
                    rm(temp.pld,par.pld);
                   }
                else if(m == "schnute2")
                   {
                    if(distr == "apnormal" | distr == "aplnormal")
                      {
                       if(length(par) != 6)
                         {
                          stop("For a two-sexes Schnute model of type 2 with adjusted profile likelihood approximation to the normal \n
                                or lognormal, par must be vector of length 6")
                         }
                       par.fem     <- par[1:3];
                       par.mal     <- par[4:6];
                       par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]));
                      }
                    if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                      {
                       if(length(par) != 8)
                         {
                          stop("For a two-sexes Schnute model of type 2 with exact normal, lognormal or gamma likelihood, \n
                                par must be vector of length 8")
                         }
                       par.fem     <- par[c(1:3,7)];
                       par.mal     <- par[c(4:6,8)];
                       par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]), 0.5*(par[7]+par[7]));
                      }
                    age.fem     <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2];
                    obslen.fem  <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3];
                    results.fem <- do.call(optimx, list(par=par.fem,
                                                        fn=.schnute2,
                                                        gr=NULL,
                                                        age=age.fem,
                                                        obslen=obslen.fem,
                                                        distr=distr,
                                                        method=method,
                                                        lower = -Inf,
                                                        upper = Inf,
                                                        control=list(),
                                                        hessian=hessian,
                                                        itnmax=itnmax,
                                                        output="estimate"));
                    age.mal     <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2];
                    obslen.mal  <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3];
                    results.mal <- do.call(optimx, list(par=par.mal,
                                                        fn=.schnute2,
                                                        gr=NULL,
                                                        age=age.mal,
                                                        obslen=obslen.mal,
                                                        distr=distr,
                                                        method=method,
                                                        lower = -Inf,
                                                        upper = Inf,
                                                        control=list(),
                                                        hessian=hessian,
                                                        itnmax=itnmax,
                                                        output="estimate"));
                    age.pld     <- c(age.fem,age.mal);
                    obslen.pld  <- c(obslen.fem,obslen.mal);
                    results.pld <- do.call(optimx, list(par=par.pld,
                                                        fn=.schnute2,
                                                        gr=NULL,
                                                        age=age.pld,
                                                        obslen=obslen.pld,
                                                        distr=distr,
                                                        method=method,
                                                        lower = -Inf,
                                                        upper = Inf,
                                                        control=list(),
                                                        hessian=hessian,
                                                        itnmax=itnmax,
                                                        output="estimate"));
                    rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                    results2                <- vector("list",3);
                    names(results2)         <- c("Females","Males","Pooled");
                    #Females
                    results2$Females        <- vector("list",length(method));
                    names(results2$Females) <- method;
                    temp.fem                <- attr(results.fem, "details")
                    for(i in 1:length(method))
                      {
                       results2$Females[[i]]$Type         <- m;
                       results2$Females[[i]]$Distr        <- distr;
                       par.names                  <- c("Linf","gamma","g1");
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          par.names <- c(par.names,"Dispersion")
                         }
                       results2$Females[[i]]$converg          <- "FAIL";
                       results2$Females[[i]]$kkt              <- NA;
                       results2$Females[[i]]$AIC              <- NA;
                       results2$Females[[i]]$bt.par           <- NA;
                       results2$Females[[i]]$num.grads        <- NA;
                       results2$Females[[i]]$bt.stdev         <- NA;
                       results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                       if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                         {
                          results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                          results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                          results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                          results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                          results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                v        <- matrix(0,length(par.names),length(par.names));
                                                if(distr == "apnormal" | distr == "aplnormal")
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                           mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                           cov=try(solve(temp.fem[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                                                else
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                           mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                           cov=try(solve(temp.fem[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                          results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                          results2$Females[[i]]$Cor              <- cor(v);
                          names(results2$Females[[i]]$num.grads) <- par.names;
                          names(results2$Females[[i]]$bt.par)    <- par.names;
                          names(results2$Females[[i]]$bt.stdev)  <- par.names;
                          colnames(results2$Females[[i]]$Cor)    <- par.names;
                          rownames(results2$Females[[i]]$Cor)    <- par.names;
                         }
                      }
                    rm(v);
                    rm(temp.fem,par.fem);
                    #Males
                    results2$Males        <- vector("list",length(method));
                    names(results2$Males)   <- method;
                    temp.mal                <- attr(results.mal, "details")
                    for(i in 1:length(method))
                      {
                       results2$Males[[i]]$Type         <- m;
                       results2$Males[[i]]$Distr        <- distr;
                       par.names                  <- c("Linf","gamma","g1");
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          par.names <- c(par.names,"Dispersion")
                         }
                       results2$Males[[i]]$converg          <- "FAIL";
                       results2$Males[[i]]$kkt              <- NA;
                       results2$Males[[i]]$AIC              <- NA;
                       results2$Males[[i]]$bt.par           <- NA;
                       results2$Males[[i]]$num.grads        <- NA;
                       results2$Males[[i]]$bt.stdev         <- NA;
                       results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                       if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                         {
                          results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                          results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                          results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                          results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                          results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                v        <- matrix(0,length(par.names),length(par.names));
                                                if(distr == "apnormal" | distr == "aplnormal")
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                           mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                           cov=try(solve(temp.mal[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                                                else
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                           mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                           cov=try(solve(temp.mal[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                          results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                          results2$Males[[i]]$Cor              <- cor(v);
                          names(results2$Males[[i]]$num.grads) <- par.names;
                          names(results2$Males[[i]]$bt.par)    <- par.names;
                          names(results2$Males[[i]]$bt.stdev)  <- par.names;
                          colnames(results2$Males[[i]]$Cor)    <- par.names;
                          rownames(results2$Males[[i]]$Cor)    <- par.names;
                         }
                      }
                    rm(v);
                    rm(temp.mal,par.mal);
                    #Pooled
                    results2$Pooled        <- vector("list",length(method));
                    names(results2$Pooled)   <- method;
                    temp.pld                <- attr(results.pld, "details")
                    for(i in 1:length(method))
                      {
                       results2$Pooled[[i]]$Type         <- m;
                       results2$Pooled[[i]]$Distr        <- distr;
                       par.names                  <- c("Linf","gamma","g1");
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          par.names <- c(par.names,"Dispersion")
                         }
                       results2$Pooled[[i]]$converg          <- "FAIL";
                       results2$Pooled[[i]]$kkt              <- NA;
                       results2$Pooled[[i]]$AIC              <- NA;
                       results2$Pooled[[i]]$bt.par           <- NA;
                       results2$Pooled[[i]]$num.grads        <- NA;
                       results2$Pooled[[i]]$bt.stdev         <- NA;
                       results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                       if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                         {
                          results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                          results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                          results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                          results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                          results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                v        <- matrix(0,length(par.names),length(par.names));
                                                if(distr == "apnormal" | distr == "aplnormal")
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                           mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                           cov=try(solve(temp.pld[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                                                else
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                           mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                           cov=try(solve(temp.pld[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                          results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                          results2$Pooled[[i]]$Cor              <- cor(v);
                          names(results2$Pooled[[i]]$num.grads) <- par.names;
                          names(results2$Pooled[[i]]$bt.par)    <- par.names;
                          names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                          colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                          rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                         }
                      }
                    rm(v);
                    rm(temp.pld,par.pld);
                   }
                else if(m == "schnute3")
                   {
                    if(distr == "apnormal" | distr == "aplnormal")
                      {
                       if(length(par) != 4)
                         {
                          stop("For a two-sexes Schnute model of type 3 with adjusted profile likelihood approximation to the normal \n
                                or lognormal, par must be vector of length 4")
                         }
                       par.fem     <- par[1:2];
                       par.mal     <- par[3:4];
                       par.pld     <- c(0.5*(par[1]+par[3]), 0.5*(par[2]+par[4]));
                      }
                    if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                      {
                       if(length(par) != 6)
                         {
                          stop("For a two-sexes Schnute model of type 3 with exact normal, lognormal or gamma likelihood, \n
                                par must be vector of length 6")
                         }
                       par.fem     <- par[c(1:2,5)];
                       par.mal     <- par[c(3:4,6)];
                       par.pld     <- c(0.5*(par[1]+par[3]), 0.5*(par[2]+par[4]), 0.5*(par[5]+par[6]));
                      }
                    age.fem     <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2];
                    obslen.fem  <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3];
                    results.fem <- do.call(optimx, list(par=par.fem,
                                                        fn=.schnute3,
                                                        gr=NULL,
                                                        age=age.fem,
                                                        obslen=obslen.fem,
                                                        distr=distr,
                                                        method=method,
                                                        lower = -Inf,
                                                        upper = Inf,
                                                        control=list(),
                                                        hessian=hessian,
                                                        itnmax=itnmax,
                                                        output="estimate"));
                    age.mal     <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2];
                    obslen.mal  <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3];
                    results.mal <- do.call(optimx, list(par=par.mal,
                                                        fn=.schnute3,
                                                        gr=NULL,
                                                        age=age.mal,
                                                        obslen=obslen.mal,
                                                        distr=distr,
                                                        method=method,
                                                        lower = -Inf,
                                                        upper = Inf,
                                                        control=list(),
                                                        hessian=hessian,
                                                        itnmax=itnmax,
                                                        output="estimate"));
                    age.pld     <- c(age.fem,age.mal);
                    obslen.pld  <- c(obslen.fem,obslen.mal);
                    results.pld <- do.call(optimx, list(par=par.pld,
                                                        fn=.schnute3,
                                                        gr=NULL,
                                                        age=age.pld,
                                                        obslen=obslen.pld,
                                                        distr=distr,
                                                        method=method,
                                                        lower = -Inf,
                                                        upper = Inf,
                                                        control=list(),
                                                        hessian=hessian,
                                                        itnmax=itnmax,
                                                        output="estimate"));
                    rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                    results2                <- vector("list",3);
                    names(results2)         <- c("Females","Males","Pooled");
                    #Females
                    results2$Females        <- vector("list",length(method));
                    names(results2$Females) <- method;
                    temp.fem                <- attr(results.fem, "details")
                    for(i in 1:length(method))
                      {
                       results2$Females[[i]]$Type         <- m;
                       results2$Females[[i]]$Distr        <- distr;
                       par.names                  <- c("g1","g2");
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          par.names <- c(par.names,"Dispersion")
                         }
                       results2$Females[[i]]$converg          <- "FAIL";
                       results2$Females[[i]]$kkt              <- NA;
                       results2$Females[[i]]$AIC              <- NA;
                       results2$Females[[i]]$bt.par           <- NA;
                       results2$Females[[i]]$num.grads        <- NA;
                       results2$Females[[i]]$bt.stdev         <- NA;
                       results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                       if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                         {
                          results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                          results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                          results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                          results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                          results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                v        <- matrix(0,length(par.names),length(par.names));
                                                if(distr == "apnormal" | distr == "aplnormal")
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2)),
                                                                           mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                           cov=try(solve(temp.fem[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                                                else
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                           mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                           cov=try(solve(temp.fem[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                          results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                          results2$Females[[i]]$Cor              <- cor(v);
                          names(results2$Females[[i]]$num.grads) <- par.names;
                          names(results2$Females[[i]]$bt.par)    <- par.names;
                          names(results2$Females[[i]]$bt.stdev)  <- par.names;
                          colnames(results2$Females[[i]]$Cor)    <- par.names;
                          rownames(results2$Females[[i]]$Cor)    <- par.names;
                         }
                      }
                    rm(v);
                    rm(temp.fem,par.fem);
                    #Males
                    results2$Males        <- vector("list",length(method));
                    names(results2$Males)   <- method;
                    temp.mal                <- attr(results.mal, "details")
                    for(i in 1:length(method))
                      {
                       results2$Males[[i]]$Type         <- m;
                       results2$Males[[i]]$Distr        <- distr;
                       par.names                  <- c("g1","g2");
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          par.names <- c(par.names,"Dispersion")
                         }
                       results2$Males[[i]]$converg          <- "FAIL";
                       results2$Males[[i]]$kkt              <- NA;
                       results2$Males[[i]]$AIC              <- NA;
                       results2$Males[[i]]$bt.par           <- NA;
                       results2$Males[[i]]$num.grads        <- NA;
                       results2$Males[[i]]$bt.stdev         <- NA;
                       results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                       if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                         {
                          results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                          results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                          results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                          results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                          results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                v        <- matrix(0,length(par.names),length(par.names));
                                                if(distr == "apnormal" | distr == "aplnormal")
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2)),
                                                                           mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                           cov=try(solve(temp.mal[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                                                else
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                           mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                           cov=try(solve(temp.mal[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                          results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                          results2$Males[[i]]$Cor              <- cor(v);
                          names(results2$Males[[i]]$num.grads) <- par.names;
                          names(results2$Males[[i]]$bt.par)    <- par.names;
                          names(results2$Males[[i]]$bt.stdev)  <- par.names;
                          colnames(results2$Males[[i]]$Cor)    <- par.names;
                          rownames(results2$Males[[i]]$Cor)    <- par.names;
                         }
                      }
                    rm(v);
                    rm(temp.mal,par.mal);
                    #Pooled
                    results2$Pooled        <- vector("list",length(method));
                    names(results2$Pooled)   <- method;
                    temp.pld                <- attr(results.pld, "details")
                    for(i in 1:length(method))
                      {
                       results2$Pooled[[i]]$Type         <- m;
                       results2$Pooled[[i]]$Distr        <- distr;
                       par.names                  <- c("g1","g2");
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          par.names <- c(par.names,"Dispersion")
                         }
                       results2$Pooled[[i]]$converg          <- "FAIL";
                       results2$Pooled[[i]]$kkt              <- NA;
                       results2$Pooled[[i]]$AIC              <- NA;
                       results2$Pooled[[i]]$bt.par           <- NA;
                       results2$Pooled[[i]]$num.grads        <- NA;
                       results2$Pooled[[i]]$bt.stdev         <- NA;
                       results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                       if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                         {
                          results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                          results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                          results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                          results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                          results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                v        <- matrix(0,length(par.names),length(par.names));
                                                if(distr == "apnormal" | distr == "aplnormal")
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2)),
                                                                           mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                           cov=try(solve(temp.pld[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                                                else
                                                  {
                                                   v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                           mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                           cov=try(solve(temp.pld[i,]$nhatend)),
                                                                           ses=FALSE)
                                                  }
                          results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                          results2$Pooled[[i]]$Cor              <- cor(v);
                          names(results2$Pooled[[i]]$num.grads) <- par.names;
                          names(results2$Pooled[[i]]$bt.par)    <- par.names;
                          names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                          colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                          rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                         }
                      }
                    rm(v);
                    rm(temp.pld,par.pld);
                   }
               }
             else if(x$Properties$Sex == "Total")
               {
                if(unsex.action == "split")
                  {
                   unsex.n       <- length(which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length)));
                   unsex.sample  <- sample(1:unsex.n,size=round(unsex.n/2));
                   unsex.fem.age <- x$Data$Unsexed[which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length)),2][unsex.sample];
                   unsex.fem.len <- x$Data$Unsexed[which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length)),3][unsex.sample];
                   unsex.mal.age <- x$Data$Unsexed[which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length)),2][-unsex.sample];
                   unsex.mal.len <- x$Data$Unsexed[which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length)),3][-unsex.sample];
                   if(m == "vonber")
                      {
                       if(distr == "apnormal" | distr == "aplnormal")
                         {
                          if(length(par) != 6)
                            {
                             stop("For a two-sexes Von Bertalanffy model with adjusted profile likelihood approximation to the normal \n
                                   or lognormal, par must be vector of length 6")
                            }
                          par.fem     <- par[1:3];
                          par.mal     <- par[4:6];
                          par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]));
                         }
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          if(length(par) != 8)
                            {
                             stop("For a two-sexes Von Bertalanffy model with exact normal, lognormal or gamma likelihood, \n
                                   par must be vector of length 8")
                            }
                          par.fem     <- par[c(1:3,7)];
                          par.mal     <- par[c(4:6,8)];
                          par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]), 0.5*(par[7]+par[8]));
                         }
                       age.fem     <- c(x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2],unsex.fem.age);
                       obslen.fem  <- c(x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3],unsex.fem.len);
                       results.fem <- do.call(optimx, list(par=par.fem,
                                                           fn=.vonBer,
                                                           gr=NULL,
                                                           age=age.fem,
                                                           obslen=obslen.fem,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.mal     <- c(x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2],unsex.mal.age);
                       obslen.mal  <- c(x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Mmales[,3])),3],unsex.mal.len);
                       results.mal <- do.call(optimx, list(par=par.mal,
                                                           fn=.vonBer,
                                                           gr=NULL,
                                                           age=age.mal,
                                                           obslen=obslen.mal,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.pld     <- c(age.fem,age.mal);
                       obslen.pld  <- c(obslen.fem,obslen.mal);
                       results.pld <- do.call(optimx, list(par=par.pld,
                                                           fn=.vonBer,
                                                           gr=NULL,
                                                           age=age.pld,
                                                           obslen=obslen.pld,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                       results2                <- vector("list",3);
                       names(results2)         <- c("Females","Males","Pooled");
                       #Females
                       results2$Females        <- vector("list",length(method));
                       names(results2$Females) <- method;
                       temp.fem                <- attr(results.fem, "details")
                       for(i in 1:length(method))
                         {
                          results2$Females[[i]]$Type         <- m;
                          results2$Females[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","L0","K");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Females[[i]]$converg          <- "FAIL";
                          results2$Females[[i]]$kkt              <- NA;
                          results2$Females[[i]]$AIC              <- NA;
                          results2$Females[[i]]$bt.par           <- NA;
                          results2$Females[[i]]$num.grads        <- NA;
                          results2$Females[[i]]$bt.stdev         <- NA;
                          results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                            {
                             results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                             results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                             results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                             results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                             results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Females[[i]]$Cor              <- cor(v);
                             names(results2$Females[[i]]$num.grads) <- par.names;
                             names(results2$Females[[i]]$bt.par)    <- par.names;
                             names(results2$Females[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Females[[i]]$Cor)    <- par.names;
                             rownames(results2$Females[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.fem,par.fem);
                       #Males
                       results2$Males        <- vector("list",length(method));
                       names(results2$Males)   <- method;
                       temp.mal                <- attr(results.mal, "details")
                       for(i in 1:length(method))
                         {
                          results2$Males[[i]]$Type         <- m;
                          results2$Males[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","L0","K");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Males[[i]]$converg          <- "FAIL";
                          results2$Males[[i]]$kkt              <- NA;
                          results2$Males[[i]]$AIC              <- NA;
                          results2$Males[[i]]$bt.par           <- NA;
                          results2$Males[[i]]$num.grads        <- NA;
                          results2$Males[[i]]$bt.stdev         <- NA;
                          results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                            {
                             results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                             results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                             results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                             results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                             results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Males[[i]]$Cor              <- cor(v);
                             names(results2$Males[[i]]$num.grads) <- par.names;
                             names(results2$Males[[i]]$bt.par)    <- par.names;
                             names(results2$Males[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Males[[i]]$Cor)    <- par.names;
                             rownames(results2$Males[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.mal,par.mal);
                       #Pooled
                       results2$Pooled        <- vector("list",length(method));
                       names(results2$Pooled)   <- method;
                       temp.pld                <- attr(results.pld, "details")
                       for(i in 1:length(method))
                         {
                          results2$Pooled[[i]]$Type         <- m;
                          results2$Pooled[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","L0","K");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Pooled[[i]]$converg          <- "FAIL";
                          results2$Pooled[[i]]$kkt              <- NA;
                          results2$Pooled[[i]]$AIC              <- NA;
                          results2$Pooled[[i]]$bt.par           <- NA;
                          results2$Pooled[[i]]$num.grads        <- NA;
                          results2$Pooled[[i]]$bt.stdev         <- NA;
                          results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                            {
                             results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                             results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                             results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                             results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                             results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Pooled[[i]]$Cor              <- cor(v);
                             names(results2$Pooled[[i]]$num.grads) <- par.names;
                             names(results2$Pooled[[i]]$bt.par)    <- par.names;
                             names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                             rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.pld,par.pld);
                      }
                   else if(m == "schnute1")
                      {
                       if(distr == "apnormal" | distr == "aplnormal")
                         {
                          if(length(par) != 10)
                            {
                             stop("For a two-sexes Schnute model of type 1 with adjusted profile likelihood approximation to the normal \n
                             or lognormal, par must be vector of length 10")
                            }
                          par.fem     <- par[1:5];
                          par.mal     <- par[6:10];
                          par.pld     <- c(0.5*(par[1]+par[6]), 0.5*(par[2]+par[7]), 0.5*(par[3]+par[8]), 0.5*(par[4]+par[9]), 0.5*(par[5]+par[10]));
                         }
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          if(length(par) != 12)
                            {
                             stop("For a two-sexes Schnute model of type 1 with exact normal, lognormal or gamma likelihood, \n
                             par must be vector of length 12")
                            }
                          par.fem     <- par[c(1:5,11)];
                          par.mal     <- par[c(6:10,12)];
                          par.pld     <- c(0.5*(par[1]+par[6]), 0.5*(par[2]+par[7]), 0.5*(par[3]+par[8]), 0.5*(par[4]+par[9]), 0.5*(par[5]+par[10]), 0.5*(par[11]+par[12]));
                         }
                       age.fem     <- c(x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2],unsex.fem.age);
                       obslen.fem  <- c(x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3],unsex.fem.len);
                       results.fem <- do.call(optimx, list(par=par.fem,
                                                           fn=.schnute1,
                                                           gr=NULL,
                                                           age=age.fem,
                                                           obslen=obslen.fem,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.mal     <- c(x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2],unsex.mal.age);
                       obslen.mal  <- c(x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3],unsex.mal.len);
                       results.mal <- do.call(optimx, list(par=par.mal,
                                                           fn=.schnute1,
                                                           gr=NULL,
                                                           age=age.mal,
                                                           obslen=obslen.mal,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.pld     <- c(age.fem,age.mal);
                       obslen.pld  <- c(obslen.fem,obslen.mal);
                       results.pld <- do.call(optimx, list(par=par.pld,
                                                           fn=.schnute1,
                                                           gr=NULL,
                                                           age=age.pld,
                                                           obslen=obslen.pld,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                       results2                <- vector("list",3);
                       names(results2)         <- c("Females","Males","Pooled");
                       #Females
                       results2$Females        <- vector("list",length(method));
                       names(results2$Females) <- method;
                       temp.fem                <- attr(results.fem, "details")
                       for(i in 1:length(method))
                         {
                          results2$Females[[i]]$Type         <- m;
                          results2$Females[[i]]$Distr        <- distr;
                          par.names                  <- c("l1","a1","mu","g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Females[[i]]$converg          <- "FAIL";
                          results2$Females[[i]]$kkt              <- NA;
                          results2$Females[[i]]$AIC              <- NA;
                          results2$Females[[i]]$bt.par           <- NA;
                          results2$Females[[i]]$num.grads        <- NA;
                          results2$Females[[i]]$bt.stdev         <- NA;
                          results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                            {
                             results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                             results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                             results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                             results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                             results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Females[[i]]$Cor              <- cor(v);
                             names(results2$Females[[i]]$num.grads) <- par.names;
                             names(results2$Females[[i]]$bt.par)    <- par.names;
                             names(results2$Females[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Females[[i]]$Cor)    <- par.names;
                             rownames(results2$Females[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.fem,par.fem);
                       #Males
                       results2$Males        <- vector("list",length(method));
                       names(results2$Males)   <- method;
                       temp.mal                <- attr(results.mal, "details")
                       for(i in 1:length(method))
                         {
                          results2$Males[[i]]$Type         <- m;
                          results2$Males[[i]]$Distr        <- distr;
                          par.names                  <- c("l1","a1","mu","g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Males[[i]]$converg          <- "FAIL";
                          results2$Males[[i]]$kkt              <- NA;
                          results2$Males[[i]]$AIC              <- NA;
                          results2$Males[[i]]$bt.par           <- NA;
                          results2$Males[[i]]$num.grads        <- NA;
                          results2$Males[[i]]$bt.stdev         <- NA;
                          results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                            {
                             results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                             results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                             results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                             results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                             results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Males[[i]]$Cor              <- cor(v);
                             names(results2$Males[[i]]$num.grads) <- par.names;
                             names(results2$Males[[i]]$bt.par)    <- par.names;
                             names(results2$Males[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Males[[i]]$Cor)    <- par.names;
                             rownames(results2$Males[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.mal,par.mal);
                       #Pooled
                       results2$Pooled        <- vector("list",length(method));
                       names(results2$Pooled)   <- method;
                       temp.pld                <- attr(results.pld, "details")
                       for(i in 1:length(method))
                         {
                          results2$Pooled[[i]]$Type         <- m;
                          results2$Pooled[[i]]$Distr        <- distr;
                          par.names                  <- c("l1","a1","mu","g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Pooled[[i]]$converg          <- "FAIL";
                          results2$Pooled[[i]]$kkt              <- NA;
                          results2$Pooled[[i]]$AIC              <- NA;
                          results2$Pooled[[i]]$bt.par           <- NA;
                          results2$Pooled[[i]]$num.grads        <- NA;
                          results2$Pooled[[i]]$bt.stdev         <- NA;
                          results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                            {
                             results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                             results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                             results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                             results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                             results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Pooled[[i]]$Cor              <- cor(v);
                             names(results2$Pooled[[i]]$num.grads) <- par.names;
                             names(results2$Pooled[[i]]$bt.par)    <- par.names;
                             names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                             rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.pld,par.pld);
                      }
                   else if(m == "schnute2")
                      {
                       if(distr == "apnormal" | distr == "aplnormal")
                         {
                          if(length(par) != 6)
                            {
                             stop("For a two-sexes Schnute model of type 2 with adjusted profile likelihood approximation to the normal \n
                                   or lognormal, par must be vector of length 6")
                            }
                          par.fem     <- par[1:3];
                          par.mal     <- par[4:6];
                          par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]));
                         }
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          if(length(par) != 8)
                            {
                             stop("For a two-sexes Schnute model of type 2 with exact normal, lognormal or gamma likelihood, \n
                                   par must be vector of length 8")
                            }
                          par.fem     <- par[c(1:3,7)];
                          par.mal     <- par[c(4:6,8)];
                          par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]), 0.5*(par[7]+par[7]));
                         }
                       age.fem     <- c(x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2],unsex.fem.age);
                       obslen.fem  <- c(x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3],unsex.fem.len);
                       results.fem <- do.call(optimx, list(par=par.fem,
                                                           fn=.schnute2,
                                                           gr=NULL,
                                                           age=age.fem,
                                                           obslen=obslen.fem,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.mal     <- c(x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2],unsex.mal.age);
                       obslen.mal  <- c(x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3],unsex.mal.len);
                       results.mal <- do.call(optimx, list(par=par.mal,
                                                           fn=.schnute2,
                                                           gr=NULL,
                                                           age=age.mal,
                                                           obslen=obslen.mal,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.pld     <- c(age.fem,age.mal);
                       obslen.pld  <- c(obslen.fem,obslen.mal);
                       results.pld <- do.call(optimx, list(par=par.pld,
                                                           fn=.schnute2,
                                                           gr=NULL,
                                                           age=age.pld,
                                                           obslen=obslen.pld,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                       results2                <- vector("list",3);
                       names(results2)         <- c("Females","Males","Pooled");
                       #Females
                       results2$Females        <- vector("list",length(method));
                       names(results2$Females) <- method;
                       temp.fem                <- attr(results.fem, "details")
                       for(i in 1:length(method))
                         {
                          results2$Females[[i]]$Type         <- m;
                          results2$Females[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","gamma","g1");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Females[[i]]$converg          <- "FAIL";
                          results2$Females[[i]]$kkt              <- NA;
                          results2$Females[[i]]$AIC              <- NA;
                          results2$Females[[i]]$bt.par           <- NA;
                          results2$Females[[i]]$num.grads        <- NA;
                          results2$Females[[i]]$bt.stdev         <- NA;
                          results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                            {
                             results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                             results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                             results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                             results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                             results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Females[[i]]$Cor              <- cor(v);
                             names(results2$Females[[i]]$num.grads) <- par.names;
                             names(results2$Females[[i]]$bt.par)    <- par.names;
                             names(results2$Females[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Females[[i]]$Cor)    <- par.names;
                             rownames(results2$Females[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.fem,par.fem);
                       #Males
                       results2$Males        <- vector("list",length(method));
                       names(results2$Males)   <- method;
                       temp.mal                <- attr(results.mal, "details")
                       for(i in 1:length(method))
                         {
                          results2$Males[[i]]$Type         <- m;
                          results2$Males[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","gamma","g1");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Males[[i]]$converg          <- "FAIL";
                          results2$Males[[i]]$kkt              <- NA;
                          results2$Males[[i]]$AIC              <- NA;
                          results2$Males[[i]]$bt.par           <- NA;
                          results2$Males[[i]]$num.grads        <- NA;
                          results2$Males[[i]]$bt.stdev         <- NA;
                          results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                            {
                             results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                             results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                             results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                             results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                             results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Males[[i]]$Cor              <- cor(v);
                             names(results2$Males[[i]]$num.grads) <- par.names;
                             names(results2$Males[[i]]$bt.par)    <- par.names;
                             names(results2$Males[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Males[[i]]$Cor)    <- par.names;
                             rownames(results2$Males[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.mal,par.mal);
                       #Pooled
                       results2$Pooled        <- vector("list",length(method));
                       names(results2$Pooled)   <- method;
                       temp.pld                <- attr(results.pld, "details")
                       for(i in 1:length(method))
                         {
                          results2$Pooled[[i]]$Type         <- m;
                          results2$Pooled[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","gamma","g1");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Pooled[[i]]$converg          <- "FAIL";
                          results2$Pooled[[i]]$kkt              <- NA;
                          results2$Pooled[[i]]$AIC              <- NA;
                          results2$Pooled[[i]]$bt.par           <- NA;
                          results2$Pooled[[i]]$num.grads        <- NA;
                          results2$Pooled[[i]]$bt.stdev         <- NA;
                          results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                            {
                             results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                             results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                             results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                             results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                             results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Pooled[[i]]$Cor              <- cor(v);
                             names(results2$Pooled[[i]]$num.grads) <- par.names;
                             names(results2$Pooled[[i]]$bt.par)    <- par.names;
                             names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                             rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.pld,par.pld);
                      }
                   else if(m == "schnute3")
                      {
                       if(distr == "apnormal" | distr == "aplnormal")
                         {
                          if(length(par) != 4)
                            {
                             stop("For a two-sexes Schnute model of type 3 with adjusted profile likelihood approximation to the normal \n
                                   or lognormal, par must be vector of length 4")
                            }
                          par.fem     <- par[1:2];
                          par.mal     <- par[3:4];
                          par.pld     <- c(0.5*(par[1]+par[3]), 0.5*(par[2]+par[4]));
                         }
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          if(length(par) != 6)
                            {
                             stop("For a two-sexes Schnute model of type 3 with exact normal, lognormal or gamma likelihood, \n
                                   par must be vector of length 6")
                            }
                          par.fem     <- par[c(1:2,5)];
                          par.mal     <- par[c(3:4,6)];
                          par.pld     <- c(0.5*(par[1]+par[3]), 0.5*(par[2]+par[4]), 0.5*(par[5]+par[6]));
                         }
                       age.fem     <- c(x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2],unsex.fem.age);
                       obslen.fem  <- c(x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3],unsex.fem.len);
                       results.fem <- do.call(optimx, list(par=par.fem,
                                                           fn=.schnute3,
                                                           gr=NULL,
                                                           age=age.fem,
                                                           obslen=obslen.fem,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.mal     <- c(x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2],unsex.mal.age);
                       obslen.mal  <- c(x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3],unsex.mal.len);
                       results.mal <- do.call(optimx, list(par=par.mal,
                                                           fn=.schnute3,
                                                           gr=NULL,
                                                           age=age.mal,
                                                           obslen=obslen.mal,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.pld     <- c(age.fem,age.mal);
                       obslen.pld  <- c(obslen.fem,obslen.mal);
                       results.pld <- do.call(optimx, list(par=par.pld,
                                                           fn=.schnute3,
                                                           gr=NULL,
                                                           age=age.pld,
                                                           obslen=obslen.pld,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                       results2                <- vector("list",3);
                       names(results2)         <- c("Females","Males","Pooled");
                       #Females
                       results2$Females        <- vector("list",length(method));
                       names(results2$Females) <- method;
                       temp.fem                <- attr(results.fem, "details")
                       for(i in 1:length(method))
                         {
                          results2$Females[[i]]$Type         <- m;
                          results2$Females[[i]]$Distr        <- distr;
                          par.names                  <- c("g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Females[[i]]$converg          <- "FAIL";
                          results2$Females[[i]]$kkt              <- NA;
                          results2$Females[[i]]$AIC              <- NA;
                          results2$Females[[i]]$bt.par           <- NA;
                          results2$Females[[i]]$num.grads        <- NA;
                          results2$Females[[i]]$bt.stdev         <- NA;
                          results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                            {
                             results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                             results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                             results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                             results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                             results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Females[[i]]$Cor              <- cor(v);
                             names(results2$Females[[i]]$num.grads) <- par.names;
                             names(results2$Females[[i]]$bt.par)    <- par.names;
                             names(results2$Females[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Females[[i]]$Cor)    <- par.names;
                             rownames(results2$Females[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.fem,par.fem);
                       #Males
                       results2$Males        <- vector("list",length(method));
                       names(results2$Males)   <- method;
                       temp.mal                <- attr(results.mal, "details")
                       for(i in 1:length(method))
                         {
                          results2$Males[[i]]$Type         <- m;
                          results2$Males[[i]]$Distr        <- distr;
                          par.names                  <- c("g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Males[[i]]$converg          <- "FAIL";
                          results2$Males[[i]]$kkt              <- NA;
                          results2$Males[[i]]$AIC              <- NA;
                          results2$Males[[i]]$bt.par           <- NA;
                          results2$Males[[i]]$num.grads        <- NA;
                          results2$Males[[i]]$bt.stdev         <- NA;
                          results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                            {
                             results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                             results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                             results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                             results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                             results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Males[[i]]$Cor              <- cor(v);
                             names(results2$Males[[i]]$num.grads) <- par.names;
                             names(results2$Males[[i]]$bt.par)    <- par.names;
                             names(results2$Males[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Males[[i]]$Cor)    <- par.names;
                             rownames(results2$Males[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.mal,par.mal);
                       #Pooled
                       results2$Pooled        <- vector("list",length(method));
                       names(results2$Pooled)   <- method;
                       temp.pld                <- attr(results.pld, "details")
                       for(i in 1:length(method))
                         {
                          results2$Pooled[[i]]$Type         <- m;
                          results2$Pooled[[i]]$Distr        <- distr;
                          par.names                  <- c("g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Pooled[[i]]$converg          <- "FAIL";
                          results2$Pooled[[i]]$kkt              <- NA;
                          results2$Pooled[[i]]$AIC              <- NA;
                          results2$Pooled[[i]]$bt.par           <- NA;
                          results2$Pooled[[i]]$num.grads        <- NA;
                          results2$Pooled[[i]]$bt.stdev         <- NA;
                          results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                            {
                             results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                             results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                             results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                             results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                             results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Pooled[[i]]$Cor              <- cor(v);
                             names(results2$Pooled[[i]]$num.grads) <- par.names;
                             names(results2$Pooled[[i]]$bt.par)    <- par.names;
                             names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                             rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.pld,par.pld);
                      }
                  }
                if(unsex.action == "males")
                  {
                   unsex.mal.age <- x$Data$Unsexed[which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length)),2];
                   unsex.mal.len <- x$Data$Unsexed[which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length)),3];
                   if(m == "vonber")
                      {
                       if(distr == "apnormal" | distr == "aplnormal")
                         {
                          if(length(par) != 6)
                            {
                             stop("For a two-sexes Von Bertalanffy model with adjusted profile likelihood approximation to the normal \n
                                   or lognormal, par must be vector of length 6")
                            }
                          par.fem     <- par[1:3];
                          par.mal     <- par[4:6];
                          par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]));
                         }
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          if(length(par) != 8)
                            {
                             stop("For a two-sexes Von Bertalanffy model with exact normal, lognormal or gamma likelihood, \n
                                   par must be vector of length 8")
                            }
                          par.fem     <- par[c(1:3,7)];
                          par.mal     <- par[c(4:6,8)];
                          par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]), 0.5*(par[7]+par[8]));
                         }
                       age.fem     <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2];
                       obslen.fem  <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3];
                       results.fem <- do.call(optimx, list(par=par.fem,
                                                           fn=.vonBer,
                                                           gr=NULL,
                                                           age=age.fem,
                                                           obslen=obslen.fem,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.mal     <- c(x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2],unsex.mal.age);
                       obslen.mal  <- c(x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3],unsex.mal.len);
                       results.mal <- do.call(optimx, list(par=par.mal,
                                                           fn=.vonBer,
                                                           gr=NULL,
                                                           age=age.mal,
                                                           obslen=obslen.mal,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.pld     <- c(age.fem,age.mal);
                       obslen.pld  <- c(obslen.fem,obslen.mal);
                       results.pld <- do.call(optimx, list(par=par.pld,
                                                           fn=.vonBer,
                                                           gr=NULL,
                                                           age=age.pld,
                                                           obslen=obslen.pld,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                       results2                <- vector("list",3);
                       names(results2)         <- c("Females","Males","Pooled");
                       #Females
                       results2$Females        <- vector("list",length(method));
                       names(results2$Females) <- method;
                       temp.fem                <- attr(results.fem, "details")
                       for(i in 1:length(method))
                         {
                          results2$Females[[i]]$Type         <- m;
                          results2$Females[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","L0","K");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Females[[i]]$converg          <- "FAIL";
                          results2$Females[[i]]$kkt              <- NA;
                          results2$Females[[i]]$AIC              <- NA;
                          results2$Females[[i]]$bt.par           <- NA;
                          results2$Females[[i]]$num.grads        <- NA;
                          results2$Females[[i]]$bt.stdev         <- NA;
                          results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                            {
                             results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                             results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                             results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                             results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                             results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Females[[i]]$Cor              <- cor(v);
                             names(results2$Females[[i]]$num.grads) <- par.names;
                             names(results2$Females[[i]]$bt.par)    <- par.names;
                             names(results2$Females[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Females[[i]]$Cor)    <- par.names;
                             rownames(results2$Females[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.fem,par.fem);
                       #Males
                       results2$Males        <- vector("list",length(method));
                       names(results2$Males)   <- method;
                       temp.mal                <- attr(results.mal, "details")
                       for(i in 1:length(method))
                         {
                          results2$Males[[i]]$Type         <- m;
                          results2$Males[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","L0","K");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Males[[i]]$converg          <- "FAIL";
                          results2$Males[[i]]$kkt              <- NA;
                          results2$Males[[i]]$AIC              <- NA;
                          results2$Males[[i]]$bt.par           <- NA;
                          results2$Males[[i]]$num.grads        <- NA;
                          results2$Males[[i]]$bt.stdev         <- NA;
                          results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                            {
                             results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                             results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                             results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                             results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                             results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Males[[i]]$Cor              <- cor(v);
                             names(results2$Males[[i]]$num.grads) <- par.names;
                             names(results2$Males[[i]]$bt.par)    <- par.names;
                             names(results2$Males[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Males[[i]]$Cor)    <- par.names;
                             rownames(results2$Males[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.mal,par.mal);
                       #Pooled
                       results2$Pooled        <- vector("list",length(method));
                       names(results2$Pooled)   <- method;
                       temp.pld                <- attr(results.pld, "details")
                       for(i in 1:length(method))
                         {
                          results2$Pooled[[i]]$Type         <- m;
                          results2$Pooled[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","L0","K");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Pooled[[i]]$converg          <- "FAIL";
                          results2$Pooled[[i]]$kkt              <- NA;
                          results2$Pooled[[i]]$AIC              <- NA;
                          results2$Pooled[[i]]$bt.par           <- NA;
                          results2$Pooled[[i]]$num.grads        <- NA;
                          results2$Pooled[[i]]$bt.stdev         <- NA;
                          results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                            {
                             results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                             results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                             results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                             results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                             results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Pooled[[i]]$Cor              <- cor(v);
                             names(results2$Pooled[[i]]$num.grads) <- par.names;
                             names(results2$Pooled[[i]]$bt.par)    <- par.names;
                             names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                             rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.pld,par.pld);
                      }
                   else if(m == "schnute1")
                      {
                       if(distr == "apnormal" | distr == "aplnormal")
                         {
                          if(length(par) != 10)
                            {
                             stop("For a two-sexes Schnute model of type 1 with adjusted profile likelihood approximation to the normal \n
                             or lognormal, par must be vector of length 10")
                            }
                          par.fem     <- par[1:5];
                          par.mal     <- par[6:10];
                          par.pld     <- c(0.5*(par[1]+par[6]), 0.5*(par[2]+par[7]), 0.5*(par[3]+par[8]), 0.5*(par[4]+par[9]), 0.5*(par[5]+par[10]));
                         }
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          if(length(par) != 12)
                            {
                             stop("For a two-sexes Schnute model of type 1 with exact normal, lognormal or gamma likelihood, \n
                             par must be vector of length 12")
                            }
                          par.fem     <- par[c(1:5,11)];
                          par.mal     <- par[c(6:10,12)];
                          par.pld     <- c(0.5*(par[1]+par[6]), 0.5*(par[2]+par[7]), 0.5*(par[3]+par[8]), 0.5*(par[4]+par[9]), 0.5*(par[5]+par[10]), 0.5*(par[11]+par[12]));
                         }
                       age.fem     <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2];
                       obslen.fem  <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3];
                       results.fem <- do.call(optimx, list(par=par.fem,
                                                           fn=.schnute1,
                                                           gr=NULL,
                                                           age=age.fem,
                                                           obslen=obslen.fem,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.mal     <- c(x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2],unsex.mal.age);
                       obslen.mal  <- c(x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3],unsex.mal.len);
                       results.mal <- do.call(optimx, list(par=par.mal,
                                                           fn=.schnute1,
                                                           gr=NULL,
                                                           age=age.mal,
                                                           obslen=obslen.mal,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.pld     <- c(age.fem,age.mal);
                       obslen.pld  <- c(obslen.fem,obslen.mal);
                       results.pld <- do.call(optimx, list(par=par.pld,
                                                           fn=.schnute1,
                                                           gr=NULL,
                                                           age=age.pld,
                                                           obslen=obslen.pld,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                       results2                <- vector("list",3);
                       names(results2)         <- c("Females","Males","Pooled");
                       #Females
                       results2$Females        <- vector("list",length(method));
                       names(results2$Females) <- method;
                       temp.fem                <- attr(results.fem, "details")
                       for(i in 1:length(method))
                         {
                          results2$Females[[i]]$Type         <- m;
                          results2$Females[[i]]$Distr        <- distr;
                          par.names                  <- c("l1","a1","mu","g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Females[[i]]$converg          <- "FAIL";
                          results2$Females[[i]]$kkt              <- NA;
                          results2$Females[[i]]$AIC              <- NA;
                          results2$Females[[i]]$bt.par           <- NA;
                          results2$Females[[i]]$num.grads        <- NA;
                          results2$Females[[i]]$bt.stdev         <- NA;
                          results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                            {
                             results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                             results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                             results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                             results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                             results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Females[[i]]$Cor              <- cor(v);
                             names(results2$Females[[i]]$num.grads) <- par.names;
                             names(results2$Females[[i]]$bt.par)    <- par.names;
                             names(results2$Females[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Females[[i]]$Cor)    <- par.names;
                             rownames(results2$Females[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.fem,par.fem);
                       #Males
                       results2$Males        <- vector("list",length(method));
                       names(results2$Males)   <- method;
                       temp.mal                <- attr(results.mal, "details")
                       for(i in 1:length(method))
                         {
                          results2$Males[[i]]$Type         <- m;
                          results2$Males[[i]]$Distr        <- distr;
                          par.names                  <- c("l1","a1","mu","g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Males[[i]]$converg          <- "FAIL";
                          results2$Males[[i]]$kkt              <- NA;
                          results2$Males[[i]]$AIC              <- NA;
                          results2$Males[[i]]$bt.par           <- NA;
                          results2$Males[[i]]$num.grads        <- NA;
                          results2$Males[[i]]$bt.stdev         <- NA;
                          results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                            {
                             results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                             results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                             results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                             results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                             results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Males[[i]]$Cor              <- cor(v);
                             names(results2$Males[[i]]$num.grads) <- par.names;
                             names(results2$Males[[i]]$bt.par)    <- par.names;
                             names(results2$Males[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Males[[i]]$Cor)    <- par.names;
                             rownames(results2$Males[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.mal,par.mal);
                       #Pooled
                       results2$Pooled        <- vector("list",length(method));
                       names(results2$Pooled)   <- method;
                       temp.pld                <- attr(results.pld, "details")
                       for(i in 1:length(method))
                         {
                          results2$Pooled[[i]]$Type         <- m;
                          results2$Pooled[[i]]$Distr        <- distr;
                          par.names                  <- c("l1","a1","mu","g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Pooled[[i]]$converg          <- "FAIL";
                          results2$Pooled[[i]]$kkt              <- NA;
                          results2$Pooled[[i]]$AIC              <- NA;
                          results2$Pooled[[i]]$bt.par           <- NA;
                          results2$Pooled[[i]]$num.grads        <- NA;
                          results2$Pooled[[i]]$bt.stdev         <- NA;
                          results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                            {
                             results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                             results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                             results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                             results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                             results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Pooled[[i]]$Cor              <- cor(v);
                             names(results2$Pooled[[i]]$num.grads) <- par.names;
                             names(results2$Pooled[[i]]$bt.par)    <- par.names;
                             names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                             rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.pld,par.pld);
                      }
                   else if(m == "schnute2")
                      {
                       if(distr == "apnormal" | distr == "aplnormal")
                         {
                          if(length(par) != 6)
                            {
                             stop("For a two-sexes Schnute model of type 2 with adjusted profile likelihood approximation to the normal \n
                                   or lognormal, par must be vector of length 6")
                            }
                          par.fem     <- par[1:3];
                          par.mal     <- par[4:6];
                          par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]));
                         }
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          if(length(par) != 8)
                            {
                             stop("For a two-sexes Schnute model of type 2 with exact normal, lognormal or gamma likelihood, \n
                                   par must be vector of length 8")
                            }
                          par.fem     <- par[c(1:3,7)];
                          par.mal     <- par[c(4:6,8)];
                          par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]), 0.5*(par[7]+par[7]));
                         }
                       age.fem     <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2];
                       obslen.fem  <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3];
                       results.fem <- do.call(optimx, list(par=par.fem,
                                                           fn=.schnute2,
                                                           gr=NULL,
                                                           age=age.fem,
                                                           obslen=obslen.fem,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.mal     <- c(x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2],unsex.mal.age);
                       obslen.mal  <- c(x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3],unsex.mal.len);
                       results.mal <- do.call(optimx, list(par=par.mal,
                                                           fn=.schnute2,
                                                           gr=NULL,
                                                           age=age.mal,
                                                           obslen=obslen.mal,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.pld     <- c(age.fem,age.mal);
                       obslen.pld  <- c(obslen.fem,obslen.mal);
                       results.pld <- do.call(optimx, list(par=par.pld,
                                                           fn=.schnute2,
                                                           gr=NULL,
                                                           age=age.pld,
                                                           obslen=obslen.pld,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                       results2                <- vector("list",3);
                       names(results2)         <- c("Females","Males","Pooled");
                       #Females
                       results2$Females        <- vector("list",length(method));
                       names(results2$Females) <- method;
                       temp.fem                <- attr(results.fem, "details")
                       for(i in 1:length(method))
                         {
                          results2$Females[[i]]$Type         <- m;
                          results2$Females[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","gamma","g1");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Females[[i]]$converg          <- "FAIL";
                          results2$Females[[i]]$kkt              <- NA;
                          results2$Females[[i]]$AIC              <- NA;
                          results2$Females[[i]]$bt.par           <- NA;
                          results2$Females[[i]]$num.grads        <- NA;
                          results2$Females[[i]]$bt.stdev         <- NA;
                          results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                            {
                             results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                             results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                             results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                             results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                             results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Females[[i]]$Cor              <- cor(v);
                             names(results2$Females[[i]]$num.grads) <- par.names;
                             names(results2$Females[[i]]$bt.par)    <- par.names;
                             names(results2$Females[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Females[[i]]$Cor)    <- par.names;
                             rownames(results2$Females[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.fem,par.fem);
                       #Males
                       results2$Males        <- vector("list",length(method));
                       names(results2$Males)   <- method;
                       temp.mal                <- attr(results.mal, "details")
                       for(i in 1:length(method))
                         {
                          results2$Males[[i]]$Type         <- m;
                          results2$Males[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","gamma","g1");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Males[[i]]$converg          <- "FAIL";
                          results2$Males[[i]]$kkt              <- NA;
                          results2$Males[[i]]$AIC              <- NA;
                          results2$Males[[i]]$bt.par           <- NA;
                          results2$Males[[i]]$num.grads        <- NA;
                          results2$Males[[i]]$bt.stdev         <- NA;
                          results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                            {
                             results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                             results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                             results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                             results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                             results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Males[[i]]$Cor              <- cor(v);
                             names(results2$Males[[i]]$num.grads) <- par.names;
                             names(results2$Males[[i]]$bt.par)    <- par.names;
                             names(results2$Males[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Males[[i]]$Cor)    <- par.names;
                             rownames(results2$Males[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.mal,par.mal);
                       #Pooled
                       results2$Pooled        <- vector("list",length(method));
                       names(results2$Pooled)   <- method;
                       temp.pld                <- attr(results.pld, "details")
                       for(i in 1:length(method))
                         {
                          results2$Pooled[[i]]$Type         <- m;
                          results2$Pooled[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","gamma","g1");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Pooled[[i]]$converg          <- "FAIL";
                          results2$Pooled[[i]]$kkt              <- NA;
                          results2$Pooled[[i]]$AIC              <- NA;
                          results2$Pooled[[i]]$bt.par           <- NA;
                          results2$Pooled[[i]]$num.grads        <- NA;
                          results2$Pooled[[i]]$bt.stdev         <- NA;
                          results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                            {
                             results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                             results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                             results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                             results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                             results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Pooled[[i]]$Cor              <- cor(v);
                             names(results2$Pooled[[i]]$num.grads) <- par.names;
                             names(results2$Pooled[[i]]$bt.par)    <- par.names;
                             names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                             rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.pld,par.pld);
                      }
                   else if(m == "schnute3")
                      {
                       if(distr == "apnormal" | distr == "aplnormal")
                         {
                          if(length(par) != 4)
                            {
                             stop("For a two-sexes Schnute model of type 3 with adjusted profile likelihood approximation to the normal \n
                                   or lognormal, par must be vector of length 4")
                            }
                          par.fem     <- par[1:2];
                          par.mal     <- par[3:4];
                          par.pld     <- c(0.5*(par[1]+par[3]), 0.5*(par[2]+par[4]));
                         }
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          if(length(par) != 6)
                            {
                             stop("For a two-sexes Schnute model of type 3 with exact normal, lognormal or gamma likelihood, \n
                                   par must be vector of length 6")
                            }
                          par.fem     <- par[c(1:2,5)];
                          par.mal     <- par[c(3:4,6)];
                          par.pld     <- c(0.5*(par[1]+par[3]), 0.5*(par[2]+par[4]), 0.5*(par[5]+par[6]));
                         }
                       age.fem     <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2];
                       obslen.fem  <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3];
                       results.fem <- do.call(optimx, list(par=par.fem,
                                                           fn=.schnute3,
                                                           gr=NULL,
                                                           age=age.fem,
                                                           obslen=obslen.fem,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.mal     <- c(x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2],unsex.mal.age);
                       obslen.mal  <- c(x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3],unsex.mal.len);
                       results.mal <- do.call(optimx, list(par=par.mal,
                                                           fn=.schnute3,
                                                           gr=NULL,
                                                           age=age.mal,
                                                           obslen=obslen.mal,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.pld     <- c(age.fem,age.mal);
                       obslen.pld  <- c(obslen.fem,obslen.mal);
                       results.pld <- do.call(optimx, list(par=par.pld,
                                                           fn=.schnute3,
                                                           gr=NULL,
                                                           age=age.pld,
                                                           obslen=obslen.pld,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                       results2                <- vector("list",3);
                       names(results2)         <- c("Females","Males","Pooled");
                       #Females
                       results2$Females        <- vector("list",length(method));
                       names(results2$Females) <- method;
                       temp.fem                <- attr(results.fem, "details")
                       for(i in 1:length(method))
                         {
                          results2$Females[[i]]$Type         <- m;
                          results2$Females[[i]]$Distr        <- distr;
                          par.names                  <- c("g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Females[[i]]$converg          <- "FAIL";
                          results2$Females[[i]]$kkt              <- NA;
                          results2$Females[[i]]$AIC              <- NA;
                          results2$Females[[i]]$bt.par           <- NA;
                          results2$Females[[i]]$num.grads        <- NA;
                          results2$Females[[i]]$bt.stdev         <- NA;
                          results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                            {
                             results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                             results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                             results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                             results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                             results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Females[[i]]$Cor              <- cor(v);
                             names(results2$Females[[i]]$num.grads) <- par.names;
                             names(results2$Females[[i]]$bt.par)    <- par.names;
                             names(results2$Females[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Females[[i]]$Cor)    <- par.names;
                             rownames(results2$Females[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.fem,par.fem);
                       #Males
                       results2$Males        <- vector("list",length(method));
                       names(results2$Males)   <- method;
                       temp.mal                <- attr(results.mal, "details")
                       for(i in 1:length(method))
                         {
                          results2$Males[[i]]$Type         <- m;
                          results2$Males[[i]]$Distr        <- distr;
                          par.names                  <- c("g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Males[[i]]$converg          <- "FAIL";
                          results2$Males[[i]]$kkt              <- NA;
                          results2$Males[[i]]$AIC              <- NA;
                          results2$Males[[i]]$bt.par           <- NA;
                          results2$Males[[i]]$num.grads        <- NA;
                          results2$Males[[i]]$bt.stdev         <- NA;
                          results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                            {
                             results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                             results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                             results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                             results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                             results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Males[[i]]$Cor              <- cor(v);
                             names(results2$Males[[i]]$num.grads) <- par.names;
                             names(results2$Males[[i]]$bt.par)    <- par.names;
                             names(results2$Males[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Males[[i]]$Cor)    <- par.names;
                             rownames(results2$Males[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.mal,par.mal);
                       #Pooled
                       results2$Pooled        <- vector("list",length(method));
                       names(results2$Pooled)   <- method;
                       temp.pld                <- attr(results.pld, "details")
                       for(i in 1:length(method))
                         {
                          results2$Pooled[[i]]$Type         <- m;
                          results2$Pooled[[i]]$Distr        <- distr;
                          par.names                  <- c("g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Pooled[[i]]$converg          <- "FAIL";
                          results2$Pooled[[i]]$kkt              <- NA;
                          results2$Pooled[[i]]$AIC              <- NA;
                          results2$Pooled[[i]]$bt.par           <- NA;
                          results2$Pooled[[i]]$num.grads        <- NA;
                          results2$Pooled[[i]]$bt.stdev         <- NA;
                          results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                            {
                             results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                             results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                             results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                             results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                             results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Pooled[[i]]$Cor              <- cor(v);
                             names(results2$Pooled[[i]]$num.grads) <- par.names;
                             names(results2$Pooled[[i]]$bt.par)    <- par.names;
                             names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                             rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.pld,par.pld);
                      }
                  }
                if(unsex.action == "females")
                  {
                   unsex.fem.age <- x$Data$Unsexed[which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length)),2];
                   unsex.fem.len <- x$Data$Unsexed[which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length)),3];
                   if(m == "vonber")
                      {
                       if(distr == "apnormal" | distr == "aplnormal")
                         {
                          if(length(par) != 6)
                            {
                             stop("For a two-sexes Von Bertalanffy model with adjusted profile likelihood approximation to the normal \n
                                   or lognormal, par must be vector of length 6")
                            }
                          par.fem     <- par[1:3];
                          par.mal     <- par[4:6];
                          par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]));
                         }
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          if(length(par) != 8)
                            {
                             stop("For a two-sexes Von Bertalanffy model with exact normal, lognormal or gamma likelihood, \n
                                   par must be vector of length 8")
                            }
                          par.fem     <- par[c(1:3,7)];
                          par.mal     <- par[c(4:6,8)];
                          par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]), 0.5*(par[7]+par[8]));
                         }
                       age.fem     <- c(x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2],unsex.fem.age);
                       obslen.fem  <- c(x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3],unsex.fem.len);
                       results.fem <- do.call(optimx, list(par=par.fem,
                                                           fn=.vonBer,
                                                           gr=NULL,
                                                           age=age.fem,
                                                           obslen=obslen.fem,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.mal     <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2];
                       obslen.mal  <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3];
                       results.mal <- do.call(optimx, list(par=par.mal,
                                                           fn=.vonBer,
                                                           gr=NULL,
                                                           age=age.mal,
                                                           obslen=obslen.mal,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.pld     <- c(age.fem,age.mal);
                       obslen.pld  <- c(obslen.fem,obslen.mal);
                       results.pld <- do.call(optimx, list(par=par.pld,
                                                           fn=.vonBer,
                                                           gr=NULL,
                                                           age=age.pld,
                                                           obslen=obslen.pld,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                       results2                <- vector("list",3);
                       names(results2)         <- c("Females","Males","Pooled");
                       #Females
                       results2$Females        <- vector("list",length(method));
                       names(results2$Females) <- method;
                       temp.fem                <- attr(results.fem, "details")
                       for(i in 1:length(method))
                         {
                          results2$Females[[i]]$Type         <- m;
                          results2$Females[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","L0","K");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Females[[i]]$converg          <- "FAIL";
                          results2$Females[[i]]$kkt              <- NA;
                          results2$Females[[i]]$AIC              <- NA;
                          results2$Females[[i]]$bt.par           <- NA;
                          results2$Females[[i]]$num.grads        <- NA;
                          results2$Females[[i]]$bt.stdev         <- NA;
                          results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                            {
                             results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                             results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                             results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                             results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                             results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Females[[i]]$Cor              <- cor(v);
                             names(results2$Females[[i]]$num.grads) <- par.names;
                             names(results2$Females[[i]]$bt.par)    <- par.names;
                             names(results2$Females[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Females[[i]]$Cor)    <- par.names;
                             rownames(results2$Females[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.fem,par.fem);
                       #Males
                       results2$Males        <- vector("list",length(method));
                       names(results2$Males)   <- method;
                       temp.mal                <- attr(results.mal, "details")
                       for(i in 1:length(method))
                         {
                          results2$Males[[i]]$Type         <- m;
                          results2$Males[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","L0","K");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Males[[i]]$converg          <- "FAIL";
                          results2$Males[[i]]$kkt              <- NA;
                          results2$Males[[i]]$AIC              <- NA;
                          results2$Males[[i]]$bt.par           <- NA;
                          results2$Males[[i]]$num.grads        <- NA;
                          results2$Males[[i]]$bt.stdev         <- NA;
                          results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                            {
                             results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                             results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                             results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                             results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                             results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Males[[i]]$Cor              <- cor(v);
                             names(results2$Males[[i]]$num.grads) <- par.names;
                             names(results2$Males[[i]]$bt.par)    <- par.names;
                             names(results2$Males[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Males[[i]]$Cor)    <- par.names;
                             rownames(results2$Males[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.mal,par.mal);
                       #Pooled
                       results2$Pooled        <- vector("list",length(method));
                       names(results2$Pooled)   <- method;
                       temp.pld                <- attr(results.pld, "details")
                       for(i in 1:length(method))
                         {
                          results2$Pooled[[i]]$Type         <- m;
                          results2$Pooled[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","L0","K");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Pooled[[i]]$converg          <- "FAIL";
                          results2$Pooled[[i]]$kkt              <- NA;
                          results2$Pooled[[i]]$AIC              <- NA;
                          results2$Pooled[[i]]$bt.par           <- NA;
                          results2$Pooled[[i]]$num.grads        <- NA;
                          results2$Pooled[[i]]$bt.stdev         <- NA;
                          results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                            {
                             results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                             results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                             results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                             results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                             results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Pooled[[i]]$Cor              <- cor(v);
                             names(results2$Pooled[[i]]$num.grads) <- par.names;
                             names(results2$Pooled[[i]]$bt.par)    <- par.names;
                             names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                             rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.pld,par.pld);
                      }
                   else if(m == "schnute1")
                      {
                       if(distr == "apnormal" | distr == "aplnormal")
                         {
                          if(length(par) != 10)
                            {
                             stop("For a two-sexes Schnute model of type 1 with adjusted profile likelihood approximation to the normal \n
                             or lognormal, par must be vector of length 10")
                            }
                          par.fem     <- par[1:5];
                          par.mal     <- par[6:10];
                          par.pld     <- c(0.5*(par[1]+par[6]), 0.5*(par[2]+par[7]), 0.5*(par[3]+par[8]), 0.5*(par[4]+par[9]), 0.5*(par[5]+par[10]));
                         }
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          if(length(par) != 12)
                            {
                             stop("For a two-sexes Schnute model of type 1 with exact normal, lognormal or gamma likelihood, \n
                             par must be vector of length 12")
                            }
                          par.fem     <- par[c(1:5,11)];
                          par.mal     <- par[c(6:10,12)];
                          par.pld     <- c(0.5*(par[1]+par[6]), 0.5*(par[2]+par[7]), 0.5*(par[3]+par[8]), 0.5*(par[4]+par[9]), 0.5*(par[5]+par[10]), 0.5*(par[11]+par[12]));
                         }
                       age.fem     <- c(x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2],unsex.fem.age);
                       obslen.fem  <- c(x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3],unsex.fem.len);
                       results.fem <- do.call(optimx, list(par=par.fem,
                                                           fn=.schnute1,
                                                           gr=NULL,
                                                           age=age.fem,
                                                           obslen=obslen.fem,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.mal     <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2];
                       obslen.mal  <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3];
                       results.mal <- do.call(optimx, list(par=par.mal,
                                                           fn=.schnute1,
                                                           gr=NULL,
                                                           age=age.mal,
                                                           obslen=obslen.mal,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.pld     <- c(age.fem,age.mal);
                       obslen.pld  <- c(obslen.fem,obslen.mal);
                       results.pld <- do.call(optimx, list(par=par.pld,
                                                           fn=.schnute1,
                                                           gr=NULL,
                                                           age=age.pld,
                                                           obslen=obslen.pld,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                       results2                <- vector("list",3);
                       names(results2)         <- c("Females","Males","Pooled");
                       #Females
                       results2$Females        <- vector("list",length(method));
                       names(results2$Females) <- method;
                       temp.fem                <- attr(results.fem, "details")
                       for(i in 1:length(method))
                         {
                          results2$Females[[i]]$Type         <- m;
                          results2$Females[[i]]$Distr        <- distr;
                          par.names                  <- c("l1","a1","mu","g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Females[[i]]$converg          <- "FAIL";
                          results2$Females[[i]]$kkt              <- NA;
                          results2$Females[[i]]$AIC              <- NA;
                          results2$Females[[i]]$bt.par           <- NA;
                          results2$Females[[i]]$num.grads        <- NA;
                          results2$Females[[i]]$bt.stdev         <- NA;
                          results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                            {
                             results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                             results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                             results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                             results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                             results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Females[[i]]$Cor              <- cor(v);
                             names(results2$Females[[i]]$num.grads) <- par.names;
                             names(results2$Females[[i]]$bt.par)    <- par.names;
                             names(results2$Females[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Females[[i]]$Cor)    <- par.names;
                             rownames(results2$Females[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.fem,par.fem);
                       #Males
                       results2$Males        <- vector("list",length(method));
                       names(results2$Males)   <- method;
                       temp.mal                <- attr(results.mal, "details")
                       for(i in 1:length(method))
                         {
                          results2$Males[[i]]$Type         <- m;
                          results2$Males[[i]]$Distr        <- distr;
                          par.names                  <- c("l1","a1","mu","g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Males[[i]]$converg          <- "FAIL";
                          results2$Males[[i]]$kkt              <- NA;
                          results2$Males[[i]]$AIC              <- NA;
                          results2$Males[[i]]$bt.par           <- NA;
                          results2$Males[[i]]$num.grads        <- NA;
                          results2$Males[[i]]$bt.stdev         <- NA;
                          results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                            {
                             results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                             results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                             results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                             results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                             results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Males[[i]]$Cor              <- cor(v);
                             names(results2$Males[[i]]$num.grads) <- par.names;
                             names(results2$Males[[i]]$bt.par)    <- par.names;
                             names(results2$Males[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Males[[i]]$Cor)    <- par.names;
                             rownames(results2$Males[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.mal,par.mal);
                       #Pooled
                       results2$Pooled        <- vector("list",length(method));
                       names(results2$Pooled)   <- method;
                       temp.pld                <- attr(results.pld, "details")
                       for(i in 1:length(method))
                         {
                          results2$Pooled[[i]]$Type         <- m;
                          results2$Pooled[[i]]$Distr        <- distr;
                          par.names                  <- c("l1","a1","mu","g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Pooled[[i]]$converg          <- "FAIL";
                          results2$Pooled[[i]]$kkt              <- NA;
                          results2$Pooled[[i]]$AIC              <- NA;
                          results2$Pooled[[i]]$bt.par           <- NA;
                          results2$Pooled[[i]]$num.grads        <- NA;
                          results2$Pooled[[i]]$bt.stdev         <- NA;
                          results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                            {
                             results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                             results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                             results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                             results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                             results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Pooled[[i]]$Cor              <- cor(v);
                             names(results2$Pooled[[i]]$num.grads) <- par.names;
                             names(results2$Pooled[[i]]$bt.par)    <- par.names;
                             names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                             rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.pld,par.pld);
                      }
                   else if(m == "schnute2")
                      {
                       if(distr == "apnormal" | distr == "aplnormal")
                         {
                          if(length(par) != 6)
                            {
                             stop("For a two-sexes Schnute model of type 2 with adjusted profile likelihood approximation to the normal \n
                                   or lognormal, par must be vector of length 6")
                            }
                          par.fem     <- par[1:3];
                          par.mal     <- par[4:6];
                          par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]));
                         }
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          if(length(par) != 8)
                            {
                             stop("For a two-sexes Schnute model of type 2 with exact normal, lognormal or gamma likelihood, \n
                                   par must be vector of length 8")
                            }
                          par.fem     <- par[c(1:3,7)];
                          par.mal     <- par[c(4:6,8)];
                          par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]), 0.5*(par[7]+par[7]));
                         }
                       age.fem     <- c(x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2],unsex.fem.age);
                       obslen.fem  <- c(x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3],unsex.fem.len);
                       results.fem <- do.call(optimx, list(par=par.fem,
                                                           fn=.schnute2,
                                                           gr=NULL,
                                                           age=age.fem,
                                                           obslen=obslen.fem,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.mal     <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2];
                       obslen.mal  <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3];
                       results.mal <- do.call(optimx, list(par=par.mal,
                                                           fn=.schnute2,
                                                           gr=NULL,
                                                           age=age.mal,
                                                           obslen=obslen.mal,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.pld     <- c(age.fem,age.mal);
                       obslen.pld  <- c(obslen.fem,obslen.mal);
                       results.pld <- do.call(optimx, list(par=par.pld,
                                                           fn=.schnute2,
                                                           gr=NULL,
                                                           age=age.pld,
                                                           obslen=obslen.pld,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                       results2                <- vector("list",3);
                       names(results2)         <- c("Females","Males","Pooled");
                       #Females
                       results2$Females        <- vector("list",length(method));
                       names(results2$Females) <- method;
                       temp.fem                <- attr(results.fem, "details")
                       for(i in 1:length(method))
                         {
                          results2$Females[[i]]$Type         <- m;
                          results2$Females[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","gamma","g1");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Females[[i]]$converg          <- "FAIL";
                          results2$Females[[i]]$kkt              <- NA;
                          results2$Females[[i]]$AIC              <- NA;
                          results2$Females[[i]]$bt.par           <- NA;
                          results2$Females[[i]]$num.grads        <- NA;
                          results2$Females[[i]]$bt.stdev         <- NA;
                          results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                            {
                             results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                             results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                             results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                             results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                             results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Females[[i]]$Cor              <- cor(v);
                             names(results2$Females[[i]]$num.grads) <- par.names;
                             names(results2$Females[[i]]$bt.par)    <- par.names;
                             names(results2$Females[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Females[[i]]$Cor)    <- par.names;
                             rownames(results2$Females[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.fem,par.fem);
                       #Males
                       results2$Males        <- vector("list",length(method));
                       names(results2$Males)   <- method;
                       temp.mal                <- attr(results.mal, "details")
                       for(i in 1:length(method))
                         {
                          results2$Males[[i]]$Type         <- m;
                          results2$Males[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","gamma","g1");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Males[[i]]$converg          <- "FAIL";
                          results2$Males[[i]]$kkt              <- NA;
                          results2$Males[[i]]$AIC              <- NA;
                          results2$Males[[i]]$bt.par           <- NA;
                          results2$Males[[i]]$num.grads        <- NA;
                          results2$Males[[i]]$bt.stdev         <- NA;
                          results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                            {
                             results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                             results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                             results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                             results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                             results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Males[[i]]$Cor              <- cor(v);
                             names(results2$Males[[i]]$num.grads) <- par.names;
                             names(results2$Males[[i]]$bt.par)    <- par.names;
                             names(results2$Males[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Males[[i]]$Cor)    <- par.names;
                             rownames(results2$Males[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.mal,par.mal);
                       #Pooled
                       results2$Pooled        <- vector("list",length(method));
                       names(results2$Pooled)   <- method;
                       temp.pld                <- attr(results.pld, "details")
                       for(i in 1:length(method))
                         {
                          results2$Pooled[[i]]$Type         <- m;
                          results2$Pooled[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","gamma","g1");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Pooled[[i]]$converg          <- "FAIL";
                          results2$Pooled[[i]]$kkt              <- NA;
                          results2$Pooled[[i]]$AIC              <- NA;
                          results2$Pooled[[i]]$bt.par           <- NA;
                          results2$Pooled[[i]]$num.grads        <- NA;
                          results2$Pooled[[i]]$bt.stdev         <- NA;
                          results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                            {
                             results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                             results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                             results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                             results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                             results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Pooled[[i]]$Cor              <- cor(v);
                             names(results2$Pooled[[i]]$num.grads) <- par.names;
                             names(results2$Pooled[[i]]$bt.par)    <- par.names;
                             names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                             rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.pld,par.pld);
                      }
                   else if(m == "schnute3")
                      {
                       if(distr == "apnormal" | distr == "aplnormal")
                         {
                          if(length(par) != 4)
                            {
                             stop("For a two-sexes Schnute model of type 3 with adjusted profile likelihood approximation to the normal \n
                                   or lognormal, par must be vector of length 4")
                            }
                          par.fem     <- par[1:2];
                          par.mal     <- par[3:4];
                          par.pld     <- c(0.5*(par[1]+par[3]), 0.5*(par[2]+par[4]));
                         }
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          if(length(par) != 6)
                            {
                             stop("For a two-sexes Schnute model of type 3 with exact normal, lognormal or gamma likelihood, \n
                                   par must be vector of length 6")
                            }
                          par.fem     <- par[c(1:2,5)];
                          par.mal     <- par[c(3:4,6)];
                          par.pld     <- c(0.5*(par[1]+par[3]), 0.5*(par[2]+par[4]), 0.5*(par[5]+par[6]));
                         }
                       age.fem     <- c(x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2],unsex.fem.age);
                       obslen.fem  <- c(x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3],unsex.fem.len);
                       results.fem <- do.call(optimx, list(par=par.fem,
                                                           fn=.schnute3,
                                                           gr=NULL,
                                                           age=age.fem,
                                                           obslen=obslen.fem,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.mal     <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2];
                       obslen.mal  <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3];
                       results.mal <- do.call(optimx, list(par=par.mal,
                                                           fn=.schnute3,
                                                           gr=NULL,
                                                           age=age.mal,
                                                           obslen=obslen.mal,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.pld     <- c(age.fem,age.mal);
                       obslen.pld  <- c(obslen.fem,obslen.mal);
                       results.pld <- do.call(optimx, list(par=par.pld,
                                                           fn=.schnute3,
                                                           gr=NULL,
                                                           age=age.pld,
                                                           obslen=obslen.pld,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                       results2                <- vector("list",3);
                       names(results2)         <- c("Females","Males","Pooled");
                       #Females
                       results2$Females        <- vector("list",length(method));
                       names(results2$Females) <- method;
                       temp.fem                <- attr(results.fem, "details")
                       for(i in 1:length(method))
                         {
                          results2$Females[[i]]$Type         <- m;
                          results2$Females[[i]]$Distr        <- distr;
                          par.names                  <- c("g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Females[[i]]$converg          <- "FAIL";
                          results2$Females[[i]]$kkt              <- NA;
                          results2$Females[[i]]$AIC              <- NA;
                          results2$Females[[i]]$bt.par           <- NA;
                          results2$Females[[i]]$num.grads        <- NA;
                          results2$Females[[i]]$bt.stdev         <- NA;
                          results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                            {
                             results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                             results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                             results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                             results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                             results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Females[[i]]$Cor              <- cor(v);
                             names(results2$Females[[i]]$num.grads) <- par.names;
                             names(results2$Females[[i]]$bt.par)    <- par.names;
                             names(results2$Females[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Females[[i]]$Cor)    <- par.names;
                             rownames(results2$Females[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.fem,par.fem);
                       #Males
                       results2$Males        <- vector("list",length(method));
                       names(results2$Males)   <- method;
                       temp.mal                <- attr(results.mal, "details")
                       for(i in 1:length(method))
                         {
                          results2$Males[[i]]$Type         <- m;
                          results2$Males[[i]]$Distr        <- distr;
                          par.names                  <- c("g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Males[[i]]$converg          <- "FAIL";
                          results2$Males[[i]]$kkt              <- NA;
                          results2$Males[[i]]$AIC              <- NA;
                          results2$Males[[i]]$bt.par           <- NA;
                          results2$Males[[i]]$num.grads        <- NA;
                          results2$Males[[i]]$bt.stdev         <- NA;
                          results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                            {
                             results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                             results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                             results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                             results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                             results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Males[[i]]$Cor              <- cor(v);
                             names(results2$Males[[i]]$num.grads) <- par.names;
                             names(results2$Males[[i]]$bt.par)    <- par.names;
                             names(results2$Males[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Males[[i]]$Cor)    <- par.names;
                             rownames(results2$Males[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.mal,par.mal);
                       #Pooled
                       results2$Pooled        <- vector("list",length(method));
                       names(results2$Pooled)   <- method;
                       temp.pld                <- attr(results.pld, "details")
                       for(i in 1:length(method))
                         {
                          results2$Pooled[[i]]$Type         <- m;
                          results2$Pooled[[i]]$Distr        <- distr;
                          par.names                  <- c("g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Pooled[[i]]$converg          <- "FAIL";
                          results2$Pooled[[i]]$kkt              <- NA;
                          results2$Pooled[[i]]$AIC              <- NA;
                          results2$Pooled[[i]]$bt.par           <- NA;
                          results2$Pooled[[i]]$num.grads        <- NA;
                          results2$Pooled[[i]]$bt.stdev         <- NA;
                          results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                            {
                             results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                             results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                             results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                             results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                             results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Pooled[[i]]$Cor              <- cor(v);
                             names(results2$Pooled[[i]]$num.grads) <- par.names;
                             names(results2$Pooled[[i]]$bt.par)    <- par.names;
                             names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                             rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.pld,par.pld);
                      }
                  }
                if(unsex.action == "ignore")
                  {
                   if(m == "vonber")
                      {
                       if(distr == "apnormal" | distr == "aplnormal")
                         {
                          if(length(par) != 6)
                            {
                             stop("For a two-sexes Von Bertalanffy model with adjusted profile likelihood approximation to the normal \n
                                   or lognormal, par must be vector of length 6")
                            }
                          par.fem     <- par[1:3];
                          par.mal     <- par[4:6];
                          par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]));
                         }
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          if(length(par) != 8)
                            {
                             stop("For a two-sexes Von Bertalanffy model with exact normal, lognormal or gamma likelihood, \n
                                   par must be vector of length 8")
                            }
                          par.fem     <- par[c(1:3,7)];
                          par.mal     <- par[c(4:6,8)];
                          par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]), 0.5*(par[7]+par[8]));
                         }
                       age.fem     <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2];
                       obslen.fem  <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3];
                       results.fem <- do.call(optimx, list(par=par.fem,
                                                           fn=.vonBer,
                                                           gr=NULL,
                                                           age=age.fem,
                                                           obslen=obslen.fem,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.mal     <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2];
                       obslen.mal  <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3];
                       results.mal <- do.call(optimx, list(par=par.mal,
                                                           fn=.vonBer,
                                                           gr=NULL,
                                                           age=age.mal,
                                                           obslen=obslen.mal,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.pld     <- c(age.fem,age.mal);
                       obslen.pld  <- c(obslen.fem,obslen.mal);
                       results.pld <- do.call(optimx, list(par=par.pld,
                                                           fn=.vonBer,
                                                           gr=NULL,
                                                           age=age.pld,
                                                           obslen=obslen.pld,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                       results2                <- vector("list",3);
                       names(results2)         <- c("Females","Males","Pooled");
                       #Females
                       results2$Females        <- vector("list",length(method));
                       names(results2$Females) <- method;
                       temp.fem                <- attr(results.fem, "details")
                       for(i in 1:length(method))
                         {
                          results2$Females[[i]]$Type         <- m;
                          results2$Females[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","L0","K");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Females[[i]]$converg          <- "FAIL";
                          results2$Females[[i]]$kkt              <- NA;
                          results2$Females[[i]]$AIC              <- NA;
                          results2$Females[[i]]$bt.par           <- NA;
                          results2$Females[[i]]$num.grads        <- NA;
                          results2$Females[[i]]$bt.stdev         <- NA;
                          results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                            {
                             results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                             results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                             results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                             results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                             results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Females[[i]]$Cor              <- cor(v);
                             names(results2$Females[[i]]$num.grads) <- par.names;
                             names(results2$Females[[i]]$bt.par)    <- par.names;
                             names(results2$Females[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Females[[i]]$Cor)    <- par.names;
                             rownames(results2$Females[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.fem,par.fem);
                       #Males
                       results2$Males        <- vector("list",length(method));
                       names(results2$Males)   <- method;
                       temp.mal                <- attr(results.mal, "details")
                       for(i in 1:length(method))
                         {
                          results2$Males[[i]]$Type         <- m;
                          results2$Males[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","L0","K");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Males[[i]]$converg          <- "FAIL";
                          results2$Males[[i]]$kkt              <- NA;
                          results2$Males[[i]]$AIC              <- NA;
                          results2$Males[[i]]$bt.par           <- NA;
                          results2$Males[[i]]$num.grads        <- NA;
                          results2$Males[[i]]$bt.stdev         <- NA;
                          results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                            {
                             results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                             results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                             results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                             results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                             results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Males[[i]]$Cor              <- cor(v);
                             names(results2$Males[[i]]$num.grads) <- par.names;
                             names(results2$Males[[i]]$bt.par)    <- par.names;
                             names(results2$Males[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Males[[i]]$Cor)    <- par.names;
                             rownames(results2$Males[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.mal,par.mal);
                       #Pooled
                       results2$Pooled        <- vector("list",length(method));
                       names(results2$Pooled)   <- method;
                       temp.pld                <- attr(results.pld, "details")
                       for(i in 1:length(method))
                         {
                          results2$Pooled[[i]]$Type         <- m;
                          results2$Pooled[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","L0","K");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Pooled[[i]]$converg          <- "FAIL";
                          results2$Pooled[[i]]$kkt              <- NA;
                          results2$Pooled[[i]]$AIC              <- NA;
                          results2$Pooled[[i]]$bt.par           <- NA;
                          results2$Pooled[[i]]$num.grads        <- NA;
                          results2$Pooled[[i]]$bt.stdev         <- NA;
                          results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                            {
                             results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                             results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                             results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                             results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                             results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Pooled[[i]]$Cor              <- cor(v);
                             names(results2$Pooled[[i]]$num.grads) <- par.names;
                             names(results2$Pooled[[i]]$bt.par)    <- par.names;
                             names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                             rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.pld,par.pld);
                      }
                   else if(m == "schnute1")
                      {
                       if(distr == "apnormal" | distr == "aplnormal")
                         {
                          if(length(par) != 10)
                            {
                             stop("For a two-sexes Schnute model of type 1 with adjusted profile likelihood approximation to the normal \n
                             or lognormal, par must be vector of length 10")
                            }
                          par.fem     <- par[1:5];
                          par.mal     <- par[6:10];
                          par.pld     <- c(0.5*(par[1]+par[6]), 0.5*(par[2]+par[7]), 0.5*(par[3]+par[8]), 0.5*(par[4]+par[9]), 0.5*(par[5]+par[10]));
                         }
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          if(length(par) != 12)
                            {
                             stop("For a two-sexes Schnute model of type 1 with exact normal, lognormal or gamma likelihood, \n
                             par must be vector of length 12")
                            }
                          par.fem     <- par[c(1:5,11)];
                          par.mal     <- par[c(6:10,12)];
                          par.pld     <- c(0.5*(par[1]+par[6]), 0.5*(par[2]+par[7]), 0.5*(par[3]+par[8]), 0.5*(par[4]+par[9]), 0.5*(par[5]+par[10]), 0.5*(par[11]+par[12]));
                         }
                       age.fem     <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2];
                       obslen.fem  <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3];
                       results.fem <- do.call(optimx, list(par=par.fem,
                                                           fn=.schnute1,
                                                           gr=NULL,
                                                           age=age.fem,
                                                           obslen=obslen.fem,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.mal     <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2];
                       obslen.mal  <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3];
                       results.mal <- do.call(optimx, list(par=par.mal,
                                                           fn=.schnute1,
                                                           gr=NULL,
                                                           age=age.mal,
                                                           obslen=obslen.mal,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.pld     <- c(age.fem,age.mal);
                       obslen.pld  <- c(obslen.fem,obslen.mal);
                       results.pld <- do.call(optimx, list(par=par.pld,
                                                           fn=.schnute1,
                                                           gr=NULL,
                                                           age=age.pld,
                                                           obslen=obslen.pld,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                       results2                <- vector("list",3);
                       names(results2)         <- c("Females","Males","Pooled");
                       #Females
                       results2$Females        <- vector("list",length(method));
                       names(results2$Females) <- method;
                       temp.fem                <- attr(results.fem, "details")
                       for(i in 1:length(method))
                         {
                          results2$Females[[i]]$Type         <- m;
                          results2$Females[[i]]$Distr        <- distr;
                          par.names                  <- c("l1","a1","mu","g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Females[[i]]$converg          <- "FAIL";
                          results2$Females[[i]]$kkt              <- NA;
                          results2$Females[[i]]$AIC              <- NA;
                          results2$Females[[i]]$bt.par           <- NA;
                          results2$Females[[i]]$num.grads        <- NA;
                          results2$Females[[i]]$bt.stdev         <- NA;
                          results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                            {
                             results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                             results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                             results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                             results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                             results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Females[[i]]$Cor              <- cor(v);
                             names(results2$Females[[i]]$num.grads) <- par.names;
                             names(results2$Females[[i]]$bt.par)    <- par.names;
                             names(results2$Females[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Females[[i]]$Cor)    <- par.names;
                             rownames(results2$Females[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.fem,par.fem);
                       #Males
                       results2$Males        <- vector("list",length(method));
                       names(results2$Males)   <- method;
                       temp.mal                <- attr(results.mal, "details")
                       for(i in 1:length(method))
                         {
                          results2$Males[[i]]$Type         <- m;
                          results2$Males[[i]]$Distr        <- distr;
                          par.names                  <- c("l1","a1","mu","g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Males[[i]]$converg          <- "FAIL";
                          results2$Males[[i]]$kkt              <- NA;
                          results2$Males[[i]]$AIC              <- NA;
                          results2$Males[[i]]$bt.par           <- NA;
                          results2$Males[[i]]$num.grads        <- NA;
                          results2$Males[[i]]$bt.stdev         <- NA;
                          results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                            {
                             results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                             results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                             results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                             results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                             results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Males[[i]]$Cor              <- cor(v);
                             names(results2$Males[[i]]$num.grads) <- par.names;
                             names(results2$Males[[i]]$bt.par)    <- par.names;
                             names(results2$Males[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Males[[i]]$Cor)    <- par.names;
                             rownames(results2$Males[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.mal,par.mal);
                       #Pooled
                       results2$Pooled        <- vector("list",length(method));
                       names(results2$Pooled)   <- method;
                       temp.pld                <- attr(results.pld, "details")
                       for(i in 1:length(method))
                         {
                          results2$Pooled[[i]]$Type         <- m;
                          results2$Pooled[[i]]$Distr        <- distr;
                          par.names                  <- c("l1","a1","mu","g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Pooled[[i]]$converg          <- "FAIL";
                          results2$Pooled[[i]]$kkt              <- NA;
                          results2$Pooled[[i]]$AIC              <- NA;
                          results2$Pooled[[i]]$bt.par           <- NA;
                          results2$Pooled[[i]]$num.grads        <- NA;
                          results2$Pooled[[i]]$bt.stdev         <- NA;
                          results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                            {
                             results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                             results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                             results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                             results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                             results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4),~exp(x5),~exp(x6)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Pooled[[i]]$Cor              <- cor(v);
                             names(results2$Pooled[[i]]$num.grads) <- par.names;
                             names(results2$Pooled[[i]]$bt.par)    <- par.names;
                             names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                             rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.pld,par.pld);
                      }
                   else if(m == "schnute2")
                      {
                       if(distr == "apnormal" | distr == "aplnormal")
                         {
                          if(length(par) != 6)
                            {
                             stop("For a two-sexes Schnute model of type 2 with adjusted profile likelihood approximation to the normal \n
                                   or lognormal, par must be vector of length 6")
                            }
                          par.fem     <- par[1:3];
                          par.mal     <- par[4:6];
                          par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]));
                         }
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          if(length(par) != 8)
                            {
                             stop("For a two-sexes Schnute model of type 2 with exact normal, lognormal or gamma likelihood, \n
                                   par must be vector of length 8")
                            }
                          par.fem     <- par[c(1:3,7)];
                          par.mal     <- par[c(4:6,8)];
                          par.pld     <- c(0.5*(par[1]+par[4]), 0.5*(par[2]+par[5]), 0.5*(par[3]+par[6]), 0.5*(par[7]+par[7]));
                         }
                       age.fem     <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2];
                       obslen.fem  <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3];
                       results.fem <- do.call(optimx, list(par=par.fem,
                                                           fn=.schnute2,
                                                           gr=NULL,
                                                           age=age.fem,
                                                           obslen=obslen.fem,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.mal     <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2];
                       obslen.mal  <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3];
                       results.mal <- do.call(optimx, list(par=par.mal,
                                                           fn=.schnute2,
                                                           gr=NULL,
                                                           age=age.mal,
                                                           obslen=obslen.mal,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.pld     <- c(age.fem,age.mal);
                       obslen.pld  <- c(obslen.fem,obslen.mal);
                       results.pld <- do.call(optimx, list(par=par.pld,
                                                           fn=.schnute2,
                                                           gr=NULL,
                                                           age=age.pld,
                                                           obslen=obslen.pld,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                       results2                <- vector("list",3);
                       names(results2)         <- c("Females","Males","Pooled");
                       #Females
                       results2$Females        <- vector("list",length(method));
                       names(results2$Females) <- method;
                       temp.fem                <- attr(results.fem, "details")
                       for(i in 1:length(method))
                         {
                          results2$Females[[i]]$Type         <- m;
                          results2$Females[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","gamma","g1");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Females[[i]]$converg          <- "FAIL";
                          results2$Females[[i]]$kkt              <- NA;
                          results2$Females[[i]]$AIC              <- NA;
                          results2$Females[[i]]$bt.par           <- NA;
                          results2$Females[[i]]$num.grads        <- NA;
                          results2$Females[[i]]$bt.stdev         <- NA;
                          results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                            {
                             results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                             results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                             results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                             results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                             results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Females[[i]]$Cor              <- cor(v);
                             names(results2$Females[[i]]$num.grads) <- par.names;
                             names(results2$Females[[i]]$bt.par)    <- par.names;
                             names(results2$Females[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Females[[i]]$Cor)    <- par.names;
                             rownames(results2$Females[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.fem,par.fem);
                       #Males
                       results2$Males        <- vector("list",length(method));
                       names(results2$Males)   <- method;
                       temp.mal                <- attr(results.mal, "details")
                       for(i in 1:length(method))
                         {
                          results2$Males[[i]]$Type         <- m;
                          results2$Males[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","gamma","g1");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Males[[i]]$converg          <- "FAIL";
                          results2$Males[[i]]$kkt              <- NA;
                          results2$Males[[i]]$AIC              <- NA;
                          results2$Males[[i]]$bt.par           <- NA;
                          results2$Males[[i]]$num.grads        <- NA;
                          results2$Males[[i]]$bt.stdev         <- NA;
                          results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                            {
                             results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                             results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                             results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                             results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                             results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Males[[i]]$Cor              <- cor(v);
                             names(results2$Males[[i]]$num.grads) <- par.names;
                             names(results2$Males[[i]]$bt.par)    <- par.names;
                             names(results2$Males[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Males[[i]]$Cor)    <- par.names;
                             rownames(results2$Males[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.mal,par.mal);
                       #Pooled
                       results2$Pooled        <- vector("list",length(method));
                       names(results2$Pooled)   <- method;
                       temp.pld                <- attr(results.pld, "details")
                       for(i in 1:length(method))
                         {
                          results2$Pooled[[i]]$Type         <- m;
                          results2$Pooled[[i]]$Distr        <- distr;
                          par.names                  <- c("Linf","gamma","g1");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Pooled[[i]]$converg          <- "FAIL";
                          results2$Pooled[[i]]$kkt              <- NA;
                          results2$Pooled[[i]]$AIC              <- NA;
                          results2$Pooled[[i]]$bt.par           <- NA;
                          results2$Pooled[[i]]$num.grads        <- NA;
                          results2$Pooled[[i]]$bt.stdev         <- NA;
                          results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                            {
                             results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                             results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                             results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                             results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                             results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3),~exp(x4)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Pooled[[i]]$Cor              <- cor(v);
                             names(results2$Pooled[[i]]$num.grads) <- par.names;
                             names(results2$Pooled[[i]]$bt.par)    <- par.names;
                             names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                             rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.pld,par.pld);
                      }
                   else if(m == "schnute3")
                      {
                       if(distr == "apnormal" | distr == "aplnormal")
                         {
                          if(length(par) != 4)
                            {
                             stop("For a two-sexes Schnute model of type 3 with adjusted profile likelihood approximation to the normal \n
                                   or lognormal, par must be vector of length 4")
                            }
                          par.fem     <- par[1:2];
                          par.mal     <- par[3:4];
                          par.pld     <- c(0.5*(par[1]+par[3]), 0.5*(par[2]+par[4]));
                         }
                       if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                         {
                          if(length(par) != 6)
                            {
                             stop("For a two-sexes Schnute model of type 3 with exact normal, lognormal or gamma likelihood, \n
                                   par must be vector of length 6")
                            }
                          par.fem     <- par[c(1:2,5)];
                          par.mal     <- par[c(3:4,6)];
                          par.pld     <- c(0.5*(par[1]+par[3]), 0.5*(par[2]+par[4]), 0.5*(par[5]+par[6]));
                         }
                       age.fem     <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),2];
                       obslen.fem  <- x$Data$Females[which(!is.na(x$Data$Females[,2]) & !is.na(x$Data$Females[,3])),3];
                       results.fem <- do.call(optimx, list(par=par.fem,
                                                           fn=.schnute3,
                                                           gr=NULL,
                                                           age=age.fem,
                                                           obslen=obslen.fem,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.mal     <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),2];
                       obslen.mal  <- x$Data$Males[which(!is.na(x$Data$Males[,2]) & !is.na(x$Data$Males[,3])),3];
                       results.mal <- do.call(optimx, list(par=par.mal,
                                                           fn=.schnute3,
                                                           gr=NULL,
                                                           age=age.mal,
                                                           obslen=obslen.mal,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       age.pld     <- c(age.fem,age.mal);
                       obslen.pld  <- c(obslen.fem,obslen.mal);
                       results.pld <- do.call(optimx, list(par=par.pld,
                                                           fn=.schnute3,
                                                           gr=NULL,
                                                           age=age.pld,
                                                           obslen=obslen.pld,
                                                           distr=distr,
                                                           method=method,
                                                           lower = -Inf,
                                                           upper = Inf,
                                                           control=list(),
                                                           hessian=hessian,
                                                           itnmax=itnmax,
                                                           output="estimate"));
                       rm(age.fem,age.mal,age.pld,obslen.fem,obslen.mal,obslen.pld);
                       results2                <- vector("list",3);
                       names(results2)         <- c("Females","Males","Pooled");
                       #Females
                       results2$Females        <- vector("list",length(method));
                       names(results2$Females) <- method;
                       temp.fem                <- attr(results.fem, "details")
                       for(i in 1:length(method))
                         {
                          results2$Females[[i]]$Type         <- m;
                          results2$Females[[i]]$Distr        <- distr;
                          par.names                  <- c("g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Females[[i]]$converg          <- "FAIL";
                          results2$Females[[i]]$kkt              <- NA;
                          results2$Females[[i]]$AIC              <- NA;
                          results2$Females[[i]]$bt.par           <- NA;
                          results2$Females[[i]]$num.grads        <- NA;
                          results2$Females[[i]]$bt.stdev         <- NA;
                          results2$Females[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.fem[i,]$ngatend) == length(par.names) & !any(is.na(temp.fem[i,]$nhatend)) & 1/kappa(temp.fem[i,]$nhatend) > 1e-15)
                            {
                             results2$Females[[i]]$converg          <- results.fem[i,length(par.fem)+5];
                             results2$Females[[i]]$kkt              <- results.fem[i,(length(par.fem)+6):(length(par.fem)+7)];
                             results2$Females[[i]]$AIC              <- 2*length(par.fem)-2*results.fem[i,length(par.fem)+1];
                             results2$Females[[i]]$bt.par           <- exp(results.fem[i,1:length(par.fem)])
                             results2$Females[[i]]$num.grads        <- temp.fem[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.fem[i,1:length(par.fem)]),
                                                                              cov=try(solve(temp.fem[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Females[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Females[[i]]$Cor              <- cor(v);
                             names(results2$Females[[i]]$num.grads) <- par.names;
                             names(results2$Females[[i]]$bt.par)    <- par.names;
                             names(results2$Females[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Females[[i]]$Cor)    <- par.names;
                             rownames(results2$Females[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.fem,par.fem);
                       #Males
                       results2$Males        <- vector("list",length(method));
                       names(results2$Males)   <- method;
                       temp.mal                <- attr(results.mal, "details")
                       for(i in 1:length(method))
                         {
                          results2$Males[[i]]$Type         <- m;
                          results2$Males[[i]]$Distr        <- distr;
                          par.names                  <- c("g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Males[[i]]$converg          <- "FAIL";
                          results2$Males[[i]]$kkt              <- NA;
                          results2$Males[[i]]$AIC              <- NA;
                          results2$Males[[i]]$bt.par           <- NA;
                          results2$Males[[i]]$num.grads        <- NA;
                          results2$Males[[i]]$bt.stdev         <- NA;
                          results2$Males[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.mal[i,]$ngatend) == length(par.names) & !any(is.na(temp.mal[i,]$nhatend)) & 1/kappa(temp.mal[i,]$nhatend) > 1e-15)
                            {
                             results2$Males[[i]]$converg          <- results.mal[i,length(par.mal)+5];
                             results2$Males[[i]]$kkt              <- results.mal[i,(length(par.mal)+6):(length(par.mal)+7)];
                             results2$Males[[i]]$AIC              <- 2*length(par.mal)-2*results.mal[i,length(par.mal)+1];
                             results2$Males[[i]]$bt.par           <- exp(results.mal[i,1:length(par.mal)])
                             results2$Males[[i]]$num.grads        <- temp.mal[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.mal[i,1:length(par.mal)]),
                                                                              cov=try(solve(temp.mal[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Males[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Males[[i]]$Cor              <- cor(v);
                             names(results2$Males[[i]]$num.grads) <- par.names;
                             names(results2$Males[[i]]$bt.par)    <- par.names;
                             names(results2$Males[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Males[[i]]$Cor)    <- par.names;
                             rownames(results2$Males[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.mal,par.mal);
                       #Pooled
                       results2$Pooled        <- vector("list",length(method));
                       names(results2$Pooled)   <- method;
                       temp.pld                <- attr(results.pld, "details")
                       for(i in 1:length(method))
                         {
                          results2$Pooled[[i]]$Type         <- m;
                          results2$Pooled[[i]]$Distr        <- distr;
                          par.names                  <- c("g1","g2");
                          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                            {
                             par.names <- c(par.names,"Dispersion")
                            }
                          results2$Pooled[[i]]$converg          <- "FAIL";
                          results2$Pooled[[i]]$kkt              <- NA;
                          results2$Pooled[[i]]$AIC              <- NA;
                          results2$Pooled[[i]]$bt.par           <- NA;
                          results2$Pooled[[i]]$num.grads        <- NA;
                          results2$Pooled[[i]]$bt.stdev         <- NA;
                          results2$Pooled[[i]]$Cor              <- matrix(NA,length(par.names),length(par.names));
                          if(length(temp.pld[i,]$ngatend) == length(par.names) & !any(is.na(temp.pld[i,]$nhatend)) & 1/kappa(temp.pld[i,]$nhatend) > 1e-15)
                            {
                             results2$Pooled[[i]]$converg          <- results.pld[i,length(par.pld)+5];
                             results2$Pooled[[i]]$kkt              <- results.pld[i,(length(par.pld)+6):(length(par.pld)+7)];
                             results2$Pooled[[i]]$AIC              <- 2*length(par.pld)-2*results.pld[i,length(par.pld)+1];
                             results2$Pooled[[i]]$bt.par           <- exp(results.pld[i,1:length(par.pld)])
                             results2$Pooled[[i]]$num.grads        <- temp.pld[i,]$ngatend;
                                                   v        <- matrix(0,length(par.names),length(par.names));
                                                   if(distr == "apnormal" | distr == "aplnormal")
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                                                   else
                                                     {
                                                      v        <- deltamethod(g=list(~exp(x1),~exp(x2),~exp(x3)),
                                                                              mean=as.numeric(results.pld[i,1:length(par.pld)]),
                                                                              cov=try(solve(temp.pld[i,]$nhatend)),
                                                                              ses=FALSE)
                                                     }
                             results2$Pooled[[i]]$bt.stdev         <- sqrt(diag(v));
                             results2$Pooled[[i]]$Cor              <- cor(v);
                             names(results2$Pooled[[i]]$num.grads) <- par.names;
                             names(results2$Pooled[[i]]$bt.par)    <- par.names;
                             names(results2$Pooled[[i]]$bt.stdev)  <- par.names;
                             colnames(results2$Pooled[[i]]$Cor)    <- par.names;
                             rownames(results2$Pooled[[i]]$Cor)    <- par.names;
                            }
                         }
                       rm(v);
                       rm(temp.pld,par.pld);
                      }
                  }
               }
            results                     <- vector("list",4);
            names(results)              <- c("Data","Initial","Methods","Model");
            results$Data                <- x;
            results$Initial             <- exp(par);
            if(x$Properties$Sex == "Pooled" | x$Properties$Sex == "Males" | x$Properties$Sex == "Females")
              {
               if(distr == "apnormal" | distr == "aplnormal")
                 {
                  if(m == "vonber")
                    {
                     par.names <- c("linf","l0","k");
                    }
                  if(m == "schnute1")
                    {
                     par.names <- c("l1","a1","mu","g1","g2");
                    }
                  if(m == "schnute2")
                    {
                     par.names <- c("linf","gam","g1");
                    }
                  if(m == "schnute3")
                    {
                     par.names <- c("g1","g2");
                    }
                 }
               if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                 {
                  if(m == "vonber")
                    {
                     par.names <- c("linf","l0","k","Dispersion");
                    }
                  if(m == "schnute1")
                    {
                     par.names <- c("l1","a1","mu","g1","g2","Dispersion");
                    }
                  if(m == "schnute2")
                    {
                     par.names <- c("linf","gam","g1","Dispersion");
                    }
                  if(m == "schnute3")
                    {
                     par.names <- c("g1","g2","Dispersion]");
                    }
                 }
              }
            if(x$Properties$Sex == "Both" | x$Properties$Sex == "Total")
              {
               if(distr == "apnormal" | distr == "aplnormal")
                 {
                  if(m == "vonber")
                    {
                     par.names <- c("linf.fem","l0.fem","k.fem","linf.mal","l0.mal","k.mal");
                    }
                  if(m == "schnute1")
                    {
                     par.names <- c("l1.fem","a1.fem","mu.fem","g1.fem","g2.fem","l1.mal","a1.mal","mu.mal","g1.mal","g2.mal");
                    }
                  if(m == "schnute2")
                    {
                     par.names <- c("linf.fem","gam.fem","g1.fem","linf.mal","gam.mal","g1.mal");
                    }
                  if(m == "schnute3")
                    {
                     par.names <- c("g1.fem","g2.fem","g1.mal","g2.mal");
                    }
                 }
               if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                 {
                  if(m == "vonber")
                    {
                     par.names <- c("linf.fem","l0.fem","k.fem","linf.mal","l0.mal","k.mal","Dispersion.fem","Dispersion.mal");
                    }
                  if(m == "schnute1")
                    {
                     par.names <- c("l1.fem","a1.fem","mu.fem","g1.fem","g2.fem","l1.mal","a1.mal","mu.mal","g1.mal","g2.mal","Dispersion.fem","Dispersion.mal");
                    }
                  if(m == "schnute2")
                    {
                     par.names <- c("linf.fem","gam.fem","g1.fem","linf.mal","gam.mal","g1.mal","Dispersion.fem","Dispersion.mal");
                    }
                  if(m == "schnute3")
                    {
                     par.names <- c("g1.fem","g2.fem","g1.mal","g2.mal","Dispersion.fem","Dispersion.mal");
                    }
                 }
              }
            names(results$Initial)      <- par.names;
            results$Methods             <- method;
            results$Model               <- results2;
            class(results)              <- "groage";
            return(results);
  }
