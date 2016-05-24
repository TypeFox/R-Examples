groageexp <-
function(x, m, unsex.action=NULL, par, distr)
  {
     distr.set  <- c("normal","apnormal","lognormal","aplnormal","gamma")
     model.set  <- c("vonBer", "schnute1", "schnute2", "schnute3")
     if(any(is.na(par))) {stop("NAs are not allowed in the parameter vector 'par'")}
     if(length(par) < 2 || length(par) > 12) {stop("For any of the model versions the number of parameters is > 1 and < 13")}
     if(class(par) != "numeric") {stop("'par' must be a numeric vector")}
     if(!distr%in%distr.set)
       {stop("distr must be a length-1 character vector with either 'normal', 'apnormal', 'lognormal', 'aplnormal', or 'gamma' ")}
     if(!m%in%model.set)
       {stop("m must be a length-1 character vector with either 'vonBer', 'schnute1', 'schnute2', or 'schnute3' ")}
     if(x$Properties$Sex == "Males" | x$Properties$Sex == "Females" | x$Properties$Sex == "Pooled")
       {
        parlist <- list(par=par,
                        age=x$Data$Age[which(!is.na(x$Data$Age) & !is.na(x$Data$Length))],
                        obslen=x$Data$Length[which(!is.na(x$Data$Age) & !is.na(x$Data$Length))],
                        distr=distr,
                        output="predict");
        if(m=="vonBer")
          {
          if(distr == "apnormal" | distr == "aplnormal")
            {
             if(length(par) != 3)
               {
                stop("For a one sex or pooled sexes von Bertalanffy model fit with adjusted profile likelihood par must be a vector of length 3")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
            {
             if(length(par) != 4)
               {
                stop("For a one sex or pooled sexes von Bertalanffy model fit with, normal, lognormal or gamma distribution, par must be a vector of length 4")
               }
            }
          results.usx <- do.call(.vonBer, parlist);
          }
        else if(m=="schnute1")
          {
          if(distr == "apnormal" | distr == "aplnormal")
            {
             if(length(par) != 5)
               {
                stop("For a one sex or pooled sexes Schnute 1 model fit with adjusted profile likelihood par must be a vector of length 5")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
            {
             if(length(par) != 6)
               {
                stop("For a one sex or pooled sexes Schnute 1 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 6")
               }
            }
          results.usx <- do.call(.schnute1, parlist);
          }
        else if(m=="schnute2")
          {
          if(distr == "apnormal" | distr == "aplnormal")
            {
             if(length(par) != 3)
               {
                stop("For a one sex or pooled sexes Schnute 2 model fit with adjusted profile likelihood par must be a vector of length 3")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
            {
             if(length(par) != 4)
               {
                stop("For a one sex or pooled sexes Schnute 2 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 4")
               }
            }
          results.usx <- do.call(.schnute2, parlist);
          }
        else if(m=="schnute3")
          {
          if(distr == "apnormal" | distr == "aplnormal")
            {
             if(length(par) != 2)
               {
                stop("For a one sex or pooled sexes Schnute 3 model fit with adjusted profile likelihood par must be a vector of length 2")
               }
            }
          if(distr == "normal" | distr == "lognormal" | distr == "gamma")
            {
             if(length(par) != 3)
               {
                stop("For a one sex or pooled sexes Schnute 3 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 3")
               }
            }
          results.usx <- do.call(.schnute3, parlist);
          }
        results                  <- vector("list",2)
        names(results)           <- c("Properties","Results")
        results[["Properties"]]  <- x$Properties
        results$Properties$Model <- m
        results$Properties$Distr <- distr
        results[["Results"]]     <- results.usx
        class(results)           <- "GroAgeExp"
        return(results)
       }
     else if(x$Properties$Sex == "Both")
       {
        if(m=="vonBer")
          {
           if(distr == "apnormal" | distr == "aplnormal")
             {
              if(length(par) != 6)
                {
                 stop("For a two sexes von Bertalanffy model fit with adjusted profile likelihood par must be a vector of length 6")
                }
              parlist     <- list(par=par[1:3],
                                  age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                  obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                  distr=distr,
                                  output="predict");
              results.fem <- do.call(.vonBer, parlist);
              parlist     <- list(par=par[4:6],
                                  age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                  obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                  distr=distr,
                                  output="predict");
              results.mal  <- do.call(.vonBer, parlist);
             }
           else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
             {
              if(length(par) != 8)
                {
                 stop("For a two sexes von Bertalanffy model fit with, normal, lognormal or gamma distribution, par must be a vector of length 8")
                }
              parlist     <- list(par=c(par[1:3],par[7]),
                                  age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                  obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                  distr=distr,
                                  output="predict");
              results.fem <- do.call(.vonBer, parlist);
              parlist     <- list(par=c(par[4:6],par[8]),
                                  age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                  obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                  distr=distr,
                                  output="predict");
              results.mal  <- do.call(.vonBer, parlist);
             }
          }
        if(m=="schnute1")
          {
           if(distr == "apnormal" | distr == "aplnormal")
             {
              if(length(par) != 10)
                {
                 stop("For a two sexes Schnute 1 model fit with adjusted profile likelihood par must be a vector of length 10")
                }
              parlist     <- list(par=par[1:5],
                                  age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                  obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                  distr=distr,
                                  output="predict");
              results.fem <- do.call(.schnute1, parlist);
              parlist     <- list(par=par[6:10],
                                  age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                  obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                  distr=distr,
                                  output="predict");
              results.mal  <- do.call(.schnute1, parlist);
             }
           else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
             {
              if(length(par) != 12)
                {
                 stop("For a two sexes Schnute 1 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 12")
                }
              parlist     <- list(par=c(par[1:5],par[11]),
                                  age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                  obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                  distr=distr,
                                  output="predict");
              results.fem <- do.call(.schnute1, parlist);
              parlist     <- list(par=c(par[6:10],par[12]),
                                  age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                  obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                  distr=distr,
                                  output="predict");
              results.mal  <- do.call(.schnute1, parlist);
             }
          }
        if(m=="schnute2")
          {
           if(distr == "apnormal" | distr == "aplnormal")
             {
              if(length(par) != 6)
                {
                 stop("For a two sexes Schnute 1 model fit with adjusted profile likelihood par must be a vector of length 6")
                }
              parlist     <- list(par=par[1:3],
                                  age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                  obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                  distr=distr,
                                  output="predict");
              results.fem <- do.call(.schnute2, parlist);
              parlist     <- list(par=par[4:6],
                                  age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                  obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                  distr=distr,
                                  output="predict");
              results.mal  <- do.call(.schnute2, parlist);
             }
           else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
             {
              if(length(par) != 8)
                {
                 stop("For a two sexes Schnute 1 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 8")
                }
              parlist     <- list(par=c(par[1:3],par[7]),
                                  age=x$Data$Females$Age,
                                  obslen=x$Data$Females$Length,
                                  distr=distr,
                                  output="predict");
              results.fem <- do.call(.schnute2, parlist);
              parlist     <- list(par=c(par[4:6],par[8]),
                                  age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                  obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                  distr=distr,
                                  output="predict");
              results.mal  <- do.call(.schnute2, parlist);
             }
          }
        if(m=="schnute3")
          {
           if(distr == "apnormal" | distr == "aplnormal")
             {
              if(length(par) != 4)
                {
                 stop("For a two sexes Schnute 1 model fit with adjusted profile likelihood par must be a vector of length 4")
                }
              parlist     <- list(par=par[1:2],
                                  age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                  obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                  distr=distr,
                                  output="predict");
              results.fem <- do.call(.schnute3, parlist);
              parlist     <- list(par=par[3:4],
                                  age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                  obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                  distr=distr,
                                  output="predict");
              results.mal  <- do.call(.schnute3, parlist);
             }
           else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
             {
              if(length(par) != 6)
                {
                 stop("For a two sexes Schnute 1 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 6")
                }
              parlist     <- list(par=c(par[1:2],par[5]),
                                  age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                  obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                  distr=distr,
                                  output="predict");
              results.fem <- do.call(.schnute3, parlist);
              parlist     <- list(par=c(par[3:4],par[6]),
                                  age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                  obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                  distr=distr,
                                  output="predict");
              results.mal  <- do.call(.schnute3, parlist);
             }
          }
        results                  <- vector("list",3)
        names(results)           <- c("Properties","Females","Males")
        results[["Properties"]]  <- x$Properties
        results$Properties$Model <- m
        results$Properties$Distr <- distr
        results[["Females"]]     <- results.fem
        results[["Males"]]       <- results.mal
        class(results)           <- "GroAgeExp"
        return(results)
       }
     else if(x$Properties$Sex == "Total")
       {
        if(is.null(unsex.action)) {stop("When Sex has been set at Total, an action as to be provided for the unsexed individuals \n See options for argument 'unsex.action'")}
        if(unsex.action == "split")
          {
           unsex.n       <- length(which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length)));
           unsex.sample  <- sample(1:unsex.n,size=round(unsex.n/2));
           unsex.fem.age <- x$Data$Unsexed[which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length)),2][unsex.sample];
           unsex.fem.len <- x$Data$Unsexed[which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length)),3][unsex.sample];
           unsex.mal.age <- x$Data$Unsexed[which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length)),2][-unsex.sample];
           unsex.mal.len <- x$Data$Unsexed[which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length)),3][-unsex.sample];
           if(m=="vonBer")
             {
              if(distr == "apnormal" | distr == "aplnormal")
                {
                 if(length(par) != 6)
                   {
                    stop("For a two sexes von Bertalanffy model fit with adjusted profile likelihood par must be a vector of length 6")
                   }
                 parlist     <- list(par=par[1:3],
                                     age=c(x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                           unsex.fem.age),
                                     obslen=c(x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                              unsex.fem.len),
                                     distr=distr,
                                     output="predict",
                                     properties=x$Properties);
                 results.fem <- do.call(.vonBer, parlist);
                 parlist     <- list(par=par[4:6],
                                     age=c(x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                           unsex.mal.age),
                                     obslen=c(x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                              unsex.mal.len),
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.vonBer, parlist);
                }
              else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                {
                 if(length(par) != 8)
                   {
                    stop("For a two sexes von Bertalanffy model fit with, normal, lognormal or gamma distribution, par must be a vector of length 8")
                   }
                 parlist     <- list(par=c(par[1:3],par[7]),
                                     age=c(x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                           unsex.fem.age),
                                     obslen=c(x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                              unsex.fem.len),
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.vonBer, parlist);
                 parlist     <- list(par=c(par[4:6],par[8]),
                                     age=c(x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                           unsex.mal.age),
                                     obslen=c(x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                              unsex.mal.len),
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.vonBer, parlist);
                }
             }
           if(m=="schnute1")
             {
              if(distr == "apnormal" | distr == "aplnormal")
                {
                 if(length(par) != 10)
                   {
                    stop("For a two sexes Schnute 1 model fit with adjusted profile likelihood par must be a vector of length 10")
                   }
                 parlist     <- list(par=par[1:5],
                                     age=c(x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                           unsex.fem.age),
                                     obslen=c(x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                              unsex.fem.len),
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute1, parlist);
                 parlist     <- list(par=par[6:10],
                                     age=c(x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                           unsex.mal.age),
                                     obslen=c(x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                              unsex.mal.len),
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute1, parlist);
                }
              else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                {
                 if(length(par) != 12)
                   {
                    stop("For a two sexes Schnute 1 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 12")
                   }
                 parlist     <- list(par=c(par[1:5],par[11]),
                                     age=c(x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                           unsex.fem.age),
                                     obslen=c(x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                              unsex.fem.len),
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute1, parlist);
                 parlist     <- list(par=c(par[6:10],par[12]),
                                     age=c(x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                           unsex.mal.age),
                                     obslen=c(x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                              unsex.mal.len),
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute1, parlist);
                }
             }
           if(m=="schnute2")
             {
              if(distr == "apnormal" | distr == "aplnormal")
                {
                 if(length(par) != 6)
                   {
                    stop("For a two sexes Schnute 1 model fit with adjusted profile likelihood par must be a vector of length 6")
                   }
                 parlist     <- list(par=par[1:3],
                                     age=c(x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                           unsex.fem.age),
                                     obslen=c(x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                              unsex.fem.len),
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute2, parlist);
                 parlist     <- list(par=par[4:6],
                                     age=c(x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                           unsex.mal.age),
                                     obslen=c(x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                              unsex.mal.len),
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute2, parlist);
                }
              else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                {
                 if(length(par) != 8)
                   {
                    stop("For a two sexes Schnute 1 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 8")
                   }
                 parlist     <- list(par=c(par[1:3],par[7]),
                                     age=c(x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                           unsex.fem.age),
                                     obslen=c(x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                              unsex.fem.len),
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute2, parlist);
                 parlist     <- list(par=c(par[4:6],par[8]),
                                     age=c(x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                           unsex.mal.age),
                                     obslen=c(x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                              unsex.mal.len),
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute2, parlist);
                }
             }
           if(m=="schnute3")
             {
              if(distr == "apnormal" | distr == "aplnormal")
                {
                 if(length(par) != 4)
                   {
                    stop("For a two sexes Schnute 1 model fit with adjusted profile likelihood par must be a vector of length 4")
                   }
                 parlist     <- list(par=par[1:2],
                                     age=c(x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                           unsex.fem.age),
                                     obslen=c(x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                              unsex.fem.len),
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute3, parlist);
                 parlist     <- list(par=par[3:4],
                                     age=c(x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                           unsex.mal.age),
                                     obslen=c(x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                              unsex.mal.len),
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute3, parlist);
                }
              else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                {
                 if(length(par) != 6)
                   {
                    stop("For a two sexes Schnute 1 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 6")
                   }
                 parlist     <- list(par=c(par[1:2],par[5]),
                                     age=c(x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                           unsex.fem.age),
                                     obslen=c(x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                              unsex.fem.len),
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute3, parlist);
                 parlist     <- list(par=c(par[3:4],par[6]),
                                     age=c(x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                           unsex.mal.age),
                                     obslen=c(x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                              unsex.mal.len),
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute3, parlist);
                }
             }
           results                  <- vector("list",3)
           names(results)           <- c("Properties","Females","Males")
           results[["Properties"]]  <- x$Properties
           results$Properties$Model <- m
           results$Properties$Distr <- distr
           results[["Females"]]     <- results.fem
           results[["Males"]]       <- results.mal
           class(results)           <- "GroAgeExp"
           return(results)
          }
        if(unsex.action == "males")
          {
           unsex.mal.age <- x$Data$Unsexed$Age[which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length))];
           unsex.mal.len <- x$Data$Unsexed$Length[which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length))];
           if(m=="vonBer")
             {
              if(distr == "apnormal" | distr == "aplnormal")
                {
                 if(length(par) != 6)
                   {
                    stop("For a two sexes von Bertalanffy model fit with adjusted profile likelihood par must be a vector of length 6")
                   }
                 parlist     <- list(par=par[1:3],
                                     age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.vonBer, parlist);
                 parlist     <- list(par=par[4:6],
                                     age=c(x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                           unsex.mal.age),
                                     obslen=c(x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                              unsex.mal.len),
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.vonBer, parlist);
                }
              else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                {
                 if(length(par) != 8)
                   {
                    stop("For a two sexes von Bertalanffy model fit with, normal, lognormal or gamma distribution, par must be a vector of length 8")
                   }
                 parlist     <- list(par=c(par[1:3],par[7]),
                                     age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.vonBer, parlist);
                 parlist     <- list(par=c(par[4:6],par[8]),
                                     age=c(x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                           unsex.mal.age),
                                     obslen=c(x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                              unsex.mal.len),
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.vonBer, parlist);
                }
             }
           if(m=="schnute1")
             {
              if(distr == "apnormal" | distr == "aplnormal")
                {
                 if(length(par) != 10)
                   {
                    stop("For a two sexes Schnute 1 model fit with adjusted profile likelihood par must be a vector of length 10")
                   }
                 parlist     <- list(par=par[1:5],
                                     age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute1, parlist);
                 parlist     <- list(par=par[6:10],
                                     age=c(x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                           unsex.mal.age),
                                     obslen=c(x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                              unsex.mal.len),
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute1, parlist);
                }
              else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                {
                 if(length(par) != 12)
                   {
                    stop("For a two sexes Schnute 1 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 12")
                   }
                 parlist     <- list(par=c(par[1:5],par[11]),
                                     age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute1, parlist);
                 parlist     <- list(par=c(par[6:10],par[12]),
                                     age=c(x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                           unsex.mal.age),
                                     obslen=c(x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                              unsex.mal.len),
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute1, parlist);
                }
             }
           if(m=="schnute2")
             {
              if(distr == "apnormal" | distr == "aplnormal")
                {
                 if(length(par) != 6)
                   {
                    stop("For a two sexes Schnute 1 model fit with adjusted profile likelihood par must be a vector of length 6")
                   }
                 parlist     <- list(par=par[1:3],
                                     age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute2, parlist);
                 parlist     <- list(par=par[4:6],
                                     age=c(x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                           unsex.mal.age),
                                     obslen=c(x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                              unsex.mal.len),
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute2, parlist);
                }
              else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                {
                 if(length(par) != 8)
                   {
                    stop("For a two sexes Schnute 1 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 8")
                   }
                 parlist     <- list(par=c(par[1:3],par[7]),
                                     age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute2, parlist);
                 parlist     <- list(par=c(par[4:6],par[8]),
                                     age=c(x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                           unsex.mal.age),
                                     obslen=c(x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                              unsex.mal.len),
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute2, parlist);
                }
             }
           if(m=="schnute3")
             {
              if(distr == "apnormal" | distr == "aplnormal")
                {
                 if(length(par) != 4)
                   {
                    stop("For a two sexes Schnute 1 model fit with adjusted profile likelihood par must be a vector of length 4")
                   }
                 parlist     <- list(par=par[1:2],
                                     age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute3, parlist);
                 parlist     <- list(par=par[3:4],
                                     age=c(x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                           unsex.mal.age),
                                     obslen=c(x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                              unsex.mal.len),
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute3, parlist);
                }
              else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                {
                 if(length(par) != 6)
                   {
                    stop("For a two sexes Schnute 1 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 6")
                   }
                 parlist     <- list(par=c(par[1:2],par[5]),
                                     age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute3, parlist);
                 parlist     <- list(par=c(par[3:4],par[6]),
                                     age=c(x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                           unsex.mal.age),
                                     obslen=c(x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                              unsex.mal.len),
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute3, parlist);
                }
             }
           results                  <- vector("list",3)
           names(results)           <- c("Properties","Females","Males")
           results[["Properties"]]  <- x$Properties
           results$Properties$Model <- m
           results$Properties$Distr <- distr
           results[["Females"]]     <- results.fem
           results[["Males"]]       <- results.mal
           class(results)           <- "GroAgeExp"
           return(results)
          }
        if(unsex.action == "females")
          {
           unsex.fem.age <- x$Data$Unsexed$Age[which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length))];
           unsex.fem.len <- x$Data$Unsexed$Length[which(!is.na(x$Data$Unsexed$Age) & !is.na(x$Data$Unsexed$Length))];
           if(m=="vonBer")
             {
              if(distr == "apnormal" | distr == "aplnormal")
                {
                 if(length(par) != 6)
                   {
                    stop("For a two sexes von Bertalanffy model fit with adjusted profile likelihood par must be a vector of length 6")
                   }
                 parlist     <- list(par=par[1:3],
                                     age=c(x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                           unsex.fem.age),
                                     obslen=c(x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                              unsex.fem.len),
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.vonBer, parlist);
                 parlist     <- list(par=par[4:6],
                                     age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.vonBer, parlist);
                }
              else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                {
                 if(length(par) != 8)
                   {
                    stop("For a two sexes von Bertalanffy model fit with, normal, lognormal or gamma distribution, par must be a vector of length 8")
                   }
                 parlist     <- list(par=c(par[1:3],par[7]),
                                     age=c(x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                           unsex.fem.age),
                                     obslen=c(x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                              unsex.fem.len),
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.vonBer, parlist);
                 parlist     <- list(par=c(par[4:6],par[8]),
                                     age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.vonBer, parlist);
                }
             }
           if(m=="schnute1")
             {
              if(distr == "apnormal" | distr == "aplnormal")
                {
                 if(length(par) != 10)
                   {
                    stop("For a two sexes Schnute 1 model fit with adjusted profile likelihood par must be a vector of length 10")
                   }
                 parlist     <- list(par=par[1:5],
                                     age=c(x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                           unsex.fem.age),
                                     obslen=c(x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                              unsex.fem.len),
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute1, parlist);
                 parlist     <- list(par=par[6:10],
                                     age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute1, parlist);
                }
              else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                {
                 if(length(par) != 12)
                   {
                    stop("For a two sexes Schnute 1 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 12")
                   }
                 parlist     <- list(par=c(par[1:5],par[11]),
                                     age=c(x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                           unsex.fem.age),
                                     obslen=c(x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                              unsex.fem.len),
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute1, parlist);
                 parlist     <- list(par=c(par[6:10],par[12]),
                                     age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute1, parlist);
                }
             }
           if(m=="schnute2")
             {
              if(distr == "apnormal" | distr == "aplnormal")
                {
                 if(length(par) != 6)
                   {
                    stop("For a two sexes Schnute 1 model fit with adjusted profile likelihood par must be a vector of length 6")
                   }
                 parlist     <- list(par=par[1:3],
                                     age=c(x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                           unsex.fem.age),
                                     obslen=c(x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                              unsex.fem.len),
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute2, parlist);
                 parlist     <- list(par=par[4:6],
                                     age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute2, parlist);
                }
              else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                {
                 if(length(par) != 8)
                   {
                    stop("For a two sexes Schnute 1 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 8")
                   }
                 parlist     <- list(par=c(par[1:3],par[7]),
                                     age=c(x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                           unsex.fem.age),
                                     obslen=c(x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                              unsex.fem.len),
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute2, parlist);
                 parlist     <- list(par=c(par[4:6],par[8]),
                                     age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute2, parlist);
                }
             }
           if(m=="schnute3")
             {
              if(distr == "apnormal" | distr == "aplnormal")
                {
                 if(length(par) != 4)
                   {
                    stop("For a two sexes Schnute 1 model fit with adjusted profile likelihood par must be a vector of length 4")
                   }
                 parlist     <- list(par=par[1:2],
                                     age=c(x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                           unsex.fem.age),
                                     obslen=c(x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                              unsex.fem.len),
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute3, parlist);
                 parlist     <- list(par=par[3:4],
                                     age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute3, parlist);
                }
              else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                {
                 if(length(par) != 6)
                   {
                    stop("For a two sexes Schnute 1 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 6")
                   }
                 parlist     <- list(par=c(par[1:2],par[5]),
                                     age=c(x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                           unsex.fem.age),
                                     obslen=c(x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                              unsex.fem.len),
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute3, parlist);
                 parlist     <- list(par=c(par[3:4],par[6]),
                                     age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute3, parlist);
                }
             }
           results                  <- vector("list",3)
           names(results)           <- c("Properties","Females","Males")
           results[["Properties"]]  <- x$Properties
           results$Properties$Model <- m
           results$Properties$Distr <- distr
           results[["Females"]]     <- results.fem
           results[["Males"]]       <- results.mal
           class(results)           <- "GroAgeExp"
           return(results)
          }
        if(unsex.action == "ignore")
          {
           if(m=="vonBer")
             {
              if(distr == "apnormal" | distr == "aplnormal")
                {
                 if(length(par) != 6)
                   {
                    stop("For a two sexes von Bertalanffy model fit with adjusted profile likelihood par must be a vector of length 6")
                   }
                 parlist     <- list(par=par[1:3],
                                     age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.vonBer, parlist);
                 parlist     <- list(par=par[4:6],
                                     age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.vonBer, parlist);
                }
              else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                {
                 if(length(par) != 8)
                   {
                    stop("For a two sexes von Bertalanffy model fit with, normal, lognormal or gamma distribution, par must be a vector of length 8")
                   }
                 parlist     <- list(par=c(par[1:3],par[7]),
                                     age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.vonBer, parlist);
                 parlist     <- list(par=c(par[4:6],par[8]),
                                     age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.vonBer, parlist);
                }
             }
           if(m=="schnute1")
             {
              if(distr == "apnormal" | distr == "aplnormal")
                {
                 if(length(par) != 10)
                   {
                    stop("For a two sexes Schnute 1 model fit with adjusted profile likelihood par must be a vector of length 10")
                   }
                 parlist     <- list(par=par[1:5],
                                     age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute1, parlist);
                 parlist     <- list(par=par[6:10],
                                     age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute1, parlist);
                }
              else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                {
                 if(length(par) != 12)
                   {
                    stop("For a two sexes Schnute 1 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 12")
                   }
                 parlist     <- list(par=c(par[1:5],par[11]),
                                     age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute1, parlist);
                 parlist     <- list(par=c(par[6:10],par[12]),
                                     age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute1, parlist);
                }
             }
           if(m=="schnute2")
             {
              if(distr == "apnormal" | distr == "aplnormal")
                {
                 if(length(par) != 6)
                   {
                    stop("For a two sexes Schnute 1 model fit with adjusted profile likelihood par must be a vector of length 6")
                   }
                 parlist     <- list(par=par[1:3],
                                     age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute2, parlist);
                 parlist     <- list(par=par[4:6],
                                     age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute2, parlist);
                }
              else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                {
                 if(length(par) != 8)
                   {
                    stop("For a two sexes Schnute 1 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 8")
                   }
                 parlist     <- list(par=c(par[1:3],par[7]),
                                     age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute2, parlist);
                 parlist     <- list(par=c(par[4:6],par[8]),
                                     age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute2, parlist);
                }
             }
           if(m=="schnute3")
             {
              if(distr == "apnormal" | distr == "aplnormal")
                {
                 if(length(par) != 4)
                   {
                    stop("For a two sexes Schnute 1 model fit with adjusted profile likelihood par must be a vector of length 4")
                   }
                 parlist     <- list(par=par[1:2],
                                     age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute3, parlist);
                 parlist     <- list(par=par[3:4],
                                     age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute3, parlist);
                }
              else if(distr == "normal" | distr == "lognormal" | distr == "gamma")
                {
                 if(length(par) != 6)
                   {
                    stop("For a two sexes Schnute 1 model fit with, normal, lognormal or gamma distribution, par must be a vector of length 6")
                   }
                 parlist     <- list(par=c(par[1:2],par[5]),
                                     age=x$Data$Females$Age[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     obslen=x$Data$Females$Length[which(!is.na(x$Data$Females$Age) & !is.na(x$Data$Females$Length))],
                                     distr=distr,
                                     output="predict");
                 results.fem <- do.call(.schnute3, parlist);
                 parlist.mal <- list(par=c(par[3:4],par[6]),
                                     age=x$Data$Males$Age[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     obslen=x$Data$Males$Length[which(!is.na(x$Data$Males$Age) & !is.na(x$Data$Males$Length))],
                                     distr=distr,
                                     output="predict");
                 results.mal  <- do.call(.schnute3, parlist.mal);
                }
             }
           results                  <- vector("list",3)
           names(results)           <- c("Properties","Females","Males")
           results[["Properties"]]  <- x$Properties
           results$Properties$Model <- m
           results$Properties$Distr <- distr
           results[["Females"]]     <- results.fem
           results[["Males"]]       <- results.mal
           class(results)           <- "GroAgeExp"
           return(results)
          }
       }
  }
