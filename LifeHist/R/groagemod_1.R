groagemod <-
function(x, method)
  {
     z <- which(x$Methods == method)
     if(x$Data$Properties$Sex == "Males" | x$Data$Properties$Sex == "Females" | x$Data$Properties$Sex == "Pooled")
       {
        if(x[[4]][[1]][[z]][3] == "FAIL")
          {
           stop("This numerical method failed. Try another method")
          }
        m     <- x[[4]][[1]][[z]][1];
        distr <- x[[4]][[1]][[z]][2];
        parlist <- list(par=log(unlist(x[[4]][[1]][[z]][6])),
                        age=x$Data$Data[[1]][,2][which(!is.na(x$Data$Data[[1]][,2]) & !is.na(x$Data$Data[[1]][,3]))],
                        obslen=x$Data$Data[[1]][,3][which(!is.na(x$Data$Data[[1]][,2]) & !is.na(x$Data$Data[[1]][,3]))],
                        output="predict",
                        distr=x[[4]][[1]][[z]][2]);
        if(m=="vonBer")
          {
          results.usx <- do.call(.vonBer, parlist);
          }
        else if(m=="schnute1")
          {
          results.usx <- do.call(.schnute1, parlist);
          }
        else if(m=="schnute2")
          {
          results.usx <- do.call(.schnute2, parlist);
          }
        else if(m=="schnute3")
          {
          results.usx <- do.call(.schnute3, parlist);
          }
        results                          <- vector("list",3)
        names(results)                   <- c("Properties","Model","Predictions")
        results[["Properties"]]          <- x$Data$Properties
        results[["Model"]]               <- vector("list",3)
        names(results[["Model"]])        <- c("Growth","Distribution","Method")
        results[["Model"]]$Growth        <- m
        results[["Model"]]$Distribution  <- distr
        results[["Model"]]$Method        <- method
        results[["Predictions"]]         <- results.usx
        class(results)                   <- "GroAgeMod"
        return(results)
       }
     else if(x$Data$Properties$Sex == "Both" | x$Data$Properties$Sex == "Total")
       {
        m     <- x[[4]][[1]][[z]][1]
        distr <- x[[4]][[1]][[z]][2]
        if(x[[4]][[1]][[z]]$converg != "FAIL")
          {
           parlist.fem     <- list(par=log(unlist(x[[4]][[1]][[z]][6])),
                                   age=x$Data$Data$Females$Age[which(!is.na(x$Data$Data$Females$Age) & !is.na(x$Data$Data$Females$Length))],
                                   obslen=x$Data$Data$Females$Length[which(!is.na(x$Data$Data$Females$Age) & !is.na(x$Data$Data$Females$Length))],
                                   output="predict",
                                   distr=x[[4]][[1]][[z]][2]);
           if(m == "vonber")
             {
              results.fem  <- do.call(.vonBer, parlist.fem);
             }
           else if(m == "schnute1")
             {
              results.fem  <- do.call(.schnute1, parlist.fem);
             }
           else if(m == "schnute2")
             {
              results.fem  <- do.call(.schnute2, parlist.fem);
             }
           else if(m == "schnute3")
             {
              results.fem  <- do.call(.schnute3, parlist.fem);
             }
          }
        if(x[[4]][[2]][[z]]$converg != "FAIL")
          {
           parlist.mal     <- list(par=log(unlist(x[[4]][[2]][[z]][6])),
                                   age=x$Data$Data$Males$Age[which(!is.na(x$Data$Data$Males$Age) & !is.na(x$Data$Data$Males$Length))],
                                   obslen=x$Data$Data$Males$Length[which(!is.na(x$Data$Data$Males$Age) & !is.na(x$Data$Data$Males$Length))],
                                   output="predict",
                                   distr=x[[4]][[2]][[z]][2]);
           if(m == "vonber")
             {
              results.mal  <- do.call(.vonBer, parlist.mal);
             }
           else if(m == "schnute1")
             {
              results.mal  <- do.call(.schnute1, parlist.mal);
             }
           else if(m == "schnute2")
             {
              results.mal  <- do.call(.schnute2, parlist.mal);
             }
           else if(m == "schnute3")
             {
              results.mal  <- do.call(.schnute3, parlist.mal);
             }
          }
        if(x[[4]][[3]][[z]]$converg != "FAIL")
          {
           parlist.pld     <- list(par=log(unlist(x[[4]][[3]][[z]][6])),
                                   age=c(x$Data$Data$Females$Age[which(!is.na(x$Data$Data$Females$Age) & !is.na(x$Data$Data$Females$Length))],
                                         x$Data$Data$Males$Age[which(!is.na(x$Data$Data$Males$Age) & !is.na(x$Data$Data$Males$Length))]),
                                   obslen=c(x$Data$Data$Females$Length[which(!is.na(x$Data$Data$Females$Age) & !is.na(x$Data$Data$Females$Length))],
                                            x$Data$Data$Males$Length[which(!is.na(x$Data$Data$Males$Age) & !is.na(x$Data$Data$Males$Length))]),
                                   output="predict",
                                   distr=x[[4]][[3]][[z]][2]);
           if(m == "vonber")
             {
              results.pld  <- do.call(.vonBer, parlist.pld);
             }
           else if(m == "schnute1")
             {
              results.pld  <- do.call(.schnute1, parlist.pld);
             }
           else if(m == "schnute2")
             {
              results.pld  <- do.call(.schnute2, parlist.pld);
             }
           else if(m == "schnute3")
             {
              results.pld  <- do.call(.schnute3, parlist.pld);
             }
          }
        results                          <- vector("list",3)
        names(results)                   <- c("Properties","Model","Predictions")
        results[["Properties"]]          <- x$Data$Properties
        results[["Model"]]               <- vector("list",3)
        names(results[["Model"]])        <- c("Growth","Distribution","Method")
        results[["Model"]]$Growth        <- m
        results[["Model"]]$Distribution  <- distr
        results[["Model"]]$Method        <- method
        if(x[[4]][[1]][[z]]$converg != "FAIL" & x[[4]][[2]][[z]]$converg != "FAIL" & x[[4]][[3]][[z]]$converg != "FAIL")
          {
           results[["Predictions"]]         <- vector("list",3)
           names(results[["Predictions"]])  <- c("Females","Males","Pooled")
           results[["Predictions"]]$Females <- results.fem
           results[["Predictions"]]$Males   <- results.mal
           results[["Predictions"]]$Pooled  <- results.pld
          }
        else if(x[[4]][[1]][[z]]$converg == "FAIL" & x[[4]][[2]][[z]]$converg != "FAIL" & x[[4]][[3]][[z]]$converg != "FAIL")
          {
           results[["Predictions"]]         <- vector("list",2)
           names(results[["Predictions"]])  <- c("Males","Pooled")
           results[["Predictions"]]$Males   <- results.mal
           results[["Predictions"]]$Pooled  <- results.pld
          }
        else if(x[[4]][[1]][[z]]$converg != "FAIL" & x[[4]][[2]][[z]]$converg == "FAIL" & x[[4]][[3]][[z]]$converg != "FAIL")
          {
           results[["Predictions"]]         <- vector("list",2)
           names(results[["Predictions"]])  <- c("Females","Pooled")
           results[["Predictions"]]$Females <- results.fem
           results[["Predictions"]]$Pooled  <- results.pld
          }
        else if(x[[4]][[1]][[z]]$converg != "FAIL" & x[[4]][[2]][[z]]$converg != "FAIL" & x[[4]][[3]][[z]]$converg == "FAIL")
          {
           results[["Predictions"]]         <- vector("list",2)
           names(results[["Predictions"]])  <- c("Females","Males")
           results[["Predictions"]]$Females <- results.fem
           results[["Predictions"]]$Males   <- results.mal
          }
        else if(x[[4]][[1]][[z]]$converg == "FAIL" & x[[4]][[2]][[z]]$converg == "FAIL" & x[[4]][[3]][[z]]$converg != "FAIL")
          {
           results[["Predictions"]]         <- vector("list",1)
           names(results[["Predictions"]])  <- c("Pooled")
           results[["Predictions"]]$Pooled  <- results.pld
          }
        else if(x[[4]][[1]][[z]]$converg != "FAIL" & x[[4]][[2]][[z]]$converg == "FAIL" & x[[4]][[3]][[z]]$converg == "FAIL")
          {
           results[["Predictions"]]         <- vector("list",1)
           names(results[["Predictions"]])  <- c("Females")
           results[["Predictions"]]$Females <- results.fem
          }
        else if(x[[4]][[1]][[z]]$converg == "FAIL" & x[[4]][[2]][[z]]$converg != "FAIL" & x[[4]][[3]][[z]]$converg == "FAIL")
          {
           results[["Predictions"]]         <- vector("list",1)
           names(results[["Predictions"]])  <- c("Males")
           results[["Predictions"]]$Males   <- results.mal
          }
        else if(x[[4]][[1]][[z]]$converg == "FAIL" & x[[4]][[2]][[z]]$converg == "FAIL" & x[[4]][[3]][[z]]$converg == "FAIL")
          {
           stop("This numerical method failed to fit models for females, males, and all data pooled. Try another method")
          }
        class(results)            <- "GroAgeMod"
        return(results)
       }
  }
