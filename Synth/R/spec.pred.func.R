spec.pred.func <-
function(
                           list.object = NULL,
                           tr.numb = NULL,
                           co.numb = NULL,
                           unit.var = NULL,
                           time.var = NULL,
                           foo.object = NULL,
                           X0.inner = NULL,
                           X1.inner = NULL
                           )
  {

    special.units.tr <-
      which(is.element(foo.object[,unit.var], tr.numb) == TRUE)

    special.units.co <-
      which(is.element(foo.object[,unit.var], co.numb) == TRUE)

    special.times <-
      which(is.element(foo.object[,time.var], list.object[[2]]) == TRUE)

    if(mode(list.object[[1]]) == "character")
      {
        list.object[[1]] <- which(names(foo.object) == list.object[[1]])
      }

    # create predictor name:
    name.predictor <- names(foo.object)[list.object[[1]]]
    if(length(list.object[[2]])>1)
      {
      name.predictor <-
      paste("special.", name.predictor, ".",
            paste(list.object[[2]][1],".",
                  list.object[[2]][length(list.object[[2]])]
                                  , sep = ""), sep = "")

      } else {

      name.predictor <-
      paste("special.", name.predictor, ".",
            paste(list.object[[2]]
            , sep = ""), sep = "")
      }

      X1.special <- as.matrix(foo.object[
                                       intersect(
                                                 special.units.tr,
                                                 special.times
                                                ),
                                       list.object[[1]]
                                       ]
                            )

    rownames(X1.special) <-
      foo.object[
                 intersect(special.units.tr, special.times),
                 time.var
                 ]
    colnames(X1.special) <- tr.numb

   # deep missing checker:
    for(i in 1:ncol(X1.special))
     {
       if(sum(is.na(X1.special[,i])) == nrow(X1.special))
        {stop(paste("\n special predictor",name.predictor,"has missing data for all periods in time.predictors.prior\n"))}

      for(j in 1:nrow(X1.special))
       {
         if(is.na(X1.special[j,i])){
         cat(paste("\n Missing data: treated unit; special predictor:",name.predictor,"; for period:",rownames(X1.special)[j],
         "\n","We ignore (na.rm = TRUE) all missing values for predictors.op.\n"))}
       }
     }

    # aggregate if only one time period is there
    if(length(list.object[[2]]) > 1)
      {
       X1.special <- apply(X1.special,
                            2,
                            paste(list.object[[3]]),
                            na.rm = TRUE
                            )
      }
    X1.special <- t(as.matrix(X1.special))
    X1.inner <- rbind(X1.inner, X1.special)
    rownames(X1.inner)[nrow(X1.inner)] <- name.predictor

   # controls
    X0.special <-
      as.matrix(foo.object[intersect(
                                     special.units.co,
                                     special.times
                                     ),
                           c(list.object[[1]])
                           ]
                )

    X0.special <-
      matrix(X0.special[,1], byrow = FALSE, ncol = length(co.numb))

    # Define row and column names

    row.names(X0.special) <-
      foo.object[intersect(special.units.tr, special.times), time.var]

    colnames(X0.special) <- co.numb

    # deep missing checker:
    for(i in 1:ncol(X0.special))
     {
      if(sum(is.na(X0.special[,i])) == nrow(X0.special))
       {stop(paste("\n control unit:",co.numb[i],"; special predictor:",name.predictor,"has missing data for all periods specified:",list.object[[2]],"\n"))}

       for(j in 1:nrow(X0.special))
       {
         if(is.na(X0.special[j,i])){
         cat(paste("\n Missing data - control unit:",co.numb[i],"; special predictor:",name.predictor,"; for period:",rownames(X0.special)[j],
         "\n","We ignore (na.rm = TRUE) all missing values for predictors.op.\n"))}
       }
     }

    # Continue with object creation
    if(length(list.object[[2]]) > 1)
     {
        X0.special <- apply(
                            X0.special,
                            2,
                            paste(list.object[[3]]),
                            na.rm = TRUE
                            )
                            
        X0.special <- t(as.matrix(X0.special))
     }

    X0.inner <- rbind(X0.inner,X0.special)
    rownames(X0.inner)[nrow(X0.inner)] <- name.predictor




    # Prepare output
    special.output <- list(X0.inner = X0.inner,
                           X1.inner = X1.inner
                           )

    return(invisible(special.output))

  }

