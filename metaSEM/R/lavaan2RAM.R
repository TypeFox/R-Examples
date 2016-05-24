lavaan2RAM <- function(model, obs.variables = NULL, A.notation="ON", S.notation="WITH",
                       auto.var = TRUE, std.lv = TRUE, ...) {
    if (!requireNamespace("lavaan", quietly=TRUE))    
        stop("\"lavaan\" package is required for this function.")

    ## fix the latent independent variables at 1
    my.model <- lavaan::lavaanify(model, fixed.x = FALSE, auto.var=auto.var, std.lv=TRUE, ...)
    
    ## replace "." with ""
    my.model$plabel <- gsub("\\.", "", my.model$plabel)
    ## remove constraints 
    ## there are constraints such as .p1. and .p2; remove them
    if (length(grep("^\\.", my.model$lhs)) >0 ) my.model <- my.model[-grep("^\\.", my.model$lhs), ]
    ## label as the variable labels
    #my.model$label <- with(my.model, ifelse(label=="", yes=plabel, no=label))
    ## set the starting values in A as 0 if NA
    my.model[(my.model$op=="=~"|my.model$op=="~")&is.na(my.model$ustart), ]$ustart <- 0
    ## set the starting values in S and free parameters as 0 if NA
    my.model[my.model$op=="~~"&is.na(my.model$ustart)&my.model$free!=0, ]$ustart <- 0
  
    ## all variables
    all.var <- sort(unique(c(my.model$lhs, my.model$rhs)))
    ## latent variables: (with indicators)
    latent <- unique(my.model[my.model$op== "=~", ]$lhs)
    ## observed variables: not latent
    observed <- all.var[!(all.var %in% latent)]
    ## check whether observed in model = observed in argument
    if (!is.null(obs.variables)) {
        if (!identical(sort(observed), sort(obs.variables))) {
            stop("Names in \"obs.variables\" do not agree with those in model.\n")
        } else {
            ## arrange the observed variables according to obs.var argument
            observed <- obs.variables
        }
    }

    ## if there are latent variables
    if (length(latent)>0) {
        ## arrange variable list as observed + latent
        all.var <- c(observed, latent)
    } else {
        all.var <- observed
    }
    
    no.lat <- length(latent)
    no.obs <- length(observed)
    no.all <- no.lat+no.obs
  
    Amatrix <- matrix(0, ncol=no.all, nrow=no.all, dimnames=list(all.var, all.var))
    Smatrix <- matrix(0, ncol=no.all, nrow=no.all, dimnames=list(all.var, all.var))
  
  ## ## latent variable
  ## lhs <- (all.var %in% my.model[my.model$op== "=~", ]$lhs) 
  ## ## observed variables
  ## rhs <- (all.var %in% my.model[my.model$op== "=~", ]$rhs)

    ## Prepare the labels
    for (i in seq_len(nrow(my.model))) {
      ## if there is no label
      if (my.model[i, ]$label=="") {
        switch(my.model[i, ]$op,
               "=~" = my.model[i, ]$label <- paste0(my.model[i, ]$rhs, A.notation, my.model[i, ]$lhs),
               "~"  = my.model[i, ]$label <- paste0(my.model[i, ]$lhs, A.notation, my.model[i, ]$rhs),
               "~~" = my.model[i, ]$label <- paste0(my.model[i, ]$lhs, S.notation, my.model[i, ]$rhs))
      }
    }
    
    ## keys in as.mxMatrix format
    key <- with(my.model, ifelse(free==0, yes=ustart, no=paste(ustart, label, sep="*")))  
    
    for (i in seq_len(nrow(my.model))) {
        my.line <- my.model[i, ]
        switch(my.line$op,
               ## lhs: IV; rhs: DV
               "=~" = Amatrix[my.line$rhs, my.line$lhs] <- key[i],
               ## lhs: DV; rhs: IV
               "~" = Amatrix[my.line$lhs, my.line$rhs] <- key[i],
               "~~" = Smatrix[my.line$lhs, my.line$rhs] <- Smatrix[my.line$rhs, my.line$lhs] <- key[i])
    }

    Fmatrix <- create.Fmatrix(c(rep(1, no.obs), rep(0, no.lat)), as.mxMatrix=FALSE)
    dimnames(Fmatrix) <- list(observed, all.var)

    list(A=Amatrix, S=Smatrix, F=Fmatrix)
}
