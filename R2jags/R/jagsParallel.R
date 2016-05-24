jags.parallel <- function (data, inits, parameters.to.save, model.file = "model.bug",
    n.chains = 2, n.iter = 2000, n.burnin = floor(n.iter/2),
    n.thin = max(1, floor((n.iter - n.burnin)/1000)), n.cluster = n.chains,
    DIC = TRUE, working.directory = NULL, jags.seed = 123, digits = 5,
    RNGname = c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper",
        "Mersenne-Twister"), jags.module = c("glm", "dic"),  export_obj_names=NULL, envir = .GlobalEnv)
{
    jags.params <- parameters.to.save
    jags.inits <- if (missing(inits)) {
        NULL
    }
    else {
        inits
    }
    jags.model <- model.file

    if( is.character( data ) && length(data) == 1
                           && regexpr( "\\.txt$", data ) > 0 ) {
    ## 1. 'data' looks like a file name [UNDOCUMENTED!]
    if ( all( basename( data ) == data )) {
      fn2 <- file.path( working.directory, data )
      if (normalizePath(fn2)!=normalizePath(data)) {  ## file.copy() to same place trashes the file
        try( file.copy(fn2 , data, overwrite = TRUE ) )
      }
    }
    if ( !file.exists( data ) ) {
      stop("File",data,"does not exist")
    }
    if (file.info(data)["size"]==0) {
      stop("Empty data file ",data)
    }
    e    <- new.env()
    eval( parse( data ), e )
    data <- as.list( e )
  } else if( is.character( data ) ||
             ( is.list( data ) && all( sapply( data,is.character ) ) ) ){
    ## 2. data is a character vector or a list of character
    dlist          <- lapply( as.list( data ), get, envir = parent.frame( 1 ) )
    names( dlist ) <- unlist( data )
    data           <- dlist
  } else if( !is.list( data ) ){
    stop( "data must be a character vector of object names, a list of object names, or a list of objects" )
  }

   list2env(data, envir=envir )




    .runjags <- function() {
        jagsfit <- jags(data = eval(expression(data)), inits = jags.inits,
            parameters.to.save = eval(expression(jags.params)),
            model.file = eval(expression(jags.model)), n.chains = 1,
            n.iter = eval(expression(n.iter)), n.burnin = eval(expression(n.burnin)),
            n.thin = eval(expression(n.thin)), DIC = eval(expression(DIC)),
            working.directory = eval(expression(working.directory)),
            jags.seed = eval(expression(jags.seed)), progress.bar = "none",
            digits = eval(expression(digits)), RNGname = eval(expression(RNGname)),
            jags.module = eval(expression(jags.module)) )
        return(jagsfit)
    }
    cl <- makeCluster(n.cluster, methods = FALSE)
    clusterExport(cl, c(names(data), "mcmc", "mcmc.list", export_obj_names ), envir = envir)
    clusterSetRNGStream(cl, jags.seed)
    tryCatch(res <- clusterCall(cl, .runjags), finally = stopCluster(cl))
    result <- NULL
    model <- NULL
    for (ch in 1:n.chains) {
        result <- abind(result, res[[ch]]$BUGSoutput$sims.array,
            along = 2)
        model[[ch]] <- res[[ch]]$model
    }
    if (is.function(model.file)) {
        model.file <- substitute(model.file)
    }
    result <- as.bugs.array2(result, model.file = model.file,
        program = "jags", DIC = DIC, n.iter = n.iter, n.burnin = n.burnin,
        n.thin = n.thin)
    out <- list(model = model, BUGSoutput = result, parameters.to.save = parameters.to.save,
        model.file = model.file, n.iter = n.iter, DIC = DIC)
    class(out) <- c("rjags.parallel", "rjags")
    return(out)
}
