jags <- function( data, inits,
                  parameters.to.save,
                  model.file  = "model.bug",
                  n.chains     = 3,
                  n.iter       = 2000,
                  n.burnin     = floor(n.iter/2),
                  n.thin       = max( 1, floor( ( n.iter - n.burnin )/1000 ) ),
                  DIC          = TRUE,
                  working.directory = NULL,
                  jags.seed    = 123,
                  refresh      = n.iter/50,
                  progress.bar = "text",
                  digits = 5,
                  RNGname = c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper", "Mersenne-Twister"),
                  jags.module = c("glm","dic")
                  )
{
  #require( rjags )
  if( !is.null( working.directory ) ){
    working.directory <- path.expand( working.directory )
    savedWD <- getwd()
    setwd( working.directory )
    on.exit( setwd( savedWD ) )
  } else {
    savedWD <- getwd()
    working.directory <- savedWD
  }
  ## jags.model() needs 'data' to be "a list or environment containing the data
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

  ## copied from R2WinBUGS:
  if (is.function(model.file)) {
    temp <- tempfile("model")
    temp <- if (is.R() || .Platform$OS.type != "windows") {
      paste(temp, "txt", sep = ".")
    }
    else {
      gsub("\\.tmp$", ".txt", temp)
    }
    write.model(model.file, con = temp, digits = digits) ## from R2WinBUGS
    model.file <- gsub("\\\\", "/", temp)
    if (!is.R())
      on.exit(file.remove(model.file), add = TRUE)
  }
  if( DIC ){
    parameters.to.save <- c( parameters.to.save, "deviance" )
    load.module( "dic", quiet = TRUE )
  }

  if( n.burnin > 0 ){
    n.adapt <- n.burnin
  } else{
    n.adapt <- 100
  }

  if (!missing(inits) && !is.function(inits) && !is.null(inits) && (length(inits) != n.chains)) {
    stop("Number of initialized chains (length(inits)) != n.chains")
  }

  RNGname <- match.arg(RNGname)
  if(RNGname %in% c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper", "Mersenne-Twister")){
    RNGname <- paste("base::",RNGname,sep="")
  }else{
    stop("The name of the RNG is not correctly provided!")
  }


  #load.module("lecuyer")
  if(!is.null(jags.module)){
    n.module <- length(jags.module)
    for(m in 1:n.module){
        load.module(jags.module[m])
    }
  }

  #load.module(jags.module)
  init.values <- vector("list", n.chains)
  if(missing(inits)){
    for (i in 1:n.chains){
        init.values[[i]]$.RNG.name <- RNGname
        init.values[[i]]$.RNG.seed <- runif(1, 0, 2^31)#abs(.Random.seed[i+1])
    }
  } else if (is.null(inits)){
      for (i in 1:n.chains){
        init.values[[i]]$.RNG.name <- RNGname
        init.values[[i]]$.RNG.seed <- runif(1, 0, 2^31)#abs(.Random.seed[i+1])
      }
  } else if (is.function(inits)){
      if (any(names(formals(inits)) == "chain")){
        for (i in 1:n.chains){
          init.values[[i]] <- inits(chain = i)
          init.values[[i]]$.RNG.name <- RNGname
          init.values[[i]]$.RNG.seed <- runif(1, 0, 2^31)#abs(.Random.seed[i+1])
        }
      } else{
          for (i in 1:n.chains){
            init.values[[i]] <- inits()
            init.values[[i]]$.RNG.name <- RNGname
            init.values[[i]]$.RNG.seed <- runif(1, 0, 2^31)#abs(.Random.seed[i+1])
          }
      }
  } else {
    #if (!is.function(inits) && !is.null(inits) && (length(inits) != n.chains)) {
    if (!is.list(inits)) {
      stop("Invalid inits")
    }
    if (length(inits) != n.chains) {
      stop("Number of initialized chains (length(inits)) != n.chains")
    }
    for (i in 1:n.chains){
      init.values[[i]] <- inits[[i]]
      init.values[[i]]$.RNG.name <- RNGname
      init.values[[i]]$.RNG.seed <- runif(1, 0, 2^31)#abs(.Random.seed[i+1])
    }
   }


#  if( is.null( inits ) ){
#    m <- jags.model( model.file,
#                     data     = data,
#                     n.chains = n.chains,
#                     n.adapt  = 0 )
#  } else{

  m <- jags.model(model.file,
                  data     = data,
                  inits    = init.values,
                  n.chains = n.chains,
                  n.adapt  = 0 )
  #}
  adapt( m,
         n.iter         = n.adapt,
         by             = refresh,
         progress.bar   = progress.bar,
         end.adaptation = TRUE )

  samples <- coda.samples( model          = m,
                           variable.names = parameters.to.save,
                           n.iter         = ( n.iter - n.burnin ),
                           thin           = n.thin,
                           by             = refresh,
                           progress.bar   = progress.bar )

  fit <- mcmc2bugs( samples,
                    model.file = model.file,
                    program    = "jags",
                    DIC        = DIC,
                    DICOutput  = NULL,
                    n.iter     = n.iter,
                    n.burnin   = n.burnin,
                    n.thin     = n.thin )

  out <- list( model              = m,
               BUGSoutput         = fit,
               parameters.to.save = parameters.to.save,
               model.file         = model.file,
               n.iter             = n.iter,
               DIC                = DIC)

  class(out) <- "rjags"
  return(out)
}
