"disfitqua" <-
 function(x, f, objfun=c("rmse", "mad"),
                init.lmr=NA, init.para=NA, type=NA, verbose=FALSE, ... ) {

   if(length(f) != length(x)) {
      warning("Lengths of arguments 'f' and 'x' are not equal, returning NULL")
      return(NULL)
   }
   if(is.na(type)) {
      warning("Distribution type can not be NA, returning NULL")
      return(NULL)
   }

   match.arg(objfun)

   n.para <- dist.list(type=type) # parameter count of the named distribution

   if(! is.na(init.lmr)) {
      if(length(init.lmr$L1) == 1) init.lmr <- lmorph(init.lmr)
      if(length(init.lmr$lambdas) < n.para) {
         warning("Length of the initial L-moments in argument 'init.lmr' is ",
                 "less than the number of parameters for the specified ",
                 "distribution type, returning NULL")
         return(NULL)
      }
   }

   # Having initial parameters will have priority over having initial L-moments
   # because the optimization is based on the parameters changing not the L-moments
   if(is.na(init.para)) {
      # The premise of the "guessing" at the parameters is based on the quantiles
      # but for many of intend applications (risk analysis) quantiles below the
      # median are not usually available.
      if(is.na(init.lmr)) {
         if(verbose) {
            message("Starting parameters not given, attempting to setup ",
                    "from the quantiles")
         }
         if(length(x) < n.para) {
            warning("Length the data values in argument 'x' is less than the ",
                    "number of parameters for the specified distribution ",
                    "type, returning NULL")
            return(NULL)
         }
         init.lmr    <- rep(NA, n.para)
         init.lmr[1] <- x[1] # tread the smallest quantile as the mean, in many
         # applications, x[1] would be the median

         # The mean() will work just fine without the 3rd item being present.
         if(n.para > 1) init.lmr[2] <- mean(x[2:3], na.rm=TRUE) - x[1]
         if(n.para > 2) init.lmr[3:n.para] <- 0

         # The L-moments are now formed into the lmomco object and inherently
         # checked for validity. The disfitqua() function assumes untrimmed
         # versions if the L-moments are hacked from the quantiles.
         try( init.lmr <- vec2lmom(c(init.lmr)) )
         if(is.na(init.lmr$lambdas[1])) {
            stop("FATAL: The initial L-moments (untrimmed version) deduced ",
                 "(assumed) from the quantiles are invalid through the ",
                 "vec2lmom() function, try providing alternative L-moments ",
                 "through the 'init.lmr' function or providing initial ",
                 "parameters through the 'init.para' function.")
            return(NULL)
         }
      }

      # The initial L-moments are used for the initial parameter estimates. The
      # L-moments could either be hacked at using the quantiles themselves or
      # from the staring L-moments.
      init.para <- NULL
      try( init.para <- lmom2par(init.lmr, type=type) )
      if(is.null(init.para)) {
         stop("FATAL: Initial parameters computing from either the deduced ",
              "(assumed) L-moments from the quantiles or given L-moments by ",
              "argument 'init.lmr' do not produce valid parameters for the ",
              "distribution 'type' that has been given.")
         return(NULL)
      }
   }

   if(verbose) message("     Initial parameters: ",
                                          paste(init.para, collapse=", "), "\n")
   para <- init.para.to.archive <- init.para$para

   fncontrol <- list(f=f, x=x, type=type)

   OBJfun <- function(hat, obs) { NULL }
   if(objfun == "rmse") {
      OBJfun <- function(hat, obs) { sqrt( sum(   (obs - hat)^2)/length(obs) ) }
   } else if(objfun == "mad") {
      OBJfun <- function(hat, obs) {       sum(abs(obs - hat)  )/length(obs)   }
   } else {
      stop("should not be here in logic")
   }

   "fn" <- function(mypar) {
      newpar <- NULL
      try( newpar <- list(para=mypar, type=fncontrol$type) )
      if(is.null(newpar)) {
         warning("FATAL: new parameters could not be estimated, ",
                 "stepping the returning error to optim() to Inf")
         return(Inf)
      }
      newXs <- par2qua(fncontrol$f, newpar)

      #print(length(fncontrol$x))
      if(length(newXs) != length(fncontrol$x)) {
         warning("NEAR FATAL?: length of new quantiles does not equal length of ",
                 "the available quantiles, stepping returning error to optim() ",
                 "to Inf")
         return(Inf)
      }
      err <- OBJfun(newXs, fncontrol$x)
      return(err)
   }

   rt <- NULL
   try( rt <- optim(para, fn, ...) )
   if(is.null(rt)) {
      message("optim() call returned NULL, try changing the initial ",
              "parameters or L-moments")
      return(NA)
   }
   para <- rt$par
   fit.para <- vec2par(para, type=type)
   if(verbose) message("     Optimized parameters: ",
                                           paste(fit.para, collapse=", "), "\n")
   fit.para$init.para <- init.para.to.archive
   fit.para$source    <- "disfitqua"
   fit.para$disfitqua <- rt # preserve the operations of optim
   return(fit.para)
}

