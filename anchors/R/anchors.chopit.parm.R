#######################################################################
##
## Function: anchors.chopit.parm()
## Author  : Jonathan Wand <wand(at)stanford.edu>
## Created :  2008-04-20
##  
## DESCRIPTION: Function for analyzing data with anchoring vignettes
##
## Extracted and refined from former chopit()  
#######################################################################
anchors.chopit.parm <- function( data, count, options ) {

  if (options$debug>2) { cat("anchors.chopit.parm: data\n"); print(data) }

  debug <- options$debug

  if (options$linear) {
    default.gamma1.constant <- options$defparm.lin.gamma.constant     
    default.gamma.constant  <- options$defparm.lin.gamma.constant 
    default.gamma           <- options$defparm.lin.gamma 
    default.se              <- options$defparm.lin.se    
  } else {
    default.gamma1.constant <- options$defparm.non.gamma.constant     
    default.gamma.constant  <- options$defparm.non.gamma.constant 
    default.gamma           <- options$defparm.non.gamma 
    default.se              <- options$defparm.non.se    
  }

  
  vnames <- list(x0 = colnames(data$x0) , ## this should be NULL if no.self=TRUE
                 y0 = colnames(data$y0) , ## this should be NULL if no.self=TRUE
                 v0 = colnames(data$v0v ),
                 v1 = colnames(data$v0v1),
                 z0 = colnames(data$z0) )
  

  if (count$n.tau.set > 1) {
    init.gamma1n<- c(paste("gamma1.cut",
                           rep(1, count$nvars.gamma1),
                           ".", vnames$v1,
                           sep=""))
    init.gamman <- c(paste("gamma",c(matrix( 1:count$n.tau.set,ncol=count$n.tau.set,
                                            nrow=(count$n.cat-2)*count$nvars.gamma,byrow=TRUE)),
                           ".cut",
                           rep(c(matrix( 2:(count$n.cat-1), ncol=(count$n.cat-2),
                                        nrow=count$nvars.gamma,byrow=TRUE)),count$n.tau.set),
                           ".", vnames$v0,
                           sep=""))
  } else {

    init.gamma1n<- c(paste("gamma.cut",
                           rep(1, count$nvars.gamma1),
                           ".", vnames$v1,
                           sep=""))
    init.gamman <- c(paste(rep("gamma", (count$n.cat-2)*count$nvars.gamma),
                           ".cut",
                           c(matrix( 2:(count$n.cat-1),
                                        ncol=(count$n.cat-2),
                                        nrow=count$nvars.gamma,byrow=TRUE)),
                           ".", vnames$v0,
                           sep=""))
    
  }
  
  if (debug > 0) {
    cat("TEST: Dumping contents of init.gamman:\n")
    print(as.matrix(init.gamma1n))
    print(as.matrix(init.gamman))
  }
  if (debug > 0)
    cat("CHECKPOINT: set initial values of normal distribution variances\n")

  ## WITH only one set of vign
  tmp0 <- vnames$z0
  s.sigma.vign <- rep(default.se,count$n.vign)
  l.sigma.vign <- paste("sigma.",tmp0,sep="")
  e.sigma.vign <- rep(TRUE ,count$n.vign)

  
  if (debug > 0) cat("CHECKPOINT: set initial values of parameters\n")
  
  parm <- list(
               ## list of DGP parameters (or starting values)
               start = list(
                 ## gamma START
                 gamma1      = rep(default.gamma,length(init.gamma1n)),
                 gamma      = rep(default.gamma,length(init.gamman)),
                 sigma.re = ifelse( options$random , default.se, default.se),
                 sigma.self  = default.se, ## rep(0,count$n.self), --> at some point flex 
                 sigma.vign  = s.sigma.vign
                 ),
               ## names for parameters... (automated)
               labels = list(
                 gamma1     = init.gamma1n,
                 gamma      = init.gamman,
                 sigma.re = "sigma.random.effect",
                 ##paste("sigma.self",1:count$n.self    ,sep=""),
                 sigma.self  = "sigma.self",
                 sigma.vign  = l.sigma.vign
                 ),
               ## which parms are fixed (FALSE) OR estimated (TRUE)
               estimated = list(
                 gamma1     = rep(TRUE,length(init.gamma1n)),
                 gamma      = rep(TRUE,length(init.gamman)),
                 sigma.re = options$random,
                 sigma.self  = FALSE, #rep(FALSE,count$n.self)
                 sigma.vign  = e.sigma.vign
                 )
               )
  ## Add theta to parm list
  if (debug > 0) 
    cat("CHECKPOINT: Adding theta to parmeter list\n")

  ## Single set of vign assumed here, but keep loop for simplicity
  ## -> z0 used to make tmp0 directly
  for (i in 1:count$n.vign.set) {
    zti <- ifelse( count$n.vign.set > 1, i, "")
    txt1 <- paste("parm$start$theta",zti," <- rep(0.0,",count$n.vign[i],")",sep="")
                                        #    tmp0 <- as.character(eval(parse(text=paste("alphad$z",i,sep=""))))
    tmp0 <- vnames$z0
    tmp1  <- paste("theta",zti,".",tmp0,sep="",collapse="\",\"")
    txt2 <- paste("parm$labels$theta",zti," <- c(\"",tmp1,"\")",sep="")
    txt3 <- paste("parm$estimated$theta",zti," <- rep(TRUE,",count$n.vign[i],")",sep="")
    if (debug > 0) {
      cat("TEST: dumping build strings for theta\n")
      print(txt1)
      print(txt2)
      print(txt3)
    }
    eval(parse(text=txt1))
    eval(parse(text=txt2))
    eval(parse(text=txt3))
  }
  ##
  if (count$n.self > 0) {
    parm$start$beta        <- rep( options$defparm.beta ,length(vnames$x0))
    parm$labels$beta       <- paste("beta.",vnames$x0,sep="")
    parm$estimated$beta    <- rep(TRUE,length(vnames$x0))
  }

  ## cannot estimate variance if
  if (count$n.self == 0) {
    if (!options$silence) cat("\nWARNING! No self-question selected\n")
    parm$estimated$sigma.vign[1]  <- FALSE # rep(FALSE ,count$n.vign.set)
    parm$estimated$theta[1] <- FALSE
  }

  #############################################################
  ## 
  ## A.3 CUSTOMIZE DATA and PARAMETERS
  ## 
  #############################################################  
  if (debug > 0) {
    cat("CHECKPOINT: customize data and parameters\n")
  }
  
  ## Find the column of ones and zero out the intercept in beta!
  idx.constants.beta <- NULL
  if (count$n.self > 0) {
    for (i in 1:length(vnames$x0)) {
      idx.constants.beta <- c(idx.constants.beta,
                              ifelse( length(unique(data$x0[,i])) == 1, i, -i)
                              )
    }
    if (all(idx.constants.beta < 0)) {
      if (options$verbose)
        cat("\nFYI! No column of data in x0 looks like a constant\n")
                                        #cat(  "WARNING! therefore no beta is set to zero\n\n")
    } else if ( sum(idx.constants.beta>0) > 1) {
      cat("\nERROR! Multiple columns are constants in x0:",idx.constants.beta,"\n")
      cat(  "ERROR! This will result in a short rank model of mu\n")
      cat(  "ERROR! Solution: omit all but one of the above columns from x0\n")
      cat(  "ERROR! Program is quitting now.\n\n")
      stop("")
    }
    
    if (!options$silence) cat("\nchopit() will be identified/normalized by:\n")
    if (options$normalize=="self") {
      if ( sum(idx.constants.beta>0) == 1) {
        ii <- idx.constants.beta[idx.constants.beta>0]
        if (!options$silence) cat("1. omitting the intercept in mu, (constant set to zero is",
                                  parm$labels$beta[ii],")\n")
        parm$start$beta[ii] <- 0
        parm$estimated$beta[ii] <- FALSE
      } else {
        if (!options$silence) cat("1. no intercept is estimated in mu\n")
      }
      if (!options$silence) cat("2. setting variance of first self question to 1\n\n")
      parm$start$sigma.self[1] <- default.se
      parm$estimated$sigma.self[1] <- FALSE 
      if (debug > 0)
        cat("actual value",parm$start$sigma.self[1],"\n")
    } else if (options$normalize=="hilo") {
      if (!options$silence) cat("1. setting first theta in first vignette to 0\n")
      if (!options$silence) cat("2. setting last theta in first vignette to 1\n\n")

      parm$start$theta[1] <- 0
      parm$start$theta[length(parm$start$theta)] <- 1
      parm$estimated$theta[1] <- FALSE
      parm$estimated$theta[length(parm$estimated$theta)] <- FALSE
      
      if ( sum(idx.constants.beta>0) == 1) {
        ii <- idx.constants.beta[idx.constants.beta>0]
        parm$estimated$beta[ii] <- TRUE
      }
      parm$estimated$sigma.self[1] <- TRUE

    } else if (options$normalize=="vign") {

      parm$start$theta[1] <- 0
      parm$start$sigma.vign[1] <- default.se
      parm$estimated$theta[1] <- FALSE
      parm$estimated$sigma.vign[1] <- FALSE

      if (!options$silence) cat("1. setting MEAN of first vignette to 0\n")
      if (!options$silence) cat("2. setting VARIANCE of first vignette to 1\n\n")
      
      if ( sum(idx.constants.beta>0) == 1) {
        ii <- idx.constants.beta[idx.constants.beta>0]
        parm$estimated$beta[ii] <- TRUE
      }
      parm$estimated$sigma.self[1] <- TRUE

    } else if (!is.null(options$estimated)) {
      if (!options$silence) cat("1. by user's custom setting (normalize option is ignored)\n\n")
      if ( sum(idx.constants.beta>0) == 1) {
        ii <- idx.constants.beta[idx.constants.beta>0]
        parm$estimated$beta[ii] <- TRUE
      }
      parm$estimated$sigma.self[1] <- TRUE

    } else {
      stop(paste("No valid 'normalize' or 'estimated' list option has been passed to chopit\n"))
    }
  }

  if (debug > 0)
    cat("CHECKPOINT: find constant in gamma\n")
  ## Find the column of ones and -1 the intercept(s) in gamma
  ## JW: currently uses vignette subset to determine constants
  idx.constants.gamma <- NULL
  for (i in 1:count$nvars.gamma) {
    idx.constants.gamma <- c(idx.constants.gamma,
                             ifelse( length(unique(data$v0v[,i])) == 1, i, -i)
                             )
  }
  if (all(idx.constants.gamma < 0)) {
    if (!options$silence) cat("\nWARNING! No column of data in v0v looks like a constant\n\n")
    no.gamma.constant <- TRUE
  } else if ( sum( idx.constants.gamma > 0) > 1 ) {
    cat("\nERROR! Multiple columns are constants in v0v:",idx.constants.gamma,"\n")
    cat(  "ERROR! This will result in a short rank model of the taus\n")
    cat(  "ERROR! Solution: omit all but one of the above columns from v0\n")
    cat(  "ERROR! Program is quitting now.\n\n")
    stop("")
  } else {
    if (debug > 0)
      cat("no.gamma.constant == FALSE\n")
    no.gamma.constant <- FALSE
    ii <- seq( idx.constants.gamma[idx.constants.gamma>0], length(init.gamman), count$nvars.gamma )
    ntmp <- sum(ii)
    tmp  <- 1:ntmp      
    parm$start$gamma[ii]  <- default.gamma.constant # (tmp-1)/ntmp
    if (debug > 0) {
      cat("SET constant gamma",default.gamma.constant,"\n")
      print(ii)
      print(parm$start$gamma)
    }
    fix.gamma.oprobit <- rep(FALSE,length(init.gamman))
    fix.gamma.oprobit[ii] <- TRUE
  }

  ## now do gamma1
  idx.constants.gamma1 <- NULL
  for (i in 1:count$nvars.gamma1) {
    idx.constants.gamma1 <- c(idx.constants.gamma1,
                             ifelse( length(unique(data$v0v1[,i])) == 1, i, -i)
                             )
  }
  if (all(idx.constants.gamma1 < 0)) {
    if (!options$silence) cat("\nWARNING! No column of data in v0v1 looks like a constant\n\n")
    no.gamma1.constant <- TRUE
  } else if ( sum( idx.constants.gamma1 > 0) > 1 ) {
    cat("\nERROR! Multiple columns are constants in v0v1:",idx.constants.gamma1,"\n")
    cat(  "ERROR! This will result in a short rank model of the taus\n")
    cat(  "ERROR! Solution: omit all but one of the above columns from v0\n")
    cat(  "ERROR! Program is quitting now.\n\n")
    stop("")
  } else {
    if (debug > 0)
      cat("no.gamma1.constant == FALSE\n")
    no.gamma1.constant <- FALSE
    ii <- seq( idx.constants.gamma1[idx.constants.gamma1>0], length(init.gamma1n), count$nvars.gamma1 )
    ntmp <- sum(ii)
    tmp  <- 1:ntmp      
    parm$start$gamma1[ii]  <- default.gamma1.constant # (tmp-1)/ntmp
    if (debug > 0) {
      cat("SET constant gamma1",default.gamma1.constant,"\n")
      print(ii)
      print(parm$start$gamma1)
    }
  }

  
  ## OVERRIDE DEFAULTS HERE!
  if (!is.null(options$start$start)) {
    if (!options$silence) cat("\n anchors.chopit.parm(): Using user specified starting values\n")
    parm$start <- replace.list(parm$start,options$start$start)
  }
  if (!is.null(options$start$labels)) {
    if (!options$silence) cat("\n anchors.chopit.parm(): Using user specified labels\n")
    parm$labels <- replace.list(parm$labels,options$start$labels)
  }
  if (!is.null(options$start$estimated)) {
    if (!options$silence) cat("\n anchors.chopit.parm(): Using user specified list of identification restrictions\n")
    parm$estimated <- replace.list(parm$estimated,options$start$estimated)
  }
  if (debug > 0)
    cat("CHECKPOINT: done overriding defaults\n")

#  ## OVERRIDE DEFAULTS HERE!
#  if (!is.null(options$start$start)) {
#    if (!options$silence) cat("\n anchors.chopit.parm(): Using user specified starting values\n")
#    parm$start <- options$start$start
#  }
#  if (!is.null(options$start$labels)) {
#    if (!options$silence) cat("\n anchors.chopit.parm(): Using user specified labels\n")
#    parm$labels <- options$start$labels
#  }
#  if (!is.null(options$start$estimated)) {
#    if (!options$silence) cat("\n anchors.chopit.parm(): Using user specified list of identification restrictions\n")
#    parm$estimated <- options$start$estimated
#  }
#  if (debug > 0)
#    cat("CHECKPOINT: done overriding defaults\n")

  
  ## NOW invoke options
  if (options$vign.var == "homo") {
    if (count$n.vign.set == 1) {
      if (!options$silence) cat("Forcing vignettes to have common variance parameter because vign.var='homo'\n")
      parm$start$sigma.vign     <- parm$start$sigma.vign[1]     
      parm$estimated$sigma.vign <- parm$estimated$sigma.vign[1]     
      parm$labels$sigma.vign    <- parm$labels$sigma.vign[1]
    } else {
      warning("Cannot force multiple sets of vignettes to have common variance parameter")
    }
  }

  if (!is.null(options$vign.map))
    parm$vign.map  <- options$vign.map
  else
    parm$vign.map  <- 1:count$n.vign.set

  ## 
  if (debug > 0) {
    cat("TEST: dumping pre-packed 'parm' list:\n")
    print(parm)
  }
  
  ## pack stuff
  parm <- packv(parm)
  if (is.null(parm)) stop("Packing of starting values failed!\n")
  parm$var.names <- vnames

  class(parm) <- "anchors.chopit.parm"

  if (debug > 0) {
    cat("TEST: dumping 'parm' list:\n")
    print(parm)
  }
  
  return(parm)
}
