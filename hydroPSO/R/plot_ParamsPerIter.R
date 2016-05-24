# File plot_ParamsPerIter.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2011-2012 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################
#                          'plot_ParamsPerIter'                                #
################################################################################
# Purpose: It plots the value of each parameter and the objective              #
# function against the Number of Model Evaluations                             #

plot_ParamsPerIter <- function(params, ...) UseMethod("plot_ParamsPerIter")
 
################################################################################
# Author : Mauricio Zambrano-Bigarini                                          #
################################################################################
# Started: 30-Nov-2010 at JRC, Ispra                                           # 
# Updates: 13-Oct-2011                                                         #
#          21-Feb-2012 ; 23-Mar-2012                                           #
################################################################################
plot_ParamsPerIter.default <- function(params, 
                                       param.names=colnames(params),
                                       main=NULL,
                                       xlab="Number of evaluations",
                                       nrows="auto",
                                       cex=0.5,
                                       cex.main=1.2,
                                       cex.axis=1.7,
                                       cex.lab=1.5,
                                       col=rainbow(ncol(params)),
                                       lty=3,
                                       verbose=TRUE,
                                       ...,                                       
                                       #### PNG options ### 
                                       do.png=FALSE,
                                       png.width=1500,
                                       png.height=900,
                                       png.res=90,
                                       png.fname="Params_ValuePerRun.png"  
                                       ) {    

    # number of parameters
    nparam <- ncol(params)
    
    # Checking 'param.names'
    if (!identical(param.names, colnames(params))) {
      if (length(param.names) != nparam) {
        stop(paste("Invalid argument: 'length(param.names) = ", length(param.names), " != ", nparam, " = nparam'", sep=""))
      } else colnames(params) <- param.names
    } # IF end

    # Number of parameter sets
    nparsets <- nrow(params)
    
    ############################################################################
    # 2)                            Plotting                                   #
    ############################################################################  
    if (verbose) message( "                                        ")  
    if (verbose) message( "[            Plotting ...              ]")  
  
    if (do.png) png(filename=png.fname, width=png.width, height=png.height, res=png.res)
    
    # Computing the number of rows for the plot 
    if (nrows == "auto") {
      if ( nparam <= 5 )                   lnr <- 1
      if ( (nparam > 5) & (nparam <= 14) ) lnr <- 2
      if ( nparam > 14 )                   lnr <- ceiling(nparam/7)
    } else lnr <- nrows      
    
    # Saving default plotting parameters
    old.par <- par(no.readonly=TRUE)
    if (!do.png) on.exit(par(old.par))
    
    # Defining the plotting window
    nr <- lnr
    nc <- ceiling(nparam/lnr)
    par(mfrow=c(nr,nc))   
    par(mar=c(5,4.5,1,2)+0.1) # Default: par(mar=c(5,4,4,2)+0.1)
    if (!is.null(main)) par(oma=c(1,1,3,0))
    
    for ( i in 1:nparam )
      plot(1:nparsets, params[,i], type="o", lty=lty, col=col[i], 
           cex.main=cex.main, cex=cex, cex.lab=cex.lab, cex.axis=cex.axis, font.lab=2, 
           xlim=c(1,nparsets), ylim=range(params[,i], na.rm=TRUE), 
           xlab=xlab, ylab=param.names[i] )
 
    # Adding a main title for the plot
    if (!is.null(main)) mtext(main, side=3, line=1, cex=cex.main, outer=TRUE)
    
    if (do.png) dev.off()
} # 'plot_ParamsPerIter.default' END


###################################################################
# Author : Mauricio Zambrano Bigarini                             #
# Started: 13-Oct-2011 at JRC, Ispra                              # 
# Updates: 13-Oct-2011                                            #
###################################################################
plot_ParamsPerIter.data.frame <- function(params, 
                                  param.names=colnames(params),
                                  main=NULL,
                                  xlab="Number of evaluations",
                                  nrows="auto",
                                  cex=0.5,
                                  cex.main=1.2,
                                  cex.axis=1.7,
                                  cex.lab=1.5,
                                  col=rainbow(ncol(params)),
                                  lty=3,                                  
                                  verbose=TRUE,
                                  ...,
                                  #### PNG options ### 
                                  do.png=FALSE,
                                  png.width=1500,
                                  png.height=900,
                                  png.res=90,
                                  png.fname="Params_ValuePerRun.png"  
                                  ) {
    
  #params <- as.matrix(params)
 
  plot_ParamsPerIter.default(params=params, 
                    param.names=param.names,
                    main=main,
                    xlab=xlab,
                    nrows=nrows,  
                    cex=cex,
                    cex.main=cex.main, 
                    cex.axis=cex.axis, 
                    cex.lab=cex.lab,
                    col=col,
                    lty=lty,
                    verbose=verbose,
                    ...,
                    #### PNG options ### 
                    do.png=do.png,
                    png.width=png.width,
                    png.height=png.height,
                    png.res=png.res,
                    png.fname=png.fname
                    )                                    
} # 'plot_ParamsPerIter.data.frame' END


###################################################################
# Author : Mauricio Zambrano Bigarini                             #
# Started: 20-Oct-2011 at JRC, Ispra                              # 
# Updates: 20-Oct-2011                                            #
###################################################################
plot_ParamsPerIter.matrix <- function(params, 
                                  param.names=colnames(params),
                                  main=NULL,
                                  xlab="Number of evaluations",
                                  nrows="auto",
                                  cex=0.5,
                                  cex.main=1.2,
                                  cex.axis=1.7,
                                  cex.lab=1.5,
                                  col=rainbow(ncol(params)),
                                  lty=3,                                                                   
                                  verbose=TRUE,
                                  ...,
                                  #### PNG options ### 
                                  do.png=FALSE,
                                  png.width=1500,
                                  png.height=900,
                                  png.res=90,
                                  png.fname="Params_ValuePerRun.png"
                                  ) {
    
  plot_ParamsPerIter.default(params=params, 
                    param.names=param.names,
                    main=main,
                    xlab=xlab,
                    nrows=nrows,  
                    cex=cex,
                    cex.main=cex.main, 
                    cex.axis=cex.axis, 
                    cex.lab=cex.lab,
                    col=col,
                    lty=lty,                    
                    verbose=verbose,
                    ...,
                    #### PNG options ### 
                    do.png=do.png,
                    png.width=png.width,
                    png.height=png.height,
                    png.res=png.res,
                    png.fname=png.fname
                    )                                    
} # 'plot_ParamsPerIter.matrix' END
