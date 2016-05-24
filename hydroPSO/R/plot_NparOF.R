# File plot_NparOF.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2010-2014 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################
#                             'plot_NparOF'                                    #
################################################################################
# Author : Mauricio Zambrano Bigiarini                                         #
# Started: Nov 30th, 2010                                                      #   
# Updates: 17-Jan-2011 ; 28-Jan-2011 ; 09-Mar-2011                             #
#          17-Feb-2012 ; 21-Feb-2012 ; 09-Mar-2012 ; 23-Mar-2012 ; 19-Nov-2012 # 
#          20-Nov-2012 ; 28-Nov-2012                                           #   
#          09-May-2013                                                         #
#          09-Abr-2013                                                         #
################################################################################
# Purpose: For 'n' user-defined parameters, it produces 'sum(1:(npar-1))'      #
#         'plot_2parOF' plots, with the  values of the objective function in   #
#         a 2D box,  where the boundaries of each parameter are used as axis.  #
#         The 'sum(1:(npar-1)) plots corresponds to all the possible           #
#         combinations of 2 parameters among all the 'n' parameters provided   #
################################################################################
# nrows  : numeric, with the amount of rows to be used in the plotting window. 
#          If \code{nrows='auto'} the number of columns is automatically computed 
#          depending on the number of parameters in \code{params}

plot_NparOF <- function(params, 
                        gofs,
                        param.names=colnames(params),
                        MinMax=c(NULL, "min", "max"),
                        beh.thr=NA, 
                        nrows="auto",
                        gof.name="GoF", 
                        main=paste(gof.name, "Surface"),
                        GOFcuts="auto",
                        colorRamp= colorRampPalette(c("darkred", "red", "orange", "yellow", "green", "darkgreen", "cyan")),
                        points.cex=0.7, 
                        alpha=0.65,                       
                        axis.rot=c(0, 0),
                        verbose=TRUE
                        ) {

    
  ##############################################################################
  # 1)                            Checkings                                    #
  ##############################################################################
  
    # Checking 'params'
    if (missing(params)) 
      stop("Missing argument: 'params' must be provided !!" )
      
    # Number of parameter sets
    n <- NROW(params)

    # Checking 'gofs'
    if (missing(gofs)) {
      stop("Missing argument: 'gofs' must be provided !!" )
    } else if (length(gofs) != n)
        stop("Invalid argument: 'length(gofs) != nrow(params)' (", length(gofs), "!=", n, ") !!" )    
        
    # Setting 'MinMax' 
    MinMax <- match.arg(MinMax)

    # Checking 'beh.thr'
    if ( !is.na(beh.thr) ) {  
      if ( is.null(MinMax) )
         stop("Missing argument: 'MinMax' has to be provided before using 'beh.thr' !!")       
      if ( is.null(gofs) )
        stop("Missing argument: 'gofs' has to be provided before using 'beh.thr' !!")
    } # IF end

    # Number of parameters that will be analysed
    npar <- length(param.names)

    # creating the variable that will store the position of the selected parameters within 'params'
    par.pos <- numeric(npar)

    # Checking 'param.names'
    for ( i in 1:npar) {
      if ( !(param.names[i] %in% colnames(params)) )
        stop("Invalid argument: the field '", param.names[i], "' does not exist in 'params'")
      
      par.pos[i] <- which(colnames(params) == param.names[i])
    } # FOR end
  
    
  ##############################################################################
  # 2)                            Computations                                 #
  ##############################################################################
  
    # Filtering out those parameter sets above/below a certain threshold
    if (!is.na(beh.thr)) {  
       # Checking 'beh.thr'
       mx <- max(gofs, na.rm=TRUE)
       if (beh.thr > mx)
         stop("Invalid argument: 'beh.thr' must be lower than ", mx ,"!!")
    
      # Computing the row index of the behavioural parameter sets
      ifelse(MinMax=="min", beh.row.index <- which(gofs <= beh.thr), 
                            beh.row.index <- which(gofs >= beh.thr) )
    
      # Removing non-behavioural parameter sets & gofs
      params <- params[beh.row.index, ]
      gofs   <- gofs[beh.row.index]
   
      # Amount of behavioural parameter sets 
      nbeh <- nrow(params)
      if (verbose) message( "[ Number of behavioural parameter sets: ", nbeh, " ]" )
    } # IF end  



    # If the user didn't provide 'GOFcuts', the 5 quantiles are used
    if (length(GOFcuts) == 1){
      if (GOFcuts=="auto") {
        if (MinMax=="min") { 
           GOFcuts <- unique( quantile( as.numeric(gofs), 
                              probs=c(0, 0.25, 0.5, 0.85, 0.9, 0.97, 1), na.rm=TRUE) )                                          
        } else if (MinMax=="max") {
            GOFcuts <- unique( quantile( as.numeric(gofs), 
                               probs=c(0, 0.03, 0.1, 0.15, 0.5, 0.75, 1), na.rm=TRUE) ) 
          } else  # MinMax==NULL
               GOFcuts <- unique( quantile( as.numeric(gofs), 
                                  probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm=TRUE) )
        
        suppressWarnings(                          
          if (verbose) message( "[ Computed GOFcuts: ", 
                       paste(as.numeric(formatC( GOFcuts, format="E", digits=3, flag=" ")), collapse=" "), " ]" ) 
        )                         
          
      } # IF end
    } # IF end
    
  ##############################################################################
  # 3)                            Plotting                                     #
  ##############################################################################  
   
    # Number of plots that will be drawn   
    nplots <- sum(1:(npar-1))

    plots <- vector("list", nplots)
    
    pos <- 1
    for ( i in 1:(npar-1) ) {
      for (j in ((i+1):npar) ) {
        if (verbose) message("[ Plotting '", param.names[i], "' vs '", param.names[j], "' ]")
        plots[[pos]] <- plot_2parOF(params=params, gofs=gofs, MinMax=MinMax, p1.name=param.names[i], 
                                     p2.name=param.names[j], gof.name=gof.name,
                                     type="sp", main="", GOFcuts=GOFcuts, colorRamp=colorRamp, 
                                     alpha=alpha, axis.rot=axis.rot, auto.key=FALSE, 
                                     points.cex=points.cex )
                                     
        pos <- pos + 1
      } # FOR j end
    } # For i end

    # Computing the number of rows for the plot 
    nplots <- nplots + 1
    if (nrows == "auto") {
      if ( nplots <= 5 )                   lnr <- 1
      if ( (nplots > 5) & (nplots <= 14) ) lnr <- 2
      if ( nplots > 14 )                   lnr <- ceiling(nplots/7)
    } else lnr <- nrows          
    
    # Defining the plotting window
    nr <- lnr
    nc <- ceiling(nplots/lnr)
    #par(oma=c(1,1,3,0))
    
    pos <- 1
    for (row in 1:nr) {
      for (col in 1:nc) {
        if (pos <= nplots) {
          if ( ( (row==nr) & (col==nc) ) | (pos==nplots) ) {
            #print(plots[[pos]], split = c(col, row, nc, nr), more = FALSE)
          } else print(plots[[pos]], split = c(col, row, nc, nr), more = TRUE)
        } # IF end
        pos <- pos + 1
      } # FOR end
    } # FOR end

    # Drawing the legend, with a dummy empty plot
    #gof.levels <- cut(gofs, GOFcuts)
    suppressWarnings(
      gof.levels <- cut(gofs, unique(as.numeric(formatC( GOFcuts, format="E", digits=4, flag=" "))))
    )
    nlevels    <- length(levels(gof.levels)) 
    
    #require(grid)
    a <- lattice::xyplot(1~1, 
                groups=gof.levels,
                type="n", xlab="", ylab="", scales=list(draw=FALSE),
                key = list(x = .5, y = .5, corner = c(0.5, 0.5),
                           title=gof.name,
                           points = list(pch=16, col=colorRamp(nlevels), cex=1.5),
                           text = list(levels(gof.levels), cex=0.8)                              # MZB: cex=0.8=leg.cex
                           #text = list(formatC( as.numeric(levels(gof.levels)), format="E", digits=2, flag=" "))                     
                           ),
                # removing outter box. From: https://stat.ethz.ch/pipermail/r-help/2007-September/140098.html
                par.settings = list(axis.line = list(col = "transparent")),
                axis = function(side, ...) {
                    lattice::axis.default(side = side, ...)
                }
             )
    col <- nc -(nc*nr-nplots)
    print(a, split = c(col, nr, nc, nr), more = FALSE)

} # 'plot_NparOF' END
