# File read.param.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2008-2011 Mauricio Zambrano-Bigiarini
# Distributed under GPL 2 or later

################################################################################
#                                'plot_params'                                 #
################################################################################
# Author : Mauricio Zambrano-Bigiarini                                         #  
################################################################################
# Started: 28-June-2010,                                                       #
# Updates: 16-May-2011  ; 24-Jun-2011 ; 01-Jul-2011 ; 01-Sep-2011              #
#          19-Jan-2012  ; 02-Feb-2012 ; 15-Feb-2012 ; 07-Mar-2012 ; 23-Mar-2012#        
################################################################################
# This function makes dotty plots of different parameter values vs the 
# corresponding objective function value (usually for plotting the 
# efficiencies of different parameter sets)

# params : data.frame, whose columns represents the behavioural parameter 
#          sets and the goodness of fit of each one of them 
# of.name: character, with the name of the column that contains the values 
#          of the objective function for each parameter set
# beh.thr: OPTIONAL, only used when 'plot=TRUE'. \cr
#          numeric, with a behavioural threshold to be used for ploting 
#          a horizontal line 
# beh.col: OPTIONAL, only used when 'plot=TRUE'. \cr
#          color to be used for plotting the horizontal line on 'beh.thr'
# beh.lty: OPTIONAL, only used when 'plot=TRUE'. \cr
#          type of line to be used for plotting the horizontal line on 'beh.thr'
# beh.lwd: OPTIONAL, only used when 'plot=TRUE'. \cr
#          line width to be used for plotting the horizontal line on 'beh.thr'
# nrows   : OPTIONAL, only used when 'plot=TRUE'. \cr
#          Number of rows to be used in the plotting window. The number 
#          of columns is automatically computed depending on the number 
#          of columns of 'params'

plot_params <- function(params, ...) UseMethod("plot_params")

plot_params.default <- function(params, 
                        gofs=NULL,
                        ptype=c("histogram", "dottyplot", "boxplot", "vioplot", "pairs"),
                        param.cols=1:ncol(params),
                        param.names=colnames(params),
                        #of.col=NULL, 
                        of.name="GoF", 
                        MinMax=NULL, 
                        beh.thr=NA, 
                        beh.col="red", 
                        beh.lty=1, 
                        beh.lwd=2, 
                        nrows="auto",
                        #col="black", 
                        col="#00000030",
                        ylab=of.name, 
                        main=NULL,
                        pch=19,  # Only used for dotty plots
                        cex=0.5, 
                        cex.main=1.5,
                        cex.axis=1.5,
                        cex.lab=1.5,
                        breaks="Scott",
                        freq=TRUE,
                        verbose=TRUE,
                        ...,
                        #### PNG options ### 
                        do.png=FALSE,
                        png.width=1500,
                        png.height=900,
                        png.res=90,
                        png.fname="Parameters.png"
                        ) {
                      
  ##############################################################################
  # 1)                            Checkings                                    #
  ##############################################################################
  
  # Setting 'ptype' 
  ptype <- match.arg(ptype)  
  
  # Checking 'vioplot' 
  if (ptype=="vioplot") {  
     if ( !require(vioplot) ) {
        warning("Package 'vioplot' is not installed =>  ptype='boxplot'")
        ptype <- "boxplot"
     } # IF end 
  } # IF end
         
  # Checking 'pairs' 
  if (ptype=="pairs") {  
     if ( !require(hydroTSM) ) {
        warning("Package 'hydroTSM' is not installed =>  ptype='histogram'")
        ptype <- "histogram"
     } # IF end 
  } # IF end
  
  # Checking 'beh.thr'
  if ( !is.na(beh.thr) ) {
     if ( is.null(MinMax) )
        stop("Missing argument: 'MinMax' has to be provided before using 'beh.thr' !!")        
      if ( is.null(gofs) )
        stop("Missing argument: 'gofs' has to be provided before using 'beh.thr' !!")
  } # IF end
         
  # Checking 'MinMax'
  if ( !is.null(MinMax) ) {
     if ( !(MinMax %in% c("min", "max")) )
        stop("Invalid argument: 'MinMax' must be in c('min', 'max')")
  } # IF end

  # Checking 'gofs' for dotty plots
  if ( ptype == "dottyplot" )
    if (is.null(gofs)) stop("Invalid argument: 'gofs' must be provided for dotty plots !!")

  # Checking gofs' length
  if (!is.null(gofs)) {
    if ( length(gofs) != nrow(params) ) 
      stop("Invalid argument: 'length(gofs) != nrow(params)' ", length(gofs), "!=", nrow(params), " !!")
  } #IF end 
                
  if (class(params) != "data.frame") params <- as.data.frame(params)

  # Checking that 'param.cols' and 'param.names' have the same length
  if ( length(param.cols) != length(param.names) )
    stop( paste("Invalid argument: length(param.cols) != length(param.names)' (", 
                length(param.cols), "!=", length(param.names), ")", sep="") ) 
                
  ##############################################################################
  # 2)                            Computations                                 #
  ##############################################################################
     
  # Keeping only the columns with parameter values
  params <- params[, param.cols] 

  # computing the number of parameters
  nparams <- NCOL(params)
  
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
  
  ##############################################################################
  # 3)                            Plotting                                     #
  ##############################################################################  
  # Plotting ALL the PARAMETER SETS
  if (verbose) message( "                                        ")  
  if (verbose) message( "[            Plotting ...              ]")  
  
  if (do.png) png(filename=png.fname, width=png.width, height=png.height, res=png.res)
  
  # Computing the number of rows for the plot 
  if (nrows == "auto") {
    if ( nparams <= 5 )                    lnr <- 1
    if ( (nparams > 5) & (nparams <= 14) ) lnr <- 2
    if ( nparams > 14 )                    lnr <- ceiling(nparams/7)
  } else lnr <- nrows 

  ncol <- ceiling( nparams/lnr )
  
  # Computing the position of the 'optimum' parameter set
  if ( !is.null(MinMax) ) 
    ifelse(MinMax=="min", best.index <- which.min(gofs), best.index <- which.max(gofs) )
    
  # Saving default plotting parameters
  old.par <- par(no.readonly=TRUE)
  if (!do.png) on.exit(par(old.par)) 

  # Plotting the distribution of all the parameter sampled with LH, not only the behavioural ones
  MinMax.colour <- "coral"
  par(mfrow=c(lnr, ncol))
  par(mar=c(5,4.5,1,2)+0.1) # Default: par(mar=c(5,4,4,2)+0.1)
  if (!is.null(main)) par(oma=c(1,1,3,0))
  if (ptype != "pairs") {
    for ( i in 1:nparams ) {
  
      if (ptype == "dottyplot") {
    
        plot(params[,i], gofs, type="p", col=col, cex=cex, cex.main=cex.main, 
             cex.axis=cex.axis, cex.lab=cex.lab, pch=pch, font.lab=2,
             xlab=colnames(params)[i], ylab=ylab )
    
        # If the user provided a behavioural threshold: 
        if ( !is.na(beh.thr) ) {      
          # Plotting an horizontal line at the behavioural threshold
          abline(h=beh.thr, col=beh.col, lwd=beh.lwd, lty=beh.lty)      
          # Showing the value of the behavioural threshold
          axis(side=4, at = beh.thr, labels=round(beh.thr,2), col=beh.col, 
               cex.axis=cex.axis, cex.lab=cex.lab )      
        } # IF end
      
        # If the user provided 'MinMax': 
        if ( !is.null(MinMax) ) {
           abline(h=gofs[best.index], col=MinMax.colour, lwd=1) 
           xpos <- c(params[best.index,i], params[best.index,i])
           ypos <- range(gofs, na.rm=TRUE)
           lines(xpos, ypos, col=MinMax.colour, lwd=1)  
        } # IF end
      
      } else if (ptype == "histogram") {
    
          h <- hist(params[,i], breaks=breaks, plot=FALSE)
          ifelse( freq==TRUE, yvals <- h$counts, yvals <- h$density )
          ylim <- c(0, max(yvals, pretty(yvals)) )
          hist(params[,i], col=col, cex=cex, cex.main=cex.main, 
               cex.axis=cex.axis, cex.lab=cex.lab, font.lab=2, 
               xlab=colnames(params)[i], main="",
               breaks=breaks, freq=freq, yaxt="n")
          axis(side=2, at = pretty(yvals), cex.axis=cex.axis, 
               cex.lab=cex.lab, ylim=ylim ) 
             
          # If the user provided 'MinMax': 
          if ( !is.null(MinMax) ) {
            xpos <- c(params[best.index,i], params[best.index,i])
            ypos <- ylim
            lines(xpos, ypos, col=MinMax.colour, lwd=0.5)
          } # IF end
      
        } else if (ptype == "boxplot") {
    
             boxplot(params[,i], col="lightblue", cex=cex, cex.main=cex.main, 
                     cex.axis=cex.axis, cex.lab=cex.lab, font.lab=2,
                     xlab=colnames(params)[i], ylab=ylab )
                   
             # If the user provided 'MinMax': 
             if ( !is.null(MinMax) ) 
                abline(h=params[best.index,i], col=MinMax.colour)
      
          } else if (ptype == "vioplot") {
          
               if ( length(unique(params[,i])) != 1 ) {
                 require(vioplot)
                 vioplot(params[,i], col="lightblue", names=colnames(params)[i])
               
                 # If the user provided 'MinMax': 
                 if ( !is.null(MinMax) ) abline(h=params[best.index,i], col=MinMax.colour)
               } else message("Vioplot is not possible: all elements of '", 
                              colnames(params)[i], "' are equal !")
        
           
            } # ELSE end  
      
    } # FOR end
     
  } else { # if (ptype == "pairs")          
      if ( !is.null(gofs) ) {
        params                         <- cbind(params, gofs) 
        colnames(params)[ncol(params)] <- of.name
        } # IF end
             
      library(hydroTSM)
      hydropairs(params)             
     }  # ELSE end   
  
  
  if (!is.null(main)) mtext(main, side=3, line=1, cex=cex.main, outer=TRUE)
  
  if (do.png) dev.off()
    
}  # 'plot_params.default' END


################################################################################
# Author : Mauricio Zambrano-Bigiarini                                         #  
################################################################################
# Started: 23-Feb-2012                                                         #
# Updates:                                                                     #
################################################################################
plot_params.data.frame <- function(params, 
                        gofs=NULL,
                        ptype=c("histogram", "dottyplot", "boxplot", "vioplot", "pairs"),
                        param.cols=1:ncol(params),
                        param.names=colnames(params),
                        #of.col=NULL, 
                        of.name="GoF", 
                        MinMax=NULL, 
                        beh.thr=NA, 
                        beh.col="red", 
                        beh.lty=1, 
                        beh.lwd=2, 
                        nrows="auto",
                        #col="black", 
                        col="#00000030",
                        ylab=of.name, 
                        main=NULL,
                        pch=19,  # Only used for dotty plots
                        cex=0.5, 
                        cex.main=1.5,
                        cex.axis=1.5,
                        cex.lab=1.5,
                        breaks="Scott",
                        freq=TRUE,
                        verbose=TRUE,
                        ...,
                        #### PNG options ### 
                        do.png=FALSE,
                        png.width=1500,
                        png.height=900,
                        png.res=90,
                        png.fname="Parameters.png"
                        ) {
                        
    plot_params.default(params=params, 
                        gofs=gofs,
                        ptype=ptype,
                        param.cols=param.cols,
                        param.names=param.names,
                        #of.col=NULL, 
                        of.name=of.name, 
                        MinMax=MinMax, 
                        beh.thr=beh.thr, 
                        beh.col=beh.col, 
                        beh.lty=beh.lty, 
                        beh.lwd=beh.lwd, 
                        nrows=nrows,
                        #col="black", 
                        col=col,
                        ylab=ylab, 
                        main=main,
                        pch=pch,  # Only used for dotty plots
                        cex=cex, 
                        cex.main=cex.main,
                        cex.axis=cex.axis,
                        cex.lab=cex.lab,
                        breaks=breaks,
                        freq=freq,
                        verbose=verbose,
                        ...,
                        #### PNG options ### 
                        do.png=do.png,
                        png.width=png.width,
                        png.height=png.height,
                        png.res=png.res,
                        png.fname=png.fname
                        )
                        
} # 'plot_params.data.frame' END


################################################################################
# Author : Mauricio Zambrano-Bigiarini                                         #  
################################################################################
# Started: 23-Feb-2012                                                         #
# Updates:                                                                     #
################################################################################
plot_params.matrix <- function(params, 
                        gofs=NULL,
                        ptype=c("histogram", "dottyplot", "boxplot", "vioplot", "pairs"),
                        param.cols=1:ncol(params),
                        param.names=colnames(params),
                        #of.col=NULL, 
                        of.name="GoF", 
                        MinMax=NULL, 
                        beh.thr=NA, 
                        beh.col="red", 
                        beh.lty=1, 
                        beh.lwd=2, 
                        nrows="auto",
                        #col="black", 
                        col="#00000030",
                        ylab=of.name, 
                        main=NULL,
                        pch=19,  # Only used for dotty plots
                        cex=0.5, 
                        cex.main=1.5,
                        cex.axis=1.5,
                        cex.lab=1.5,
                        breaks="Scott",
                        freq=TRUE,
                        verbose=TRUE,
                        ...,
                        #### PNG options ### 
                        do.png=FALSE,
                        png.width=1500,
                        png.height=900,
                        png.res=90,
                        png.fname="Parameters.png"
                        ) {
     
    params <- as.data.frame(params)
                     
    plot_params.data.frame(params=params, 
                        gofs=gofs,
                        ptype=ptype,
                        param.cols=param.cols,
                        param.names=param.names,
                        #of.col=NULL, 
                        of.name=of.name, 
                        MinMax=MinMax, 
                        beh.thr=beh.thr, 
                        beh.col=beh.col, 
                        beh.lty=beh.lty, 
                        beh.lwd=beh.lwd, 
                        nrows=nrows,
                        #col="black", 
                        col=col,
                        ylab=ylab, 
                        main=main,
                        pch=pch,  # Only used for dotty plots
                        cex=cex, 
                        cex.main=cex.main,
                        cex.axis=cex.axis,
                        cex.lab=cex.lab,
                        breaks=breaks,
                        freq=freq,
                        verbose=verbose,
                        ...,
                        #### PNG options ### 
                        do.png=do.png,
                        png.width=png.width,
                        png.height=png.height,
                        png.res=png.res,
                        png.fname=png.fname
                        )
                        
} # 'plot_params.matrix' END
