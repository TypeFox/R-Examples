# File params2ecdf.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2008-2011 Mauricio Zambrano-Bigiarini
# Distributed under GPL 2 or later

################################################################################
#                                'params2ecdf'                                 #
################################################################################
# Purpose:      This function computes (weighted) empirical CDFs               #
#               (ECDFs) for each calibrated parameter, by using the parameter  #
#               values obtained during the optimisation with PSO, with optional#
#               plot                                                           #
################################################################################                                                                              #  

# params      : matrix with the behavioural parameter values, where each row represent a different parameter set, and each column represent the value 
#               of a different model's parameter.
# param.names : character vector, which meaningful names to be used for each model's parameter in \code{params} (by default column namaes)
# weights     : numeric vector, with the values of the weights to be used for computing the empirical CDFs. \cr
#               Omitting the \code{weights} argument or specifying \code{NULL} or a zero-length vector will result in the usual unweighted estimates.
# byrow       : logical, indicating if the computations have to be made for each column or for each row of \code{params}. \cr
#               When the parameter sets are stored in rows, i.e., values for different model's parameter are stored in columns, \code{byrow} must be \kbd{FALSE}. \cr
#               When the parameter sets are stored in columns, i.e., values for different model's parameter are stored in rows, \code{byrow} must be \kbd{TRUE}.
# plot        : logical, indicating if a plot with the Empirical CDFs for each model's parameter has to be produced or not.
# nrows       : Number of rows to be used in the plotting window. The number 
#               of columns is automatically computed depending on the number 
#               of columns of 'params'

params2ecdf <- function(params, ...) UseMethod("params2ecdf")


################################################################################                                                                                
# Author : Mauricio Zambrano-Bigiarini                                         #
# Started: 12-Oct-2011                                                         #        
# Updates: 15-Feb-2012 ; 21-Feb-2012 ; 19-Nov-2012                             #
#          09-Abr-2014                                                         #
################################################################################
params2ecdf.default <- function(params, 
                                param.names=NULL,
                                gofs=NULL,
                                MinMax=NULL, 
                                beh.thr=NA, 
                                weights=NULL,                                                  
                                byrow=FALSE, 
                                plot=TRUE,
                                obs=NULL,
                                main=NULL,
                                nrows="auto",  
                                ylab="Probability",
                                col="blue",
                                leg.cex=1.2,
                                leg.pos="topleft",
                                cex.axis=1.2, 
                                cex.main=1.2, 
                                cex.lab=1.2,
                                verbose=TRUE, 
                                ...,
                                #### PNG options ### 
                                do.png=FALSE,
                                png.width=1500,
                                png.height=900,
                                png.res=90,
                                png.fname="Params_ECDFs.png" 
                                ) {
                       
 if (is.null(param.names)) param.names <- deparse(substitute(params))
          
 # number of parameters
 nparam <- NCOL(params)

 # Number of parameter sets
 n <- NROW(params) 
    
 # Checking 'param.names'
 if (length(param.names) != nparam)
   stop("Invalid argument: 'length(param.names) = ", length(param.names), " != ", nparam, " = nparam'")

 # Checking 'beh.thr'
 if ( !is.na(beh.thr) ) {
   if ( is.null(MinMax) )
     stop("Missing argument: 'MinMax' has to be provided before using 'beh.thr' !!")        
   if ( is.null(gofs) ) {
     stop("Missing argument: 'gofs' has to be provided before using 'beh.thr' !!")
   } else if (length(gofs) != n)
       stop("Invalid argument: 'length(gofs) != nrow(params)' (", length(gofs), "!=", n, ") !!" ) 
 } # IF end
         
 # Checking 'MinMax'
 if ( !is.null(MinMax) ) {
   if ( !(MinMax %in% c("min", "max")) )
     stop("Invalid argument: 'MinMax' must be in c('min', 'max')")
 } # IF end
      
 # checking that the user provided 1 weight for each behavioural parameter set
 if ( !is.null(weights) ) {
   if (length(weights) != n )
     stop("Invalid argument: 'length(w) != nrow(params)' (", length(weights), "!=", n, ")" )
 } # IF end
    
 # creating the final output, a list with the ECDFs 
 ecdf <- vector("list", nparam)  
    
 # Checking 'do.png' and 'plot'
 if (do.png==TRUE & plot==FALSE)
   stop("Invalid argument: 'plot=FALSE' & 'do.png=TRUE' is not possible !!")

 if (nparam==1) params <- matrix(params, ncol=1)

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
    
 ########################     Plotting Preliminars   ########################
 # If there are too many parameters to plot,more than 1 plot is produced
 nfigs <- ceiling(nparam/21)
    
 # Creating a backup of the original 'params'
 params.bak <- params
 nparam.bak <- nparam
    
 # plot/parameter number
 p <- 1
    
 for (f in 1:nfigs) {
    
   # Setting the filename of the PNG file
   ifelse(f==1, fig.png.fname <- png.fname, 
                fig.png.fname <- paste(substr(png.fname, 1,nchar(png.fname)-4), f, ".png", sep="") )
       
   # PNG ?            
   if (do.png) {
     png(filename=fig.png.fname, width=png.width, height=png.height, res=png.res)
   } else if (f >1) dev.new()
         
   # Subsetting   
   if (nparam.bak <= 21) {
     L <- nparam.bak
   } else ifelse(f*21 <= nparam.bak, L <- f*21, L <- nparam.bak) 
   params <- params.bak[ ,((f-1)*21+1):L]
   nparam <- NCOL(params)
   if (nparam==1) params <- matrix(params, ncol=1)
    
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
    
   if (!is.null(main)) par(oma=c(1,1,3,0))
    
   # Computting the weighted ECDF for each parameter
   ncharmax <- max(nchar(param.names))
   for ( i in 1:nparam ) {
    
       if ( !is.null(weights) ) {
           char1 <- "weighted ECDF"
           char2 <- "wECDF"            
         } else {
           char1 <- "ECDF"
           char2 <- "ECDF"
         } # ELSE end
    
       if (verbose) message("[ Computing the ", char1, " for '", 
                             format(param.names[p], width=ncharmax), 
                             "' , ", p, "/", nparam.bak, " => ", 
                             format(round(100*i/nparam.bak, 2), width=5, 
                             nsmall=2, justify="left"),
                             "% ]" )
    
       # Weighted ECDF for the "i-th" desired quantile, where the unweighted 
       # 'i-th' quantile of each behavioural parameter set is now weighted by the 
       # weights given by 'w' (usually, the normalized less-formal likelihood)
       p.ecdf <- Hmisc::wtd.Ecdf(params[, i], weights = weights, normwt=TRUE)
        
       ecdf[[p]]      <- p.ecdf
       names(ecdf)[p] <- param.names[p]  
       
       ################### PLOTTING ###########################################
       if (plot) {
        
         if ( !is.null(obs) & is.numeric(obs) ) { 
            if (is.na(match(class(obs), c("zoo", "numeric", "integer") ) ) )
              stop("Invalid argument: 'class(obs)' must be in c('zoo', 'numeric', 'integer')") 
               
             # Observed value
             quantile.obs <- as.numeric(obs[p])
         } # IF end      
        
         # plot label
         main.loc <- paste(char2, "of", param.names[p], sep=" ")
        
         # Drawing the plotting the area, but without Y-axis
         # cex fix the point size in the ecdf
         plot(p.ecdf$x, p.ecdf$ecdf, xlab= param.names[p], main=main.loc, 
              ylab=ylab, col=col, yaxt = "n", type="b", cex=0.2, 
              cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab, font.lab=2, ... )
        
         # Drawing the labels in the 'y' axis
         Axis(side = 2, at = seq(0.0, 1, by=0.05), labels = FALSE, 
              cex.axis=cex.axis, cex.lab=cex.lab)
         Axis(side = 2, at = seq(0.0, 1, by=0.1), labels = seq(0.0, 1, by=0.1), 
              cex.axis=cex.axis, cex.lab=cex.lab, font.lab=2 ) 
               
         # Drawing a vertical line on the observed quantile Q5, Q50, Q95
         if ( !is.null(obs) & is.numeric(obs) )
           abline(v=quantile.obs, lty=3, col="black", lwd=2)   
        
         # Drawing an horizontal line on Probability = 0.5
         abline(h=0.5, lty=2, col="grey", lwd=2)        
        
         # Computing a function that give the 'x' value that corresponds to a 
         # given value of cdf, by using linear interpolation
         f <- approxfun(p.ecdf$ecdf, p.ecdf$x)
        
         # Quantile corresponding to a ecdf=0.5
         quantile.sim <- f(0.5)
          
         # Drawing a vertical line on the simulated quantile Q5, Q50, Q95
         abline(v=quantile.sim, lty=3, col="grey", lwd=2)
        
         # Drawing a legend
         if ( !is.null(obs) & is.numeric(obs) ) { 
          
           # Bias of the simulated streamflows, in percentage [%]          
           bias <- 100 * (quantile.sim - quantile.obs) / quantile.obs  
           if (bias == 0) txt.col <- "green"   
           if (bias < 0) txt.col <- "red"
           if (bias > 0) txt.col <- "blue"
        
           # Presenting the value of the observed Q5 as a legend
           leg.txt <- c( paste("Obs    : ", round(quantile.obs,3) ) ,
                         paste("Q50 sim: ", round(quantile.sim,3) ),
                         paste("Bias   : ", round(bias,1), "[%]", sep=" " ) ) 
           legend(leg.pos, legend=leg.txt,  inset=0.02, bty="n", 
                  cex =leg.cex, text.col=c("black", "black", txt.col)) 
                              
         } else {
             leg.txt <- paste("Q50 sim: ", round(quantile.sim,3) ) 
             legend(leg.pos, legend=leg.txt, bty="n", cex=leg.cex)
           } # ELSE end      
        
       } # IF end
          
     p <- p+1
    
   } # FOR i end
    
   if (!is.null(main)) mtext(main, side=3, line=1, cex=cex.main, outer=TRUE)
      
   if (do.png) dev.off()
      
 } # FOR f end
    
 return(ecdf)
    
} # END 'params2ecdf.default'


################################################################################  
# Author : Mauricio Zambrano-Bigiarini                                         #
# Started: 12-Oct-2011                                                         #        
# Updates: 12-Oct-2011 ; 19-Nov-2012                                           #
################################################################################
params2ecdf.matrix <- function(params, 
                               param.names=colnames(params),
                               gofs=NULL,
                               MinMax=NULL, 
                               beh.thr=NA, 
                               weights=NULL,                                                  
                               byrow=FALSE, 
                               plot=TRUE,
                               obs=NULL,
                               main=NULL,
                               nrows="auto",  
                               ylab="Probability",
                               col="blue",
                               leg.cex=1.2,
                               leg.pos="topleft",
                               cex.axis=1.2, 
                               cex.main=1.2, 
                               cex.lab=1.2,
                               verbose=TRUE, 
                               ...,
                               #### PNG options ### 
                               do.png=FALSE,
                               png.width=1500,
                               png.height=900,
                               png.res=90,
                               png.fname="Params_ECDFs.png" 
                               ) {
 params2ecdf.default(params=params, 
                     param.names=param.names,
                     gofs=gofs,
                     MinMax=MinMax, 
                     beh.thr=beh.thr, 
                     weights=weights,                                                  
                     byrow=byrow, 
                     plot=plot,
                     obs=obs,
                     main=main,
                     nrows=nrows,  
                     ylab=ylab,
                     col=col,
                     leg.cex=leg.cex,
                     leg.pos=leg.pos,
                     cex.axis=cex.axis, 
                     cex.main=cex.main, 
                     cex.lab=cex.lab,
                     verbose=verbose, 
                     ...,
                     # PNG options
                     do.png=do.png,
                     png.width=png.width,
                     png.height=png.height,
                     png.res=png.res,
                     png.fname=png.fname 
                     )                         
} # END 'params2ecdf.data.frame'

################################################################################  
# Author : Mauricio Zambrano-Bigiarini                                         #
# Started: 12-Oct-2011                                                         #        
# Updates: 12-Oct-2011 ; 19-Nov-2012                                           #
################################################################################
params2ecdf.data.frame <- function(params, 
                                   param.names=colnames(params),
                                   gofs=NULL,
                                   MinMax=NULL, 
                                   beh.thr=NA, 
                                   weights=NULL,                                                  
                                   byrow=FALSE, 
                                   plot=TRUE,
                                   obs=NULL,
                                   main=NULL,
                                   nrows="auto",  
                                   ylab="Probability",
                                   col="blue",
                                   leg.cex=1.2,
                                   leg.pos="topleft",
                                   cex.axis=1.2, 
                                   cex.main=1.2, 
                                   cex.lab=1.2,
                                   verbose=TRUE, 
                                   ...,
                                   #### PNG options ### 
                                   do.png=FALSE,
                                   png.width=1500,
                                   png.height=900,
                                   png.res=90,
                                   png.fname="Params_ECDFs.png"  
                                   ) {
 params <- as.matrix(params)
 
 params2ecdf.default(params=params, 
                     param.names=param.names,
                     gofs=gofs,
                     MinMax=MinMax, 
                     beh.thr=beh.thr, 
                     weights=weights,                                                  
                     byrow=byrow, 
                     plot=plot,
                     obs=obs,
                     main=main,
                     nrows=nrows,  
                     ylab=ylab,
                     col=col,
                     leg.cex=leg.cex,
                     leg.pos=leg.pos,
                     cex.axis=cex.axis, 
                     cex.main=cex.main, 
                     cex.lab=cex.lab,
                     verbose=verbose, 
                     ...,
                     # PNG options
                     do.png=do.png,
                     png.width=png.width,
                     png.height=png.height,
                     png.res=png.res,
                     png.fname=png.fname 
                     )                         
} # END 'params2ecdf.data.frame'
