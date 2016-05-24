# File quant2ecdf.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2010-2012 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################                                                                               
#                                 'quant2ecdf'                                 #
################################################################################                                                                                
# Purpose: This function computes ECDFs for user-defined quantiles of the      #
#          streamflows simulated with different behavioural parameter sets,    # 
#          with optional plot                                                  #  
################################################################################    
# Author : Mauricio Zambrano-Bigiarini                                         #
# Started: June 28th, 2010                                                     #        
# Updates: 01-Dec-2010                                                         #
#          28-Jan-2011 ; 13-Feb-2011 ; 27-Apr-2011  ; 12-Oct-2011              #
#          15-Feb-2012 ; 29-Nov-2012                                           #
################################################################################                                                                                
# Steps:
# 1) it computes un-weighted quantiles (e.g., Q5, Q50, Q95) for the 
#    streamflows simulated with EACH behavioural parameter set
# 2) it computes ECDFs for each desired quantile, by weighting the quantile
#    of each behavioural parameter set by its corresponding less-formal likelihood

# sim    : matrix with the behavioural model ouputs
# obs    : vector with the observed values
# weights: numeric vector, with the values of the weights to be used for computing the quantiles. \cr
#          Omitting the \code{weights} argument or specifying \code{NULL} or a zero-length vector will result in the usual unweighted estimates.
# byrow  : logical, indicating if the computations have to be made for each column or for each row of \code{x}.
#          When the simulated values obtained with different behavioural parameter sets are stored in columns, \code{byrow} must be \kbd{TRUE}. \cr
#          When the simulated values obtained with different behavioural parameter sets are stored in rows, \code{byrow} must be \kbd{FALSE}.

quant2ecdf <-function(sim, ...) UseMethod("quant2ecdf")

quant2ecdf.default <- function(sim, 
                       weights=NULL,
                       byrow=TRUE, 
                       quantiles.desired= c(0.05, 0.5, 0.95),
                       plot=TRUE,
                       obs=NULL, 
                       quantiles.labels= c("Q5", "Q50", "Q95"),
                       main=NULL,
                       ylab="Probability",
                       col="blue",
                       leg.cex=1.2,
                       leg.pos="bottomright",
                       cex.axis=1.2, 
                       cex.main=1.2, 
                       cex.lab=1.2,
                       verbose=TRUE, 
                       ...
                       ) {
                       
    if (is.na(match(class(sim), c("matrix", "data.frame") ) ) )
          stop("Invalid argument: 'class(sim)' must be in c('matrix', 'data.frame')") 
    
    # Unweighted quantiles (usually, Q5, Q50, Q95).
    # Due to the fact that 'weights' is missing, and its default value 
    # is NULL, the computation returns un-weighted quantiles 
    ifelse(byrow==TRUE, stg <- "ROW", stg <- "COLUMN")
    if (verbose) message( "[ Computing un-weighted quantiles for each ", stg, " of 'sim' ... ]" )
    q5.50.95 <- wquantile(sim, byrow=byrow, probs=quantiles.desired) 
    
    # number of desired quantiles, usually 3
    nquantiles <- length(quantiles.desired) 
    
    # creating the final output, a list with the ECDFs 
    ecdf <- vector("list", nquantiles)
    
    if (plot) {
      # Saving default plotting parameters
      old.par <- par(no.readonly=TRUE)
      #if (!do.png) on.exit(par(old.par))
    
      par(mfrow=c(1,nquantiles)) 
      if (!is.null(main)) par(oma=c(1,1,3,0))   
    } # IF end
    
    ###########################   MAIN LOOP   ##################################
    ifelse(is.null(weights), stg <- "the UN-weighted", stg <- "the weighted")
    for ( i in 1:nquantiles ) {
    
        if (verbose) message("[ Computing ", stg, " ECDF for quantile ", 
                                  format(quantiles.desired[i], width=5, justify="left"), 
                                  " , ", i, "/", nquantiles, "=>", 
                                  format(round(100*i/nquantiles, 2), width=7, justify="left"),
                                  "% ... ]" )
        
        # Weighted ECDF for the "i-th" desired quantile, where the unweighted 
        # 'i-th' quantile of each behavioural parameter set is now weighted by the 
        # weights given by 'w' (usually, the normalized less-formal likelihood)
        x.ecdf <- Hmisc::wtd.Ecdf(q5.50.95[, i], weights = weights, normwt=TRUE)
        
        ecdf[[i]]      <- x.ecdf
        names(ecdf)[i] <- quantiles.labels[i]
        
        if (plot) {        
          if ( !is.null(obs) & is.numeric(obs) ) {      
             if (is.na(match(class(obs), c("zoo", "numeric", "integer") ) ) )
               stop("Invalid argument: 'class(obs)' must be in c('zoo', 'numeric', 'integer')") 
               
              # Observed quantile (Q5, Q50, Q95) during the calibration period
              quantile.obs <- quantile(as.numeric(obs), probs=quantiles.desired[i], na.rm=TRUE) 
          } # IF end
    
          # plot label for each ECDF
          main.loc <- paste("Empirical CDF of", quantiles.labels[i], sep=" ")
         
          # Drawing the plotting the area, but without Y-axis
          plot(x.ecdf$x, x.ecdf$ecdf, xlab= paste(quantiles.labels[i], "[m3/s]", sep=" "), 
               col=col, yaxt = "n", type="b", cex=0.2, main=main.loc, ylab=ylab, 
               cex.axis=cex.axis, cex.main=cex.main, cex.lab=cex.lab, font.lab=2, ... )
        
          # Drawing the labels in the 'y' axis
          Axis(side = 2, at = seq(0.0, 1, by=0.05), labels = FALSE, 
               cex.axis=cex.axis, cex.lab=cex.lab)
          Axis(side = 2, at = seq(0.0, 1, by=0.1), labels = seq(0.0, 1, by=0.1), 
               cex.axis=cex.axis, cex.lab=cex.lab, font.lab=2)    
        
          # Drawing a vertical line on the observed quantile Q5, Q50, Q95
          if ( !is.null(obs) & is.numeric(obs) )
            abline(v=quantile.obs, lty=3, col="black", lwd=2)
        
          # Drawing an horizontal line on Probability = 0.5
          abline(h=0.5, lty=2, col="grey", lwd=2)
        
          # Computing a function that give the 'x' value that corresponds to a 
          # given value of cdf, by using linear interpolation
          f <- approxfun(x.ecdf$ecdf, x.ecdf$x)
        
          # Quantile corresponding to a cdf=0.5
          quantile.sim <- f(0.5)
        
          # Drawing a vertical line on the simulated quantile Q5, Q50, Q95
          abline(v=quantile.sim, lty=3, col="grey", lwd=2)
          
          # Bias of the simulated streamflows, in percentage [%]
          if ( !is.null(obs) & is.numeric(obs) ) { 
            bias <- 100 * (quantile.sim - quantile.obs) / quantile.obs  
            if (bias == 0) txt.col <- "green"   
            if (bias < 0) txt.col <- "red"
            if (bias > 0) txt.col <- "blue"
        
            # Presenting the value of the observed Q5 as a legend
            leg.txt <- c( paste(quantiles.labels[i], "obs= ", round(quantile.obs,2) ) ,
                          paste(quantiles.labels[i], "sim= ", round(quantile.sim,2) ),
                          paste("Bias=", round(bias,1), "[%]", sep=" " ) ) 
            legend(leg.pos, legend=leg.txt,  inset=0.02, bty="n", 
                   cex =leg.cex, text.col=c("black", "black", txt.col))
          } else {
              leg.txt <- paste(quantiles.labels[i], "sim= ", round(quantile.sim,3) )
              legend(leg.pos, legend=leg.txt, bty="n", cex=leg.cex)
            } # ELSE end  
          
        } # IF end
    
    } # FOR end
    
    if (plot) {
      # Adding a main title for the plot
      if (!is.null(main)) mtext(main, side=3, line=1, cex=cex.main, outer=TRUE)
    } # IF end
    
    return(ecdf)
    
} # END 'quant2ecdf.default'


################################################################################
# Author : Mauricio Zambrano-Bigiarini                                         #
# Started: 12-Oct-2011                                                         #        
# Updates: 12-Oct-2011                                                         #
################################################################################
quant2ecdf.matrix <- function(sim,  
                       weights=NULL,
                       byrow=TRUE, 
                       quantiles.desired= c(0.05, 0.5, 0.95),
                       plot=TRUE,
                       obs=NULL,
                       quantiles.labels= c("Q5", "Q50", "Q95"),
                       main=NULL,
                       ylab="Probability",
                       col="blue",
                       leg.cex=1.2,
                       leg.pos="bottomright",
                       cex.axis=1.2, 
                       cex.main=1.2, 
                       cex.lab=1.2,
                       verbose=TRUE, 
                       ...
                       ) {
    
  quant2ecdf.default(sim=sim, 
                 weights=weights,
                 byrow=byrow, 
                 quantiles.desired=quantiles.desired,
                 plot=plot,
                 obs=obs, 
                 quantiles.labels= quantiles.labels,
                 main=main,
                 ylab=ylab,
                 col=col,
                 leg.cex=leg.cex,
                 leg.pos=leg.pos,
                 cex.axis=cex.axis, 
                 cex.main=cex.main, 
                 cex.lab=cex.lab,
                 verbose=verbose, 
                  ...    
                    )            
                       
} # 'quant2ecdf.matrix' END


################################################################################
# Author : Mauricio Zambrano-Bigiarini                                         #
# Started: 12-Oct-2011                                                         #        
# Updates: 12-Oct-2011                                                         #
###############################################################################
quant2ecdf.data.frame <- function(sim,  
                       weights=NULL,
                       byrow=TRUE, 
                       quantiles.desired= c(0.05, 0.5, 0.95),
                       plot=TRUE,
                       obs=NULL,
                       quantiles.labels= c("Q5", "Q50", "Q95"),
                       main=NULL,
                       ylab="Probability",
                       col="blue",
                       leg.cex=1.2,
                       leg.pos="bottomright",
                       cex.axis=1.2, 
                       cex.main=1.2, 
                       cex.lab=1.2,
                       verbose=TRUE, 
                       ...
                       ) {
                       
  sim <- as.matrix(sim) 
    
  quant2ecdf.matrix(sim=sim, 
                 weights=weights,
                 byrow=byrow, 
                 quantiles.desired=quantiles.desired,
                 plot=plot,
                 obs=obs, 
                 quantiles.labels= quantiles.labels,
                 main=main,
                 ylab=ylab,
                 col=col,
                 leg.cex=leg.cex,
                 leg.pos=leg.pos,
                 cex.axis=cex.axis, 
                 cex.main=cex.main, 
                 cex.lab=cex.lab,
                 verbose=verbose, 
                  ...    
                    )            
                       
} # 'quant2ecdf.data.frame' END
