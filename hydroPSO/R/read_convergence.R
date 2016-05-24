# File read_convergence.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2008-2011 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################
#                         'read_convergence'                                   # 
################################################################################
# Author : Mauricio Zambrano-Bigiarini & Rodrigo Rojas                         #  
# Started: 08-Nov-2011,                                                        #
# Updates: 13-Ene-2012 ; 29-Ene-2012 ; 01-Feb-2012 ; 14-Nov-2012 ; 19-Nov-2012 #        
################################################################################
# Purpose: To read the output file 'ConvergenceMeasures.txt'                   #
################################################################################
# Output:                                                                      #
# a data.frame with the following columns for each partile and each iteration: #
# 1) Iter       : iteration number                                             #
# 2) Gbest      : global best for the iteration corresponding to 'iter'        #
# 3) GbestRate  : ratio between the global best of the current iteration and   #
#                 the previous one                                             #
# 4) IterBestFit: best fitness value  for the iteration corresponding to 'iter'#
# 5) normSwarmRadius: rormalised swarm radious for the iteration corresponding #
#                     to 'iter'                                                #
# 6) GbestPbestRatio: ratio between the global best and the mean of the        #
#                     personal best of all the particles for the iteration     #
#                     corresponding to 'iter'                                  #
################################################################################
 
                      
read_convergence <- function(file="ConvergenceMeasures.txt",
                             MinMax=NULL, 
                             beh.thr=NA, 
                             verbose=TRUE, 
                             # Plotting arguments #
                             plot=TRUE,
                             col=c("black", "darkolivegreen"),
                             lty=c(1,3), 
                             lwd=c(2,2), 
                             main="Global Optimum & Normalized Swarm Radius vs Iteration Number", 
                             xlab="Iteration Number", 
                             ylab=c("Global Optimum", expression(delta[norm]) ), 
                             pch=c(15, 18), 
                             cex=1, 
                             cex.main=1.4, 
                             cex.axis=1.2, 
                             cex.lab=1.2,     
                             legend.pos="topright",  
                             ...,
                             #### PNG options ### 
                             do.png=FALSE,
                             png.width=1500,
                             png.height=900,
                             png.res=90,
                             png.fname="ConvergenceMeasures.png"                  
                             ) {
                         
                         
  # Checking that 'file' exists
  if ( !file.exists(file) )
     stop( "Invalid argument value: The file '", basename(file), "' does not exist !!")

   # Checking 'beh.thr'
  if ( !is.na(beh.thr) ) {
     if ( is.null(MinMax) )
        stop("Missing argument: 'MinMax' has to be provided before using 'beh.thr' !!")  
  } # IF end
  
  # Checking 'MinMax'
  if ( !is.null(MinMax) ) {
     if ( !(MinMax %in% c("min", "max")) )
        stop("Invalid argument: 'MinMax' must be in c('min', 'max')")
  } # IF end

  # Reading the file
  if (verbose) message( "                                                     ")  
  if (verbose) message( "[ Reading the file '", basename(file), "' ... ]" )  
  conv       <- read.table(file=file, header=FALSE, skip=1) 
  conv.names <- c("Iter", "Gbest", "GbestRate", "IterBestFit", "normSwarmRadius", "GbestPbestRatio" )
  
  # Giving the right name to the column with the conveters
  colnames(conv) <- conv.names

  # computing the number of convergence measures
  nconv <- nrow(conv)
  if (verbose) message( "[ Total number of iterations: ", nconv, " ]" )  

  # Filtering out those parameter sets above/below a certain threshold
  if (!is.na(beh.thr)) {  
    # Getting the goodness-of-fit of each particle
    gbests <- conv[, "Gbest"]

    # Checking 'beh.thr'
    mx <- max(gbests, na.rm=TRUE)
    if (beh.thr > mx)
      stop("Invalid argument: 'beh.thr' must be lower than ", mx ,"!!")
    
    # Computing the row index of the behavioural parameter sets
    if (MinMax=="min") {
       beh.row.index <- which(gbests <= beh.thr)
       beh.symb <- "<="
    } else {
             beh.row.index <- which(gbests >= beh.thr)
             beh.symb <- ">="
           } # ELSE end
    
    # Removing non-behavioural iterations
    conv   <- conv[beh.row.index, ]
    gbests <- gbests[beh.row.index]
   
    # Amount of behavioural parameter sets 
    nbeh <- nrow(conv)
    if (verbose) message( "[ Number of iterations with Gbest ", beh.symb, " ", beh.thr, 
                          ": ", nbeh, " ]" )
    
    # To avoid problems with 'plot_params'
    if (plot) beh.thr <- NA
  } # IF end
  
  if (plot) {
     plot_convergence(x=conv,
                      col=col,
                      lty=lty, 
                      lwd=lwd, 
                      main=main,
                      xlab=xlab, 
                      ylab=ylab, 
                      cex=cex, 
                      pch=pch, 
                      cex.main=cex.main, 
                      cex.axis=cex.axis,
                      cex.lab=cex.lab, 
                      legend.pos=legend.pos,      
                      do.png=do.png,
                      ...,
                      png.width=png.width,
                      png.height=png.height,
                      png.res=png.res,
                      png.fname=png.fname                    
                      )
  } # IF end
  
  return(conv)
  
}  # 'read_convergence' END



