################################################################################
#                          'plot_GofPerParticle'                               #
################################################################################
# Purpose: It plots the value of each parameter and the objective              #
#          function against the iteration number.                              #
################################################################################
# Author  : Mauricio Zambrano-Bigarini                                         #
# Started : 20-Dec-2010 at JRC, Ispra                                          #
# Modified: 11-Nov-2011 ; 21-Feb-2012                                          #
################################################################################
 
plot_GofPerParticle <- function(x, # data.frame with the GoF for each particle per iteration. ncol=number of particles ;  nrow = number of iterations
                                ptype="one", # character, representing the type of plot. Valid values are: in c("one", "many"), for plotting all the particles in the smae figure or in one windows per particle, respectively
                                nrows="auto",
                                main=NULL,
                                xlab="Number of Iterations",
                                cex=0.4,
                                cex.main=1.5,
                                cex.axis=1.7,
                                cex.lab=1.5,                             
                                col=rainbow(ncol(x)),
                                lty=3,
                                ylim=NULL,
                                verbose=TRUE,
                                ...,
                                #### PNG options ### 
                                do.png=FALSE,
                                png.width=1500,
                                png.height=900,
                                png.res=90,
                                png.fname="Particles_GofPerIter.png"                                     
                                ) {    
                                
    # Checking the user provide a valid value for 'ptype'
    if (is.na(match(ptype, c("one", "many"))))
        stop("Invalid argument: 'ptype' must be in c('one', 'many')")

    # Number of iterations that will be analysed
    niter <- nrow(x)
    iters <- 1:niter

    # Number of particles
    npart <- ncol(x)
 
    # Computing the number of rows for the plot 
    if (nrows == "auto") {
      if ( npart <= 5 )                  lnr <- 1
      if ( (npart > 5) & (npart <= 14) ) lnr <- 2
      if ( npart > 14 )                  lnr <- ceiling(npart/7)
    } else lnr <- nrows 
    
    # Defining the plotting window
    nr <- lnr
    nc <- ceiling(npart/lnr) 

    # Range of the objective function, considering ALL the particles
    if (is.null(ylim)) ylim <- range( x, na.rm=TRUE )
    
    ##############################################################################
    # 2)                             Plotting                                    #
    ##############################################################################
    msg <- "[ Plotting GoF for each particle vs Number of Model Evaluations"
    if (do.png) msg <- paste(msg, " into '", basename(png.fname), sep="")
    msg <- paste(msg, "' ... ]", sep="")
    message(msg) 

    if (do.png) png(filename=png.fname, width=png.width, height=png.height, res=png.res)   

    # Saving default plotting parameters
    old.par <- par(no.readonly=TRUE)
    if (!do.png) on.exit(par(old.par))        
    
    # Doing the plots
    if (ptype=="many") {

    par(mfrow=c(nr,nc))   
    if (!is.null(main)) par(oma=c(1,1,3,0)) 
    for ( i in 1:npart ) {
      plot(iters, 
           x[,i], 
           type="o", lty=lty, col=col[i], 
           cex=cex, cex.lab=cex.lab, cex.axis=1.5, 
           main=paste("Particle", i, sep=" "), 
           xlim=c(1,niter), ylim=ylim, 
           xlab=xlab, 
           ylab=paste("Particle", i, sep=" ")
          )
    } # FOR end
    
    # Adding a main title for the plot
    if (!is.null(main)) mtext(main, side=3, line=1, cex=cex.main, outer=TRUE)

   } else if (ptype=="one") {
   
     if (is.null(main)) main <- "GoF per Particle"
     # Plotting the first particle 
     plot(iters, 
           x[,2], 
           type="o", lty=lty, col=col[1], 
           cex=cex, cex.lab=cex.lab, cex.axis=1.5, 
           main=main, 
           xlim=c(1,niter), ylim=ylim, 
           xlab=xlab, 
           ylab="GoF"
          )

     if (npart > 1) {
       for ( i in 2:npart ) {
         lines(iters, 
               x[,i], 
               type="o", lty=lty, col=col[i], 
               cex=cex, cex.lab=cex.lab, cex.axis=1.5, 
               xlim=c(1,niter), ylim=ylim, 
               xlab=xlab, 
               ylab="GoF"
              )
      } # FOR end
     } # IF end
    } # ELSE end
    
    if (do.png) dev.off()
   
} # 'plot_GofPerParticle' END

#plot_GofPerParticle(ptype="many", nrows=3 )
