################################################################################
#                           'read_GofPerParticle'                              #
################################################################################
# Purpose: This function reads the output file 'Particles_GofPerIter.txt' and  #
#          plots the value of each parameter and the objective function against#
#          the iteration number.                                               #
################################################################################
# Author  : Mauricio Zambrano Bigarini                                         #
# Started : 20-Dec-2010 at JRC, Ispra                                          #
# Modified: 11-Nov-2011 ; 21-Feb-2012                                          #
################################################################################
 
read_GofPerParticle <- function(file="Particles_GofPerIter.txt", 
                                na.strings="NA",
                                plot=TRUE,
                                ### Plotting arguments ###
                                ptype="one", # Valid values are: in c("one", "many")
                                nrows="auto",
                                main=NULL,
                                xlab="Number of Iterations",
                                cex=0.4,
                                cex.main=1.5,
                                cex.axis=1.7,
                                cex.lab=1.5,                             
                                col,
                                lty=3,
                                ylim=NULL,
                                verbose=TRUE,
                                #### PNG options ### 
                                do.png=FALSE,
                                png.width=1500,
                                png.height=900,
                                png.res=90,
                                png.fname="Particles_GofPerIter.png"   
                                ) {    

    # Checking that 'file' exists
    if ( !file.exists(file) )
       stop( "Invalid argument value: The file '", basename(file), "' doesn't exist")

    # Reading ALL the PARAMETER SETS
    if (verbose) message( "                                                     ")  
    if (verbose) message( "[ Reading the file '", basename(file), "' ... ]" )  
    Particles.GofPerIter <- read.table(file=file, header=TRUE, skip=0, na.strings=na.strings) 

    # Removing the column with the iteration number
    Particles.GofPerIter <- Particles.GofPerIter[, 2:ncol(Particles.GofPerIter)]

    # Number of iterations that will be analysed
    niter <- nrow(Particles.GofPerIter)

    # Number of particles
    npart <- ncol(Particles.GofPerIter)
 
    # Printing the number of particles read
    if (verbose) message( "[ Number of particles : ", npart, " ]" )

    # Printing the number of iterations read
    if (verbose) message( "[ Number of iterations: ", niter, " ]" )
    
    # Setting the colours
    if (missing(col)) col <- rainbow(ncol(Particles.GofPerIter))
    
    # Plotting
    if (plot) {
       plot_GofPerParticle(x=Particles.GofPerIter,
                           ptype=ptype, 
                           nrows=nrows,
                           main=main,
                           xlab=xlab,
                           cex=cex,
                           cex.main=cex.main,
                           cex.axis=cex.axis,
                           cex.lab=cex.lab,                             
                           col=col,
                           lty=lty,
                           ylim=ylim,
                           verbose=verbose,
                           do.png=do.png,
                           png.width=png.width,
                           png.height=png.height,
                           png.res=png.res,
                           png.fname=png.fname   
                           )
    } # IF end       
                           
    return(Particles.GofPerIter)
   
} # 'read_GofPerParticle' END

#plotParticlesGof(ptype="many", nrows=3 )
