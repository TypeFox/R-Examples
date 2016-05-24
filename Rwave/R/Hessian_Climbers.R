#########################################################################
#      $Log: Pca_Climbers.S,v $
#
#               (c) Copyright  1997
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University
#                  All right reserved                           
#########################################################################





hescrc <- function(tfrep, tfspec = numeric(dim(tfrep)[2]), grida = 10,
                   gridb = 20, bstep = 3, iteration = 10000, rate = .001,
                   seed = -7, nbclimb = 10, flag.int = TRUE, chain= TRUE,
                   flag.temp = FALSE, lineflag = FALSE)
#########################################################################
#  hescrc:   Time-frequency multiple ridge estimation (hessian climbers)
#  ------
# 	use the hessian climber algorithm to evaluate ridges of
#    	   continuous Time-Frequency transform
#
#      input:
#      ------
# 	tfrep: wavelet or Gabor transform
#	tfspec: additional potential (coming from learning the noise)
#       grida : window of a
#       gridb : window of b
#       pct : the percent number of points in window 2*grida and
#             2*gridb to select histogram
#       count : the maximal number of repetitive selection, if location
#               was chosen before
#       iteration: number of iterations
#       rate: initial value of the temperature
#       seed: initialization for random numbers generator
#       nbclimb: number of crazy climbers
#       bstep: step size for the climber walk in the time direction
#	flag.int: if set to TRUE, computes the integral on the ridge.
#	chain: if set to TRUE, chains the ridges.
#	flag.temp: if set to TRUE, keeps a constant temperature.
#       linfflag: if set to TRUE, the line segments oriented to the
#                 where the climber will move freely at the block is
#                 shown in image   	
#       
#      output:
#      -------
#       beemap: 2D array containing the (weighted or unweighted)
#               occupation measure (integrated with respect to time)
#       hesmap: 2D array containing number from 1 to 4, denoting the
#               direction a climber will go at some position; where
#               1 moving freely along b, 3 freely along a, 2 along 
#               line a = -b, and 4 along the line a = b. The restricted
#               move is perdendicular to the free move.
#
#########################################################################
{

   tfspectrum <- tfspec 

   d <- dim(tfrep)
   sigsize <- d[1]
   nscale <- d[2]

   sqmodulus <- Re(tfrep*Conj(tfrep))
   for (k in 1:nscale)
     sqmodulus[,k] <- sqmodulus[,k] - tfspectrum[k]
   image(sqmodulus)

   ## number of grid size
   nbx <- as.integer(sigsize/gridb)
   nby <- as.integer(nscale/grida)
   if((sigsize/gridb - nbx) > 0) nbx <- nbx + 1
   if((nscale/grida - nby) > 0) nby <- nby + 1

   nbblock <- nbx * nby

   ## first two locations at each block stores the lower-left corner
   ## of the block followed by the locations of up-right corner and by
   ## the negative values of hessian matrix from xx, xy, yx, and yy at
   ## the last.

   tst <- matrix(0, 8, nbblock)
   dim(tst) <- c(nbblock * 8,1)

   dim(sqmodulus) <- c(sigsize * nscale, 1)

   z <- .C("Shessianmap",
           as.double(sqmodulus),
           as.integer(sigsize),
           as.integer(nscale),
           ublock = as.integer(nbblock),
           as.integer(gridb),
           as.integer(grida),
           tst =as.double(tst),
           PACKAGE="Rwave")
  
   tst <- z$tst
   dim(tst) <- c(8, nbblock)

   ## actual number of block 	
   ublock <- z$ublock

   ## principle component analysis and calculate the direction of each block
   pcamap <- matrix(1,sigsize,nscale)

   ## first eigenvector of a block
   eigv1 <- matrix(0,ublock,2)

   ## to draw eigenvector in a block 
   if(lineflag) {	
     lng <- min(grida, gridb)
     oldLng <- lng
     linex <- numeric(oldLng)
     liney <- numeric(oldLng)
   }


   for(j in 1:ublock) {
        
	left <- max(1,as.integer(tst[1,j]))
	down <- max(1,as.integer(tst[2,j]))
	right <- min(as.integer(tst[3,j]),sigsize)
	up <- min(as.integer(tst[4,j]),nscale)

	ctst <- tst[5:8,j]
	dim(ctst) <- c(2,2)
        ctst <- t(ctst)
        eig<- eigen(ctst)

	
        ## cat("first eig = ",eig$values[1]," ratio in eigenvalues =
        ## ",abs(eig$values[1]/eig$values[2])," \n")

        ## eigen vector 1
        eigv1[j,] <- eig$vector[,1]


        ## shifted center position of a block
	if(lineflag) {
          centerx <- as.integer((left + right)/2)
          centery <- as.integer((down + up)/2)
	}

        ## theta of the eigen vector 1; -pi < theta <= pi

        theta <- atan2((eigv1[j,])[2],(eigv1[j,])[1])

        ## along x : denote 1 in the hesmap

        if(((theta <= pi/8) && (theta > -pi/8)) ||
	   (theta > 7*pi/8 || theta <= -7*pi/8)) {

	  pcamap[left:(right-1),down:(up-1)] <- 1

	  if(lineflag) {
             Lng <- min(lng,(sigsize-centerx))
             if(Lng != oldLng) {
               linex <- numeric(Lng)
               liney <- numeric(Lng)
             }
             llng <- as.integer(Lng/2)
	     x1 <- max(centerx-llng,1)
             x2 <- min(sigsize,x1 + Lng-1)
             linex[] <- x1:x2
             liney[] <- centery
            lines(linex,liney)
	}
	}

        ## along x=y : denoted as 4 in hesmap

        if(((theta > pi/8) && (theta <= 3*pi/8)) ||
	   ((theta > -7*pi/8) && (theta <= -5*pi/8))) {

	  pcamap[left:(right-1),down:(up-1)] <- 4

	  if(lineflag) {
            Lng <- min(lng,sigsize-centerx)
            Lng <- min(Lng,nscale-centery)
            if(Lng != oldLng) {
              linex <- numeric(Lng)
              liney <- numeric(Lng)
            }
            llng <- as.integer(Lng/2)
            x1 <- max(centerx-llng,1)
            x2 <- min(sigsize,x1 + Lng-1)
	    y1 <- max(centery-llng,1)
            y2 <- min(sigsize,y1 + Lng-1)

            linex[] <- x1:x2
            liney[] <- y1:y2
            lines(linex,liney)
           }
	}


        ## along y : as 3 in hesmap

        if(((theta > 3*pi/8) && (theta <= 5*pi/8)) ||
	   ((theta > -5*pi/8) && (theta <= -3*pi/8))) {

	  pcamap[left:(right-1),down:(up-1)] <- 3

	  if(lineflag) {
             Lng <- min(lng, nscale-centery)
             if(Lng != oldLng) {
               linex <- numeric(Lng)
               liney <- numeric(Lng)
             }
            llng <- as.integer(Lng/2)
	    y1 <- max(centery-llng,1)
            y2 <- min(sigsize,y1 + Lng-1)

            linex[] <- centerx
            liney[] <- y1:y2
            lines(linex,liney)
          }
	}


        ## along x=-y : as 2 in hesmap

        if(((theta > 5*pi/8) && (theta <= 7*pi/8)) ||
	   ((theta > -3*pi/8) && (theta <= -pi/8))) {

	  pcamap[left:(right-1),down:(up-1)] <- 2

	  if(lineflag) {
            Lng <- min(lng,nscale-centery)
            Lng <- min(Lng,sigsize-centerx)
            if(Lng != oldLng) {
               linex <- numeric(Lng)
               liney <- numeric(Lng)
            }
            llng <- as.integer(Lng/2)

	    x1 <- max(centerx-llng,1)
            x2 <- min(sigsize,x1 + Lng-1)
	    y1 <- max(centery-llng,1)
            y2 <- min(sigsize,y1 + Lng-1)

            linex[] <- x1:x2
            liney[] <- y2:y1
            lines(linex,liney)
          }
	}

	  if(lineflag) oldLng <- Lng
   }

   ##   if(lineflag==F) image(pcamap)	

   ## From the following on, it is similar to the process of crazy climbers...

   beemap <- matrix(0,sigsize,nscale) 
   dim(beemap) <- c(nscale * sigsize, 1)
   dim(pcamap) <- c(nscale * sigsize, 1)
     
   z <- .C("Spca_annealing",
           as.double(sqmodulus),
           beemap= as.double(beemap),
           as.integer(pcamap),
           as.double(rate),
           as.integer(sigsize),
           as.integer(nscale),
           as.integer(iteration),
           as.integer(seed),
           as.integer(bstep),
           as.integer(nbclimb),
           as.integer(flag.int),
           as.integer(chain),
           as.integer(flag.temp),
           PACKAGE="Rwave")

   beemap <- z$beemap
   dim(beemap) <- c(sigsize, nscale)
   dim(pcamap) <- c(sigsize, nscale)
   image(beemap)
   list(beemap = beemap, pcamap = pcamap)
}
