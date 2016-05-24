#########################################################################
#      $Log: Pca_Climbers.S,v $
#
#               (c) Copyright  1997
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University 
#                  All right reserved                           
#########################################################################





pcacrc <- function(tfrep, tfspec = numeric(dim(tfrep)[2]), grida = 10,
	gridb = 20, pct = 0.1, count = 10000, bstep = 3, iteration = 10000,
	rate = .001, seed = -7, nbclimb = 10, flag.int = TRUE, chain= TRUE,
	flag.temp = FALSE, lineflag = FALSE, flag.cwt = FALSE, nvoice = 0)
#########################################################################
#  pcacrc:   Time-frequency multiple ridge estimation (pca climbers)
#  ------
# 	use the pca climber algorithm to evaluate ridges of
#    	   continuous Time-Frequency transform
#
#      input:
#      ------
# 	tfrep: wavelet or Gabor transform
#	tfspec: additional potential (coming from learning the noise)
#       grida : window of a
#       gridb : window of b
#       pct : the percent number of points in window 2*grida and
#	  2*gridb to select histogram
#       count : the maximal number of repetitive selection, if
#	  location was chosen before
#       iteration: number of iterations
#       rate: initial value of the temperature
#       seed: initialization for random numbers generator
#       nbclimb: number of crazy climbers
#       bstep: step size for the climber walk in the time direction
#	flag.int: if set to TRUE, computes the integral on the ridge.
#	chain: if set to TRUE, chains the ridges.
#	flag.temp: if set to TRUE, keeps a constant temperature.
#	flag.cwt: if set to TRUE, the bstep as well as the block size
#	 of a and b are adjusted according to the scale of wavelet
#	transform(not implemented yet).
#       nvoice: number of voice per octave (cwt)
#       linfflag: if set to TRUE, the line segments oriented to the
#                 where the climber will move freely at the block is shown in 
#                 image   	
#       
#      output:
#      -------
#       beemap: 2D array containing the (weighted or unweighted)
#               occupation measure (integrated with respect to time)
#       pcamap: 2D array containing number from 1 to 4, denoting the
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

  # percentage of points to be selected in a block  
  # considering overlapping 1/2 on each coordinate with neighborhood 
   nbpoint <- as.integer(2* grida * 2 * gridb * pct)
#   nbpoint <- as.integer(grida * gridb * pct)

  # number of grid size
   nbx <- as.integer(sigsize/gridb)
   nby <- as.integer(nscale/grida)
   if((sigsize/gridb - nbx) > 0) nbx <- nbx + 1
   if((nscale/grida - nby) > 0) nby <- nby + 1

   nbblock <- nbx * nby

  # first two locations at each block stores the lower-left corner of the block
  # followed by the locations of up-right corner and by the coordinates of the 
  # selected points in the block. The locations are in the order of (x, y)

   tstsize <- 2 * nbpoint + 4

   tst <- matrix(0, tstsize, nbblock)
   dim(tst) <- c(nbblock * tstsize,1)

  # the map of points selected in all the block
   pointmap <- matrix(0,sigsize,nscale)
   dim(pointmap) <- c(sigsize * nscale, 1)
   dim(sqmodulus) <- c(sigsize * nscale, 1)

   z <- .C("Spointmap",
           as.double(sqmodulus),
           as.integer(sigsize),
           as.integer(nscale),
           as.integer(gridb),
           as.integer(grida),
           as.integer(nbblock),
           as.integer(nbpoint),
           pointmap = as.integer(pointmap),
           tst = as.double(tst),
           as.integer(tstsize),
           as.integer(count),
           as.integer(seed),
           as.integer(flag.cwt),
           as.integer(nvoice),
           PACKAGE="Rwave")
  
   pointmap <- z$pointmap
   tst <- z$tst
   dim(pointmap) <- c(sigsize,nscale)
   dim(tst) <- c(tstsize, nbblock)

   

  # principle component analysis and calculate the direction at each block


  pcamap <- matrix(1,sigsize,nscale)

  # first eigenvector of a block
   eigv1 <- matrix(0,nbblock,2)

   # to draw eigenvector in a block 
   if(lineflag) {	
     lng <- min(grida, gridb)
     oldLng <- lng
     linex <- numeric(oldLng)
     liney <- numeric(oldLng)
   }


   for(j in 1:nbblock) {
        
	left <- max(1,as.integer(tst[1,j]))
	down <- max(1,as.integer(tst[2,j]))
	right <- min(as.integer(tst[3,j]),sigsize)
	up <- min(as.integer(tst[4,j]),nscale)

	input <- tst[5:tstsize,j]
	dim(input) <- c(2,nbpoint)
        input <- t(input)
        ctst <- var(input)
        eig<- eigen(ctst)

	# cat(" ratio in eigenvalues =  ",abs(eig$values[1]/eig$values[2])," \n")

        # eigen vector 1
        eigv1[j,] <- eig$vector[,1]


        # shifted center position of a block
	if(lineflag) {
          centerx <- as.integer((left + right)/2)
          centery <- as.integer((down + up)/2)
	}

        # theta of the eigen vector 1; -pi < theta <= pi

        theta <- atan2((eigv1[j,])[2],(eigv1[j,])[1])

        # along x : denote 1 in the pcamap

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

        # along x=y : denoted as 2 in pcamap

        if(((theta > pi/8) && (theta <= 3*pi/8)) ||
	   ((theta > -7*pi/8) && (theta <= -5*pi/8))) {

	  pcamap[left:(right-1),down:(up-1)] <- 2

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


        # along y : as 3 in pcamap

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


        # along x=-y : as 4 in pcamap

        if(((theta > 5*pi/8) && (theta <= 7*pi/8)) ||
	   ((theta > -3*pi/8) && (theta <= -pi/8))) {

	  pcamap[left:(right-1),down:(up-1)] <- 4

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

# From the following on, it is similar to the process of crazy climbers ....

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

   


pcafamily <-function(pcaridge,orientmap,
                     maxchnlng=as.numeric(dim(pcaridge)[1])+10,
                     bstep = 1, nbchain = 100, ptile = 0.05)
#########################################################################
#     pcafamily:
#     ---------
#     chain the ridges obtained by pca climber.
#
#      input:
#      ------
#       pcaridge: ridge found by pca climber
#       orientmap : the first eigen vector direction at each position
#       bstep: maximal length for a gap in a ridge
#       nbchain: maximum number of chains
#       ptile: relative threshold for the ridges 
#       maxchnlng: maximal chain length
#
#      output:
#      ------
#       ordered: ordered map: image containing the ridges
#                (displayed with different colors)
#       chain: 2D array containing the chained ridges, according
#		to the chain data structure:
#		chain[,1]: first point of the ridge
#		chain[,2]: length of the chain
#		chain[,3:(chain[,2]+2)]: values of the ridge
#	nbchain: number of chains.
#
#########################################################################
{
   d <- dim(pcaridge)
   sigsize <- d[1]
   nscale <- d[2]
   threshold <- range(pcaridge)[2] * ptile
   sz <- 2 * (maxchnlng);
   chain <- matrix(-1,nbchain,sz)
   orderedmap <- matrix(0,sigsize,nscale)	
   dim(chain) <- c(nbchain * sz, 1)
   dim(pcaridge) <- c(nscale * sigsize, 1)
   dim(orderedmap) <- c(nscale * sigsize, 1)
   dim(orientmap) <- c(nscale * sigsize, 1)
	
   
   z <- .C("Spca_family",
           as.double(pcaridge),
           as.integer(orientmap),
           orderedmap = as.double(orderedmap),
           chain = as.integer(chain),
           chainnb = as.integer(nbchain),
           as.integer(sigsize),
           as.integer(nscale),
           as.integer(bstep),
           as.double(threshold),
           as.integer(maxchnlng),
           PACKAGE="Rwave")
	
   orderedmap <- z$orderedmap
   chain <- z$chain

   dim(orderedmap) <- c(sigsize,nscale)   
   dim(chain) <- c(nbchain,sz)   
   nbchain <- z$chainnb
   chain <- chain +1
   chain[,1] <- chain[,1]-1

   list(ordered = orderedmap, chain = chain, nbchain=nbchain)
}   



pcamaxima <- function(beemap,orientmap)
{
   d <- dim(beemap)
   sigsize <- d[1]
   nscale <- d[2]  
   tfmaxima <- matrix(0,sigsize,nscale)
   dim(orientmap) <- c(nscale * sigsize, 1)
   dim(beemap) <- c(nscale * sigsize, 1)   
   dim(tfmaxima) <- c(nscale * sigsize, 1)   

   z <- .C("Stf_pcaridge",
           as.double(beemap),
           tfmaxima = as.double(tfmaxima),
           as.integer(sigsize),
           as.integer(nscale),
           as.integer(orientmap),
           PACKAGE="Rwave")
	
   tfmaxima <- z$tfmaxima
   dim(tfmaxima) <- c(sigsize,nscale)   
   tfmaxima
}


simplepcarec <- function(siginput, tfinput, beemap, orientmap, bstep = 5,
                         ptile = .01, plot = 2)
#########################################################################
#     simplepcarec:
#     -------------
#      Simple reconstruction of a real valued signal from ridges found by 
#      pca climbers.
#
#      input:
#      ------
#      siginput: input signal
#      tfinput: continuous time-frequency transform (output of cwt or cgt)
#      beemap: output of pca climber algorithm
#      orientmap: pca dirction of climber
#      bstep: used for the chaining
#      ptile: 
#      plot: if set to 1, displays the signal, the components, and
#            signal and reconstruction one after another. If set to
#            2, displays the signal, the components, and the
#            reconstruction on the same page. Else, no plot.
#
#      output:
#      -------
#      rec: reconstructed signal
#      ordered: image of the ridges (with different colors)
#      comp: 2D array containing the signals reconstructed from ridges
#
#########################################################################
{
   tmp <- pcafamily(beemap,orientmap,ptile=ptile)
   chain <- tmp$chain
   nbchain <- tmp$nbchain
   rec <- numeric(length(siginput))

   if(plot != FALSE){
     par(mfrow=c(2,1))
     plot.ts(siginput)
     title("Original signal")
     image(tmp$ordered)
     title("Chained Ridges")
   }

   npl(1)

   sigsize <- dim(tfinput)[1]
   nscale <- dim(tfinput)[2]
   rec <- numeric(sigsize)

   tmp3 <- matrix(0+0i,nbchain,sigsize)

   for (j in 1:nbchain){
      for(i in 1:chain[j,1]) {
        cat("The (b,a) at chain", j,"is (",chain[j,2*i+1],chain[j,2*i],")\n")
        tmp3[j,chain[j,2*i+1]] <- tmp3[j,chain[j,2*i+1]]+
                                  2*tfinput[chain[j,2*i+1],chain[j,2*i]]
      }
      rec <- rec + Re(tmp3[j,])
   }

   if(plot == 1){
      par(mfrow=c(2,1))
      par(cex=1.1)
      plot.ts(siginput)
      title("Original signal")
      plot.ts(rec)
      title("Reconstructed signal")
   }
   else if (plot == 2){
      par(mfrow=c(nbchain+2,1))
      par(mar=c(2,4,4,4))
      par(cex=1.1)
      par(err=-1)
      plot.ts(siginput)
      title("Original signal")

      for (j in 1:nbchain)
         plot.ts(Re(tmp3[j,]))
        
      plot.ts(rec)
      title("Reconstructed signal")
   }
   list(rec=rec, ordered=tmp$ordered,chain = tmp$chain, comp=tmp3)
}





