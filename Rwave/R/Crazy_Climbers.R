#########################################################################
#      $Log: Crazy_Climbers.S,v $
#########################################################################
#
#               (c) Copyright  1997
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University
#                  All right reserved                           
#########################################################################

crc <- function(tfrep, tfspec = numeric(dim(tfrep)[2]), bstep = 3,
                iteration = 10000, rate = .001, seed = -7, nbclimb=10,
                flag.int = TRUE, chain = TRUE, flag.temp = FALSE)
#########################################################################
#  crc:   Time-frequency multiple ridge estimation (crazy climbers)
#  ----
# 	use the crazy climber algorithm to evaluate ridges of
#    	   continuous Time-Frequency transform
#
#      input:
#      ------
# 	tfrep: wavelet or Gabor transform
#	tfspec: additional potential (coming from learning the noise)
#       iteration: number of iterations
#       rate: initial value of the temperature
#       seed: initialization for random numbers generator
#       nbclimb: number of crazy climbers
#       bstep: step size for the climber walk in the time direction
#	flag.int: if set to TRUE, computes the integral on the ridge.
#	chain: if set to TRUE, chains the ridges.
#	flag.temp: if set to TRUE, keeps a constant temperature.
#
#      output:
#      -------
#       beemap: 2D array containing the (weighted or unweighted)
#               occupation measure (integrated with respect to time)
#
#########################################################################
{
  tfspectrum <- tfspec 
  
  d <- dim(tfrep)
  sigsize <- d[1]
  nscale <- d[2]
  beemap <- matrix(0,sigsize,nscale)
  sqmodulus <- Re(tfrep*Conj(tfrep))
  for (k in 1:nscale)
    sqmodulus[,k] <- sqmodulus[,k] - tfspectrum[k]
  dim(beemap) <- c(nscale * sigsize, 1)
  dim(sqmodulus) <- c(nscale * sigsize, 1)

  z <- .C("Sbee_annealing",
          as.double(sqmodulus),
          beemap= as.double(beemap),
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
  dim(beemap) <- c(sigsize,nscale)   
  if(dev.set() != 1) image(beemap)
  beemap
}

cfamily <- function(ccridge, bstep = 1, nbchain = 100, ptile = 0.05)
#########################################################################
#     cfamily:
#     --------
#     chain the ridges obtained by crazy climber.
#
#      input:
#      ------
#       ccridge: unchained ridge (output of Bee_Annealing)
#       bstep: maximal length for a gap in a ridge
#       nbchain: maximum number of chains
#       ptile: relative threshold for the ridges 
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
  d <- dim(ccridge)
  sigsize <- d[1]
  nscale <- d[2]
  threshold <- range(ccridge)[2] * ptile
  sz <- sigsize + 2   	
  chain <- matrix(-1,nbchain,sz)
  orderedmap <- matrix(0,sigsize,nscale)	
  
  dim(chain) <- c(nbchain * sz, 1)
  dim(ccridge) <- c(nscale * sigsize, 1)
  dim(orderedmap) <- c(nscale * sigsize, 1)
   
  z <- .C("Scrazy_family",
          as.double(ccridge),
          orderedmap = as.double(orderedmap),
          chain = as.integer(chain),
          chainnb = as.integer(nbchain),
          as.integer(sigsize),
          as.integer(nscale),
          as.integer(bstep),
          as.double(threshold),
          PACKAGE="Rwave")
	
  orderedmap <- z$orderedmap
  chain <- z$chain

  dim(orderedmap) <- c(sigsize,nscale)   
  dim(chain) <- c(nbchain,sz)   
  nbchain <- z$chainnb
  chain <- chain +1
  chain[,2] <- chain[,2]-1
  
  list(ordered = orderedmap, chain = chain, nbchain=nbchain)
}   

crfview <- function(beemap, twod=TRUE)
#########################################################################
#     crfview:
#     --------
#     displays a family of chained ridges
#
#########################################################################
{
  if (twod)
    image(beemap$ordered > 0)
  else {
    Dim <- dim(beemap$ordered)
    matmp <- matrix(0,Dim[1],Dim[2])
    matmp[1,1] <- 1
    image(matmp)
    for (k in 1:beemap$nbchain){
      start <- beemap$chain[k,1]
      end <- start + beemap$chain[k,2]  
      lines(start:(end-1), beemap$chain[k,3:(end-start+2)], type="l")
    }
  }
  title(" Chained Ridges")
  cat("")
}
