matrixFDM <- function(S, rf, q, vol, rho) {

  # check inputs
  if (nargs() != 5) {
    stop("incorrect number of arguments")
  }
  if(is.list(S) == FALSE) {
    stop(paste("S must be a list containing the vectors of the spatial",
               "grid points of the underlying assets"))
  } else {
    # determine number of underlying assets
    nAsset <- length(S)
  }
  if(!is.numeric(rf) | !is.numeric(q) | !is.numeric(vol) | 
       !is.numeric(rho)) {
    stop("rf, q, vol, and rho must be numeric")
  }
  if (length(rf) != 1) {
    stop("length of rf must be 1")
  }
  if (length(S) != nAsset | length(q) != nAsset | length(vol) != nAsset) {
    stop("length of S, q, and vol must equal the number of underlying assets")
  }
  if (sum(dim(rho) != nAsset) != 0) {
    stop("rho must be a square correlation matrix")
  }
  
	# calculate number of non-zero diagonals in FDM matrix and find main 
  # diagonal column number
	nCol <- 2 * nAsset^2 + 1
	ctr <- nCol / 2 + 0.5

	# find dimension of underlying spatial nodes, create difference vectors, 
  # create interior and boundary vectors
	naDiff <- 1
  dimS <- dSminus <- dSplus <- intS <- lBoundS <- rBoundS <- c()
	for (i in 1:nAsset) {
	  dimS <- c(dimS, length(S[[i]]))
  	dSminus <- c(dSminus, list(c(naDiff, diff(S[[i]]))))
  	dSplus <- c(dSplus, list(c(diff(S[[i]]), naDiff)))
  	intS <- c(intS, list(c(0, S[[i]][2:(length(S[[i]])-1)], 0)))
  	lBoundS <- c(lBoundS, list(c(S[[i]][1], rep(0,times=
                  length(S[[i]])-1))))
  	rBoundS <- c(rBoundS, list(c(rep(0,times=length(S[[i]])-1), 
                  S[[i]][length(S[[i]])])))
	}

	# initialize component matrices of mDiags
	mInterior <- mLBound <- mRBound <- mLL <- mIntL <- matrix(0, prod(dimS), nCol)
  mRL <- mLInt <- mLR <- mRInt <- mIntR <- mRR <- matrix(0, prod(dimS), nCol)

	# populate component matricies
	for (i in 1:nAsset) {
    
    # first and second spatial difference terms for given underlying asset i
		mInterior[,(ctr-(i*(i-1)+1))] <- mInterior[,(ctr-(i*(i-1)+1))] + 
      (((vol[i]^2 * expand.grid(intS)[i]^2) / (expand.grid(dSminus)[i] * 
      (expand.grid(dSminus)[i] + expand.grid(dSplus)[i]))) - (((rf - q[i]) * 
      expand.grid(intS)[i]) / (expand.grid(dSminus)[i] + expand.grid(dSplus)
      [i])))[,1]
		mInterior[,ctr] <- mInterior[,ctr] - (((vol[i]^2 * expand.grid(intS)[i]^2) / 
      (expand.grid(dSminus)[i] * expand.grid(dSplus)[i])))[,1]
		mInterior[,(ctr+(i*(i-1)+1))] <- mInterior[,(ctr+(i*(i-1)+1))] + 
      (((vol[i]^2 * expand.grid(intS)[i]^2)  / (expand.grid(dSplus)[i] * 
      (expand.grid(dSminus)[i] + expand.grid(dSplus)[i]))) + (((rf - q[i]) * 
      expand.grid(intS)[i]) / (expand.grid(dSminus)[i] + expand.grid(dSplus)
      [i])))[,1]
		
		mLBound[,ctr] <- mLBound[,ctr] + (-((rf - q[i]) * expand.grid(lBoundS)[i]) / 
      expand.grid(dSplus)[i])[,1]
		mLBound[,(ctr+(i*(i-1)+1))] <- mLBound[,(ctr+(i*(i-1)+1))] + (((rf - q[i]) * 
      expand.grid(lBoundS)[i]) / expand.grid(dSplus)[i])[,1]

		mRBound[,(ctr-(i*(i-1)+1))] <- mRBound[,(ctr-(i*(i-1)+1))] + (-((rf - q[i]) *
      expand.grid(rBoundS)[i]) / expand.grid(dSminus)[i])[,1]
		mRBound[,ctr] <- mRBound[,ctr] + (((rf - q[i]) * expand.grid(rBoundS)[i]) /
      expand.grid(dSminus)[i])[,1]

		# pairwise correlation terms for given underlying asset i
		if (i > 1) {
			for (j in 1:(i-1)) {

				mInterior[,(ctr-(i*(i-1)+1)-j)] <- mInterior[,(ctr-(i*(i-1)+1)-j)] + 
          ((rho[i,j] * vol[j] * vol[i] * expand.grid(intS)[j] * expand.grid
          (intS)[i]) / ((expand.grid(dSminus)[j] + expand.grid(dSplus)[j]) *
          (expand.grid(dSminus)[i] + expand.grid(dSplus)[i])))[,1]
				mInterior[,(ctr-(i*(i-1)+1)+j)] <- mInterior[,(ctr-(i*(i-1)+1)+j)] -
          ((rho[i,j] * vol[j] * vol[i] * expand.grid(intS)[j] * expand.grid
          (intS)[i]) / ((expand.grid(dSminus)[j] + expand.grid(dSplus)[j]) *
          (expand.grid(dSminus)[i] + expand.grid(dSplus)[i])))[,1]
				mInterior[,(ctr+(i*(i-1)+1)-j)] <- mInterior[,(ctr+(i*(i-1)+1)-j)] -
          ((rho[i,j] * vol[j] * vol[i] * expand.grid(intS)[j] * expand.grid
          (intS)[i]) / ((expand.grid(dSminus)[j] + expand.grid(dSplus)[j]) *
          (expand.grid(dSminus)[i] + expand.grid(dSplus)[i])))[,1]
				mInterior[,(ctr+(i*(i-1)+1)+j)] <- mInterior[,(ctr+(i*(i-1)+1)+j)] +
          ((rho[i,j] * vol[j] * vol[i] * expand.grid(intS)[j] * expand.grid
          (intS)[i]) / ((expand.grid(dSminus)[j] + expand.grid(dSplus)[j]) *
          (expand.grid(dSminus)[i] + expand.grid(dSplus)[i])))[,1]
				
				mLL[,ctr] <- mLL[,ctr]  + ((rho[i,j] * vol[j] * vol[i] * expand.grid
          (lBoundS)[j] * expand.grid(lBoundS)[i]) / ((expand.grid(dSplus)[j]) 
          * (expand.grid(dSplus)[i])))[,1] 
				mLL[,(ctr+(j*(j-1)+1))] <- mLL[,(ctr+(j*(j-1)+1))] - ((rho[i,j] * 
          vol[j] * vol[i] * expand.grid(lBoundS)[j] * expand.grid(lBoundS)[i]) 
          / ((expand.grid(dSplus)[j]) * (expand.grid(dSplus)[i])))[,1] 
				mLL[,(ctr+(i*(i-1)+1))] <- mLL[,(ctr+(i*(i-1)+1))] - ((rho[i,j] * 
          vol[j] * vol[i] * expand.grid(lBoundS)[j] * expand.grid(lBoundS)[i]) 
          / ((expand.grid(dSplus)[j]) * (expand.grid(dSplus)[i])))[,1]   
				mLL[,(ctr+(i*(i-1)+1)+j)] <- mLL[,(ctr+(i*(i-1)+1)+j)] + ((rho[i,j] * 
          vol[j] * vol[i] * expand.grid(lBoundS)[j] * expand.grid(lBoundS)[i]) 
          / ((expand.grid(dSplus)[j]) * (expand.grid(dSplus)[i])))[,1]  

				mIntL[,(ctr-(j*(j-1)+1))] <- mIntL[,(ctr-(j*(j-1)+1))] + ((rho[i,j] * 
          vol[j] * vol[i] * expand.grid(intS)[j] * expand.grid(lBoundS)[i]) / 
          ((expand.grid(dSminus)[j] + expand.grid(dSplus)[j]) * (expand.grid
          (dSplus)[i])))[,1]
				mIntL[,(ctr+(j*(j-1)+1))] <- mIntL[,(ctr+(j*(j-1)+1))] - ((rho[i,j] * 
          vol[j] * vol[i] * expand.grid(intS)[j] * expand.grid(lBoundS)[i]) / 
          ((expand.grid(dSminus)[j] + expand.grid(dSplus)[j]) * (expand.grid(
          dSplus)[i])))[,1]
				mIntL[,(ctr+(i*(i-1)+1)-j)] <- mIntL[,(ctr+(i*(i-1)+1)-j)] - ((rho[i,j] 
          * vol[j] * vol[i] * expand.grid(intS)[j] * expand.grid(lBoundS)[i]) / 
          ((expand.grid(dSminus)[j] + expand.grid(dSplus)[j]) * (expand.grid
          (dSplus)[i])))[,1]
				mIntL[,(ctr+(i*(i-1)+1)+j)] <- mIntL[,(ctr+(i*(i-1)+1)+j)] + ((rho[i,j] 
          * vol[j] * vol[i] * expand.grid(intS)[j] * expand.grid(lBoundS)[i]) / 
          ((expand.grid(dSminus)[j] + expand.grid(dSplus)[j]) * (expand.grid
          (dSplus)[i])))[,1]

				mRL[,(ctr-(j*(j-1)+1))] <- mRL[,(ctr-(j*(j-1)+1))] + ((rho[i,j] * vol[j] 
          * vol[i] * expand.grid(rBoundS)[j] * expand.grid(lBoundS)[i]) / 
          ((expand.grid(dSminus)[j]) * (expand.grid(dSplus)[i])))[,1]
				mRL[,ctr] <- mRL[,ctr] - ((rho[i,j] * vol[j] * vol[i] * expand.grid
          (rBoundS)[j] * expand.grid(lBoundS)[i]) / ((expand.grid(dSminus)[j]) * 
          (expand.grid(dSplus)[i])))[,1]
				mRL[,(ctr+(i*(i-1)+1)-j)] <- mRL[,(ctr+(i*(i-1)+1)-j)] - ((rho[i,j] * 
          vol[j] * vol[i] * expand.grid(rBoundS)[j] * expand.grid(lBoundS)[i]) / 
          ((expand.grid(dSminus)[j]) * (expand.grid(dSplus)[i])))[,1]
				mRL[,(ctr+(i*(i-1)+1))] <- mRL[,(ctr+(i*(i-1)+1))] + ((rho[i,j] * vol[j] 
          * vol[i] * expand.grid(rBoundS)[j] * expand.grid(lBoundS)[i]) / 
          ((expand.grid(dSminus)[j]) * (expand.grid(dSplus)[i])))[,1]

				mLInt[,(ctr-(i*(i-1)+1))] <- mLInt[,(ctr-(i*(i-1)+1))] + ((rho[i,j] * 
          vol[i] * vol[j] * expand.grid(lBoundS)[i] * expand.grid(intS)[j]) / 
          ((expand.grid(dSplus)[j]) * (expand.grid(dSminus)[i] + expand.grid
          (dSplus)[i])))[,1]
				mLInt[,(ctr-(i*(i-1)+1)+j)] <- mLInt[,(ctr-(i*(i-1)+1)+j)] - ((rho[i,j] 
          * vol[i] * vol[j] * expand.grid(lBoundS)[i] * expand.grid(intS)[j]) / 
          ((expand.grid(dSplus)[j]) * (expand.grid(dSminus)[i] + expand.grid
          (dSplus)[i])))[,1]
				mLInt[,(ctr+(i*(i-1)+1))] <- mLInt[,(ctr+(i*(i-1)+1))] - ((rho[i,j] * 
          vol[i] * vol[j] * expand.grid(lBoundS)[i] * expand.grid(intS)[j]) / 
          ((expand.grid(dSplus)[j]) * (expand.grid(dSminus)[i] + expand.grid
          (dSplus)[i])))[,1]
				mLInt[,(ctr+(i*(i-1)+1)+j)] <- mLInt[,(ctr+(i*(i-1)+1)+j)] + ((rho[i,j] 
          * vol[i] * vol[j] * expand.grid(lBoundS)[i] * expand.grid(intS)[j]) / 
          ((expand.grid(dSplus)[j]) * (expand.grid(dSminus)[i] + expand.grid
          (dSplus)[i])))[,1]

				mLR[,(ctr-(i*(i-1)+1))] <- mLR[,(ctr-(i*(i-1)+1))] + ((rho[i,j] * 
          vol[j] * vol[i] * expand.grid(lBoundS)[j] * expand.grid(rBoundS)[i]) / 
          ((expand.grid(dSplus)[j]) * (expand.grid(dSminus)[i])))[,1]
				mLR[,(ctr-(i*(i-1)+1)+j)] <- mLR[,(ctr-(i*(i-1)+1)+j)] - ((rho[i,j] * 
          vol[j] * vol[i] * expand.grid(lBoundS)[j] * expand.grid(rBoundS)[i]) / 
          ((expand.grid(dSplus)[j]) * (expand.grid(dSminus)[i])))[,1]
				mLR[,ctr] <- mLR[,ctr] - ((rho[i,j] * vol[j] * vol[i] * expand.grid
          (lBoundS)[j] * expand.grid(rBoundS)[i]) / ((expand.grid(dSplus)[j]) * 
          (expand.grid(dSminus)[i])))[,1]
				mLR[,(ctr+(j*(j-1)+1))] <- mLR[,(ctr+(j*(j-1)+1))] + ((rho[i,j] * vol[j] 
          * vol[i] * expand.grid(lBoundS)[j] * expand.grid(rBoundS)[i]) / 
          ((expand.grid(dSplus)[j]) * (expand.grid(dSminus)[i])))[,1]

				mRInt[,(ctr-(i*(i-1)+1)-j)] <- mRInt[,(ctr-(i*(i-1)+1)-j)] + ((rho[i,j] 
          * vol[j] * vol[i] * expand.grid(rBoundS)[j] * expand.grid(intS)[i]) / 
          ((expand.grid(dSminus)[j]) * (expand.grid(dSminus)[i] + expand.grid
          (dSplus)[i])))[,1]
				mRInt[,(ctr-(i*(i-1)+1))] <- mRInt[,(ctr-(i*(i-1)+1))] - ((rho[i,j] * 
          vol[j] * vol[i] * expand.grid(rBoundS)[j] * expand.grid(intS)[i]) / 
          ((expand.grid(dSminus)[j]) * (expand.grid(dSminus)[i] + expand.grid
          (dSplus)[i])))[,1]
				mRInt[,(ctr+(i*(i-1)+1)-j)] <- mRInt[,(ctr+(i*(i-1)+1)-j)] - ((rho[i,j] 
          * vol[j] * vol[i] * expand.grid(rBoundS)[j] * expand.grid(intS)[i]) / 
          ((expand.grid(dSminus)[j]) * (expand.grid(dSminus)[i] + expand.grid
          (dSplus)[i])))[,1]
				mRInt[,(ctr+(i*(i-1)+1))] <- mRInt[,(ctr+(i*(i-1)+1))] + ((rho[i,j] * 
          vol[j] * vol[i] * expand.grid(rBoundS)[j] * expand.grid(intS)[i]) / 
          ((expand.grid(dSminus)[j]) * (expand.grid(dSminus)[i] + expand.grid
          (dSplus)[i])))[,1]

				mIntR[,(ctr-(i*(i-1)+1)-j)] <- mIntR[,(ctr-(i*(i-1)+1)-j)] + ((rho[i,j] 
          * vol[j] * vol[i] * expand.grid(intS)[j] * expand.grid(rBoundS)[i]) / 
          ((expand.grid(dSminus)[j] + expand.grid(dSplus)[j]) * (expand.grid
          (dSminus)[i])))[,1]
				mIntR[,(ctr-(i*(i-1)+1)+j)] <- mIntR[,(ctr-(i*(i-1)+1)+j)] - ((rho[i,j] 
          * vol[j] * vol[i] * expand.grid(intS)[j] * expand.grid(rBoundS)[i]) / 
          ((expand.grid(dSminus)[j] + expand.grid(dSplus)[j]) * (expand.grid
          (dSminus)[i])))[,1]
				mIntR[,(ctr-(j*(j-1)+1))] <- mIntR[,(ctr-(j*(j-1)+1))] - ((rho[i,j] * 
          vol[j] * vol[i] * expand.grid(intS)[j] * expand.grid(rBoundS)[i]) / 
          ((expand.grid(dSminus)[j] + expand.grid(dSplus)[j]) * (expand.grid
          (dSminus)[i])))[,1]
				mIntR[,(ctr+(j*(j-1)+1))] <- mIntR[,(ctr+(j*(j-1)+1))] + ((rho[i,j] * 
          vol[j] * vol[i] * expand.grid(intS)[j] * expand.grid(rBoundS)[i]) / 
          ((expand.grid(dSminus)[j] + expand.grid(dSplus)[j]) * (expand.grid
          (dSminus)[i])))[,1]

				mRR[,(ctr-(i*(i-1)+1)-j)] <- mRR[,(ctr-(i*(i-1)+1)-j)] + ((rho[i,j] * 
          vol[j] * vol[i] * expand.grid(rBoundS)[j] * expand.grid(rBoundS)[i]) 
          / ((expand.grid(dSminus)[j]) * (expand.grid(dSminus)[i])))[,1] 
				mRR[,(ctr-(i*(i-1)+1))] <- mRR[,(ctr-(i*(i-1)+1))] - ((rho[i,j] * 
          vol[j] * vol[i] * expand.grid(rBoundS)[j] * expand.grid(rBoundS)[i]) 
          / ((expand.grid(dSminus)[j]) * (expand.grid(dSminus)[i])))[,1] 
				mRR[,(ctr-(j*(j-1)+1))] <- mRR[,(ctr-(j*(j-1)+1))] - ((rho[i,j] * vol[j] 
          * vol[i] * expand.grid(rBoundS)[j] * expand.grid(rBoundS)[i]) / 
          ((expand.grid(dSminus)[j]) * (expand.grid(dSminus)[i])))[,1] 
				mRR[,ctr] <- mRR[,ctr]  + ((rho[i,j] * vol[j] * vol[i] * expand.grid
          (rBoundS)[j] * expand.grid(rBoundS)[i]) / ((expand.grid(dSminus)[j]) 
          * (expand.grid(dSminus)[i])))[,1] 

			}
		}
	}

	# combine component matrices
	mDiags <- mInterior + mLBound + mRBound + mLL + mIntL + mRL + mLInt + 
              mLR + mRInt + mIntR + mRR

	# subtract discount factor
	mDiags[,ctr] <- mDiags[,ctr] - rf
  
  # adjust columns of mDiags to comply with bandSparse function "diagonals"
  # argument
  for (i in 1:nAsset) {
   	if (i == 1) {
    			mDiags2 <- c(mDiags[2:prod(dimS),(ctr-1)],0)
  		} else {
    			for (j in 1:(i-1)) {
      			if (j == 1) {
        				mDiags2_temp <- c(mDiags[(prod(dimS[1:(i-1)])+2):prod(dimS),
                                  (ctr-(i*(i-1)+1)-1)], rep(0, times=prod
                                  (dimS[1:(i-1)])+1))
       					mDiags2_temp <- cbind(mDiags2_temp, c(mDiags[(prod(dimS[1:
                                  (i-1)])+1):prod(dimS),(ctr-(i*(i-1)+1)+0)],
                                  rep(0, times=prod(dimS[1:(i-1)])+0)))
        				mDiags2_temp <- cbind(mDiags2_temp, c(mDiags[(prod(dimS[1:
                                  (i-1)])+0):prod(dimS),(ctr-(i*(i-1)+1)+1)],
                                  rep(0, times=prod(dimS[1:(i-1)])-1)))
      			} else {
        				mDiags2_temp <- cbind(c(mDiags[(prod(dimS[1:(i-1)])+
                                  dimS[j-1]+1):prod(dimS),(ctr-(i*(i-1)
                                  +1)-j)], rep(0, times=prod(dimS[1:(i-1)])
                                  +dimS[j-1])), mDiags2_temp, c(mDiags
                                  [(prod(dimS[1:(i-1)])-dimS[j-1]+1):
                                  prod(dimS),(ctr-(i*(i-1)+1)+j)], rep(0, 
                                  times=prod(dimS[1:(i-1)])-dimS[j-1])))
      			}
    			}
    			mDiags2 <- cbind(mDiags2_temp, mDiags2)
  		}
  }
  mDiags2 <- cbind(mDiags2, mDiags[,ctr:nCol]) 

	# calculate which diagonals of sparse banded matrix to populate
	for (i in 1:nAsset) {
		if (i == 1) {
  			cDiags <- 1
		} else {
  			for (j in 1:(i-1)) {
    				if (j == 1) {
      				cDiags_temp <- prod(dimS[1:(i-1)]) + 1
    				} else {
      				cDiags_temp <- c(cDiags_temp, prod(dimS[1:(i-1)]) 
                               + dimS[j-1])
    				}
  			}
  			cDiags <- c(cDiags, 2 * (prod(dimS[1:(i-1)])) - rev(cDiags_temp), 
                    prod(dimS[1:(i-1)]), cDiags_temp)
		}
	}
	cDiags <- c(-rev(cDiags), 0, cDiags)
  
  # create sparse banded matrix mO from columns of mDiags2
 	mO <- bandSparse(prod(dimS), prod(dimS), cDiags, mDiags2)
  
	# return FDM matrix
	mO

}