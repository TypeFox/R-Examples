drawMahal <- 
		function(x,x.mean,x.cov,whichlines=c(0.975,0.90,0.75),m=360,plot=TRUE,
				xlab="",ylab="",cex.lab=1.2, ...){
# draw ellipse corresponding to Mahalanobis distances
# PF, 06.12.2006
#
# x ... data
# x.mean ... mean of data
# x.cov ... cov-matrix of data
# whichlines: lines corresponding to quantiles of chi2
# m: number of data points on the ellipse
# plot ... if plot=TRUE the lines are plotted, otherwise not
#
	
# OUTPUT
# objects "mdX" and "mdY" with m rows and length(whichlines) columns
# (the columns refer to X- and Y-coordinates of one ellipse)
	
	mdX=matrix(NA,nrow=m,ncol=length(whichlines))
	mdY=matrix(NA,nrow=m,ncol=length(whichlines))
	
	cov.svd <- svd(x.cov, nv = 0)
	r <- cov.svd[["u"]] %*% diag(sqrt(cov.svd[["d"]]))
	
	alphamd <- sqrt(qchisq(whichlines,2))
	
	for (j in 1:length(whichlines)){
		e1md <- cos(c(0:(m-1))/m * 2 * pi) * alphamd[j]
		e2md <- sin(c(0:(m-1))/m * 2 * pi) * alphamd[j]
		emd <- cbind(e1md, e2md)
		ttmd <- t(r %*% t(emd)) + rep(1, m) %o% x.mean
		
		mdX[,j] <- ttmd[,1]
		mdY[,j] <- ttmd[,2]
		
		if (plot) lines(mdX[,j],mdY[,j], type = "l", ...)
	}
	
	list(mdX=mdX,mdY=mdY)
}

impAll <-
		function(x) {
	maxLimits = apply(x, 2, min, na.rm = TRUE)
	maxLimits = -maxLimits
	maxLimits[maxLimits < 0] = 0
	
	if(any(is.na(x))) {
		temp <- x
		for(i in 1:length(maxLimits)) {
			temp[na.omit(temp[i]) < 0, i] = maxLimits[i] * 2 / 3
		}
		
		res <- adjust(impCoda(temp))
		isna = is.na(x)
		x[isna] = res$xImp[isna]
	}
	if(any(x < 0)) {
		temp <- x
		temp[temp < 0] = 0
		
		res <- impRZilr(temp, dl=maxLimits)
		x = res$xImp
	}
	
	return(x)
}

genTest <- function(){
  wind <- matrix(c(F, F, T, F, F, T, F, 
                   F, T, F, F, F, T, F,
                   T, T, F, F, F, F, F), ncol=7, byrow=TRUE)
  orig <- matrix(c(3, 5, NA, 8, 10, NA, 6,
                   34,NA,99, 4, 32, NA, 77,
                   NA,NA,3,4,5,6,7), ncol=7, byrow=TRUE)   
  xOrig <- matrix(c(3, 5, 6.5, 8, 10, 7, 6,
                   34,44.3,99,4,32,11.3,77,
                   0.75,1.25,3,4,5,6,7), ncol=7, byrow=TRUE)  
  xImp <- matrix(c(3, 5, 7, 8, 10, 7.77, 6,
                  36.2,44.8,100,4.5,33,12,78.8,
                  0.11,0.2,0.4,0.6,0.8,0.95,1.2), ncol=7, byrow=TRUE)
  result <- matrix(c(3, 5, 7, 8, 10, 7.77, 6,
                     36.2,44.8,100,4.5,33,12,78.8,
                     34,44.8*246/252.5,99,4,32,12*246/252.5,77), ncol=7, byrow=TRUE)
  list(wind=wind, orig=orig, xOrig=xOrig, xImp=xImp, result=result)
}
# 
# adjustImputed <- function(xImp, xOrig, wind){
#   ## aim: 
#   ## (1) ratios must be preserved
#   ## (2) do not change original values
#   ## (3) adapt imputations
#   xneu  <- xImp
#   s1 <- rowSums(xOrig, na.rm = TRUE)
#   ## per row: consider rowsums of imputed data
#   sumPrevious <- sumAfter <- numeric(nrow(xImp))
#   for (i in 1:nrow(xImp)) {
#     if(any(wind[i, ]) & !all(wind[i,])){
#       sumPrevious[i] <- sum(xOrig[i, !wind[i, ]]) 
#       sumAfter[i] <- sum(xImp[i, !wind[i,]])
#     } else{ 
#       sumPrevious[i] <- sumAfter[i] <- 1
#     }
#   }
#   # how much is rowsum increased by imputation:
#   fac <- sumPrevious/sumAfter
#   
#   #    # decrese rowsums of orig.
#   #    s1[i] <- s1[i]/fac
#   #  }
#   ## for non-zeros overwrite them:
#   xneu[!wind] <- xOrig[!wind]
#   ## adjust zeros:
#   for(i in 1:nrow(xImp)){
#     if(any(wind[i,])){
#       xneu[i,wind[i,]] <- fac[i]*xneu[i,wind[i,]]
#     }
#   }
#   return(xneu)
# }




adjustImps <- function (xImp, xOrig, wind){
  ## aim: 
  ## (1) ratios must be preserved
  ## (2) do not change original values
  ## (3) adapt imputations
  xneu  <- xImp
  s1 <- rowSums(xOrig, na.rm = TRUE)
  ## per row: consider rowsums of imputed data
  ## example: 
  ## wind: wind <- matrix(c(F, F, T, F, F, T, F), ncol=7)
  ## ganz orig:  orig <- matrix(c(3, 5, NA, 8, 10, NA, 6), ncol=7)   (sum=26)
  ## orig(init): xOrig <- matrix(c(3, 5, 6.5, 8, 10, 7, 6), ncol=7)  (sum=32.5)
  ## imp:        xImp <- matrix(c(3, 5, 7, 8, 10, 7.77, 6), ncol=7)    (sum=33)
  ## s: 26
  ## s2: 7
  ## fac: 26/(26+7)
  ## s1: 32.5/(26/(26+7)) =  41.25
  for (i in 1:nrow(xImp)) {
    if(any(wind[i,])) s <- sum(xImp[i, !wind[i, ]]) else s <- 1
    if(any(wind[i,])) s2 <- sum(xImp[i, wind[i, ]]) else s2 <- 0
    # how much is rowsum increased by imputation:
    fac <- s/(s + s2)
    # decrese rowsums of orig.
    s1[i] <- s1[i]/fac
  }
  ## impS: 41.25/33
  impS <- s1/rowSums(xImp)
  for (i in 1:ncol(xImp)) {
    xneu[, i] <- xImp[, i] * impS
  }
  xImp <- xneu
  return(xImp)
}


