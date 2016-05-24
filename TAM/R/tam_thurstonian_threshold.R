####################################################################
# Thurstonian thresholds (gammas)
#...................................................................
tam.threshold <- function (tamobj, prob.lvl=0.5)
{
  resp <- tamobj$resp
  nitems <- tamobj$nitems
  maxK <- tamobj$maxK
  AXsi <- tamobj$AXsi
  xsi <- ( tamobj$xsi )[,1]
  A <- tamobj$A
  B <- tamobj$B
#  if ( dim(B)[3] > 1){
#	stop("Thurstonian thresholds are only calculated for unidimensional models.")
#					}
  
  maxKi <- apply( resp , 2 , max , na.rm=TRUE )
  
  threshold <- matrix(-99, nitems, maxK-1)
  for (i in 1:nitems) {
    mc <- maxKi[i]  # max category value for item i
    initlow <- rep(-12, mc)
    inithigh <- rep(12, mc)
    lowhigh <- matrix(c(initlow, inithigh), nrow=mc, ncol=2)
    thresh <- rowMeans(lowhigh)
    oldthresh <- rep(-99, length(thresh))
    cprobs <- rep(0,mc)                 
	B0 <- B
	if ( dim(B)[[3]] > 1){
		D <- dim(B)[[3]]
		ind <- which( B[i,2,] > 0 )[1]
		sel2 <- c( ind , setdiff( 1:D , ind ) )
		B0 <- B[,,sel2]	
		}
	
    while (max(abs(oldthresh-thresh)) > 0.0001) {
      res.p <- calc_prob.v5( iIndex=i:i , A=A , AXsi=AXsi , B=B0, 
                             xsi=xsi , theta=matrix(thresh,nrow=mc,ncol=1) , 
							 nnodes=mc, maxK=maxK)        	
      rprobs <- res.p[["rprobs"]]
      for (k in 1:mc) {
        cprobs[k] <- sum(rprobs[1,(k+1):(mc+1),k], na.rm=TRUE)
      }
      lowhigh[cprobs<prob.lvl,1] <- thresh[cprobs<prob.lvl]
      lowhigh[cprobs>prob.lvl,2] <- thresh[cprobs>prob.lvl]
      oldthresh <- thresh
      thresh <- rowMeans(lowhigh)
      #      lowhigh[abs(thresh-initlow)<0.005,1] <- lowhigh[abs(thresh-initlow)<0.005,1] - 12
      #      initlow[abs(thresh-initlow)<0.005] <- initlow[abs(thresh-initlow)<0.005] - 12
      #      lowhigh[abs(thresh-inithigh)<0.005,2] <- lowhigh[abs(thresh-inithigh)<0.005,2] + 12
      #      inithigh[abs(thresh-inithigh)<0.005] <- inithigh[abs(thresh-inithigh)<0.005] + 12
    } # end of while loop
    
    threshold[i, 1:mc] <- thresh
  } #end of all items
  threshold[threshold == -99] <- NA
  rownames(threshold) <- colnames(tamobj$resp)
  colnames(threshold) <- paste0("Cat" , 1:ncol(threshold))
  return(threshold)
}
