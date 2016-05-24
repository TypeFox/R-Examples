rSPDistance <- function(x, from, to, theta, totalNet="net", method=1)
{
	if(theta < 0 | theta > 20 ) {stop("theta value out of range (between 0 and 20)")}
  if(method != 1 & method != 2) {stop("method should be either 1 or 2")}
	
	cellnri <- cellFromXY(x, from)
	cellnrj <- cellFromXY(x, to)
	transition <- .transitionSolidify(x)
	tc <- transitionCells(x)

	ci <- match(cellnri,tc)
	cj <- match(cellnrj,tc)
		
	tr <- transitionMatrix(x, inflate=FALSE)
	
	.rSPDist(tr, ci, cj, theta, totalNet, method)
}

.rSPDist <- function(tr, ci, cj, theta, totalNet, method)
{
	trR <- tr
	trR@x <- 1 / trR@x 
	nr <- dim(tr)[1] 
	Id <- Diagonal(nr) 
	P <- .normalize(tr, "row")
  

  
  if(method == 1)
  {

    W <- trR
	  W@x <- exp(-theta * trR@x) #zero values are not relevant because of next step exp(-theta * trR@x) 
    W <- W * P
    
  }	else	{
    
    adj <- adjacencyFromTransition(tr)
	  W <- trR
	  W[adj] <- exp(-theta * -log(P[adj])) #if the value is 1 then you get a natural random walk

  }
   
  #if (any(rowSums(W) < 1)) warning("one or more row sums of W are < 1")
  
	D <- matrix(0, nrow=length(ci), ncol=length(cj))
	
	for(j in 1:length(cj))
	{
	  Ij <- Diagonal(nr)
		Ij[cj[j],cj[j]] <- 0
		Wj <- Ij %*% W
		IdMinusWj <- as((Id - Wj), "dgCMatrix")		
		ej <- rep(0,times=nr)
		ej[cj[j]] <- 1
		zcj <- solve(IdMinusWj, ej)
    
		for(i in 1:length(ci))
		{
			ei <- rep(0,times=nr)
			ei[ci[i]] <- 1
			zci <- solve(t(IdMinusWj),ei)
			zcij <- sum(ei*zcj)

			N <- (Diagonal(nr, as.vector(zci)) %*% Wj %*% Diagonal(nr, as.vector(zcj))) / zcij
      
      if(totalNet == "net")
      {
        N <- skewpart(N) * 2 #N is here the NET number of passages, like McRae-random walk
        N@x[N@x<0] <- 0
      }

      # Computation of the cost dij between node i and node j
			D[i,j] <-  sum(trR * N)
			
    }
	}
	
	return(D)

}





  
  