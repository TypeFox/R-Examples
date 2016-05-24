##' getCounts function
##'
##' This function is used to count the number of observations falling inside grid cells, the output
##' is used in the function \link{lgcpPredict}.
##'
##' @param xyt stppp or ppp data object
##' @param subset Logical vector. Subset of data of interest, by default this is all data.
##' @param M number of centroids in x-direction 
##' @param N number of cnetroids in y-direction
##' @param ext how far to extend the grid eg (M,N) to (ext*M,ext*N)
##' @return The number of observations in each grid cell returned on a grid suitable for use in the extended FFT space.
##' @seealso \link{lgcpPredict}
##' @examples
##' require(spatstat)
##' xyt <- stppp(ppp(runif(100),runif(100)),t=1:100,tlim=c(1,100))
##' cts <- getCounts(xyt,M=64,N=64,ext=2) # gives an output grid of size 128 by 128
##' ctssub <- cts[1:64,1:64] # returns the cell counts in the observation
##'                          # window of interest
##' @export
getCounts <- function(xyt,subset=rep(TRUE,xyt$n),M,N,ext){
	if(M<5 | N<5){stop("M and/or N too small")}
	
	test1 <- FALSE
	test2 <- FALSE
	for (i in 1:20){
		if (2^i==M){test1 <- TRUE}
		if (2^i==N){test2 <- TRUE}
	}
	if ((!test1)&(!test2)){
		stop("Both M and N must be a power of 2 and both M,N <= 2^20")
	}
	
	xran <- xyt$window$xrange
	yran <- xyt$window$yrange
	bbox <- rbind(xran,yran)
	xwd <- diff(bbox[1,])/M
	ywd	<- diff(bbox[2,])/N
	sg <- SpatialGrid(GridTopology(c(bbox[1,1]+xwd/2,bbox[2,1]+ywd/2),c(xwd,ywd),c(M,N)))
	sp <- SpatialPoints(matrix(cbind(xyt$x,xyt$y)[subset,],sum(subset),2))
	#EJP: ol <- overlay(sg,sp)
	ol <- over(sp, sg)
	tol <- table(ol)
	idx <- as.numeric(rownames(tol))	
	smat <- matrix(0,M,N)
	smat[idx] <- tol
	smat <- smat[,N:1] # change orientation 
	
	nis <- matrix(0,ext*M,ext*N)
	nis[1:M,1:N] <- smat 
 
    return(nis)
}


##' getCellCounts function
##'
##' This function is used to count the number of observations falling inside grid cells.
##'
##' @param x x-coordinates of events
##' @param y y-coordinates of events
##' @param xgrid x-coordinates of grid centroids
##' @param ygrid y-coordinates of grid centroids
##' @return The number of observations in each grid cell.
##' @export
getCellCounts <- function(x,y,xgrid,ygrid){
	
	M <- length(xgrid)
	N <- length(ygrid)
	xwd <- diff(xgrid[1:2])
    ywd <- diff(ygrid[1:2])
	sg <- SpatialGrid(GridTopology(c(xgrid[1],ygrid[1]),c(xwd,ywd),c(M,N)))
	sp <- SpatialPoints(cbind(x,y))
	#EJP: ol <- overlay(sg,sp)
	ol <- over(sp,sg)
	tol <- table(ol)
	idx <- as.numeric(rownames(tol))	
	smat <- matrix(0,M,N)
	smat[idx] <- tol
	smat <- smat[,N:1] # change orientation 
	
	nis <- matrix(0,M,N)
	nis[1:M,1:N] <- smat 
 
    return(nis)
}
