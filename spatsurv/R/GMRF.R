##' getbb function
##'
##' A function to get the bounding box of a Spatial object
##'
##' @param obj a spatial object e.g. a SpatialPolygonsDataFrame, SpatialPolygons, etc ... anything with a bounding box that can be computed with bbox(obj)
##' @return a SpatialPolygons object: the bounding box
##' @export

getbb <- function(obj){
    bb <- bbox(obj)
    pg <- matrix(NA,5,2)
    pg[1,] <- c(bb[1,1],bb[2,1])
    pg[2,] <- c(bb[1,1],bb[2,2])
    pg[3,] <- c(bb[1,2],bb[2,2])
    pg[4,] <- c(bb[1,2],bb[2,1])
    pg[5,] <- pg[1,]
    bound <- Polygon(pg)
    bound <- Polygons(list(poly1=bound),ID=1)
    bound <- SpatialPolygons(list(poly=bound))
    proj4string(bound) <- CRS(proj4string(obj))
    return(bound)
}

##' getgrd function
##'
##' A function to create a regular grid over an observation window in order to model the spatial randome effects as a Gaussian
##' Markov random field.
##'
##' @param shape an object of class SpatialPolygons or SpatialPolygonsDataFrame 
##' @param cellwidth a scalar, the width of the grid cells
##' @return a SpatialPolygons object: the grid on which prediction of the spatial effects will occur
##' @references 
##' \enumerate{
##'     \item Benjamin M. Taylor. Auxiliary Variable Markov Chain Monte Carlo for Spatial Survival and Geostatistical Models. Benjamin M. Taylor. Submitted. \url{http://arxiv.org/abs/1501.01665}
##'     \item Finn Lindgren, Havard Rue, Johan Lindstrom. An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. Journal of the Royal Statistical Society: Series B 73(4)
##' }
##' @export

getgrd <- function(shape,cellwidth){
	shape <- gBuffer(shape,width=cellwidth)
	bb <- bbox(shape)
	xwid <- diff(bb[1,])
	ywid <- diff(bb[2,])

	delx <- xwid/cellwidth
	dely <- ywid/cellwidth

	M <- ceiling(delx)
	N <- ceiling(dely)

	xg <- bb[1,1] - (M*cellwidth-xwid)/2 + cellwidth*(0:M)
	yg <- bb[2,1] - (N*cellwidth-ywid)/2 + cellwidth*(0:N)

	spts <- SpatialPoints(expand.grid(xg,yg),proj4string=CRS(proj4string(shape)))
	ov <- over(spts,geometry(shape))
	spts <- spts[!is.na(ov),]

	return(as(SpatialPixels(spts),"SpatialPolygons"))

}


##' neighLocs function
##'
##' A function used in the computation of neighbours on non-rectangular grids. Not intended for general use.
##'
##' @param coord coordinate of interest
##' @param cellwidth a scalar, the width of the grid cells
##' @param order the order of the SPDE approximation: see Lindgren et al 2011 for details
##' @return coordinates of centroids of neighbours
##' @references 
##' \enumerate{
##'     \item Benjamin M. Taylor. Auxiliary Variable Markov Chain Monte Carlo for Spatial Survival and Geostatistical Models. Benjamin M. Taylor. Submitted. \url{http://arxiv.org/abs/1501.01665}
##'     \item Finn Lindgren, Havard Rue, Johan Lindstrom. An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. Journal of the Royal Statistical Society: Series B 73(4)
##' }
##' @export

neighLocs <- function(coord,cellwidth,order){
	M <- outer(abs(-order:order),abs(-order:order),"+")
	M[M>order] <- NA
	idx <- which(!is.na(M),arr.ind=TRUE)
	xv <- matrix(coord[1]+cellwidth*(-order:order),nrow=2*order+1,ncol=2*order+1,byrow=TRUE)
	yv <- matrix(coord[2]+cellwidth*(order:(-order)),nrow=2*order+1,ncol=2*order+1)
	return(cbind(xv[idx],yv[idx]))
}


##' neighOrder function
##'
##' A function to compute the order of a set of neighbours. Not intended for general use.
##'
##' @param neighlocs an object created by the function neighLocs
##' @return the neighbour orders
##' @references 
##' \enumerate{
##'     \item Benjamin M. Taylor. Auxiliary Variable Markov Chain Monte Carlo for Spatial Survival and Geostatistical Models. Benjamin M. Taylor. Submitted. \url{http://arxiv.org/abs/1501.01665}
##'     \item Finn Lindgren, Havard Rue, Johan Lindstrom. An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. Journal of the Royal Statistical Society: Series B 73(4)
##' }
##' @export

neighOrder <- function(neighlocs){
	mid <- neighlocs[(nrow(neighlocs)-1)/2+1,]
	del <- t(apply(neighlocs,1,function(x){x-mid}))

	md <- min(del[del>0]) # these three lines attempt to deal with floating point arithmetic issues with rank function
	del <- round(del/md)
	del <- del*md

	dis <- apply(del,1,function(x){x[1]^2+x[2]^2}) # no need to sqrt here
	rk <- rank(dis)
	tb <- table(rk)
	ds <- as.numeric(names(tb))
	ord <- 0:(length(tb)-1)
	
	return(sapply(rk,function(x){ord[which(ds==x)]}))
}


##' setupPrecMatStruct function
##'
##' A function to set up the computational grid and precision matrix structure for SPDE models.
##'
##' @param shape an object of class SpatialPolygons or SpatialPolygonsDataFrame 
##' @param cellwidth a scalar, the width of the grid cells 
##' @param no the order of the SPDE approximation: see Lindgren et al 2011 for details 
##' @return the computational grid and a function for constructing the precision matrix
##' @references 
##' \enumerate{
##'     \item Benjamin M. Taylor. Auxiliary Variable Markov Chain Monte Carlo for Spatial Survival and Geostatistical Models. Benjamin M. Taylor. Submitted. \url{http://arxiv.org/abs/1501.01665}
##'     \item Finn Lindgren, Havard Rue, Johan Lindstrom. An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. Journal of the Royal Statistical Society: Series B 73(4)
##' }
##' @export

setupPrecMatStruct <- function(shape,cellwidth,no){

	gr <- getgrd(shape,cellwidth)
	p4s <- proj4string(shape)

	ng <- neighLocs(coordinates(gr)[1,],cellwidth,no)
	nord <- neighOrder(ng)
	maxord <- max(nord)
	nneigh <- nrow(ng)
	npoly <- length(gr)

	index <- c()
	ng <- c()
	cat("Setting up computational grid ...\n")
	pb <- txtProgressBar(1,length(gr))
	for (i in 1:npoly){
		ng <- rbind(ng,neighLocs(coordinates(gr)[i,],cellwidth,no))
		setTxtProgressBar(pb, value=i)
	}
	cat("Done.\n")
	close(pb)
	idx <- over(SpatialPoints(ng,CRS(proj4string(gr))),geometry(gr))
	ind <- which(!is.na(idx))
	
	index <- cbind(rep(1:npoly,each=nneigh),idx,rep(nord,npoly))
	index <- index[ind,]
	index <- index[-which(index[,2]>index[,1]),]
	
	ord <- order(index[,3])
	index <- index[ord,]
	idxls <- lapply(0:maxord,function(x){which(index[,3]==x)})

	f <- function(fun){
		entries <- c()
		for(i in 0:maxord){
			entries[idxls[[i+1]]] <- fun(i)
		}
		return(sparseMatrix(i=index[,1],j=index[,2],x=entries,symmetric=TRUE))
	}
	attr(f,"order") <- no

	ans <- list()
	ans$f <- f
	ans$grid <- gr

	return(ans)

}



## GMRFprec function : PROBLEM: THIS CAN LEAD TO NON-POSITIVE DEFINITE MATRICES	
##
## A function to 
##
## @param par X 
## @return ...
## @export
# GMRFprec <- function(par){
# 	f <- function(i){
# 		return(par[i+1])
# 	}
# 	return(f)
# }



##' SPDEprec function
##'
##' A function to used in entering elements into the precision matrix of an SPDE model. Not intended for general use.
##'
##' @param a parameter a, see Lindgren et al 2011.
##' @param ord the order of the SPDE model, see Lindgren et al 2011.
##' @return a function used for creating the precision matrix
##' @references 
##' \enumerate{
##'     \item Benjamin M. Taylor. Auxiliary Variable Markov Chain Monte Carlo for Spatial Survival and Geostatistical Models. Benjamin M. Taylor. Submitted. \url{http://arxiv.org/abs/1501.01665}
##'     \item Finn Lindgren, Havard Rue, Johan Lindstrom. An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. Journal of the Royal Statistical Society: Series B 73(4)
##' }
##' @export

SPDEprec <- function(a,ord){
	if(ord==1){
		f <- function(i){
			if(i==0){
				return(a)
			}
			else if(i==1){
				return(-1)
			}
			else{
				stop("error in function SPDEprec")
			}			
		}
	}
	else if(ord==2){
		f <- function(i){
			if(i==0){
				return(4+a^2)
			}
			else if(i==1){
				return(-2*a)
			}
			else if(i==2){
				return(2)
			}
			else if(i==3){
				return(1)
			}
			else{
				stop("error in function SPDEprec")
			}
		}	
	}
	else if(ord==3){
		f <- function(i){
			if(i==0){
				return(a*(a^2+12))
			}
			else if(i==1){
				return(-3*(a^2+3))
			}
			else if(i==2){
				return(6*a)
			}
			else if(i==3){
				return(3*a)
			}
			else if(i==4){
				return(-3)
			}
			else if(i==5){
				return(-1)
			}
			else{
				stop("error in function SPDEprec")
			}
		}
	}
	else{
		stop("Higher order neighbourhood structures not currently supported.")
	}

	return(f)
}

##' YFromGamma_SPDE function
##'
##' A function to go from Gamma to Y
##'
##' @param gamma Gamma
##' @param U upper Cholesky matrix
##' @param mu the mean
##' @return the value of Y for the given Gamma
##' @references 
##' \enumerate{
##'     \item Benjamin M. Taylor. Auxiliary Variable Markov Chain Monte Carlo for Spatial Survival and Geostatistical Models. Benjamin M. Taylor. Submitted. \url{http://arxiv.org/abs/1501.01665}
##'     \item Finn Lindgren, Havard Rue, Johan Lindstrom. An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. Journal of the Royal Statistical Society: Series B 73(4)
##' }
##' @export

YFromGamma_SPDE <- function(gamma,U,mu){ # U= L^T
	return(mu+as.numeric(Matrix::solve(U,gamma)))
}

##' GammaFromY_SPDE function
##'
##' A function to go from Y to Gamma
##'
##' @param Y Y
##' @param U upper Cholesky matrix 
##' @param mu the mean 
##' @return the value of Gamma for the given Y
##' @references 
##' \enumerate{
##'     \item Benjamin M. Taylor. Auxiliary Variable Markov Chain Monte Carlo for Spatial Survival and Geostatistical Models. Benjamin M. Taylor. Submitted. \url{http://arxiv.org/abs/1501.01665}
##'     \item Finn Lindgren, Havard Rue, Johan Lindstrom. An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. Journal of the Royal Statistical Society: Series B 73(4)
##' }
##' @export

GammaFromY_SPDE <- function(Y,U,mu){ # U= L^T
	return(as.numeric(U%*%(Y-mu)))
}

