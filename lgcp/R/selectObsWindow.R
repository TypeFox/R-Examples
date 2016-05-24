##' selectObsWindow function
##'
##'
##' See ?selectObsWindow.stppp for further details on usage. This is a generic function for the purpose of selecting an observation window 
##' (or more precisely a bounding box) to contain the extended FFT grid.
##'
##' @param xyt an object  
##' @param ... additional arguments  
##' @return method selectObsWindow
##' @seealso \link{selectObsWindow.default}, \link{selectObsWindow.stppp}
##' @export

selectObsWindow <- function(xyt,...){
    UseMethod("selectObsWindow")
}



##' selectObsWindow.default function
##'
##' Default method, note at present, there is only an implementation for stppp objects.
##'
##' !!NOTE!! that this function also returns the grid ($xvals and $yvals) on which the FFT (and hence MALA) will be performed. It is useful to
##' define spatialAtRiskobjects on this grid to prevent loss of information from the bilinear interpolation that takes place as part of the fitting 
##' algorithm.
##'
##' @method selectObsWindow default
##' @param xyt an object
##' @param cellwidth size of the grid spacing in chosen units (equivalent to the cell width argument in \link{lgcpPredict})
##' @param ... additional arguments      
##' @return this is the same as selectObsWindow.stppp
##' @seealso \link{spatialAtRisk} \link{selectObsWindow.stppp}
##' @export

selectObsWindow.default <- function(xyt,cellwidth,...){
    return(selectObsWindow.stppp(xyt,cellwidth,...))
}



##' selectObsWindow.stppp function
##'
##' This function computes an appropriate observation window on which to perform prediction. Since the FFT grid
##' must have dimension 2^M by 2^N for some M and N, the window \code{xyt$window}, is extended to allow this to be fit in for a given cell width.
##'
##' !!NOTE!! that this function also returns the grid ($xvals and $yvals) on which the FFT (and hence MALA) will be performed. It is useful to
##' define spatialAtRiskobjects on this grid to prevent loss of information from the bilinear interpolation that takes place as part of the fitting 
##' algorithm.
##'
##' @method selectObsWindow stppp
##' @param xyt an object of class stppp
##' @param cellwidth size of the grid spacing in chosen units (equivalent to the cell width argument in \link{lgcpPredict})
##' @param ... additional arguments  
##' @return a resized stppp object together with grid sizes M and N ready for FFT, together with the FFT grid locations, can be useful for estimating lambda(s)
##' @seealso \link{spatialAtRisk}
##' @export

selectObsWindow.stppp <- function(xyt,cellwidth,...){
 
    ###
    # compute M and N
    ###

    M <- 2^ceiling(log(diff(xyt$window$xrange)/cellwidth,base=2)) # minimum gridsize in x-direction
    N <- 2^ceiling(log(diff(xyt$window$yrange)/cellwidth,base=2)) # minimum gridsize in x-direction  
    
    ###
    # Lastly, adjust size of observation window to suit
    ###
    
    xdim <- M * cellwidth
    ydim <- N * cellwidth
    xadd <- (xdim - diff(xyt$window$xrange))/2
    yadd <- (ydim - diff(xyt$window$yrange))/2    

    xyt$window$xrange <- xyt$window$xrange + c(-xadd,xadd)
    xyt$window$yrange <- xyt$window$yrange + c(-yadd,yadd) 

    xdiff <- diff(xyt$window$xrange)/M
    ydiff <- diff(xyt$window$yrange)/N
    xvals <- xyt$window$xrange[1] + xdiff/2 + xdiff*(0:(M-1))
    yvals <- xyt$window$yrange[1] + ydiff/2 + ydiff*(0:(N-1))    
    
    obj <- list(xyt=xyt,M=M,N=N,xvals=xvals,yvals=yvals) 
    
    return(obj)

}



##' getRotation function
##'
##' Generic function for the computation of rotation matrices.
##'
##' @param xyt an object
##' @param ... additional arguments
##' @return method getRotation
##' @seealso \link{getRotation.stppp}
##' @export

getRotation <- function(xyt,...){
    UseMethod("getRotation")
}



##' getRotation.default function
##'
##' Presently there is no default method, see ?getRotation.stppp
##'
##' @method getRotation default
##' @param xyt an object
##' @param ... additional arguments
##' @return currently no default implementation
##' @seealso \link{getRotation.stppp}
##' @export

getRotation.default <- function(xyt,...){
    stop("Method only implemented for objects of class stppp, see ?getRotation.stppp")
}



##' getRotation.stppp function
##'
##' Compute  rotation matrix if observation window is a polygonal boundary
##'
##' @method getRotation stppp
##' @param xyt an object of class stppp
##' @param ... additional arguments
##' @return the optimal rotation matrix and rotated data and observation window. Note it may or may not be advantageous to rotate the window, this information is displayed prior to the MALA routine when using lgcpPredict
##' @export


getRotation.stppp <- function(xyt,...){

    if (xyt$window$type=="polygonal"){
        csrotmat <- function(cs){
            return(matrix(c(cs[1],cs[2],-cs[2],cs[1]),2,2))
        }
        opfun <- function(cs){
            rwin <- affine(win,mat=csrotmat(cs),rescue=FALSE)
            return(diff(range(rwin$bdry[[1]]$x))*diff(range(rwin$bdry[[1]]$y)))
        }
        win <- xyt$window
        ch <- convexhull(win)
        chx <- ch$bdry[[1]]$x
        chy <- ch$bdry[[1]]$y
        chx <- c(chx,chx[1])
        chy <- c(chy,chy[1])
        cstheta <- matrix(NA,length(chx)-1,2)
        for (i in 1:(length(chx)-1)){
            v <- c(chx[i+1]-chx[i],chy[i+1]-chy[i])
            cstheta[i,] <- v / sqrt(v[1]^2+v[2]^2)
        }      
        areas <- apply(cstheta,1,opfun)
        rotation <- csrotmat(cstheta[which(areas==min(areas))[1],])
        return(list(xyt=affine(X=xyt,mat=rotation),rotation=rotation))
    }
    else{
        cat("Note: window rotation only implemented for polygonal windows\n")
    }          
}


##' roteffgain function
##'
##' Compute whether there might be any advantage in rotating the observation window in the object xyt for a proposed cell width.
##'
##' @param xyt an object of class stppp
##' @param cellwidth size of grid on which to do MALA
##' @return whether or not there woud be any efficiency gain in the MALA by rotating window
##' @seealso \link{getRotation.stppp}
##' @export

roteffgain <- function(xyt,cellwidth){
    ow <- selectObsWindow(xyt,cellwidth)
	xyt <- ow$xyt
	M <- ow$M
	N <- ow$N  
    
    rotinf <- getRotation(xyt)
    rot <- selectObsWindow(rotinf$xyt,cellwidth)
    rotws <- 4*rot$M*rot$N
    notrotws <- 4*M*N
    if (round(100*rotws/notrotws)<100){
        cat(paste("By rotating observation window, the efficiency gain would be: ",round(100*notrotws/rotws),"%, see ?getRotation.stppp\n",sep=""))
        cat("NOTE: efficiency gain is measured as the percentage increase in FFT grid cells from not rotating compared with rotating\n")
        return(TRUE)
    }
    else{
        return(FALSE)
    }
}


##' rotmat function
##'
##' This function returns a rotation matrix corresponding to an anticlockwise rotation of theta radians about the origin
##'
##' @param theta an angle in radians
##' @return the transformation matrix corresponding to an anticlockwise rotation of theta radians about the origin 
##' @export

rotmat <- function(theta){
    return(matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),2,2))
}



##' affine.stppp function
##'
##' An affine transformation of an object of class \code{stppp}
##'
##' @method affine stppp
##' @param X an object of class stppp
##' @param mat matrix of affine transformation
##' @param ... additional arguments
##' @return the object acted on by the transformation matrix 
##' @export


affine.stppp <- function(X,mat,...){
    return(stppp.ppp(P=affine.ppp(X,mat=mat,rescue=FALSE),t=X$t,tlim=X$tlim)) 
}



##' affine.fromFunction function
##'
##' An affine transformation of an object of class \code{fromFunction}
##'
##' @method affine fromFunction
##' @param X an object of class fromFunction
##' @param mat matrix of affine transformation
##' @param ... additional arguments
##' @return the object acted on by the transformation matrix 
##' @export


affine.fromFunction <- function(X,mat,...){
    matinv <- solve(mat)
    f2 <- function(xdash,ydash){
        xy <- as.vector(matinv%*%c(xdash,ydash))
        return(X$f(xy[1],xy[2]))
    }
    return(spatialAtRisk(X=f2))       
}



##' affine.fromXYZ function
##'
##' An affine transformation of an object of class \code{fromXYZ}. Nearest Neighbour interpolation
##'
##' @method affine fromXYZ
##' @param X an object of class fromFunction
##' @param mat matrix of affine transformation
##' @param ... additional arguments
##' @return the object acted on by the transformation matrix 
##' @export


affine.fromXYZ <- function(X,mat,...){
    xdiv <- X$X[2] - X$X[1]
    ydiv <- X$Y[2] - X$Y[1]
    pts <- matrix(c(rep(X$X,length(X$Y)),rep(X$Y,each=length(X$X))),length(X$X)*length(X$Y),2) # grid in original frame
    trpts <- t(mat%*%t(pts)) # transformed grid
    xrg <- range(trpts[,1])
    yrg <- range(trpts[,2])
    # now create new grid in transformed space
    gsx <- ceiling(diff(xrg)/xdiv)
    gsy <- ceiling(diff(yrg)/ydiv)
    newX <- seq(xrg[1],xrg[2],length.out=gsx)
    newY <- seq(yrg[1],yrg[2],length.out=gsy)
    newpts <- matrix(c(rep(newX,gsy),rep(newY,each=gsx)),gsx*gsy,2) # new grid in transformed frame
    invtrpts <- t(solve(mat)%*%t(newpts)) # inverse transformed grid in original frame
    # now do the interpolation in the original frame using interp.im
    newZm <- matrix(interp.im(im(t(X$Zm),xcol=X$X,yrow=X$Y),invtrpts[,1],invtrpts[,2]),gsx,gsy)
    return(spatialAtRisk(list(X=newX,Y=newY,Zm=newZm)))
}



##' affine.SpatialPolygonsDataFrame function
##'
##' An affine transformation of an object of class \code{SpatialPolygonsDataFrame}
##'
##' @method affine SpatialPolygonsDataFrame
##' @param X an object of class fromFunction
##' @param mat matrix of affine transformation
##' @param ... additional arguments
##' @return the object acted on by the transformation matrix 
##' @export

affine.SpatialPolygonsDataFrame <- function(X,mat,...){
    df <- as(X, "data.frame")
    polys <- slot(as(X,"SpatialPolygons"),"polygons")
    aff <- function(x,mat){
        return(Polygon(t(apply(slot(x, "coords"),1,function(co){mat%*%co}))))
    }
    getpolys <- function(x){
        pgs <- slot(x, "Polygons")
        return(Polygons(lapply(pgs,aff,mat=mat),ID = slot(x, "ID")))
    }
    newpolys <- lapply(polys,getpolys)
    rn <- row.names(df)
    sapply(1:length(newpolys),function(i){slot(newpolys[[i]], "ID") <<- rn[i]})
    return(SpatialPolygonsDataFrame(SpatialPolygons(newpolys), data = df))
}



##' affine.fromSPDF function
##'
##' An affine transformation of an object of class \code{fromSPDF}
##'
##' @method affine fromSPDF
##' @param X an object of class fromSPDF
##' @param mat matrix of affine transformation
##' @param ... additional arguments
##' @return the object acted on by the transformation matrix 
##' @export


affine.fromSPDF <- function(X,mat,...){
    obj <- list()
    obj$spdf <- affine(X$spdf,mat=mat,...)
    class(obj) <- c("fromSPDF","spatialAtRisk","SpatialPolygonsDataFrame")
    return(obj)
}


