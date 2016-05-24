make.SpatialVx <- function(X, Xhat, thresholds=NULL, loc=NULL, projection=FALSE, subset=NULL, time.vals = NULL,
			    reg.grid = TRUE, map = FALSE, loc.byrow = FALSE, field.type="", units="", data.name=c("X","Xhat"),
			    q=c(0, 0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9, 0.95), qs=NULL) {

    xdim <- dim(X)

    if(is.list(Xhat)) {

	nforecast <- length(Xhat)
	ydim <- dim(Xhat[[1]])
	ydims <- lapply(Xhat, dim)
	f <- function(x, d) all(x == d)
	if(!all(unlist(lapply(ydims, f, d=ydim)))) stop("make.SpatialVx: invalid Xhat argument.")

    } else {

	nforecast <- 1
	ydim <- dim( Xhat )

    }

    if(all(xdim != ydim)) {

	stopmsg <- paste("make.SpatialVx: dim of X (", xdim[ 1 ], " by ", xdim[ 2 ],
		") must be the same as dim of (each component of) Xhat (", ydim[ 1 ], " by ",
		ydim[ 2 ], ")", sep = "" )

	stop( stopmsg )
	# stop("make.SpatialVx: dim of X must be the same as dim of (each component of) Xhat")

    } # end of stop bc dims are not equal stmt.

    out <- list(X, Xhat)

    if(is.null(thresholds)) {

	thresholds <- cbind(quantile(c(X), probs=q, na.rm=TRUE), quantile(c(unlist(Xhat)), probs=q, na.rm=TRUE))
        qs <- as.character(q)

    } else if(!is.matrix(thresholds)) {

	qs <- as.character(thresholds)
	thresholds <- cbind(thresholds, thresholds)

    } else qs <- paste("threshold", 1:dim(thresholds)[1])

    udim <- dim(thresholds)
    if( udim[2] == 2) colnames(thresholds) <- c("X", "Xhat")
    else if( udim[2] > 2) colnames(thresholds) <- c("X", paste("Xhat", 1:(udim[2] - 1), sep=""))

    rownames(thresholds) <- qs

    if(is.null(loc)) loc <- cbind(rep(1:xdim[1],xdim[2]), rep(1:xdim[2], each=xdim[1]))

    if(is.null(time.vals)) {

	if(length(xdim) == 3) time.vals <- 1:xdim[3]
	else time.vals <- 1

    }

    if(length(data.name) == nforecast + 2) msg <- paste("", "\n", data.name[1], sep="")
    else msg <- ""

    if(field.type != "" && units != "") msg <- paste(msg, "\n", field.type, " (", units, ")", sep="")
    else if(field.type != "") msg <- paste(msg, "\n", field.type, sep="")
    else if(units != "") msg <- paste(msg, "\n", " (", units, ")", sep="")


    class(out) <- "SpatialVx"
    attr(out, "xdim") <- xdim
    attr(out, "time") <- time.vals
    attr(out, "thresholds") <- thresholds
    attr(out, "udim") <- udim
    attr(out, "loc") <- loc
    attr(out, "loc.byrow") <- loc.byrow
    attr(out, "subset") <- subset
    attr(out, "data.name") <- data.name
    attr(out, "nforecast") <- nforecast
    attr(out, "field.type") <- field.type
    attr(out, "units") <- units
    attr(out, "projection") <- projection
    attr(out, "reg.grid") <- reg.grid
    attr(out, "map") <- map
    attr(out, "qs") <- qs
    attr(out, "msg") <- msg

    return(out)
} # end of 'make.SpatialVx' function.

print.SpatialVx <- function(x, ...) {
    tmp <- attributes(x)
    cat("SpatialVx data object: ", tmp$data.name, "\n")
    cat(tmp$msg, "\n")
    # if(tmp$field.type != "") cat("Field type: ", tmp$field.type)
    # if(tmp$units != "") cat(" (", tmp$units, ")", "\n")

    cat("Field dimensions: ", tmp$xdim, "\n")

    if(tmp$nforecast > 1) cat(tmp$nforecast, " forecasts to be evaluated/verified.\n")
    else cat("1 forecast to be evaluated/verified.\n")

    u <- tmp$thresholds
    f <- function(x) {
	    if(all(x == x[1])) return(TRUE)
	    else return(FALSE)
	} # end of internal 'f' function.
    if(all(apply(u,1,f))) {
	cat("thresholds of interest (same for all fields):\n")
	print(tmp$qs)
    } else {
	cat("Different thresholds between verification and forecast(s)\n")
	print(u)
    }

    if(!is.null(tmp$levels)) {
        cat("\n", "Neighborhood levels to be applied for neighborhood methods:\n")
        print(tmp$levels)
        cat("\n", "Smoothing function to be applied for neighborhood methods:\n")
        print(tmp$smooth.fun)
        if(!is.null(tmp$smooth.params)) {
            cat("\n", "Smoothing function parameters:\n")
            print(tmp$smooth.params)
        }
    }

    invisible()
} # end of 'print.SpatialVx' function.

summary.SpatialVx <- function(object, ...) {

    out <- list()

    x <- object
    print(x)
    tmp <- attributes(x)

    print( out$X <- summary(x[[1]]))

    if(tmp$nforecast > 1) print(out$Xhat <- lapply(x[[2]], summary))
    else print(out$Xhat <- summary(x[[2]]))

    invisible(out)
} # end of 'summary.SpatialVx' function.

plot.SpatialVx <- function(x, ..., set.pw=FALSE, time.point=1, model=1, col, zlim, horizontal=TRUE) {

    tmp <- attributes(x)
    loc.byrow <- tmp$loc.byrow
    nf <- tmp$nforecast
    Nt <- length(tmp$time)

    xd <- tmp$xdim
    loc <- tmp$loc

    if(tmp$map) ax <- list(x=pretty(round(loc[,1], digits=2)), y=pretty(round(loc[,2], digits=2)))

    u <- tmp$thresholds

    if(!is.logical(set.pw)) {
	if(is.vector(set.pw)) par(mfrow=set.pw, oma=c(0,0,2,0))
    } else if(set.pw) {
	if(is.na(time.point) && is.na(model)) {
	    Npl <- nf + 2 * Nt - 1
	    par(mfrow=c(Nt, nf + 1), oma=c(0,0,2,0))
	} else if(is.function(time.point) || is.character(time.point)) {
	    if(length(model) == 1) par(mfrow=c(1,2), oma=c(0,0,2,0))
	    else par(mfrow=c(2,length(model)), oma=c(0,0,2,0))
	} else par(mfrow=c(1,2), oma=c(0,0,2,0))
    }

    if(length(tmp$data.name) == nf + 1) ptitle <- tmp$data.name
    else if(length(tmp$data.name) == nf + 2) ptitle <- tmp$data.name[-1]
    else ptitle <- NULL

    msg <- tmp$msg

    if(missing(zlim)) zlim <- range(c(unlist(x)), finite=TRUE)
    if(missing(col)) {
	if(min(c(unlist(x)), na.rm=TRUE) == 0) col <- c("gray", tim.colors(64))
	else col <- tim.colors(64)
    }

    Vx <- x[[1]]
    Fcst <- x[[2]]

    if(!is.na(model)) {
	dn <- tmp$data.name
        if(length(dn) == nf + 2) {
	    obs.name <- dn[2]
	    mod.names <- dn[-(1:2)]
        } else {
	    obs.name <- dn[1]
	    mod.names <- dn[-1]
	}
        if(!is.numeric(model)) model <- (1:nf)[mod.names == model]
	if(nf > 1) Fcst <- Fcst[[model]]
	ptitle <- c(obs.name, mod.names[model])
    }

    if(Nt > 1 && !is.na(time.point)) {

        if(is.function(time.point) || is.character(time.point)) {

	    afun <- match.fun(time.point)
	    X <- apply(Vx, c(1,2), afun, ...)
	    Xhat <- apply(Fcst, c(1,2), afun, ...)
	    if(is.character(time.point)) fname <- time.point
	    else fname <- deparse(substitute(time.point))
	    ptitle <- paste(ptitle, "\n(aggregated over time by ", fname, ")", sep="")

        } else {

	    X <- Vx[,,time.point]
	    ptitle <- paste(ptitle, "\n(t = ", tmp$time[time.point], ")", sep="")

	}

    } # end of if more than one time point, and !is.na(time.point) stmts.

    if(is.na(time.point) && is.na(model)) {

	tiden <- tmp$time
        tdim <- dim(tiden)
        if(!is.null(tdim)) {
            if(length(tdim) == 2) {
                if(nf==1) tiden <- c(tiden)
                else tiden <- c(tiden[,1], rep(tiden[,2], nf))
            }
        } else tiden <- rep(tiden, nf + 1)
        if(is.numeric(tiden)) tiden <- paste("t = ", tiden, sep="")

        if(is.null(ptitle)) ptitle <- tiden
        else ptitle <- c(paste(ptitle[1], tiden[1:Nt], sep="\n"), paste(rep(ptitle[-1], each=Nt), tiden[-(1:Nt)], sep="\n"))

        for(tiid in 1:Nt) {
    	    if(Nt == 1) X <- Vx
    	    else X <- Vx[,,tiid]
    
            if(tmp$map) {

        	r <- apply(tmp$loc, 2, range, finite=TRUE)
        	map(xlim=r[,1], ylim=r[,2], type="n")
		axis(1, at=ax$x, labels=ax$x)
		axis(2, at=ax$y, labels=ax$y)

        	if(tmp$projection && tmp$reg.grid) {

        	        poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow),
			    matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
			    X, add=TRUE, col=col, zlim=zlim)

        	} else image(as.image(X, x=loc, nx=xd[1], ny=xd[2]), add=TRUE, col=col, zlim=zlim)
        	
        	map(add=TRUE, lwd=2)
        	map(database="state", add=TRUE, lwd=1.5)
    	        # contour(X, levels=u[,"X"], col="white", add=TRUE)
    
        	if(!is.null(ptitle)) title(ptitle[tiid])
      
        	for(i in 1:nf) {
        	    if(nf > 1) Xhat <- Fcst[[i]]
    	            else Xhat <- Fcst
    
    	            if(Nt > 1) Xhat <- Xhat[,,tiid]
    
        	    map(xlim=r[,1], ylim=r[,2], lwd=2, type="n")
		    axis(1, at=ax$x, labels=ax$x)
		    axis(2, at=ax$y, labels=ax$y)

        	    if(tmp$projection && tmp$reg.grid) {

        	       poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow),
				matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
				Xhat, add=TRUE, col=col, zlim=zlim)

        	    } else image(as.image(Xhat, x=loc, nx=xd[1], ny=xd[2]), add=TRUE, col=col, zlim=zlim)

    		    map(add=TRUE, lwd=2)
        	    map(database="state", add=TRUE, lwd=1.5)

    		    # if(dim(u)[2] > 2) contour(Xhat, levels=u[,i+1], col="white", add=TRUE)
    		    # else contour(Xhat, levels=u[,"Xhat"], col="white", add=TRUE)
    
        	    if(!is.null(ptitle)) title(ptitle[Nt * i + tiid])

        	} # end of for 'i' loop.

            } else {

        	if(!is.null(ptitle)) image(X, axes=FALSE, col=col, zlim=zlim, main=ptitle[1])
                else image(X, axes=FALSE, col=col, zlim=zlim)
    	        # contour(X, levels=u[,"X"], col="white", add=TRUE)
    
    	        for(i in 1:nf) {
                    if(nf > 1) Xhat <- Fcst[[i]]
                    else Xhat <- Fcst
    
                    if(Nt > 1) Xhat <- Xhat[,,tiid]
    
                    if(tmp$projection && tmp$reg.grid) {

                        poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow),
				matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
				Xhat, col=col, zlim=zlim)

                    } else image(Xhat, axes=FALSE, col=col, zlim=zlim)
    
                    if(!is.null(ptitle)) title(ptitle[Nt * i + tiid])

                } # end of for 'i' loop.
    
            } # end of if else 'map' stmts.

        } # end of for 'tiid' loop.  

        image.plot(x[[1]], legend.only=TRUE, col=col, zlim=zlim, horizontal=horizontal)

    } else if(is.na(time.point)) {
	
	for(tiid in 1:Nt) {
            if(Nt == 1) X <- Vx
            else X <- Vx[,,tiid]

            if(tmp$map) {
                r <- apply(tmp$loc, 2, range, finite=TRUE)
                map(xlim=r[,1], ylim=r[,2], type="n")
	 	axis(1, at=ax$x, labels=ax$x)
                axis(2, at=ax$y, labels=ax$y)

                if(tmp$projection && tmp$reg.grid) {

                        poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow),
                            matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
                            X, add=TRUE, col=col, zlim=zlim)

                } else image(as.image(X, x=loc, nx=xd[1], ny=xd[2]), add=TRUE, col=col, zlim=zlim)

                map(add=TRUE, lwd=2)
                map(database="state", add=TRUE, lwd=1.5)
                # contour(X, levels=u[,"X"], col="white", add=TRUE)

                if(!is.null(ptitle)) title(ptitle[tiid])

                if(nf > 1) Xhat <- Fcst[[model]]
                else Xhat <- Fcst

                if(Nt > 1) Xhat <- Xhat[,,tiid]

                map(xlim=r[,1], ylim=r[,2], type="n")
	 	axis(1, at=ax$x, labels=ax$x)
                axis(2, at=ax$y, labels=ax$y)
                if(tmp$projection & tmp$reg.grid) {
                       poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow),
                                matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
                                Xhat, add=TRUE, col=col, zlim=zlim)
                } else image(as.image(Xhat, x=loc, nx=xd[1], ny=xd[2]), add=TRUE, col=col, zlim=zlim)
                map(add=TRUE, lwd=2)
                map(database="state", add=TRUE, lwd=1.5)
                # if(dim(u)[2] > 2) contour(Xhat, levels=u[,model+1], col="white", add=TRUE)
                # else contour(Xhat, levels=u[,"Xhat"], col="white", add=TRUE)

                if(!is.null(ptitle)) title(ptitle[Nt + tiid])

            } else {
                if(!is.null(ptitle)) image(X, axes=FALSE, col=col, zlim=zlim, main=ptitle[1])
                else image(X, axes=FALSE, col=col, zlim=zlim)
                # contour(X, levels=u[,"X"], col="white", add=TRUE)

                if(nf > 1) Xhat <- Fcst[[model]]
                else Xhat <- Fcst

                if(Nt > 1) Xhat <- Xhat[,,tiid]

                if(tmp$projection && tmp$reg.grid) {

                    poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow),
                                matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
                                Xhat, col=col, zlim=zlim)

                } else image(Xhat, axes=FALSE, col=col, zlim=zlim)
                # if(dim(u)[2] > 2) contour(Xhat, levels=u[,model+1], col="white", add=TRUE)
                # else contour(Xhat, levels=u[,"Xhat"], col="white", add=TRUE)

                if(!is.null(ptitle)) title(ptitle[Nt + tiid])
   
            } # end of if else 'map' stmts.
        } # end of for 'tiid' loop.

	image.plot(X, legend.only=TRUE, col=col, zlim=zlim, horizontal=horizontal)

    } else if(is.na(model)) {
  
        ptitle <- tmp$data.name
	if(length(ptitle) == nf + 2) ptitle <- ptitle[-1]
	ptitle <- paste(ptitle, "\n(t = ", tmp$time[time.point], ")", sep="")
	
        if(tmp$map) {
             r <- apply(tmp$loc, 2, range, finite=TRUE)
             map(xlim=r[,1], ylim=r[,2], type="n")
	     axis(1, at=ax$x, labels=ax$x)
             axis(2, at=ax$y, labels=ax$y)

             if(tmp$projection && tmp$reg.grid) {

                 poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow),
                            matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
                            X, add=TRUE, col=col, zlim=zlim)

             } else image(as.image(X, x=loc, nx=xd[1], ny=xd[2]), add=TRUE, col=col, zlim=zlim)

             map(add=TRUE, lwd=2)
             map(database="state", add=TRUE, lwd=1.5)
             # contour(X, levels=u[,"X"], col="white", add=TRUE)
   
             if(!is.null(ptitle)) title(ptitle[1])

             for(i in 1:nf) {
                 if(nf > 1) Xhat <- Fcst[[i]]
                 else Xhat <- Fcst

                 if(Nt > 1) Xhat <- Xhat[,,time.point]

                 map(xlim=r[,1], ylim=r[,2], type="n")
	 	 axis(1, at=ax$x, labels=ax$x)
                 axis(2, at=ax$y, labels=ax$y)

                 if(tmp$projection && tmp$reg.grid) {

                       poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow),
                                matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
                                Xhat, add=TRUE, col=col, zlim=zlim)

                 } else image(as.image(Xhat, x=loc, nx=xd[1], ny=xd[2]), add=TRUE, col=col, zlim=zlim)

                    map(add=TRUE, lwd=2)
                    map(database="state", add=TRUE, lwd=1.5)
                    # if(dim(u)[2] > 2) contour(Xhat, levels=u[,i+1], col="white", add=TRUE)
                    # else contour(Xhat, levels=u[,"Xhat"], col="white", add=TRUE)

                    if(!is.null(ptitle)) title(ptitle[i + 1])
                } # end of for 'i' loop.

            } else {

                if(!is.null(ptitle)) image(X, axes=FALSE, col=col, zlim=zlim, main=ptitle[1])
                else image(X, axes=FALSE, col=col, zlim=zlim)
                # contour(X, levels=u[,"X"], col="white", add=TRUE)

                for(i in 1:nf) {

                    if(nf > 1) Xhat <- Fcst[[i]]
                    else Xhat <- Fcst

                    if(Nt > 1) Xhat <- Xhat[,,time.point]

                    if(tmp$projection && tmp$reg.grid) {

                        poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow),
                                matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
                                Xhat, col=col, zlim=zlim)

                    } else image(Xhat, axes=FALSE, col=col, zlim=zlim)
                    # if(dim(u)[2] > 2) contour(Xhat, levels=u[,i+1], col="white", add=TRUE)
                    # else contour(Xhat, levels=u[,"Xhat"], col="white", add=TRUE)

                    if(!is.null(ptitle)) title(ptitle[i + 1])
                } # end of for 'i' loop.

            } # end of if else 'map' stmts.

        image.plot(X, legend.only=TRUE, col=col, zlim=zlim, horizontal=horizontal)

    } else {

	X <- Vx
	# if(nf > 1) Xhat <- Fcst[[model]]
	# else
	Xhat <- Fcst
	if(Nt > 1) Xhat <- Xhat[,,time.point]

	if(tmp$map) {

            r <- apply(tmp$loc, 2, range, finite=TRUE)
            map(xlim=r[,1], ylim=r[,2], type="n")
	    axis(1, at=ax$x, labels=ax$x)
            axis(2, at=ax$y, labels=ax$y)

            if(tmp$projection && tmp$reg.grid) {

                poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow), matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
                    X, add=TRUE, col=col, zlim=zlim)

            } else image(as.image(X, x=loc, nx=xd[1], ny=xd[2]), add=TRUE, col=col, zlim=zlim)

            map(add=TRUE, lwd=2)
            map(database="state", add=TRUE, lwd=1.5)
            # contour(X, levels=u[,"X"], col="white", add=TRUE)

            if(!is.null(ptitle)) title(ptitle[1])
 
            map(xlim=r[,1], ylim=r[,2], type="n")
	    axis(1, at=ax$x, labels=ax$x)
            axis(2, at=ax$y, labels=ax$y)

            if(tmp$projection && tmp$reg.grid) {

                poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow), matrix(loc[,2], xd[1], xd[2],
                        byrow=loc.byrow), Xhat, add=TRUE, col=col, zlim=zlim)

            } else image(as.image(Xhat, x=loc, nx=xd[1], ny=xd[2]), add=TRUE, col=col, zlim=zlim)

            map(add=TRUE, lwd=2)
            map(database="state", add=TRUE, lwd=1.5)
            # if(dim(u)[2] > 2) contour(Xhat, levels=u[,model+1], col="white", add=TRUE)
            # else contour(Xhat, levels=u[,"Xhat"], col="white", add=TRUE)

            if(!is.null(ptitle)) title(ptitle[2])

        } else {

            if(!is.null(ptitle)) image(X, axes=FALSE, col=col, zlim=zlim, main=ptitle[1])
            else image(X, axes=FALSE, col=col, zlim=zlim)
            # contour(X, levels=u[,"X"], col="white", add=TRUE)

            if(tmp$projection && tmp$reg.grid) {

                    poly.image(matrix(loc[,1], xd[1], xd[2], byrow=loc.byrow),
			matrix(loc[,2], xd[1], xd[2], byrow=loc.byrow),
			Xhat, col=col, zlim=zlim)

            } else image(Xhat, axes=FALSE, col=col, zlim=zlim)
            # if(dim(u)[2] > 2) contour(Xhat, levels=u[,model+1], col="white", add=TRUE)
            # else contour(Xhat, levels=u[,"Xhat"], col="white", add=TRUE)

            if(!is.null(ptitle)) title(ptitle[2])

        } # end of if else 'map' stmts.
	image.plot(X, legend.only=TRUE, col=col, zlim=zlim, horizontal=horizontal)
    } # end of which plots to make stmts.
    if(set.pw) {
	title("")
	mtext(msg, line=0.05, outer=TRUE)
    }
    invisible()
} # end of 'plot.SpatialVx' function.

# datagrabber <- function(x, ...) {
#     UseMethod("datagrabber", x)
# } # end of 'datagrabber' function.

datagrabber.SpatialVx <- function(x, ..., time.point=1, model=1) {

    tmp <- attributes(x)

    if(!missing(time.point)) {
        if(length(time.point) > 1) stop("datagrabber: length of time.point must be one.")
        tiid <- tmp$time
        Nt <- length(tiid)
        if(!is.numeric(time.point)) time.point <- (1:Nt)[tiid == time.point]
        if(is.na(time.point)) stop("datagrabber: invalid time.point argument.")
   }

   if(!missing(model)) {
        if(length(model) > 1) stop("datagrabber: length of model argument must be one.")
        if(!is.numeric(model)) {
            nf <- tmp$nforecast
            dn <- tmp$data.name
            if(length(dn) == nf + 2) mod.names <- dn[-(1:2)]
            else mod.names <- dn[-1]
            model <- (1:nf)[dn == model]
            if(is.na(model)) stop("datagrabber: invalid model argument.")
        }
   }    

    xdim <- tmp$xdim
    nf <- tmp$nforecast

    Vx <- x[[1]]
    Fcst <- x[[2]]
    if(nf > 1) Fcst <- Fcst[[model]]

    if(length(xdim) == 3) {
	Vx <- Vx[,,time.point]
	Fcst <- Fcst[,,time.point]
    }
    out <- list(X=Vx, Xhat=Fcst)
    return(out) 
} # end of 'datagrabber.SpatialVx' function.

datagrabber.features <- function(x, ...) {

    return(list(X = x$X, Xhat = x$Xhat))

} # end of 'data.grabber.features' function.

datagrabber.matched <- function(x, ...) {

    return(list(X = x$X, Xhat = x$Xhat))

} # end of 'data.grabber.matched' function.

hist.SpatialVx <- function(x, ..., time.point=1, model=1) {
   tmp <- attributes(x)

   if(!missing(time.point) && !missing(model)) dat <- datagrabber(x, time.point=time.point, model=model)
   else if(!missing(time.point)) dat <- datagrabber(x, time.point=time.point)
   else if(!missing(model)) dat <- datagrabber(x, model=model)
   else dat <- datagrabber(x)

   X <- dat$X
   Y <- dat$Xhat

   u <- tmp$thresholds
   udim <- dim(u)
   udim[1] <- udim[1]+1

   dn <- tmp$data.name
   nf <- tmp$nforecast
   if(length(dn) == nf + 2) {
	vxname <- dn[2]
	fcname <- dn[-(1:2)]
   } else {
	vxname <- dn[1]
	fcname <- dn[-1]
   }

   if(!missing(model)) {
	if(length(model) > 1) stop("hist.SpatialVx: length of model argument must be one.")
        if(!is.numeric(model)) {
            model <- (1:nf)[dn == model]
            if(is.na(model)) stop("hist.SpatialVx: invalid model argument.")
        }
   }
   fcname <- fcname[model]

   par(mfrow=udim)
   for(i in 1:udim[1]) {
        if(i==1) {
           m1 <- paste(vxname, "\n", "All ", tmp$field.type, " values", sep="")
           m2 <- paste(fcname, "\n", "All ", tmp$field.type, " values", sep="")
           tmpX <- c(X)
           tmpY <- c(Y)
        } else {
           m1 <- paste(vxname, "\n", "Only values >= ", u[i-1,2], " ", tmp$units, sep="")
           m2 <- paste(fcname, "\n", "Only values >= ", u[i-1,1], " ", tmp$units, sep="")
           tmpX <- c(X[X>=u[i-1,2]])
           tmpY <- c(Y[Y>=u[i-1,1]])
        }
        if(i==udim[1]) x1 <- paste(tmp$field.type, " (", tmp$units, ")", sep="")
        else x1 <- ""
        hist(tmpX, main=m1, xlab=x1, ...)
        hist(tmpY, main=m2, xlab=x1, ...)
   } # end of for 'i' loop.
   invisible()
} # end of 'hist.SpatialVx' function.
