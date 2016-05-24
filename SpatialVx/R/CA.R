clusterer <- function(X, Y=NULL, ...) {
    UseMethod("clusterer", X)
} # end of 'clusterer' function.

clusterer.SpatialVx <- function(X, Y=NULL, ..., time.point=1, model=1, xyp=TRUE, threshold=1e-8,
    linkage.method="complete", stand=TRUE, trans="identity", verbose=FALSE) {

    if(!missing(Y) && !is.null(Y)) warning("clusterer.SpatialVx: Y argument is ignored.")

    a <- attributes(X)

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(X, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(X, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(X, model=model)
    else dat <- datagrabber(X)

    X <- dat$X
    Xhat <- dat$Xhat
    ## End: Get the data sets

    out <- clusterer.default(X=X, Y=Xhat, ..., xloc=a$loc, xyp=xyp, threshold=threshold,
				linkage.method=linkage.method, stand=stand, trans=trans, a=a, verbose=verbose)

    attr(out, "time.point") <- time.point
    attr(out, "model") <- model

    if(length(a$data.name) == a$nforecast + 2) {
        dn <- a$data.name[-(1:2)]
        vxname <- a$data.name[1:2]
    } else {
        dn <- a$data.name[-1]
        vxname <- a$data.name[1]
    }
    if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model

    attr(out, "data.name") <- c(vxname, dn[model.num])

    return(out)
} # end of 'clusterer.SpatialVx' function.

clusterer.default <- function(X, Y=NULL, ..., xloc=NULL, xyp=TRUE, threshold=1e-8, linkage.method="complete",
    stand=TRUE, trans="identity", a=NULL, verbose=FALSE) {

    out <- list()
    if(!is.null(a)) {

	if(!is.null(a$names)) a$names <- NULL
	attributes(out) <- a
	xdim <- a$xdim

    } else {

	data.name <- c(as.character(substitute(X)), as.character(substitute(Y)))
        names(data.name) <- c("verification","forecast")
        attr(out, "data.name") <- data.name
	xdim <- dim(X)
	attr(out, "xdim") <- xdim

    }

    out$data <- list(X=X, Xhat=Y)

    out$linkage.method <- linkage.method
    out$trans <- trans
    N <- prod(xdim)
    out$N <- N

    if(is.null(xloc)) {
	xloc <- cbind(rep(1:xdim[1],xdim[2]),rep(1:xdim[2],each=xdim[1]))
	goodloc <- !logical(N)
    } else goodloc <- !is.na(xloc[,1]) & !is.na(xloc[,2])
    if(is.null(a$loc)) attr(out, "loc") <- xloc

    if(length(threshold)==1) u <- c(threshold, threshold)
    else u <- threshold 
    attr(out, "threshold") <- u

    idX <- !is.na(c(X)) & c(X) >= u[1] & goodloc
    idY <- !is.na(c(Y)) & c(Y) >= u[2] & goodloc

    if(xyp) {
	tmpX <- cbind(xloc[idX,], c(X)[idX])
   	tmpY <- cbind(xloc[idY,], c(Y)[idY])
	tmpX[,3] <- c(do.call(trans, list(x=tmpX[,3])))
	tmpY[,3] <- c(do.call(trans, list(x=tmpY[,3])))
    } else {
	tmpX <- xloc[idX,]
	tmpY <- xloc[idY,]
    }
    NCx  <- dim(tmpX)[1]
    NCy  <- dim(tmpY)[1]
    out$NCo <- NCx:1
    out$NCf <- NCy:1
    if(verbose) {
	cat("\n", NCx, " initial clusters in verification field.\n")
	cat("\n", NCy, " initial clusters in forecast field.\n")
	cat("\n", "Performing cluster analysis on verification field.\n")
    }
    if(NCx < 2 | NCy < 2) stop("clusterer: Less than two initial clusters in one or both fields.\n")
    if(stand) {
        tmpX <- (tmpX - colMeans(tmpX, na.rm=TRUE))/sqrt(apply(tmpX,2,var,na.rm=TRUE))
        tmpY <- (tmpY - colMeans(tmpY, na.rm=TRUE))/sqrt(apply(tmpY,2,var,na.rm=TRUE))
    }

    args <- list(...)
    if(is.null(args$method)) metric <- "euclidean"
    else metric <- args$method

    if(is.element(linkage.method,c("ward","centroid","median")) & metric != "euclidean") stop(paste("clusterer: method must be eulcidean with ", linkage.method, " linkage.", sep=""))
 
    dX <- dist(tmpX, ...)
    XclustObj <- hclust(dX, method=linkage.method)
    XclustObj2 <- MakeClusterList(XclustObj)
    if(verbose) cat("\n", "Performing cluster analysis on forecast field.\n")
    dY <- dist(tmpY, ...)
    YclustObj <- hclust(dY, method=linkage.method)
    YclustObj2 <- MakeClusterList(YclustObj)
    out$cluster.identifiers <- list(X=XclustObj2, Y=YclustObj2)

    out$idX <- idX
    out$idY <- idY
    out$cluster.objects <- list()
    out$cluster.objects$X <- XclustObj
    out$cluster.objects$Y <- YclustObj
 
    zDiss <- matrix(NA, NCx, NCy)
    if(is.element(linkage.method,c("ward","median","centroid"))) {
	XclustObjent <- list()
	XclustObjent[[1]] <- NULL
	YclustObjent <- XclustObjent
    }
    if(verbose) cat("\n", "Finding distances between verification and forecast clusters and sorting them out.\n")
    Indy <- 1:NCx
    for(i in 1:NCx) {
	if(verbose) cat(i, " ")
	if(metric == "euclidean") {
	   if(is.matrix(tmpY)) zDiss[i,] <- sqrt(rowSums((matrix(tmpX[i,], nrow=nrow(tmpY), ncol=ncol(tmpX), byrow=TRUE) - tmpY)^2, na.rm=TRUE))
	   else zDiss[i,] <- sqrt(rowSums((matrix(tmpX[i,], nrow=nrow(tmpY), ncol=ncol(tmpX), byrow=TRUE) - matrix(tmpY, ncol(tmpY)))^2, na.rm=TRUE))
	} else if(metric == "manhattan") {
	   if(is.matrix(tmpY)) zDiss[i,] <- rowSums(abs(matrix(tmpX[i,], nrow=nrow(tmpY), ncol=ncol(tmpX), byrow=TRUE) - tmpY), na.rm=TRUE)
	   else zDiss[i,] <- rowSums(abs(matrix(tmpX[i,], nrow=nrow(tmpY), ncol=ncol(tmpX), byrow=TRUE) - matrix(tmpY, ncol(tmpY))), na.rm=TRUE)
	} else stop("clusterer: metric can only be euclidean or manhattan.")
	if((i>1) & is.element(linkage.method,c("ward","median","centroid"))) {
	   mm <- matrix(XclustObj$merge[1:(i-1),], ncol=2)
	   if(linkage.method=="median") {
		if(mm[i-1,1] < 0 & mm[i-1,2] < 0) XclustObjent[[i]] <- matrix(colMeans(tmpX[XclustObj2[[i]][[1]],],na.rm=TRUE), nrow=1)
		else if(mm[i-1,1] < 0 & mm[i-1,2] > 0) XclustObjent[[i]] <- (tmpX[abs(mm[i-1,1]),] + XclustObjent[[mm[i-1,2]+1]])/2
		else if(mm[i-1,1] > 0 & mm[i-1,2] < 0) XclustObjent[[i]] <- (tmpX[abs(mm[i-1,2]),] + XclustObjent[[mm[i-1,1]+1]])/2
		else XclustObjent[[i]] <- (XclustObjent[[mm[i-1,1]+1]] + XclustObjent[[mm[i-1,2]+1]])/2
	   } else XclustObjent[[i]] <- matrix(colMeans(tmpX[XclustObj2[[i]][[1]],],na.rm=TRUE), nrow=1)
	} # end of if 'i>1' stmt.
    } # end of for 'i' loop.
    if(verbose) cat("\n", "Inter-cluster distances found.  Verification clusters sorted out.  Sorting out forecast clusters.\n")
    Indy <- 1:NCy
    for(i in 1:NCy) {
	if(verbose) cat(i, " ")
	if(i == 1) YclustObj2[[1]] <- as.list(1:NCy)
        else if(i < NCy)  {
	   mm <- matrix(YclustObj$merge[1:(i-1),], ncol=2)
	   if((i>1) & is.element(linkage.method,c("ward","median","centroid"))) {
              if(linkage.method=="median") {
                if(mm[i-1,1] < 0 & mm[i-1,2] < 0) YclustObjent[[i]] <- matrix(colMeans(tmpY[YclustObj2[[i]][[1]],],na.rm=TRUE), nrow=1)
                else if(mm[i-1,1] < 0 & mm[i-1,2] > 0) YclustObjent[[i]] <- (tmpY[abs(mm[i-1,1]),] + YclustObjent[[mm[i-1,2]+1]])/2
                else if(mm[i-1,1] > 0 & mm[i-1,2] < 0) YclustObjent[[i]] <- (tmpY[abs(mm[i-1,2]),] + YclustObjent[[mm[i-1,1]+1]])/2
                else YclustObjent[[i]] <- (YclustObjent[[mm[i-1,1]+1]] + YclustObjent[[mm[i-1,2]+1]])/2
              } else YclustObjent[[i]] <- matrix(colMeans(tmpY[YclustObj2[[i]][[1]],],na.rm=TRUE), nrow=1)
           } # end of if 'i>1' stmt.
	} # end of if 'i < NCy' stmts.
    } # end of for 'i' loop.

    if(verbose) cat("\n", "Finding inter-cluster distances (NCf X NCo double loop within a double loop).\n")
    ICdists <- list()
    ICdists[[1]] <- list()
    ICdists[[1]][[1]] <- zDiss
    minICdists <- numeric(0)
    if(linkage.method=="mcquitty") {
	warning("clusterer: McQuitty linkage method not available for comparing clusters across fields.  Using average instead.")
	out$linkage.method <- c(linkage.method, "average")
	linkage.method <- "average"
    }
    for(i in 1:NCx) {
	if(verbose) cat("\n", i, ":\n")
	lookX <- XclustObj2[[i]]
	nx <- length(lookX)
	for(j in 1:NCy) {
	   if(i==1 & j==1) next
	   if(verbose) cat(j, " ")
	      lookY <- YclustObj2[[j]]
	      ny <- length(lookY)
	      look <- matrix(NA, nx, ny)
	      for(kx in 1:nx) {
		nx2 <- length(lookX[[kx]])
		if(is.element(linkage.method,c("ward","median","centroid"))) xcen <- matrix(colMeans(matrix(tmpX[lookX[[kx]],],ncol=ncol(tmpX),byrow=TRUE),na.rm=TRUE), ncol=ncol(tmpX), byrow=TRUE)
		if(linkage.method=="ward") Fx <- sum(rowSums((matrix(tmpX[lookX[[kx]],],ncol=ncol(tmpX),byrow=TRUE) - matrix(rep(c(xcen),nx2),ncol=ncol(tmpX),byrow=TRUE))^2,na.rm=TRUE),na.rm=TRUE)
		for(ky in 1:ny) {
		ny2 <- length(lookY[[ky]])
		if(is.element(linkage.method,c("ward","median","centroid"))) ycen <- matrix(colMeans(matrix(tmpY[lookY[[ky]],],ncol=ncol(tmpY),byrow=TRUE),na.rm=TRUE), ncol=ncol(tmpY), byrow=TRUE)
		if(linkage.method == "ward") {
		   Fy <- sum(rowSums((matrix(tmpY[lookY[[ky]],],ncol=ncol(tmpY),byrow=TRUE) - matrix(rep(c(ycen),ny2),ncol=ncol(tmpY),byrow=TRUE))^2,na.rm=TRUE),na.rm=TRUE)
		   xycen <- matrix(colMeans(rbind(tmpX[lookX[[kx]],],tmpY[lookY[[ky]],]),na.rm=TRUE), ncol=ncol(tmpX), byrow=TRUE)
		   Fxy <- sum(rowSums((rbind(tmpX[lookX[[kx]],],tmpY[lookY[[ky]],]) - matrix(rep(c(xycen),nx2+ny2),ncol=ncol(tmpX),byrow=TRUE))^2,na.rm=TRUE),na.rm=TRUE)
		   look[kx,ky] <- Fxy - Fx - Fy
		} else if(is.element(linkage.method,c("median","centroid"))) {
		   if(linkage.method=="median") {
		      if(i>1 & kx == 1) xcen <- XclustObjent[[i]]
		      if(j>1 & ky == 1) ycen <- YclustObjent[[i]]
		   }
		   look[kx,ky] <- sqrt(sum((xcen - ycen)^2,na.rm=TRUE))
		} else if(is.element(linkage.method,c("single","complete","average"))) {
		   ind <- cbind(rep(1:nx2,ny2),rep(1:ny2,each=nx2))
		   adist <- zDiss[cbind(lookX[[kx]][ind[,1]],lookY[[ky]][ind[,2]])]
		   if(linkage.method=="single") look[kx,ky] <- min(adist,na.rm=TRUE)
		   else if(linkage.method=="complete") look[kx,ky] <- max(adist,na.rm=TRUE)
		   else if(linkage.method=="average") look[kx,ky] <- mean(adist,na.rm=TRUE)
		}
		} # end of inner-inner for 'ky' double loop.
	      } # end of for 'kx' and 'ky' double loop.
	      if((i>1) & j==1) ICdists[[i]] <- list()
	      ICdists[[i]][[j]] <- look
	      minICdists <- c(minICdists, min(look,na.rm=TRUE))
	} # end of inner for 'j' loop.
    } # end of outer for 'i' loop.
   
    if(verbose) cat("\n")
    out$inter.cluster.dist <- ICdists
    out$min.intercluster.dists <- minICdists

    class(out) <- "clusterer"
    return(out)
} # end of 'clusterer.default' function.

print.clusterer <- function(x, ...) {

    a <- attributes(x)
    if(!is.null(a$msg)) cat(a$msg, "\n")
    print(a$data.name)
    cat("\n", "Cluster Analysis performed.\n")

    invisible()
} # end of 'print.clusterer' function.

summary.clusterer <- function(object, ...) {
   out <- object
   args <- list(...)
   if(is.null(args$z)) z <- 1
   else z <- args$z
   if(is.null(args$sigma)) sigma <- sqrt(var(object$min.intercluster.dists,na.rm=TRUE))
   else sigma <- args$sigma
   if(is.null(args$silent)) silent <- FALSE
   else silent <- TRUE
   medio <- quantile(object$min.intercluster.dists,probs=0.5)
   u <- medio + z*sigma
   if(!silent) {
	a <- attributes(object)
        if(!is.null(a$msg)) cat(a$msg, "\n")
        print(a$data.name)
        cat("\n")
	cat("\n", "Matched clusters determined by clusters whose inter-cluster distance is < ", u, "\n")
    }
   out$cutoff <- u
   NCf <- object$NCf
   NCo <- object$NCo
   if(!silent) {
	cat("\n", "Number of forecast objects at each iteration of CA:\n")
	print(NCf)
	cat("\n", "Number of verification objects at each iteration of CA:\n")
	print(NCo)
   }
   NCf <- max(NCf)
   NCo <- max(NCo)
   HMFtab <- array(NA, dim=c(NCo,NCf,3))
   AvgErr <- matrix(NA, NCo, NCf)
   ICDs <- object$inter.cluster.dist
   for(i in 1:NCo) {
	for(j in 1:NCf) {
	   hold <- ICDs[[i]][[j]]
	   hdim <- dim(hold)
	   if(all(hold > u)) {
		HMFtab[i,j,1] <- 0
		HMFtab[i,j,2] <- hdim[1]
		HMFtab[i,j,3] <- hdim[2]
	   } else {
	      hit <- 0
	      xhold <- cbind(rep(1:hdim[1],hdim[2]),rep(1:hdim[2],each=hdim[1]))
	      hold <- cbind(xhold,c(hold))
	      look0 <- numeric(0)
	      done <- FALSE
	      while(!done) {
		if(is.null(dim(hold))) done <- TRUE
		else if(all(hold[,3]>u)) {
		   HMFtab[i,j,2] <- length(unique(hold[,1]))
		   HMFtab[i,j,3] <- length(unique(hold[,2]))
		   done <- TRUE
		} else {
		   hit <- hit+1
		   look0 <- c(look0, c(hold[hold[,3]==min(hold[,3]),3])[1])
		   hold2 <- (1:dim(hold)[1])[hold[,3]==min(hold[,3])]
		   hold <- hold[-hold2,]
		}
	      } # end of while '!done' loop.
	      HMFtab[i,j,1] <- hit
	      AvgErr[i,j] <- mean(look0, na.rm=TRUE)
	   } # end of if else any small enough distances stmts.
	} # end of for 'j' loop.
   } # end of for 'i' loop.
   csifun <- function(x) x[1]/sum(x)
   csi <- apply(HMFtab,1:2,csifun)
   if(!silent) {
	cat("\n", "Critical Success Index for each number of verification (row) and forecast (column) objects:\n")
	print(csi)
   }
   out$csi <- csi
   out$HMF <- HMFtab
   out$AvgErr <- AvgErr
   class(out) <- "summary.clusterer"
   invisible(out)
} # end of 'summary.clusterer' function.

plot.clusterer <- function(x, ..., set.pw=FALSE, icol=c("gray", tim.colors(64)), horizontal=FALSE) {

    a <- attributes(x)
    loc.byrow <- a$loc.byrow

    if(!is.logical(set.pw)) {
	if(!is.numeric(set.pw) || length(set.pw) != 2) stop("plot.clusterer: invalid set.pw argument.")
	par(mfrow=set.pw, oma=c(0,0,2,0))
    } else if(set.pw)  par(mfrow=c(4,2), mar=c(4.1,4.1,4.1,4.1), oma=c(0,0,2,0))
    else if(!is.null(a$msg)) par(oma=c(0,0,2,0))

    X <- x$data$X
    Y <- x$data$Xhat

    X[X<x$threshold[1]] <- 0
    Y[Y<x$threshold[2]] <- 0

    zl <- range(c(c(X),c(Y)),finite=TRUE)

    if(!is.null(a$projection)) proj <- a$projection
    else proj <- FALSE

    if(!is.null(a$map)) domap <- a$map
    else domap <- FALSE

    xd <- a$xdim

    if(proj) {
	loc <- list(x=matrix(a$loc[,1], xd[1], xd[2], byrow=loc.byrow),
		    y=matrix(a$loc[,2], xd[1], xd[2], byrow=loc.byrow))
    }

    if(domap) {
	locr <- apply(a$loc, 2, range, finite=TRUE)
	ax <- list(x=pretty(round(a$loc[,1], digits=2)), y=pretty(round(a$loc[,2], digits=2)))
    }

    if(length(a$data.name) == 3) dn <- a$data.name[-1]
    else dn <- a$data.name

    mainX <- paste(dn[1], "\n(Threshold = ", a$threshold[1], ")", sep="")
    mainY <- paste(dn[2], "\n(Threshold = ", a$threshold[2], ")", sep="")

    if(domap) {

	if(proj) {
	    map(xlim=locr[,1], ylim=locr[,2], type="n")
	    axis(1, at=ax$x, labels=ax$x)
	    axis(2, at=ax$y, labels=ax$y)
	    poly.image(loc$x, loc$y, X, add=TRUE, col=icol, zlim=zl, main=mainX)
	    map(add=TRUE, lwd=1.5)
	    map(add=TRUE, database="state")

	    map(xlim=locr[,1], ylim=locr[,2], type="n")
	    axis(1, at=ax$x, labels=ax$x)
	    axis(2, at=ax$y, labels=ax$y)
	    poly.image(loc$x, loc$y, Y, add=TRUE, col=icol, zlim=zl, main=mainY)
	    map(add=TRUE, lwd=1.5)
            map(add=TRUE, database="state")
	} else {
	    map(xlim=locr[,1], ylim=locr[,2], type="n")
	    axis(1, at=ax$x, labels=ax$x)
	    axis(2, at=ax$y, labels=ax$y)
	    image(as.image(X, nx=xd[1], ny=xd[2], x=a$loc, na.rm=TRUE), add=TRUE, col=icol, zlim=zl, main=mainX)
	    map(add=TRUE, lwd=1.5)
            map(add=TRUE, database="state")

	    map(xlim=locr[,1], ylim=locr[,2], type="n")
	    axis(1, at=ax$x, labels=ax$x)
	    axis(2, at=ax$y, labels=ax$y)
            image(as.image(Y, nx=xd[1], ny=xd[2], x=a$loc, na.rm=TRUE), add=TRUE, col=icol, zlim=zl, main=mainY)
            map(add=TRUE, lwd=1.5)
            map(add=TRUE, database="state")
	}

    } else {
	if(proj) {

	    poly.image(loc$x, loc$y, X, col=icol, zlim=zl, main=mainX)
	    poly.image(loc$x, loc$y, Y, col=icol, zlim=zl, main=mainY)

	} else {

	    image(X, col=icol, zlim=zl, main=mainX)
      	    image(Y, col=icol, zlim=zl, main=mainY)

	}
    } # end of if else 'domap' stmts.  
   image.plot(Y, col=icol, zlim=zl, legend.only=TRUE, horizontal=horizontal)

   plot(x$cluster.objects$X, xlab="")
   plot(x$cluster.objects$Y, xlab="")
   hist(x$min.intercluster.dists, col="darkblue", xlab="Minimum inter-cluster distances \n(between fields)", breaks="FD", main="")
   args <- list(...)
   if(is.null(args$silent)) res <- summary(x,silent=TRUE,...) 
   else res <- summary(x, ...)
   res$par.set <- TRUE
   plot(res)

    if(!is.null(a$msg)) {
	title("")
	mtext(a$msg, line=0.05, outer=TRUE)
    }

   invisible()
} # end of 'plot.clusterer' function.

plot.summary.clusterer <- function(x, ...) {
   NCo <- max(x$NCo)
   NCf <- max(x$NCf)
   MF <- cbind(c(x$HMF[,,1]), c(x$HMF[,,2]),c(x$HMF[,,3]))
   colnames(MF) <- c("Hit", "Miss", "False \nAlarm")
   if(!is.null(x$trans)) {
	if(x$trans=="identity") m1 <- "Average Matched Object Error"
	else m1 <- paste("Average Error for matched objects \n(", x$trans, " transformed)", sep="")
   } else m1 <- "Average Matched Object Error"
   if(is.null(x$par.set)) par(mfrow=c(1,3))
   boxplot(MF, col="darkblue", notch=TRUE, main="Across Scales")
   image(x$csi, col=c("gray",tim.colors(64)), xlab="Number of Verification Objects", ylab="Number of Forecast objects", main="CSI", zlim=c(0,1), axes=FALSE) 
   axis(1, at=seq(0,1,,NCo), labels=1:NCo)
   axis(2, at=seq(0,1,,NCf), labels=1:NCf)
   image.plot(x$csi, col=c("gray",tim.colors(64)), legend.only=TRUE, zlim=c(0,1))
   image(x$AvgErr, col=c("gray",tim.colors(64)), xlab="Number of Verification Objects", ylab="Number of Forecast objects", main=m1, axes=FALSE)
   axis(1, at=seq(0,1,,NCo), labels=1:NCo)
   axis(2, at=seq(0,1,,NCf), labels=1:NCf)
   image.plot(x$AvgErr, col=c("gray",tim.colors(64)), legend.only=TRUE)
   invisible()
} # end of 'plot.summary.clusterer' function.

MakeClusterList <- function(x) {
   n <- dim(x$merge)[1]+1
   out <- list()
   for(i in 1:n) {
	look <- as.numeric(cutree(x,k=n-i+1))
	m <- length(unique(look))
	look <- cbind(1:n,look)
	out[[i]] <- list()
	for(j in 1:m) out[[i]][[j]] <- c(look[look[,2]==j,1])
   } # end of for 'i' loop.
   return(out)
} # end of 'MakeClusterList' function.

CSIsamples <- function(x, ...) {
    UseMethod("CSIsamples", x)
} # end of 'CSIsamples' function.

CSIsamples.SpatialVx <- function(x, ..., time.point=1, model=1, nbr.csi.samples = 100,
                   threshold = 20, k = 100, width = 25, stand=TRUE,
                   z.mult = 0, hit.threshold = 0.1, max.csi.clust = 100,
                   diss.metric="euclidean", linkage.method="average", verbose=FALSE) {

    a <- attributes(x)
    TheCall <- match.call()

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(x, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(x, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(x, model=model)
    else dat <- datagrabber(x)

    X <- dat$X
    Xhat <- dat$Xhat
    ## End: Get the data sets

    out <- CSIsamples.default(x=X, ..., xhat=Xhat, nbr.csi.samples=nbr.csi.samples,
				threshold=threshold, k=k, width=width, stand=stand,
				z.mult=z.mult, hit.threshold=hit.threshold,
				max.csi.clust=max.csi.clust, diss.metric=diss.metric,
				linkage.method=linkage.method, verbose=verbose)

    attr(out, "time.point") <- time.point
    attr(out, "model") <- model

    if(length(a$data.name) == a$nforecast + 2) {
        dn <- a$data.name[-(1:2)]
        vxname <- a$data.name[1:2]
    } else {
        dn <- a$data.name[-1]
        vxname <- a$data.name[1]
    }
    if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model

    attr(out, "data.name") <- c(vxname, dn[model.num])

    attr(out, "xdim") <- a$xdim
    attr(out, "msg") <- a$msg
    out$call <- TheCall
    return(out)
} # end of 'CSIsamples.SpatialVx' function.

CSIsamples.default <- function(x, ..., xhat, nbr.csi.samples = 100,
                   threshold = 20, k = 100, width = 25, stand=TRUE,
                   z.mult = 0, hit.threshold = 0.1, max.csi.clust = 100,
                   diss.metric="euclidean", linkage.method="average", verbose=FALSE) {

  ##
  ## Internal functions
  ##
    X <- x
    Y <- xhat
  out <- list()
  data.name <- c(as.character(substitute(X)),as.character(substitute(Y))) 
  names(data.name) <- c("verification","forecast")
  out$data.name <- data.name
  out$call <- match.call()

  ## convert and threshold the image
  convert.image <- function(raw.file, threshold = 0) {
        mat <- as.matrix(raw.file)
        dims <- dim(mat)
        res <- matrix(0, ncol = 3, nrow = dims[1]*dims[2])
        res[,1] <- rep(1:dims[2], rep(dims[1], dims[2]))
        res[,2] <- rep(1:dims[1], dims[2])
        res[,3] <- as.vector(mat)
        res <- res[res[,3] > threshold,]
        if(!is.matrix(res)) res <- matrix(res, ncol = 3)
        return(res)
   } # end of internal 'convert.image' function.

 ## standardize the two files, separately for each dimension
 ## return a list object with the revised data matrices
 stdize.xyz <- function(raw.file1, raw.file2, threshold = 0, stand=TRUE) {
        conv1 <- convert.image(raw.file1, threshold)
        conv2 <- convert.image(raw.file2, threshold)
        # nrows <- nrow(conv1)+nrow(conv2)
	if(stand) {
           stdze.mean <- colMeans(rbind(conv1, conv2))
           stdze.sd <- apply(rbind(conv1, conv2), 2, sd)
           for (i in 1:ncol(conv1)){
	      conv1[,i] <- (conv1[,i] - stdze.mean[i])/stdze.sd[i]
	      conv2[,i] <- (conv2[,i] - stdze.mean[i])/stdze.sd[i]
           }
	} # end of if 'stand' stmt.
        return(list(std.table1 = conv1, std.table2 = conv2))
   } # end of internal 'stdize.xyz' function.

   ## Sample points from a single dataframe with cluster assignments
   sample.pts <- function(orig.pts, clust.assign, width) {
        n <- nrow(orig.pts)
        resample <- function(x, size, ...) {
           if(length(x) == 1) out <- rep(x, size)
           else out <- sample(x, size, ...)
           return(out)
        } # end of internal-internal 'resample' function.
        #Likely faster to use sapply, but less understandable
        new.data <- matrix(0, nrow = max(clust.assign), ncol = dim(orig.pts)[2]*width)
        for (i in 1:(max(clust.assign))){
           row.sample <- resample((1:n)[clust.assign == i], size = width, rep = T)
           new.data[i,] <- as.vector(orig.pts[row.sample,])
        } # end of for 'i' loop.
        return(new.data)
   } # end of internal 'sample.pts' function.

  generate.sample.csi <- function(nbr.clusts, orig.pts, clust.assign, width, source.lab, toss.out, diss.metric, linkage.method, verbose) {
        point.sample <- sample.pts(orig.pts, clust.assign, width)
        h <- dist(point.sample, method=diss.metric)
        sampled.cluster <- hclust(h, method=linkage.method)
        temp.cutree <- cutree(sampled.cluster, 1:nbr.clusts)
        point.count <- tapply(source.lab, list(clust.assign, source.lab), length)
        if (ncol(point.count) == 1) {
           if (source.lab[1] == 0) point.count <- cbind(point.count, 0)
           else point.count <- cbind(rep(0, nrow(point.count)), point.count)
        } else {
           point.count[is.na(point.count[,1]),1] <- 0
           point.count[is.na(point.count[,2]),2] <- 0
        }
        csi.scores <- NULL
        for (i in 1:nbr.clusts){
           if(verbose) cat(i, " ")
           temp.sums <- t(sapply(1:i, function(m) if(sum(temp.cutree[,i] == m) == 1) point.count[temp.cutree[,i] == m,] else colSums(point.count[temp.cutree[,i] == m,]), simplify=TRUE))
           temp.raw.scores <- temp.sums[,2]/(temp.sums[,2]+temp.sums[,1])
           csi.scores <- c(csi.scores, mean(abs(temp.raw.scores-.5) < .5 - toss.out))
        } # end of for 'i' loop.
        if(verbose) cat("\n")
        return(csi.scores)
   } # end of internal 'generate.sample.csi' internal function.

  #rough check of inputs
  stopifnot(is.data.frame(Y) || is.matrix(Y),
    z.mult >= 0, (hit.threshold >=0 || hit.threshold <= 0.5),
    max.csi.clust >= 1, max.csi.clust <= k)

  #check of dimension of the tables
  stopifnot(nrow(Y) == nrow(X) || ncol(Y) == ncol(X))

  #create standardized matrices; list object of length 2
  std.data <- stdize.xyz(Y, X, threshold, stand=stand)
  min.pts <- min(c(nrow(std.data[[1]]), nrow(std.data[[2]])))

  #check if there are at least enough rows for kmeans step
  if(min.pts <= k) stop(paste("CSIsamples: There may not be enough valid post-threshold points (", min.pts, ") compared to k (", k, ")",  sep=""))

  #warn if data guarantees misses for at least some clusters.
  if(verbose) {
  if (min.pts < max.csi.clust){
        if(verbose) cat("Note: the minimum number of points from the two files is less \nthan the max CSI cluster; this guarantees a zero CSI for some clusters.\n")
        }
  } # end of if 'verbose' stmts.

  #calculation steps

  #create labels to identify the forecast and observation tables
  labs.fore <- rep(0, nrow(std.data[[1]]))
  labs.obs <- rep(1, nrow(std.data[[2]]))

  source.lab <- c(labs.fore, labs.obs)
  all.data <- rbind(std.data[[1]], std.data[[2]])

  if(verbose) cat("Applying scaling factor.\n")
  all.data[,3] <- all.data[,3]*z.mult

  if(verbose) cat("Finding K-means cluster assignments\n")
  init.assign <- kmeans(all.data, k)$cluster

  if(verbose) cat("Computing the appropriate number of CSI values.\n")
  results <- replicate(nbr.csi.samples, generate.sample.csi(max.csi.clust, all.data, init.assign,
                        width, source.lab, hit.threshold, diss.metric=diss.metric, linkage.method=linkage.method, verbose=verbose))
  results <- as.data.frame(results, row.names = 1:max.csi.clust)
  names(results) <- paste("sample", 1:nbr.csi.samples, sep = "")
  out$results <- results
  class(out) <- "CSIsamples"
  return(out)
} # end of 'CSIsamples.default' function.

summary.CSIsamples <- function(object, ...) {
   out <- object
   args <- list(...)
   if(is.null(args$silent)) silent <- FALSE
   else silent <- args$silent
   csi <- rowMeans(object$results, na.rm=TRUE)
   out$csi <- csi
   if(!silent) {
	cat("\n", "Sample average CSI by number of clusters:\n")
 	tmp <- matrix(csi, nrow=1)
	colnames(tmp) <- 1:length(csi)
	print(tmp)
   }
   class(out) <- "summary.CSIsamples"
   invisible(out)
} # end of 'summary.CSIsamples' function.

plot.CSIsamples <- function(x, ...) {
    y <- summary(x, silent=TRUE)
    plot(y, ...)
    invisible(y)
} # end of 'plot.CSIsamples' function.

plot.summary.CSIsamples <- function(x, ...) {
    a <- attributes(x)
    if(!is.null(a$data.name) && length(a$data.name) == 2) m1 <- paste(a$data.name[1], " vs ", a$data.name[2], sep="")
    else if(!is.null(a$data.name) && length(a$data.name) == 3) m1 <- paste(a$data.name[2], " vs ", a$data.name[3], sep="")
    else if(!is.null(x$data.name)) m1 <- paste(x$data.name[1], " vs ", x$data.name[2], sep="")
    plot(x$csi, type="l", ylim=c(0,1), main=m1, xlab="Number of Clusters", ylab="CSI", lwd=1.5)
    invisible()
} # end of 'plot.summary.CSIsamples' function.

print.CSIsamples <- function(x, ...) {
    a <- attributes(x)
    if(!is.null(a$msg)) {
	cat("\n", a$msg, "\n")
	print(a$data.name)
    }
    print(x$call)
    d <- dim(x$results)
    if(!is.null(d)) cat("\n", d[2], " samples of ", d[1], "\n")
    invisible()
} # end of 'print.CSIsamples' function.
