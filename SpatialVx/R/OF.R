OF <- function(x, ...) {
    UseMethod("OF", x)
} # end of 'OF' function.

OF.SpatialVx <- function(x, ..., time.point=1, model=1, W=5, grads.diff=1, center=TRUE, cutoffpar=4, verbose=FALSE) {

    a <- attributes(x)

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(x, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(x, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(x, model=model)
    else dat <- datagrabber(x)
   
    X <- dat$X
    Xhat <- dat$Xhat
    ## End: Get the data sets

    out <- OF.default(x=X, xhat=Xhat, W=W, grads.diff=grads.diff, center=center, verbose=verbose)

    if(length(a$data.name) == a$nforecast + 2) {
        dn <- a$data.name[-(1:2)]
        vxname <- a$data.name[2]
    } else {
        dn <- a$data.name[-1]
        vxname <- a$data.name[1]
    }
    if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model

    attr(out, "time.point") <- time.point
    attr(out, "model") <- model

    attr(out, "msg") <- a$msg
    attr(out, "data.name") <- c(vxname, dn[model.num])

    attr(out, "map") <- a$map
    attr(out, "projection") <- a$projection
    attr(out, "loc") <- a$loc

    return(out)
} # end of 'OF.SpatialVx' function.

OF.default <- function(x, ..., xhat, W=5, grads.diff=1, center=TRUE, cutoffpar=4, verbose=FALSE) {
   out <- list()
   data.name <- c(as.character(substitute(x)), as.character(substitute(xhat)))
   names(data.name) <- c("verification","forecast")
   out$data.name <- data.name
   out$call <- match.call()
   if(verbose) begin.tiid <- Sys.time()

    out$data <- list(x=x, xhat=xhat)

   ##
   ## Internal functions.
   ##
   degtorad <- function(degree) return(2*pi*degree/360)
   radtodeg <- function(radian) return(360*radian/(2*pi))

   myatan <- function(u,v){
	# This one gives angles between 0 and 360 .
	if(u>0 && v>=0) return( radtodeg (atan( v/u )) )
	else if(u<0 && v>=0) return( 180 + radtodeg ( atan( v/u )) )
	else if(u<0 && v<=0) return( 180 + radtodeg ( atan( v/u )) )
	else if(u>0 && v<=0) return( 360 + radtodeg ( atan( v/u )) )
	else if(u==0 && v>0) return(90)
	else if(u==0 && v<0) return(270)
	else if(u==0 && v==0) return(0)
    } # end of 'myatan' internal function.

   cutoff <- function(field,k) {
	# sets all elements of field to NaN, if
	# they are outside of k std dev of median.
	temp1 <- median(as.vector(field),na.rm=TRUE)
	temp2 <- sqrt(var(as.vector(field),na.rm=TRUE))
	field[ field<temp1-k*temp2 ]=NaN # sqrt(-1)
	field[ field>temp1+k*temp2 ]=NaN # sqrt(-1)
	return(field)
   } # end of internal 'cutoff' function.

   ##
   ## End of internal functions.
   ##

   if(center) {
	if(verbose) cat("\n", "Calculating mean field.\n")
	N <- sum(!is.na(xhat), na.rm=TRUE)
	mean.field <- sum(colSums(xhat, na.rm=TRUE),na.rm=TRUE)/N
   } # end of if 'center' stmt.

   nr <- nrow(x)
   nc <- ncol(x)
   rows <- seq( 2, nr-2 , 1 )  # Start with 2 to assure the smallest window is big
   cols <- seq( 2 ,nc-2 , 1 )  # enough to allow for derivatives in optflow().
   out$rows <- rows
   out$cols <- cols

   err.add.lin <- err.mag.lin <- err.ang.lin <-
   err.add.nlin <- err.mag.nlin <- err.ang.nlin <-
   err.vc.lin <- err.vr.lin <- err.vc.nlin <- err.vr.nlin <- matrix(NaN, nrow=nr,ncol=nc,byrow=T)

   if(verbose) cat("\n", "Finding optical flow for each window of size ", W, "\n")
   for(i in rows ) {   # along y axis
	if(verbose) cat(i, " ")
	for(j in cols ) {   # along x axis
	   # if(verbose) cat(j, " ")
	   initial.crop <- xhat[ pmax(i-0.5*(W-1),1) : pmin(i+0.5*(W-1),nr) , pmax(j-0.5*(W-1),1) : pmin(j+0.5*(W-1),nc) ]
	   final.crop <- x[ pmax(i-0.5*(W-1),1) : pmin(i+0.5*(W-1),nr) , pmax(j-0.5*(W-1),1) : pmin(j+0.5*(W-1),nc) ]

	   if(!center) of <- optflow(initial.crop, final.crop, grads.diff=grads.diff, ...)   # from forecast to obs.  
	   if(center) {
		if(grads.diff==1) grads <- grads1(xhat)
		else if(grads.diff==2) grads <- grads2(xhat)
		else stop("OF: grads.diff must be 1 or 2.")
		mean.grads <- apply(grads,2,mean,na.rm=TRUE)

		of <- optflow(initial.crop, final.crop, grads.diff=grads.diff, mean.field=mean.field, ...)
		err.add.nlin[i,j] <- of[1]-mean.field-of[2]*mean.grads[1]-of[3]*mean.grads[2]-
					0.5*of[2]^2*mean.grads[3]- 0.5*of[3]^2*mean.grads[4]- 0.5*of[2]*of[3]*mean.grads[5]
		err.add.lin[i,j]=of[4]-mean.field-of[5]*mean.grads[1]-of[6]*mean.grads[2]
	   } else {
		err.add.nlin[i,j] <- of[1]
		err.add.lin[i,j] <- of[4]
	   } # end of if else 'center' stmts.

	   err.vr.nlin[i,j] <- -of[2]
	   err.vc.nlin[i,j] <- -of[3]

	   err.vr.lin[i,j] <- -of[5]
	   err.vc.lin[i,j] <- -of[6]

	   # error in magnitude and angle of OF vectors:
	   err.mag.nlin[i,j] <- sqrt( err.vc.nlin[i,j]^2 + err.vr.nlin[i,j]^2 )
	   err.ang.nlin[i,j] <- myatan(err.vc.nlin[i,j],err.vr.nlin[i,j])
	   err.mag.lin[i,j] <- sqrt( err.vc.lin[i,j]^2 + err.vr.lin[i,j]^2 )
	   err.ang.lin[i,j] <- myatan(err.vc.lin[i,j],err.vr.lin[i,j])
	} # end of for 'j' loop.
   } # end of for 'i' loop.
   if(verbose) cat("\n")
   out$err.add.lin <- err.add.lin
   out$err.mag.lin <- err.mag.lin
   out$err.ang.lin <- err.ang.lin
   out$err.add.nlin <- err.add.nlin
   out$err.mag.nlin <- err.mag.nlin
   out$err.ang.nlin <- err.ang.nlin
   out$err.vc.lin <- err.vc.lin
   out$err.vr.lin <- err.vr.lin
   out$err.vc.nlin <- err.vc.nlin
   out$err.vr.nlin <- err.vr.nlin
   if(verbose) {
	cat("\n", "Total run time was:")
	print(Sys.time() - begin.tiid)
	cat("\n")
   } # end of if 'verbose' stmt.
   class(out) <- "OF"
   return(out)
} # end of 'OF.default' function.

optflow <- function(initial, final, grads.diff=1, mean.field=NULL, ...) {
   args <- list(...)
   nnr <- nrow(initial)
   nnc <- ncol(initial)

   eps <- 0.0001    # for collinearity.

   if(grads.diff==1) grads <- grads1(initial)
   else if(grads.diff==2) grads <- grads2(initial)
   else stop("OF: grads.diff must be 1 or 2.")

   Ir <- grads[,1]
   Ic <- grads[,2]
   Irr <- grads[,3]
   Icc <- grads[,4]
   Irc <- grads[,5]
   
   y <- as.vector(final)
   x1 <- as.vector(initial)
   x2 <- Ir
   x3 <- Ic 
   x4 <- Irr 
   x5 <- Icc
   x6 <- Irc

   if(!is.null(mean.field)) {
	# center all the covariates (not the response), in order to remove
	# the correlation between the estimates.
	mean.grads <- apply(grads,2,mean,na.rm=TRUE)
	x1 <- x1 - mean.field
	x2 <- x2 - mean.grads[1]
	x3 <- x3 - mean.grads[2]
	x4 <- x4 - mean.grads[3]
	x5 <- x5 - mean.grads[4]
	x6 <- x6 - mean.grads[5]
   } # end of if 'center' stmts.

   lm.1 <- lm( y ~ x2 + x3, data=data.frame(y=y-x1, x2=x2, x3=x3))
   beta <- lm.1$coef[1:3]        # for use in initializing optim(), below.

  ##
  ## Internal functions
  ##

   fr2 <- function(a) mean( (y - x1 - a[1] - a[2]*x2 - a[3]*x3 - 0.5*(a[2]^2)*x4 - 0.5*(a[3]^2)*x5 - 0.5*a[2]*a[3]*x6 )^2 ,na.rm=TRUE)

   grr2 <- function(a){
	temp <- y - x1 - a[1] - a[2]*x2 - a[3]*x3 - 0.5*(a[2]^2)*x4 - 0.5*(a[3]^2)*x5 - 0.5*a[2]*a[3]*x6
	return(c(-2*mean(temp,na.rm=T), -2*mean(temp*(x2 + x4*a[2] + 0.5*x6*a[3] ),na.rm=T), -2*mean(temp*(x3 + x5*a[3] + 0.5*x6*a[2] ),na.rm=T)))
   } # end of 'grr2' internal function.

   if(is.null(args) | length(args)==0) om.1 <- optim(beta,fr2,grr2,method="BFGS",control=list(maxit=20))
   else om.1 <- optim(beta,fr2,grr2,method="BFGS",...)

   return(c(om.1$par[1], om.1$par[2], om.1$par[3], lm.1$coefficients[1], lm.1$coefficients[2], lm.1$coefficients[3]))
} # end of 'optflow' function.

grads1 <- function(initial) {
   # first deriv, with first difference
   nnr <- nrow(initial)   # points along x direction *in image*.
   nnc <- ncol(initial)   # points along y direction *in image*.

   Ir <- Ic <- Irr <- Icc <- Irc <- matrix(nrow=nnr,ncol=nnc )
   for(i in c(1:(nnr-2)) ){
      for(j in c(1:(nnc-2)) ){
           Ir[i,j] <- initial[i+1,j]- initial[i,j]
           Ic[i,j] <- initial[i,j+1]- initial[i,j]
           Irr[i,j] <- initial[i+2,j]- 2*initial[i+1,j] + initial[i,j]
           Icc[i,j] <- initial[i,j+2]- 2*initial[i,j+1] + initial[i,j]
           Irc[i,j] <- initial[i+1,j+1]-initial[i,j+1]-initial[i+1,j]+initial[i,j]
        } # end of for 'j' loop.
    } # end of for 'i' loop.
    return( cbind(as.vector(Ir), as.vector(Ic), as.vector(Irr), as.vector(Icc), as.vector(Irc) ) )
} # end of internal 'grads1' function.

grads2 <- function(initial) {
    # first deriv, but second difference.
    nnr <- nrow(initial)
    nnc <- ncol(initial)

    Ir <- Ic <- Irr <- Icc <- Irc <- matrix(nrow=nnr,ncol=nnc )
    for(i in c(3:(nnr-2)) ){
       for(j in c(3:(nnc-2)) ){
          Ir[i,j] <- (-initial[i+2,j]+8*initial[i+1,j]-8*initial[i-1,j]+initial[i-2,j])/12
          Ic[i,j] <- (-initial[i,j+2]+8*initial[i,j+1]-8*initial[i,j-1]+initial[i,j-2])/12
          Irr[i,j] <- initial[i+2,j]- 2*initial[i+1,j] + initial[i,j]
          Icc[i,j] <- initial[i,j+2]- 2*initial[i,j+1] + initial[i,j]
          Irc[i,j] <- initial[i+1,j+1]-initial[i,j+1]-initial[i+1,j]+initial[i,j]
       } # end of for 'j' loop.
    } # end of for 'i' loop.
    return( cbind(as.vector(Ir), as.vector(Ic), as.vector(Irr), as.vector(Icc), as.vector(Irc) ) )
} # end of internal 'grads2' function.

plot.OF <- function(x, ...) {

   initial <- x$data$xhat
   final   <- x$data$x
   nr <- nrow(initial)
   nc <- ncol(initial)

   levels <- round(seq(min(initial,final),max(initial,final),40),2)
   args <- list(...)
   if(is.null(args$full)) full <- FALSE
   else full <- args$full

   if(is.null(args$scale)) scale <- 1
   else scale <- args$scale

   if(is.null(args$of.scale)) of.scale <- 1
   else of.scale <- args$of.scale

   if(is.null(args$of.step)) of.step <- 4
   else of.step <- args$of.step

   if(is.null(args$prop)) prop <- 2
   else prop <- args$prop

   ang <- rad(c(x$err.ang.nlin))
   ang <- ang[!is.na(ang)]

   if(is.null(args$nbins)) nbins <- 40
   else nbins <- args$nbins

   m1 <- paste("Optical Flow (", x$data.name[2], " vs ", x$data.name[1], ")", sep="")

   if(!full) {
	layout(matrix(c(1,1,1,1,2,3,4,5,6,7),5,2,byrow=TRUE), c(2,2), c(1,1), TRUE)
	par(mar=c(1,1,1,1))
	contour(t(final)/scale,axes=FALSE,col=4,main=m1, levels=levels,method="edge",labcex=0.75,vfont=c("serif","bold") )
   	contour(t(initial)/scale,axes=FALSE,add=TRUE,col=2,levels=levels,method="flattest",labcex=0.75,vfont=c("serif","bold"))
	box()

	for(i in x$rows ) {   # along y axis
	   for(j in x$cols ) {   # along x axis
		if((i+1)%%of.step==0 && (j+1)%%of.step==0)
		arrows((j-1)/(nc-1),(i-1)/(nr-1), (j-1+of.scale*x$err.vc.nlin[i,j])/(nc-1), (i-1+of.scale*x$err.vr.nlin[i,j])/(nr-1),length=0.03,col=1)
	   } # end of for 'j' loop.
	} # end of for 'i' loop.

	par(mar=c(2,1,1,1))
	image(t(x$err.add.nlin)/scale, col=gray((500:0)/500),axes=F,main="Intensity Error")
	contour(t(x$err.add.nlin)/scale,nlevels=10,add=TRUE,col=2,vfont=c("serif","bold"),method="simple"); box()

	par(mar=c(2,1,1,1))
	hist(x$err.add.nlin/scale,breaks=200,xlab="",ylab="",main="")
	box()

	par(mar=c(2,1,1,1))
	image(t(x$err.mag.nlin), col=gray((500:0)/500),axes=F,main="Displacement Error")
	contour(t(x$err.mag.nlin),nlevels=10,add=TRUE,col=2,vfont=c("serif","bold")); box()

	par(mar=c(2,1,1,1))
	hist(x$err.mag.nlin,breaks=500,xlab="",ylab="",main="")
	box()

	par(mar=c(2,1,1,1))
	image(t(x$err.ang.nlin), col=gray((500:0)/500),axes=F,main="Angular Error")
	contour(t(x$err.ang.nlin),nlevels=8,add=TRUE,col=2,vfont=c("serif","bold"))
	box()

	par(mar=c(2,1,1,1))
	rose.diag(ang, bins = nbins, prop=prop, main="")
   } else {
	layout(1)
	par(mfrow=c(4,4))
	par(mar = c(2, 2, 1, 1))
	contour(t(final)/scale,axes=TRUE,col=4,main=m1,levels=levels,method="edge",labcex=0.75,vfont=c("serif","bold") )
	contour(t(initial)/scale,axes=FALSE,add=TRUE,col=2,levels=levels,method="flattest",labcex=0.75,vfont=c("serif","bold"))

	for(i in x$rows ) {   # along y axis
	   for(j in x$cols ) {   # along x axis
		if((i+1)%%of.step==0 && (j+1)%%of.step==0)
		arrows((j-1)/(nc-1),(i-1)/(nr-1), (j-1+of.scale*x$err.vc.nlin[i,j])/(nc-1), (i-1+of.scale*x$err.vr.nlin[i,j])/(nr-1),length=0.03,col=1)
	   } # end of for 'j' loop.
	} # end of for 'i' loop.  
	box()

	plot(c(1,1),axes=FALSE,type="n")

	image(t(initial-final), col=gray((500:0)/500),axes=FALSE)   # min=white, max=black
	contour(t(initial-final),nlevels=10,add=T,col=2); box()
	hist(initial-final,breaks=100,main="")

	image(t(x$err.add.lin), col=gray((500:0)/500),axes=FALSE)
	contour(t(x$err.add.lin),nlevels=10,add=T,col=2); box()
	hist(x$err.add.lin,breaks=100,main="")

	image(t(x$err.add.nlin), col=gray((500:0)/500),axes=FALSE)
	contour(t(x$err.add.nlin),nlevels=10,add=T,col=2); box()
	hist(x$err.add.nlin,breaks=100,main="")

	image(t(x$err.mag.lin), col=gray((500:0)/500),axes=FALSE)
	contour(t(x$err.mag.lin),nlevels=5,add=T,col=2); box()
	hist(x$err.mag.lin[x$err.mag.lin!=0],breaks=500,main="")

	image(t(x$err.mag.nlin), col=gray((500:0)/500),axes=FALSE)
	contour(t(x$err.mag.nlin),nlevels=5,add=T,col=2)
	box()
	hist(x$err.mag.nlin,breaks=500,main="")

	image(t(x$err.ang.lin), col=gray((500:0)/500),axes=FALSE)
	contour(t(x$err.ang.lin),nlevels=5,add=T,col=2); box()
	ang.lin <- rad(c(x$err.ang.lin))
	ang.lin <- ang.lin[!is.na(ang.lin)]
	rose.diag(ang.lin, bins=nbins, prop=prop, main="")

	image(t(x$err.ang.nlin), col=gray((500:0)/500),axes=FALSE)
	contour(t(x$err.ang.nlin),nlevels=5,add=T,col=2)
	box()
	rose.diag(ang, bins=nbins, prop=prop, main="")

   } # end of if else 'full' stmts.
   invisible()
} # end of 'plot.OF' function.

hist.OF <- function(x, ...) {
   out <- x
   out$hist.call <- match.call
   m1 <- paste(x$data.name[2], " vs ", x$data.name[1], sep="")

   args <- list(...)
   if(is.null(args$xmin)) xmin <- 0
   else xmin <- args$xmin

   if(is.null(args$xmax)) xmax <- 360
   else xmax <- args$xmax

   if(is.null(args$ymin)) ymin <- 0
   else ymin <- args$ymin

   if(is.null(args$ymax)) ymax <- 4
   else ymax <- args$ymax

   if(is.null(args$nbreaks)) nbreaks <- 100
   else nbreaks <- args$nbreaks

   y <- x$err.mag.nlin
   x <- x$err.ang.nlin

   xbreaks <- seq(xmin,xmax,length=nbreaks)
   ybreaks <- seq(ymin,ymax,length=nbreaks)
   out$breaks <- list(x=xbreaks, y=ybreaks)

   xb <- numeric(length(xbreaks)-1)
   yb <- numeric(length(ybreaks)-1)
   nb <- matrix(0, nrow=length(xbreaks)-1,ncol=length(ybreaks)-1)
   for(i in 1:(length(xbreaks)-1)){
      for(j in 1:(length(ybreaks)-1)){
   	nb[i,j] <- length( x[x >= xbreaks[i] & x < xbreaks[i+1] & y >= ybreaks[j] & y < ybreaks[j+1] ] )
      } # end of for 'j' loop.
   } # end of for 'i' loop.
   out$hist.vals <- list(xb=xb, yb=yb, nb=nb)

   image(nb,col=c("#FFFFFFFF", gray((500:0)/500)),main=m1, axes=FALSE,xlim=c(0,1),ylim=c(0,1.0),xlab="Angle",ylab="Magnitude")
   axis(1, labels=round(ybreaks,1), at=seq(0,1,length=nbreaks), cex.axis=1)
   axis(2, labels = round(ybreaks,1), at=seq(0,1,length=nbreaks), cex.axis = 1)
   box()
   contour(nb,add=TRUE,col=0)
   invisible(out)
} # end of 'hist.OF' function.

summary.OF <- function(object, ...) {
    a <- attributes(object)
    if(is.null(a$data.name)) cat("Optical Flow Verification for ", object$data.name[2], " into ", object$data.name[1], "\n")
    else if(length(a$data.name == 2)) cat("Optical Flow Verification for ", a$data.name[2], " into ", a$data.name[1], "\n")
    else if(length(a$data.name == 3)) cat("Optical Flow Verification for ", a$data.name[3], " into ", a$data.name[2], "\n")
    cat("Additive Errors (linear):\n")
    print(stats(c(object$err.add.lin)))
    cat("Magnitude (Displacement) Errors (linear):\n")
    print(stats(c(object$err.mag.lin)))
    cat("Angular Errors (linear):\n")
    x <- rad(c(object$err.ang.lin))
    x <- x[!is.na(x)]
    print(circ.summary(x))
    # print(stats(c(object$err.ang.lin)))
    cat("Additive Errors (nonlinear):\n")
    print(stats(c(object$err.add.nlin)))
    cat("Magnitude (Displacement) Errors (nonlinear):\n")
    print(stats(c(object$err.mag.nlin)))
    cat("Angular Errors (nonlinear):\n")
    x <- rad(c(object$err.ang.nlin))
    x <- x[!is.na(x)]
    print(circ.summary(x))
    # print(stats(c(object$err.ang.nlin)))
    invisible()
} # end of 'summary.OF' function.

print.OF <- function(x, ...) {
    print(summary(x))
    invisible()
} # end of 'print.OF' function.
