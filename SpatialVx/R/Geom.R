Cindex <- function(x, thresh=NULL, connect.method="C", ...) {

    UseMethod("Cindex", x)

} # end of 'Cindex' function.

Cindex.SpatialVx <- function(x, thresh=NULL, connect.method="C", ..., time.point=1, model=1) {

    a <- attributes(x)

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(x, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(x, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(x, model=model)
    else dat <- datagrabber(x)

    X <- dat$X
    Xhat <- dat$Xhat
    ## End: Get the data sets

    res1 <- Cindex.default(x=X, thresh=thresh, connect.method=connect.method, ...)
    res2 <- Cindex.default(x=Xhat, thresh=thresh, connect.method=connect.method, ...)

    res <- c(res1, res2)

    if(length(a$data.name) == a$nforecast + 2) {

        dn <- a$data.name[-(1:2)]
        vxname <- a$data.name[2]

    } else {

        dn <- a$data.name[-1]
        vxname <- a$data.name[1]

    }

    if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model

    names(res) <- c(vxname, dn[model.num])

    return(res)

} # end of 'Cindex.SpatialVx' function.

Cindex.default <- function(x, thresh=NULL, connect.method="C", ...) {

   if(!is.null(thresh)) x[x < thresh] <- 0

   NP <- sum(colSums(x==0,na.rm=TRUE),na.rm=TRUE)

   x[x==0] <- NA
   x <- as.im(x)
   x <- connected(x, method = connect.method)

   NC <- max(as.numeric(x$v), na.rm = TRUE)

   return(1 - (NC - 1)/(sqrt(NP) + NC))

} # end of 'Cindex.default' function.

Sindex <- function(x, thresh=NULL, ...) {

    UseMethod("Sindex", x)

} # end of 'Sindex' function.

Sindex.SpatialVx <- function(x, thresh=NULL, ..., time.point=1, model=1) {

    a <- attributes(x)

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(x, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(x, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(x, model=model)
    else dat <- datagrabber(x)

    X <- dat$X
    Xhat <- dat$Xhat
    ## End: Get the data sets

    res1 <- Sindex.default(x=X, thresh=thresh, ..., loc=a$loc)
    res2 <- Sindex.default(x=Xhat, thresh=thresh, ..., loc=a$loc)

    res <- rbind(res1, res2)

    if(length(a$data.name) == a$nforecast + 2) {
        dn <- a$data.name[-(1:2)]
        vxname <- a$data.name[2]
    } else {
        dn <- a$data.name[-1]
        vxname <- a$data.name[1]
    }
    if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model

    rownames(res) <- c(vxname, dn[model.num])

    return(res)

} # end of 'Sindex.SpatialVx' function.

Sindex.default <- function(x, thresh=NULL, ..., loc=NULL) {

    if(!is.null(thresh)) x[x < thresh] <- 0
    n <- sum(colSums(x>0,na.rm=TRUE),na.rm=TRUE)
    n2 <- sqrt(n)

    if( floor(n2) == n2 ) Pmin <- 4 * n2
    else Pmin <- 2 * ( floor( 2 * n2 ) + 1 )

    if(is.null(loc)) {

	xdim <- dim(x)
	loc <- cbind(rep(1:xdim[1], xdim[2]), rep(1:xdim[2], each = xdim[1]))

    }

    id <- c(x) != 0
    corners <- apply(loc[id,], 2, range, finite = TRUE)
    P <- 2 * ( diff( corners[,1] ) + 1 ) + 2 * ( diff( corners[,2] ) + 1 )

    res <- c(Pmin / P, Pmin, P)
    names(res) <- c("Sindex", "Pmin", "P")

    return(res)

} # end of 'Sindex.default' function.

Aindex <- function(x, thresh=NULL, dx=1, dy=1, ...) {

    UseMethod("Aindex", x)

} # end of 'Aindex' function.

Aindex.SpatialVx <- function(x, thresh=NULL, dx=1, dy=1, ..., time.point=1, model=1) {
    a <- attributes(x)

    ## Begin: Get the data sets
    if(!missing(time.point) && !missing(model)) dat <- datagrabber(x, time.point=time.point, model=model)
    else if(!missing(time.point)) dat <- datagrabber(x, time.point=time.point)
    else if(!missing(model)) dat <- datagrabber(x, model=model)
    else dat <- datagrabber(x)
   
    X <- dat$X
    Xhat <- dat$Xhat
    ## End: Get the data sets

    res1 <- Aindex.default(x=X, thresh=thresh, dx=dx, dy=dy, ...)
    res2 <- Aindex.default(x=Xhat, thresh=thresh, dx=dx, dy=dy, ...)

    res <- rbind(res1, res2)

    if(length(a$data.name) == a$nforecast + 2) {
        dn <- a$data.name[-(1:2)]
        vxname <- a$data.name[2]
    } else {
        dn <- a$data.name[-1]
        vxname <- a$data.name[1]
    }
    if(!is.numeric(model)) model.num <- (1:a$nforecast)[dn == model]
    else model.num <- model

    rownames(res) <- c(vxname, dn[model.num])

    return(res)
} # end of 'Aindex.SpatialVx' function.

Aindex.default <- function(x, thresh=NULL, dx=1, dy=1, ...) {
    if(is.null(thresh)) thresh <- 1e-8 
    x <- as.im(x)
    x <- solutionset(x>=thresh)
    ch <- convexhull(x)
    A <- sum(colSums(x$m,na.rm=TRUE),na.rm=TRUE)*dx*dy
    Aconvex <- area.owin(ch)*dx*dy
    Aindex <- A/Aconvex

    res <- c(Aindex, A, Aconvex, dx, dy)
    names(res) <- c("Aindex", "A", "Aconvex", "dx", "dy")

    return(res)
} # end of 'Aindex.default' function.
