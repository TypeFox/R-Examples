#
# clipper.R
#
# Interface to Clipper C++ code
#
#  $Revision: 1.13 $ $Date: 2016/03/24 00:52:57 $
#

validxy <- function(P) {
  is.list(P) && all(c("x","y") %in% names(P)) &&
  is.vector(P$x) && is.vector(P$y) && length(P$x)==length(P$y)
}

validpoly <- function(P) {
  is.list(P) && all(unlist(lapply(P, validxy)))
}

xrange <- function(z) { range(z$x) }
yrange <- function(z) { range(z$y) }

ensurexydouble <- function(P) lapply(P[c("x", "y")],
                                     "storage.mode<-", value="double")

ensuredouble <- function(A) lapply(A, ensurexydouble)

aspolygonlist <- function(A) lapply(A, "names<-", value=c("x", "y"))

polysimplify <-
  function(A,
           ...,
           eps, x0, y0,
           filltype=c("evenodd", "nonzero", "positive", "negative")
           ) {
    # validate parameters and convert to integer codes
    filltype <- match.arg(filltype)
    pft <- match(filltype, c("evenodd", "nonzero", "positive", "negative"))
    # validate polygon
    if(!validpoly(A)) {
      if(validxy(A)) A <- list(A) else
      stop("Argument A should be a list of lists, each containing vectors x,y")
    }
    # determine value of 'eps' if missing
    if(missing(eps) || missing(x0) || missing(y0)) {
      xr <- range(range(unlist(lapply(A, xrange))))
      yr <- range(range(unlist(lapply(A, yrange))))
      if(missing(eps)) eps <- max(diff(xr), diff(yr))/1e9
      if(missing(x0)) x0 <- mean(xr)
      if(missing(y0)) y0 <- mean(yr)
    } 
    # call clipper library on each component path
    result <- list()
    A <- ensuredouble(A)
    storage.mode(pft) <- "integer"
    storage.mode(x0) <- storage.mode(y0) <- storage.mode(eps) <- "double"
    result <- .Call("Csimplify",
                    A, pft, x0, y0, eps)
    return(aspolygonlist(result))
  }

polyclip <-
  function(A, B, 
           op=c("intersection", "union", "minus", "xor"),
           ...,
           eps, x0, y0,
           fillA=c("evenodd", "nonzero", "positive", "negative"),
           fillB=c("evenodd", "nonzero", "positive", "negative")
           ) {
    # validate parameters and convert to integer codes
    op <- match.arg(op)
    fillA <- match.arg(fillA)
    fillB <- match.arg(fillB)
    ct <- match(op, c("intersection", "union", "minus", "xor"))
    pftA <- match(fillA, c("evenodd", "nonzero", "positive", "negative"))
    pftB <- match(fillB, c("evenodd", "nonzero", "positive", "negative"))
    # validate polygons and rescale
    if(!validpoly(A)) {
      if(validxy(A)) A <- list(A) else
      stop("Argument A should be a list of lists, each containing vectors x,y")
    }
    if(!validpoly(B)) {
      if(validxy(B)) B <- list(B) else
      stop("Argument B should be a list of lists, each containing vectors x,y")
    }
    # determine value of 'eps' if missing
    if(missing(eps) || missing(x0) || missing(y0)) {
      xr <- range(range(unlist(lapply(A, xrange))),
                  range(unlist(lapply(B, xrange))))
      yr <- range(range(unlist(lapply(A, yrange))),
                  range(unlist(lapply(B, yrange))))
      if(missing(eps)) eps <- max(diff(xr), diff(yr))/1e9
      if(missing(x0)) x0 <- mean(xr)
      if(missing(y0)) y0 <- mean(yr)
    } 
    # call clipper library
    A <- ensuredouble(A)
    B <- ensuredouble(B)
    storage.mode(ct) <- storage.mode(pftA) <- storage.mode(pftB) <- "integer"
    storage.mode(x0) <- storage.mode(y0) <- storage.mode(eps) <- "double"
    ans <- .Call("Cclipbool",
                 A, B, pftA, pftB, ct,
                 x0, y0, eps)
    return(aspolygonlist(ans))
  }

polyoffset <-
  function(A, delta, 
           ...,
           eps, x0, y0,
           miterlim=2, arctol=abs(delta)/100,
           jointype = c("square", "round", "miter")
           ) {
    # validate parameters and convert to integer codes
    jointype <- match.arg(jointype)
    jt <- match(jointype, c("square", "round", "miter")) 
    # validate polygons and rescale
    if(!validpoly(A)) {
      if(validxy(A)) A <- list(A) else
      stop("Argument A should be a list of lists, each containing vectors x,y")
    }
    # determine value of 'eps' if missing
    if(missing(eps) || missing(x0) || missing(y0)) {
      xr <- range(unlist(lapply(A, xrange)))
      yr <- range(unlist(lapply(A, yrange)))
      if(missing(eps)) eps <- max(diff(xr), diff(yr))/1e9
      if(missing(x0)) x0 <- mean(xr)
      if(missing(y0)) y0 <- mean(yr)
    }
    # arc tolerance
    arctol <- max(eps/4, arctol)
    # call clipper library
    A <- ensuredouble(A)
    storage.mode(jt) <- "integer"
    storage.mode(delta) <-
      storage.mode(miterlim) <- storage.mode(arctol) <- "double"
    storage.mode(x0) <- storage.mode(y0) <- storage.mode(eps) <- "double"
    ans <- .Call("Cpolyoffset", A, delta, jt,
                 miterlim, arctol, x0, y0, eps)
    return(aspolygonlist(ans))
  }


polylineoffset <-
  function(A, delta, 
           ...,
           eps, x0, y0,
           miterlim=2, arctol=abs(delta)/100,
           jointype = c("square", "round", "miter"),
           endtype = c("closedpolygon", "closedline",
             "openbutt", "opensquare", "openround",
             "closed", "butt", "square", "round")
           ) {
    ## validate parameters and convert to integer codes
    jointype <- match.arg(jointype)
    jt <- match(jointype, c("square", "round", "miter"))

    endtype <- match.arg(endtype)
    if(endtype == "closed") endtype <- "closedpolygon"
    if(endtype %in% c("butt", "square", "round"))
      endtype <- paste0("open", endtype)
    et <- match(endtype, c("closedpolygon", "closedline",
                           "openbutt", "opensquare", "openround"))
    
    ## validate polygons and rescale
    if(!validpoly(A)) {
      if(validxy(A)) A <- list(A) else
      stop("Argument A should be a list of lists, each containing vectors x,y")
    }
    ## determine value of 'eps' if missing
    if(missing(eps) || missing(x0) || missing(y0)) {
      xr <- range(unlist(lapply(A, xrange)))
      yr <- range(unlist(lapply(A, yrange)))
      if(missing(eps)) eps <- max(diff(xr), diff(yr))/1e9
      if(missing(x0)) x0 <- mean(xr)
      if(missing(y0)) y0 <- mean(yr)
    }
    # arc tolerance
    arctol <- max(eps/4, arctol)
    # call clipper library
    A <- ensuredouble(A)
    storage.mode(jt) <- storage.mode(et) <- "integer"
    storage.mode(delta) <- storage.mode(miterlim) <-
      storage.mode(arctol) <- "double"
    storage.mode(x0) <- storage.mode(y0) <- storage.mode(eps) <- "double"
    ans <- .Call("Clineoffset", A, delta, jt, et,
                 miterlim, arctol, x0, y0, eps)
    return(aspolygonlist(ans))
  }

polyminkowski <-
  function(A, B, 
           ...,
           eps, x0, y0,
           closed=TRUE
           ) {
    # validate parameters and convert to integer codes
    closed <- as.logical(closed)
    # validate polygons/paths
    if(!validpoly(A)) {
      if(validxy(A)) A <- list(A) else
      stop("Argument A should be a list of lists, each containing vectors x,y")
    }
    if(length(A) > 1)
      stop("Not implemented when A consists of more than one polygon")
    if(!validpoly(B)) {
      if(validxy(B)) B <- list(B) else
      stop("Argument B should be a list of lists, each containing vectors x,y")
    }
    # determine value of 'eps' if missing
    if(missing(eps) || missing(x0) || missing(y0)) {
      xr <- range(range(unlist(lapply(A, xrange))))
      yr <- range(range(unlist(lapply(A, yrange))))
      xr <- range(xr, range(unlist(lapply(B, xrange))))
      yr <- range(yr, range(unlist(lapply(B, yrange))))
      if(missing(eps)) eps <- max(diff(xr), diff(yr))/1e9
      if(missing(x0)) x0 <- xr[1]
      if(missing(y0)) y0 <- yr[1] - diff(yr)/16
    }
    # call clipper library on each component path
    A <- ensuredouble(A)
    B <- ensuredouble(B)
    storage.mode(x0) <- storage.mode(y0) <- storage.mode(eps) <- "double"
    storage.mode(closed) <- "logical"
    result <- .Call("Cminksum",
                    A, B, closed, x0, y0, eps)
    return(aspolygonlist(result))
  }
