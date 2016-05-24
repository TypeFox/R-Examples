# Methods for classes *SplineBasis defined in AllClass.R 


###################################################################################
#######          createurs
###################################################################################


MSplineBasis<-function(knots, degree=3, keep.duplicates=FALSE, log=FALSE) {
  # for M-splines,
  # knots are 1 boundary min knots, interior knots, 1 boundary max knots
  # boundary knots are not tested to be duplicated
  # code from orthogonalsplinebasis
        order<-degree+1
	n<-length(knots)
        if( n ==2 ){
# no interior knots
          knots<-c(rep(knots[1], order), rep(knots[n], order))
        }
        else {
          if(any(table(knots[2:(n-1)])>1)&&!keep.duplicates){
            warning("Duplicate interior knots. Removing duplicates.\n    (use keep.duplicates=TRUE to keep duplicates)")
                # modif MGk 06/06/2011 to keep the (order-1) first and last knots 
            oknots <- knots
            iknots<-unique(oknots[2:(n-1)])
            knots<-c(rep(oknots[1], order), iknots, rep(oknots[n], order))
          }
          else {
            # just duplicates first and last knots
            oknots <- knots
            iknots<-oknots[2:(n-1)]
            knots<-c(rep(oknots[1], order), iknots, rep(oknots[n], order))
          }            
        }
        # recompute n number of knots
        n <- length(knots)
	q <- n-order

        SB <- orthogonalsplinebasis::SplineBasis(knots=knots, order=order, keep.duplicates=keep.duplicates)
	M<-SB@Matrices
	new("MSplineBasis", knots=knots, min=knots[1], max=knots[length(knots)],
            degree=as.integer(degree), nbases=as.integer(q), Matrices=M, SplineBasis=SB, log=log)
}

MSplineBasis2<-function(allknots, degree=3, keep.duplicates=FALSE, log=FALSE) {
  # for B-splines, allknots are degree+1 boundary min knots, interior knots, degere+1 boundary max knots
  # the 2x4 boundary knots are equal
  # code from orthogonalsplinebasis
        order<-degree+1
	n<-length(allknots)
        if( n ==2*order ){
# no interior knots
          knots<-allknots  # !! all the boundary knots are given
        }
        else {
          if(any(table(allknots[order:(n-order+1)])>1)&&!keep.duplicates){
		warning("Duplicate interior knots. Removing duplicates.\n    (use keep.duplicates=TRUE to keep duplicates)")
                # modif MGk 06/06/2011 to keep the (order-1) first and last knots 
                oknots <- allknots
		iknots<-unique(oknots[(order+1):(n-order)])
		knots<-c(rep(oknots[order], order), iknots, rep(oknots[n-order+1], order))
              }
        }
	q<-n-order
        SB <- orthogonalsplinebasis::SplineBasis(knots=knots, order=order, keep.duplicates=keep.duplicates)
	M<-SB@Matrices
	new("BSplineBasis", knots=knots, min=knots[1], max=knots[length(knots)],
            degree=as.integer(degree), nbases=as.integer(q), Matrices=M, SplineBasis=SB, log=log)
}

BSplineBasis<-function(knots, degree=3, keep.duplicates=FALSE, log=FALSE) {
  # for B-splines, knots are degree+1 boundary min knots, interior knots, degree+1 boundary max knots
  # boundary knots are not tested to be duplicated
  # code from orthogonalsplinebasis
        order<-degree+1
	n<-length(knots)
        if( n ==2*order ){
# no interior knots
          knots<-knots  # !! all the boundary knots are given
        }
        else {
          if(any(table(knots[order:(n-order+1)])>1)&&!keep.duplicates){
            warning("Duplicate interior knots. Removing duplicates.\n    (use keep.duplicates=TRUE to keep duplicates)")
                # modif MGk 06/06/2011 to keep the (order-1) first and last knots 
            oknots <- knots
            iknots<-unique(oknots[order:(n-order+1)])
            knots<-c(oknots[1:(order-1)], iknots, oknots[(n-order+2):length(oknots)])
          }
        }
	q<-n-order
        SB <- orthogonalsplinebasis::SplineBasis(knots=knots, order=order, keep.duplicates=keep.duplicates)
	M<-SB@Matrices
	new("BSplineBasis", knots=knots, degree=as.integer(degree), nbases=as.integer(q), Matrices=M, SplineBasis=SB, log=log)
}

maketpdegrees <- function(knots, order){
  order - unlist(lapply(table(knots), function(x) x:1))
}


TPSplineBasis<-function(knots, degree=3, min, max, type=c("standard", "increasing"), coef=NULL, log=FALSE, keep.duplicates=FALSE) {
  # knots are interior knots
  type  <- match.arg(type)       # type of tp-spline
  order<-degree+1
  n <- length(knots)
  if( n ==0 ){
  # no interior knots
    knots <- as.numeric(knots)
  }
  else {
    knots <- sort(knots)
    if (!keep.duplicates & any(table(knots) > 1) ) {
      warning("Duplicate interior knots. Removing duplicates.\n")
      knots <- unique(knots)
      degrees <- rep(degree, length(knots))
    }
    else {
      degrees <- maketpdegrees(knots, order)
    }
  }
  n <- length(knots)
  if( is.null(coef)) {
    coef <- rep(1, n+order)
  }
  new("TPSplineBasis", knots=knots, degree=as.integer(degree), nbases=as.integer(n+order),
      min=min, max=max , coef=coef, degrees=as.integer(degrees),
      log=log, type=type)
}

TPRSplineBasis<-function(knots, degree=3, ref=0, min, max, type=c("standard", "increasing"), coef=NULL, log=FALSE, keep.duplicates=FALSE) {
  # knots are interior knots
  type  <- match.arg(type)       # type of tp-spline
  order<-degree+1
  n <- length(knots)
  if( n ==0 ){
  # no interior knots
    knots <- as.numeric(knots)
  }
  else {
    knots <- sort(knots)
    if (!keep.duplicates & any(table(knots) > 1) ) {
      warning("Duplicate interior knots. Removing duplicates.\n")
      knots <- unique(knots)
      degrees <- rep(degree, length(knots))
    }
    else {
      degrees <- maketpdegrees(knots, order)
    }
  }
  n <- length(knots)
  if( is.null(coef)) {
    coef <- rep(1, n+order)
  }
  new("TPSplineBasis", knots=knots, degree=as.integer(degree), nbases=as.integer(n+order),
      min=min, max=max , coef=coef, degrees=as.integer(degrees), ref=ref,
      log=log, type=type)

  new("TPSplineBasis", knots=knots, degree=as.integer(degree), nbases=as.integer(n+degree+1), ref=ref, min=min, max=max, log=log )
}

C0BSplineBasis<-function(knots, degree=3, keep.duplicates=FALSE, log=FALSE) {
  # B-splines, for constraints optimisation : sum(beta_i)_(i<= degree)b_i(0) = 1
  # for B-splines, knots are degree+1 boundary min knots, interior knots, degree+1 boundary max knots
  # code from orthogonalsplinebasis
        order<-degree+1
	n<-length(knots)
	if(any(table(knots[order:(n-order+1)])>1)&&!keep.duplicates){
		warning("Duplicate interior knots. Removing duplicates.\n    (use keep.duplicates=TRUE to keep duplicates)")
                # modif MGk 06/06/2011 to keep the (order-1) first and last knots 
                oknots <- knots
		iknots<-unique(oknots[order:(n-order+1)])
		knots<-c(oknots[1:(order-1)], iknots, oknots[(n-order+2):length(oknots)])
	}
	q<-n-order
        SB <- orthogonalsplinebasis::SplineBasis(knots=knots, order=degree+1, keep.duplicates=keep.duplicates)
	M<-SB@Matrices
	new("C0BSplineBasis", knots=knots, degree=as.integer(degree), nbases=as.integer(q) , Matrices=M, log=log)
}



C0TPSplineBasis<-function(knots, degree=3, log=FALSE) {
  # duplicated notes are suppressed
  # for natural splines, knots are interior knots
  order<-degree+1
    n <- length(knots)
    if (any(table(knots[2:(n - 1)]) > 1) ) {
        warning("Duplicate interior knots. Removing duplicates.\n")
        knots <- unique(knots)
    }
  new("TPSplineBasis", knots=knots, degree=as.integer(degree), log=log)
}

###################################################################################
####### fin          createurs
###################################################################################


setGeneric("getKnots",function(object)standardGeneric("getKnots"))
setMethod("getKnots",signature(".SplineBasis"),function(object)object@knots)

setGeneric("getDegree",function(object)standardGeneric("getDegree"))
setMethod("getDegree",signature(".SplineBasis"),function(object)object@degree)

setGeneric("getOrder",function(object)standardGeneric("getOrder"))
setMethod("getOrder",signature(".SplineBasis"),function(object)object@degree+1)

setGeneric("getNBases",function(object)standardGeneric("getNBases"))
setMethod("getNBases",signature(".SplineBasis"),function(object)object@nbases)

setGeneric("getLog",function(object)standardGeneric("getLog"))
setMethod("getLog",signature(".SplineBasis"),function(object)object@log)

setGeneric("getSplineMatrix",function(object)standardGeneric("getSplineMatrix"))
setMethod("getSplineMatrix",signature("BSplineBasis"),function(object)object@Matrix)

setGeneric("getRef",function(object)standardGeneric("getRef"))
setMethod("getRef",signature("TPSplineBasis"),function(object)object@ref)

setGeneric("getMin",function(object)standardGeneric("getMin"))
setMethod("getMin",signature("TPSplineBasis"),function(object)object@min)

setGeneric("getMax",function(object)standardGeneric("getMax"))
setMethod("getMax",signature("TPSplineBasis"),function(object)object@max)



###################################################################################
#######         GETEUR
###################################################################################


###################################################################################
####### fin          GETEUR
###################################################################################


######################################################################
EvaluateBBasis<-function(object,x,intercept=TRUE, xname=NULL,  ...) { 
  stopifnot(is.numeric(x))
  dots<-list(...)
  nx <- names(x)
  nax <- is.na(x)
  if (nas <- any(nax)){ 
    x <- x[!nax]
  }
  Aknots<-object@knots
  degree<-object@degree
  ol <- x < Aknots[degree+1]
  or <- x > Aknots[length(Aknots)-degree]
  outside <- ( or | ol)
  
  if (any(outside)) {
    basis <- array(, dim= c(length(x), length(Aknots) - degree - 1L))
    if (any(inside <- !outside)){ 
      basis[inside, ] <- spline.des(knots=Aknots, x=x[inside], ord=degree+1)$design
    }
  }
  else {
    basis <- spline.des(knots=Aknots, x=x, ord=degree+1)$design
 }
  if (!intercept) {
    basis <- basis[, -1L, drop = FALSE]
  }

  if (object@log) {
    basis <- cbind(basis, log(x))
  }


  
#add na x
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }

    #add dimnames
    if(!is.null(xname)){
      if (intercept) {
        dbs <- paste("BS-", xname, 1:(getNBases(object)+getLog(object)), sep="")
      }
      else {
        dbs <- paste("BS-", xname, 2:(getNBases(object)+getLog(object)), sep="")
      }
      dimnames(basis) <- list(nx, dbs)
    }



  a <- list(degree = degree, knots =  Aknots, 
            intercept = intercept, log=getLog(object))
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("bs", "basis", "matrix")
  basis
}

# fast evaluate, no dimnames
FEvaluateMBasisold<-function(object, x, intercept=TRUE, outer.ok=TRUE, ...) {

  stopifnot(is.numeric(x))
  nax <- is.na(x)
  if (nas <- any(nax)){ 
    x <- x[!nax]
  }

  Aknots<-object@knots
  degree<-object@degree
  nbases<-object@nbases
  ol <- x < Aknots[degree+1]
  or <- x > Aknots[length(Aknots)-degree]
  outside <- ( or | ol)
  
  if (any(outside)) {
    if(outer.ok) {
      basis <- array(0, dim= c(length(x), nbases))
    }
    else {
      basis <- array(NA, dim= c(length(x), nbases))
    }
    if (any(inside <- !outside)){ 
      basis[inside, ] <- spline.des(knots=Aknots, x=x[inside], ord=degree+1)$design
    }
  }
  else {
    basis <- spline.des(knots=Aknots, x=x, ord=degree+1)$design
  }
  if (!intercept) {
    basis <- basis[, -1L, drop = FALSE]
  }

  if (object@log) {
    basis <- cbind(basis, log(x))
  }

#add na x
  if (nas) {
    nabasis <- matrix(NA, length(nax), nbases)
    nabasis[!nax, ] <- basis
    nabasis
  }
  else {
    basis
  }
}

FEvaluateMBasis<-function(object, x, intercept=TRUE, outer.ok=TRUE, ...) {
  stopifnot(is.numeric(x))
  dots<-list(...)
  M<-object@Matrices
  knots<-object@knots
  order<-object@degree+1
  basis <- .Call("eval_spline_basis", as.double(knots), as.integer(order), M, 
                 as.integer(intercept), as.double(x), as.integer(outer.ok),
                 PACKAGE="flexrsurv")

  if (object@log) {
    return(cbind(basis, log(x)))
  }
 else {
   return(basis)
 }

}


EvaluateMBasis<-function(object,x,intercept=TRUE, xname=NULL, outer.ok=TRUE, namespline= "B", ...) {

  basis <- FEvaluateMBasis(object=object, x=x,
                           intercept=intercept, 
                           xname=xname,
                           outer.ok=TRUE, ...) 

  Aknots<-object@knots
  degree<-object@degree
  nbases<-object@nbases

  
    #add dimnames
  if(!is.null(xname)){
    if (intercept) {
      dbs <- paste(namespline, "-", xname, ":", 1:(getNBases(object)+getLog(object)), sep="")
    }
    else {
      dbs <- paste(namespline, "-", xname, ":", 2:(getNBases(object)+getLog(object)), sep="")
    }
    dimnames(basis)[[2]] <- dbs
  }


#  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = degree, knots =  Aknots, 
            intercept = intercept, log=getLog(object))
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("bs", "basis", "matrix")
  basis
}

# evaluate tp Basis, no ndimnames
FEvaluateTPBasisPosOld<-function(object,x,intercept=TRUE, ref=NULL, xmin=-Inf, xmax=+Inf, outer.ok=TRUE, ...) { 
  stopifnot(is.numeric(x))
  nax <- is.na(x)
  if (nas <- any(nax)){ 
    x <- x[!nax]
  }

  Aknots<-object@knots
  degree<-object@degree
  if (xmin != -Inf){
    ol <- x < xmin
    if (xmax != +Inf){
      or <- x > xmax
      outside <- ( or | ol)
    }
    else {
      outside <- ol
    }
  }
  else {
    if (xmax != +Inf){
      outside <- x > xmax
    }
    else {
      outside <- FALSE
    }
  }

  Aknotref <- Aknots 
  if( !is.null(ref)){
    x <- x - ref
    if(length(Aknots)>0){
      Aknotref <- Aknots - ref
    }
  }

  
  if (intercept) {
    nb <- degree+ 1 + length(Aknots) 
    last <- degree+1
  }
  else{
    nb <- degree  + length(Aknots) 
    last <- degree
  }    
   if( outer.ok) {
    outervalue <- 0.0
  }
  else {
    outervalue <- NA
  }

  basis <- matrix(1, nrow=length(x), ncol=nb)
  if (any(!outside)){ 
    if (intercept) {
      for(i in 1:degree){
        basis[,i+1]<-x^i
      }
    }
    else{
      for(i in 1:degree){
        basis[,i]<-x^i
      }
    }
    if(length(Aknots)>0){
      for (k in 1:(length(Aknots))){
        basis[,last+k] <- ifelse(x>Aknotref[k], (x-Aknotref[k])^degree, 0)
      }
    }
    if (object@log) {
      basis <- cbind(basis, log(x))
    }
  }    

#add outside values
  if(any(outside)){
    basis[outside,]<-outervalue
  }

#add na x
  if (nas) {
    nabasis <- matrix(NA, nrow=length(nax), ncol=nb)
    nabasis[!nax, ] <- basis
    nabasis
  }
  else {
    basis
  }
}

FEvaluateTPBasisPos<-function(object,x,intercept=TRUE, ref=NULL, xmin=-Inf, xmax=+Inf, outer.ok=TRUE, ...) { 
  stopifnot(is.numeric(x))

  min<-object@min
  max<-object@max
  allknots<-object@knots
  order<-object@degree+1
  knots <- unique(allknots)

  if( !is.null(ref)){
    x <- x - ref
    min <- min- ref
    max <- max - ref
    if(length(knots)>0){
      knot <- knots - ref
    }
  }
  
  replicates <- table(allknots)
  degrees <- object@degrees
  coef <- object@coef
  basis <- .Call("eval_trunc_power_basis", as.double(knots), as.double(replicates),
                 as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
                 as.double(x), as.integer(outer.ok),
                 PACKAGE="flexrsurv")

  return(basis)
  
}

# evaluate TP basis, add dimnames
EvaluateTPBasisPos<-function(object,x,intercept=TRUE, ref=NULL, xname=NULL, xmin=-Inf, xmax=+Inf, outer.ok=TRUE, ...) {

  basis <- FEvaluateTPBasisPos(object=object, x=x,
                             intercept=intercept, ref=ref,
                             xname=xname,
                             outer.ok=TRUE, ...) 


  Aknots<-object@knots
  degree<-object@degree
  
  #add dimnames
  if(!is.null(xname)){
    if(degree>2 ){
      if( !is.null(ref)){
        dnva <- c(paste("(",xname, "-", ref, ")", sep=""),
                  paste("(",xname, "-", ref, ")^", 2:degree, sep=""))
      }
      else {
        dnva <- c(xname, paste(xname, "^", 2:degree, sep=""))
      }
    }
    else {
      if( !is.null(ref)){
        dnva <- c(paste("(",xname, "-", ref, ")", sep=""))
      }
      else {
        dnva <- c(xname)
      }
    }
    
    if (intercept) {
      dnva <- c("Intercept", dnva)
    }
    
    if(length(Aknots)>0){
      dnvaplus <- paste("(",xname, "-", Aknots[1:(length(Aknots))], ")_+^", degree, sep="")
    }
    if (getLog(object)) {
      dlog  <- paste("log(",xname, ")", sep="")
    }
    else {
      dlog  <- NULL
    }
    dimnames(basis)[[2]]<-c(dnva, dnvaplus, dlog)
  }

  a <- list(degree = degree, knots =  Aknots, 
            intercept = intercept, ref=ref, log=getLog(object))
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("tp", "basis", "matrix")
  basis
}

  
# similar to FEvaluateTPBasisPos but if knot < ref, -(k-x)+^d if x<ref
FEvaluateTPBasisOld<-function(object,x,intercept=TRUE, ref=NULL, xmin=-Inf, xmax=+Inf, outer.ok=TRUE, ...) { 
  stopifnot(is.numeric(x))
  nax <- is.na(x)
  if (nas <- any(nax)){ 
    x <- x[!nax]
  }

  Aknots<-object@knots
  degree<-object@degree
  if (xmin != -Inf){
    ol <- x < xmin
    if (xmax != +Inf){
      or <- x > xmax
      outside <- ( or | ol)
    }
    else {
      outside <- ol
    }
  }
  else {
    if (xmax != +Inf){
      outside <- x > xmax
    }
    else {
      outside <- FALSE
    }
  }

}

FEvaluateTPBasisstd<-function(object,x,intercept=TRUE, ref=NULL, xmin=-Inf, xmax=+Inf, outer.ok=TRUE, ...) { 
  stopifnot(is.numeric(x))

  min<-object@min
  max<-object@max
  allknots<-object@knots
  order<-object@degree+1
  knots <- unique(allknots)

  if( !is.null(ref)){
    x <- x - ref
    min <- min- ref
    max <- max - ref
    if(length(knots)>0){
      knot <- knots - ref
    }
  }
  
  replicates <- table(allknots)
  degrees <- object@degrees
  coef <- object@coef
  basis <- .Call("eval_trunc_power_increasing_basis", as.double(knots), as.double(replicates),
                 as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
                 as.double(x), as.integer(outer.ok),
                 PACKAGE="flexrsurv")

  return(basis)
}


# truncated powezr basis
FEvaluateTPBasis<-function(object,x,intercept=TRUE, ref=NULL, xmin=-Inf, xmax=+Inf, outer.ok=TRUE, ...) { 
  stopifnot(is.numeric(x))

  min<-object@min
  max<-object@max
  allknots<-object@knots
  order<-object@degree+1
  knots <- unique(allknots)

  if( !is.null(ref)){
    x <- x - ref
    min <- min- ref
    max <- max - ref
    if(length(knots)>0){
      knot <- knots - ref
    }
  }
  
  replicates <- table(allknots)
  degrees <- object@degrees
  coef <- object@coef

  if(object@type == "increasing"){
    basis <- .Call("eval_trunc_power_increasing_basis", as.double(knots), as.double(replicates),
                   as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
                   as.double(x), as.integer(outer.ok),
                   PACKAGE="flexrsurv")
  }
  else {
  basis <- .Call("eval_trunc_power_basis", as.double(knots), as.double(replicates),
                 as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
                 as.double(x), as.integer(outer.ok),
                 PACKAGE="flexrsurv")
  }
  if (object@log) {
    basis <- cbind(basis, log(x))
  }
  return(basis)
}



# same as FEvaluateTPBasis, but with dimnames, attribures, ...
EvaluateTPBasis<-function(object,x,intercept=TRUE, ref=NULL, xname=NULL, 
                           outer.ok=TRUE, ...) { 


  basis <- FEvaluateTPBasis(object=object, x=x,
                             intercept=intercept, ref=ref,
                             xname=xname,
                             outer.ok= outer.ok, ...)  

  Aknots<-object@knots
  degree<-object@degree

  #add dimnames
    if(!is.null(xname)){
      if(degree>2 ){
        if( !is.null(ref)){
          dnva <- c(paste("(",xname, "-", ref, ")", sep=""),
                    paste("(",xname, "-", ref, ")^", 2:degree, sep=""))
        }
        else {
          dnva <- c(xname, paste(xname, "^", 2:degree, sep=""))
        }
      }
      else {
        if( !is.null(ref)){
          dnva <- c(paste("(",xname, "-", ref, ")", sep=""))
        }
        else {
          dnva <- c(xname)
        }
      }
      
      if (intercept) {
        dnva <- c("Intercept", dnva)
      }
      
      if(length(Aknots)>0){
        dnvaplus <- paste("(",xname, "-", Aknots[1:(length(Aknots))], ")_+^", degree, sep="")
      }
      if (getLog(object)) {
        dlog  <- paste("log(",xname, ")", sep="")
      }
      else {
        dlog  <- NULL
      }
      dimnames(basis)[[2]]<-c(dnva, dnvaplus, dlog)
    }

  a <- list(degree = degree, knots =  Aknots, 
            intercept = intercept, ref=ref, log=getLog(object))
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("tp", "basis", "matrix")
  basis
}




######################################################################
EvaluateC0SBasis<-function(object,x,intercept=TRUE, ...) { 
# output : b_1(t) , b_i(t) - b_1(t) if i in 2:(degree), b_i(t) if i > degree
  stopifnot(is.numeric(x))
  dots<-list(...)
  nx <- names(x)
  nax <- is.na(x)
  if (nas <- any(nax)){ 
    x <- x[!nax]
  }

  Aknots<-object@knots
  degree<-object@degree
  ol <- x < Aknots[degree+1]
  or <- x > Aknots[length(Aknots)-degree]
  outside <- ( or | ol)

  
  if (any(outside)) {
    basis <- array(NA, dim= c(length(x), length(Aknots) - degree - 1L))
    if (any(inside <- !outside)){ 
      thebasis <- spline.des(knots=Aknots, x=x[inside], ord=degree+1)$design
      for( i in 2:degree){
      thebasis[,i]<-thebasis[,i] - thebasis[,1]
      }
      basis[inside, ] <- thebasis
    }
  }
  else {
    basis <- spline.des(knots=Aknots, x=x, ord=degree+1)$design
    for( i in 2:degree){
      basis[,i]<-basis[,i] - basis[,1]
    }
  }
  if (!intercept) {
    basis <- basis[, -1L, drop = FALSE]
  }
  if (object@log) {
    basis <- cbind(basis, log(x))
  }

#add na x
  if (nas) {
    n.col <- ncol(basis)
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = degree, knots =  Aknots, 
            intercept = intercept, log=getLog(object))
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("bs", "basis", "matrix")
  basis
}

# here, evaluate for tp-spline is EvaluateTPBasis
# thus, evaluate(tpspline, 0 , intercept = TRUE) == c(1, rep(0, nn)) 

setGeneric("evaluate",function(object, x,...)standardGeneric("evaluate"))
setMethod("evaluate",signature("BSplineBasis","numeric"),function(object, x, ...)EvaluateBBasis(object=object, x=x, ...))
setMethod("evaluate",signature("MSplineBasis","numeric"),function(object, x, ...)EvaluateMBasis(object=object, x=x, ...))
setMethod("evaluate",signature("TPSplineBasis","numeric"),function(object, x, ...)EvaluateTPBasis(object=object, x=x, ...))
setMethod("evaluate",signature("TPRSplineBasis","numeric"),function(object, x, ...)EvaluateTPBasis(object=object, x=x, ref=object@ref, ...))
setMethod("evaluate",signature("C0BSplineBasis","numeric"),function(object, x, ...)EvaluateC0SBasis(object=object, x=x, ...))

setGeneric("fevaluate",function(object, x,...)standardGeneric("fevaluate"))
setMethod("fevaluate",signature("BSplineBasis","numeric"),function(object, x, ...)EvaluateBBasis(object=object, x=x, ...))
setMethod("fevaluate",signature("MSplineBasis","numeric"),function(object, x, ...)FEvaluateMBasis(object=object, x=x, ...))
setMethod("fevaluate",signature("TPSplineBasis","numeric"),function(object, x, ...)FEvaluateTPBasis(object=object, x=x, ...))
setMethod("fevaluate",signature("TPRSplineBasis","numeric"),function(object, x, ...)EvaluateTPBasis(object=object, x=x, ref=object@ref, ...))
setMethod("fevaluate",signature("C0BSplineBasis","numeric"),function(object, x, ...)EvaluateC0SBasis(object=object, x=x, ...))
setMethod("fevaluate",signature("NULL","numeric"),function(object, x, ...) NULL)


######################################################################
#  evaluate linear combination of spline basis

# M-spline
EvaluateLCMBasis<-function(object, x, beta, intercept=TRUE, outer.ok=TRUE, ...) {
  stopifnot(is.numeric(x))
  dots<-list(...)
  M<-object@Matrices
  knots<-object@knots
  order<-object@degree+1
  cl <- .Call("eval_lc_spline_basis", as.double(knots), as.integer(order), M,
                 as.integer(intercept), as.double(x), as.double(beta), as.integer(outer.ok),
                 PACKAGE="flexrsurv")

  if (object@log) {
    return(cl + log(x) * beta[length(beta)])
  }
 else {
   return(cl)
 }

}

EvaluateLCTPBasis<-function(object, x, beta, intercept=TRUE, ref=NULL, outer.ok=TRUE, ...) { 
  stopifnot(is.numeric(x))

  min<-object@min
  max<-object@max
  allknots<-object@knots
  order<-object@degree+1
  knots <- unique(allknots)

  if( !is.null(ref)){
    x <- x - ref
    min <- min- ref
    max <- max - ref
    if(length(knots)>0){
      knot <- knots - ref
    }
  }
  
  replicates <- table(allknots)
  degrees <- object@degrees
  coef <- object@coef
  if(object@type == "increasing"){
    cl <- .Call("eval_lc_trunc_power_increasing_basis", as.double(knots), as.double(replicates),
                   as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
                   as.double(x), as.double(beta), as.integer(outer.ok),
                   PACKAGE="flexrsurv")
  }
  else {
  cl <- .Call("eval_lc_trunc_power_basis", as.double(knots), as.double(replicates),
                 as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
                 as.double(x), as.double(beta), as.integer(outer.ok),
                 PACKAGE="flexrsurv")
  }
  if (object@log) {
    return(cl + log(x) * beta[length(beta)])
  }
 else {
   return(cl)
 }
}





setGeneric("evaluatelc",function(object, x, beta,...)standardGeneric("evaluatelc"))
setMethod("evaluatelc",signature("MSplineBasis","numeric","numeric"),function(object, x, beta, ...)EvaluateLCMBasis(object=object, x=x, beta=beta, ...))
setMethod("evaluatelc",signature("TPSplineBasis","numeric","numeric"),function(object, x, beta, ...)EvaluateLCTPBasis(object=object, x=x, beta=beta, ...))


######################################################################
#  Method predic for SplineParam

# MSplineParam
# computes f(x) = sum_i beta_i b_i(x)
PredictMBasisBeta <- function(object=object, x=x, beta=beta, intercept=TRUE, outer.ok=TRUE, ...){
  if(intercept){
    predict(object * beta, x, intercept=TRUE, ...)
  }
  else {
    predict(object * c(0, beta), x,  intercept=FALSE, ...)
  }
}

# computes f(x) = sum_i B_i(x)
# assuming B_i = beta_i * beta_i
predict.MSplineBasis <- function(object=object, x=x, intercept=TRUE, outer.ok=TRUE, ...){
  stopifnot(is.numeric(x))

#  if(!is.null(beta)){
#    if(intercept){
#      object <- object * beta
#    }
#    else {
#      object <- object * c(0, beta)
#    }
#  }

  dots<-list(...)
  M<-object@Matrices
  knots<-object@knots
  order<-object@degree+1
  cl <- .Call("predict_spline_basis", as.double(knots), as.integer(order), M,
                 as.integer(intercept), as.double(x), as.integer(outer.ok),
                 PACKAGE="flexrsurv")

  if (object@log) {
    return(cl + log(x) * beta[length(beta)])
  }
 else {
   return(cl)
 }

}


# TPSplineParam
# computes f(x) = sum_i beta_i b_i(x)
PredictTPBasisBeta <- function(object=object, x=x, beta=beta, intercept=TRUE, ref=NULL, outer.ok=TRUE, ...){
  if(intercept){
    predict(object * beta, x, intercept=TRUE, ...)
  }
  else {
    predict(object * c(0, beta), x,  intercept=FALSE, ...)
  }
}

# computes f(x) = sum_i beta_i * B_i(x)
# if beta == NULL, assuming beta_i = 1
predict.TPSplineBasis <- function(object=object, x=x, beta=NULL, intercept=TRUE, ref=NULL, outer.ok=TRUE, ...){
  stopifnot(is.numeric(x))

  if(!is.null(beta)){
    if(intercept){
      object <- object * beta
    }
    else {
      object <- object * c(0, beta)
    }
  }
  min<-object@min
  max<-object@max
  allknots<-object@knots
  order<-object@degree+1
  knots <- unique(allknots)

  if( !is.null(ref)){
    x <- x - ref
    min <- min- ref
    max <- max - ref
    if(length(knots)>0){
      knot <- knots - ref
    }
  }
  
  replicates <- table(allknots)
  degrees <- object@degrees
  coef <- object@coef
  if(object@type == "increasing"){
    cl <- .Call("predict_trunc_power_increasing_basis", as.double(knots), as.double(replicates),
                   as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
                   as.double(x), as.integer(outer.ok),
                   PACKAGE="flexrsurv")
  }
  else {
  cl <- .Call("predict_trunc_power_basis", as.double(knots), as.double(replicates),
                 as.double(min), as.double(max), as.integer(order), as.double(coef), as.double(degrees), as.integer(intercept),
                 as.double(x), as.integer(outer.ok),
                 PACKAGE="flexrsurv")
  }
  if (object@log) {
    return(cl + log(x) * beta[length(beta)])
  }
 else {
   return(cl)
 }

}


setGeneric("predictSpline",function(object, x, beta,...)standardGeneric("predictSpline"))
setMethod("predictSpline",signature(object="MSplineBasis",x="numeric", beta="missing"),function(object, x, ...)predict.MSplineBasis(object=object, x=x,  ...))
setMethod("predictSpline",signature(object="TPSplineBasis",x="numeric", beta="missing"),function(object, x, ...)predict.TPSplineBasis(object=object, x=x,  ...))

#setMethod("predict",signature("MSplineBasis","numeric","numeric"),function(object, x, beta, ...)PredictMBasisBeta(object=object, x=x, beta=beta, ...))
#setMethod("predict",signature("TPSplineBasis","numeric","numeric"),function(object, x, beta, ...)PredictTPBasisBeta(object=object, x=x, beta=beta, ...))


######################################################################
#  integrate 

# define parameters for integrated Spline Basis
integrate.TPSplineBasis<-function(object){

  min<-object@min
  max<-object@max
  allknots<-object@knots
  order<-object@degree+1
  knots <- unique(allknots)
  replicates <- table(allknots)



  idegrees <- object@degrees + 1
  coef <- object@coef
  nbases <- length(coef)

  if(object@type == "standard"){
    icoef <- c(0, coef[1:order] / (1:order), coef[(order+1):nbases] / idegrees)
  }
  else {
    icoef <- c(0, coef[1:order] / (1:order), coef[(order+1):nbases] / ifelse(allknots<0, -  idegrees,   idegrees))
  }
  new("TPSplineBasis", knots=object@knots, min=object@min, max=object@max, 
      degree=object@degree+1, nbases=nbases, coef=icoef, degrees=idegrees, log=FALSE, type=object@type)
}

# define parameters for integrated Spline Basis
integrate.MSplineBasis<-function(object){
  SB <- orthogonalsplinebasis::integrate(object@SplineBasis)
  new("MSplineBasis", knots=SB@knots, min=object@min, max=object@max,
      degree=as.integer(SB@order-1), nbases=as.integer(dim(SB@Matrices)[2]), Matrices=SB@Matrices, SplineBasis=SB, log=FALSE)
}

# compute values of integrated basis
IntegrateMBasis<-function(object,x,intercept=TRUE, xname=NULL, outer.ok=TRUE, ...) {
  stopifnot(is.numeric(x))
  if (object@log) {
    stop("no method 'integrate' for MSplineBasis with additional log basis" )
  }
  else {
    evaluate(integrate(object), x=x, intercept=intercept, outer.ok=outer.ok, ...)
  }
}



# compute values of integrated basis
IntegrateMBasisOld<-function(object,x,intercept=TRUE, xname=NULL, outer.ok=TRUE, ...) {
  stopifnot(is.numeric(x))
  nax <- is.na(x)
  if (nas <- any(nax)){ 
    x <- x[!nax]
  }

  Aknots<-object@knots
  degree<-object@degree
  nbases<-object@nbases
  ol <- x < Aknots[degree+1]
  or <- x > Aknots[length(Aknots)-degree]
  outside <- ( or | ol)
  
  if (any(outside)) {
    if(outer.ok) {
      basis <- array(0, dim= c(length(x), nbases))
    }
    else {
      basis <- array(NA, dim= c(length(x), nbases))
    }
    if (any(inside <- !outside)){ 
      basis[inside, ] <- orthogonalsplinebasis::evaluate(orthogonalsplinebasis::integrate(object@SplineBasis), x[inside] )
    }
  }
  else {
    basis <- orthogonalsplinebasis::evaluate(orthogonalsplinebasis::integrate(object@SplineBasis), x )
  }
  if (!intercept) {
    basis <- basis[, -1L, drop = FALSE]
  }

  if (object@log) {
    basis <- cbind(basis, 1/x)
  }

#add na x
  if (nas) {
    nabasis <- matrix(NA, length(nax), nbases)
    nabasis[!nax, ] <- basis
    nabasis
  }
  else {
    basis
  }


}

# compute values of integrated basis
# truncated powezr basis
integrateTPBasis<-function(object,x,intercept=TRUE, ref=NULL, xmin=-Inf, xmax=+Inf, outer.ok=TRUE, ...) { 
  stopifnot(is.numeric(x))
  if (object@log) {
    stop("no method 'integrate' for MSplineBasis with additional log basis" )
  }
  else {
    evaluate(integrate(object), x=x, intercept=intercept, ref=ref, xmin=xmin, xmax=xmax, outer.ok=outer.ok, ...)
  }
}

FIntegrateTPBasis0<-function(object,x,intercept=TRUE, ref=NULL, xmin=-Inf, xmax=+Inf, ...) { 
# similar to FEvaluateTPBasis but if knots<ref, non 0 if x<0
# thus, active bases around ref are the full monomials (x-ref)^k
  stopifnot(is.numeric(x))
  nax <- is.na(x)
  if (nas <- any(nax)){ 
    x <- x[!nax]
  }

  Aknots<-object@knots
  degree<-object@degree
  degreep1 <- degree+1
  if (xmin != -Inf){
    ol <- x < xmin
    if (xmax != +Inf){
      or <- x > xmax
      outside <- ( or | ol)
    }
    else {
      outside <- ol
    }
  }
  else {
    if (xmax != +Inf){
      outside <- x > xmax
    }
    else {
      outside <- FALSE
    }
  }

  Aknotref <- Aknots 
  if( !is.null(ref)){
    x <- x - ref
    if(length(Aknots)>0){
      Aknotref <- Aknots - ref
    }
  }

  
  if (intercept) {
    nb <- degree+ 1 + length(Aknots) 
    last <- degree+1
  }
  else{
    nb <- degree  + length(Aknots) 
    last <- degree
  }    
  basis <- matrix(1, nrow=length(x), ncol=nb)



  if (any(!outside)){ 
    if (intercept) {
      for(i in 1:(degree+1)){
        basis[,i]<-x^i/i
      }
    }
    else{
      for(i in 2:(degree+1)){
        basis[,i+1]<-x^i/i
      }
    }
    if(length(Aknots)>0){
      for (k in 1:(length(Aknots))){
        if(Aknots[k]>0){
          basis[,last+k] <- ifelse(x>Aknotref[k], (x-Aknotref[k])^degreep1/degreep1, 0)
        } else {
          basis[,last+k] <- ifelse(x>Aknotref[k], 0, (x-Aknotref[k])^degreep1/degreep1)
        }
      }
    }
    if (object@log) {
      basis <- cbind(basis, log(x))
    }
  }    
  if (any(outside)) {
    basis[outside,]<-rep(NA, dim(basis)[2])
    if (object@log) {
      basis <- cbind(basis, rep(NA, dim(basis)[1]))
    }
  }

#add na x
  if (nas) {
    nabasis <- matrix(NA, nrow=length(nax), ncol=nb)
    nabasis[!nax, ] <- basis
    nabasis
  }
  else {
    basis
  }
}




FIntegrateTPBasis<-function(object,x,intercept=TRUE, ref=NULL, xmin=-Inf, xmax=+Inf, ...) { 
  stopifnot(is.numeric(x))
  nax <- is.na(x)
  if (nas <- any(nax)){ 
    x <- x[!nax]
  }
  Aknots<-object@knots
  degree<-object@degree
  degreep1 <- degree+1
  if (xmin != -Inf){
    ol <- x < xmin
    if (xmax != +Inf){
      or <- x > xmax
      outside <- ( or | ol)
    }
    else {
      outside <- ol
    }
  }
  else {
    if (xmax != +Inf){
      outside <- x > xmax
    }
    else {
      outside <- FALSE
    }
  }

  Aknotref <- Aknots 
  if( !is.null(ref)){
    x <- x - ref
    if(length(Aknots)>0){
      Aknotref <- Aknots - ref
    }
  }

  
  if (intercept) {
    nb <- degree+ 1 + length(Aknots) 
    last <- degree+1
  }
  else{
    nb <- degree  + length(Aknots) 
    last <- degree
  }    
  basis <- matrix(1, nrow=length(x), ncol=nb)



  if (any(!outside)){ 
    if (intercept) {
      for(i in 1:(degree+1)){
        basis[,i]<-x^i/i
      }
    }
    else{
      for(i in 2:(degree+1)){
        basis[,i+1]<-x^i/i
      }
    }
    if(length(Aknots)>0){
      for (k in 1:(length(Aknots))){
        basis[,last+k] <- ifelse(x>Aknotref[k], (x-Aknotref[k])^degreep1/degreep1, 0)
      }
    }
    if (object@log) {
      basis <- cbind(basis, 1/x)
    }
  }    
  if (any(outside)) {
    basis[outside,]<-rep(NA, dim(basis)[2])
    if (object@log) {
      basis <- cbind(basis, rep(NA, dim(basis)[1]))
    }
  }

#add na x
  if (nas) {
    nabasis <- matrix(NA, nrow=length(nax), ncol=nb)
    nabasis[!nax, ] <- basis
    nabasis
  }
  else {
    basis
  }
}



setGeneric("integrate",function(object, x,...)standardGeneric("integrate"))
setMethod("integrate",signature("MSplineBasis", "missing"),integrate.MSplineBasis)
setMethod("integrate",signature("TPSplineBasis", "missing"),integrate.TPSplineBasis)

setMethod("integrate",signature("MSplineBasis","numeric"),function(object, x, ...)IntegrateMBasis(object=object, x=x, ...))
setMethod("integrate",signature("TPSplineBasis","numeric"),function(object, x, ...)FIntegrateTPBasis(object=object, x=x, ...))
#setMethod("integrate",signature("TPRSplineBasis","numeric"),function(object, x, ...)IntegrateTPBasis(object=object, x=x, ref=object@ref, ...))
#setMethod("integrate",signature("C0BSplineBasis","numeric"),function(object, x, ...)IntegrateC0SBasis(object=object, x=x, ...))



InitCoefSBasis<-function(object,ncol, init=1, intercept=TRUE, xname=NULL, ...) { 
  # knots are all knots (first and last replicated
  # output matrix with all "init" 
  stopifnot(is.integer(ncol))
  nb <- object@nbases + 1 - intercept + object@log
  matrix(init, ncol=ncol , nrow=nb)

}

InitCoefMBasis<-function(object,ncol, init=1, intercept=TRUE, xname=NULL, ...) { 
  # output matrix with all "init" 
  # knots are all knots (first and last replicated
  stopifnot(is.integer(ncol))
  nb <- object@nbases + 1 - intercept + object@log

  matrix(init, ncol=ncol , nrow=nb)

}


InitCoefTPBasis<-function(object,ncol, init=1, intercept=TRUE, xname=NULL, ...) { 
  stopifnot(is.integer(ncol))
  # knots are ordered  interior knots
  # init vectoe to (1, 0, 0, ...)
  nb <- object@nbases + 1 - intercept + object@log

  vv <- rep(0,nb)
  vv[1] <- 1

  matrix(vv, ncol=ncol , nrow=nb)

}

setGeneric("initcoef",function(object, ncol,...)standardGeneric("initcoef"))
setMethod("initcoef",
          signature("BSplineBasis","integer"),
          function(object, ncol, ...)InitCoefSBasis(object=object, ncol=ncol, ...)
          )
setMethod("initcoef",
          signature("MSplineBasis","integer"),
          function(object, ncol, ...)InitCoefMBasis(object=object, ncol=ncol, ...)
          )
setMethod("initcoef",
          signature("TPSplineBasis","integer"),
          function(object, ncol, ...)InitCoefTPBasis(object=object, ncol=ncol, ...)
          )


# idem but for constraints parameters
# for B-spline, all coefs are 1
InitCoefCSBasis<-function(object,ncol,intercept=TRUE, xname=NULL, ...) { 
  # knots are all knots (first and last replicated
  stopifnot(is.integer(ncol))
  nb <- length(object@knots) - object@degree - 3 + intercept  + object@log

  matrix(1, ncol=ncol , nrow=nb)
}



# for natural spline, all coefs are null for non intercvept term
InitCoefCTPBasis<-function(object,ncol,intercept=TRUE, xname=NULL, ...) {
  # knots are knot_min and interior knots
  stopifnot(is.integer(ncol))
  nb <- object@degree+ length(object@knots) -3 + intercept + object@log

  vv <- rep(0,nb)
  matrix(vv, ncol=ncol , nrow=nb)

}

setGeneric("initcoefC",function(object, ncol, ...)standardGeneric("initcoefC"))
setMethod("initcoefC",
          signature("BSplineBasis","numeric"),
          function(object, ncol, ...)InitCoefCSBasis(object=object, ncol=ncol, ...)
          )
  
setMethod("initcoefC",
          signature("TPSplineBasis","numeric"),
          function(object, ncol, ...)InitCoefCTPBasis(object=object, ncol=ncol, ...))


ExpandMCoefSBasis0<-function(object,ncol,coef, intercept=TRUE, xname=NULL, ...) { 
# add a first row of nrow - sum(matcoef) 
  stopifnot(is.integer(ncol))
  nb <- dim(coef)[1]
  
  rbind(intercept - rep(1, nb) %*% coef , coef)
}

ExpandVCoefSBasis0<-function(object,ncol,coef, intercept=TRUE, xname=NULL, ...) { 
# add a first row of nrow - sum(matcoef) 
  stopifnot(is.integer(ncol))
  if(!is.matrix(coef)){
    coef <- matrix(coef, ncol=ncol)
  }
  
  nb <- dim(coef)[1]+1
  
  rbind(intercept - rep(1, nb-1) %*% coef , coef)
}



######################################################################
# operators


# SplineBasis
Prod.S.n <- function(e1, e2) {
    d <- dim(e1@Matrices)
    le <- length(e2)
    if(le == 1 ) { # matching dim
        e1@Matrices <- e1@Matrices * e2
        e1
    } else if (le >= d[2]){
# on multiplie chaque Matrices[,,i] par diag(e2)
      for(i in 1:d[3]){
        e1@Matrices[,,i] <- e1@Matrices[,,i] %*%diag(e2[1:d[2]])
      }
      e1
    }
    else stop ("length of args does not")
}
setMethod("*", signature(e1 = "SplineBasis", e2 = "numeric"), Prod.S.n)
Prod.n.S <- function(e1, e2) {
  e2 * e1
}
setMethod("*", signature(e1 ="numeric" , e2 = "SplineBasis"), Prod.n.S)


Sum.S.S <- function(e1, e2) {
  if(e1@knots== e2@knots && e1@order == e2@order){
    e1@Matrices <- e1@Matrices + e2@Matrices
    e1
  }
  else stop("args cannot be added")
}
setMethod("+", signature(e1 = "SplineBasis", e2 = "SplineBasis"), Sum.S.S)


Sum.S.n <- function(e1, e2) {
 e1 + e2 * SplineBasis(knots=e1@knots, order=e1@order, keep.duplicates=TRUE)
}
setMethod("+", signature(e1 = "SplineBasis", e2 = "numeric"), Sum.S.n)

Sum.n.S <- function(e1, e2) {
  e2 + e1
}
setMethod("+", signature(e1 ="numeric" , e2 = "SplineBasis"), Sum.n.S)



Dif.S.S <- function(e1, e2) {
  if(e1@knots== e2@knots && e1@order == e2@order){
    e1@Matrices <- e1@Matrices - e2@Matrices
    e1
  }
  else stop("args cannot be added")
}
setMethod("-", signature(e1 = "SplineBasis", e2 = "SplineBasis"), Dif.S.S)

Dif.S.n <- function(e1, e2) {
  e1 + ((-1) * e2)
}
setMethod("-", signature(e1 = "SplineBasis", e2 = "numeric"), Dif.S.n)

Dif.n.S <- function(e1, e2) {
  ( e2 * (-1)) + e1
}
setMethod("-", signature(e1 = "numeric", e2 = "SplineBasis"), Dif.n.S)

######################################################################
# MSplineBasis
Prod.MS.n <- function(e1, e2) {
    e1@SplineBasis <- e1@SplineBasis * e2
    e1@Matrices <- e1@SplineBasis@Matrices
    e1
}
setMethod("*", signature(e1 = "MSplineBasis", e2 = "numeric"), Prod.MS.n)

Prod.n.MS <- function(e1, e2) {
  e2 * e1
}
setMethod("*", signature(e1 ="numeric" , e2 = "MSplineBasis"), Prod.n.MS)


Sum.MS.MS <- function(e1, e2) {
    e1@SplineBasis <- e1@SplineBasis + e2@SplineBasis 
    e1@Matrices <- e1@SplineBasis@Matrices
    e1
}
setMethod("+", signature(e1 = "MSplineBasis", e2 = "MSplineBasis"), Sum.MS.MS)


Sum.MS.n <- function(e1, e2) {
    e1@SplineBasis <- e1@SplineBasis + e2 
    e1@Matrices <- e1@SplineBasis@Matrices
    e1
}
setMethod("+", signature(e1 = "MSplineBasis", e2 = "numeric"), Sum.MS.n)

Sum.n.MS <- function(e1, e2) {
  e2 + e1
}
setMethod("+", signature(e1 ="numeric" , e2 = "MSplineBasis"), Sum.n.MS)



Dif.MS.MS <- function(e1, e2) {
    e1@SplineBasis <- e1@SplineBasis - e2@SplineBasis 
    e1@Matrices <- e1@SplineBasis@Matrices
    e1
}
setMethod("-", signature(e1 = "MSplineBasis", e2 = "MSplineBasis"), Dif.MS.MS)

Dif.MS.n <- function(e1, e2) {
  e1 + ((-1) * e2)
}
setMethod("-", signature(e1 = "MSplineBasis", e2 = "numeric"), Dif.MS.n)

Dif.n.MS <- function(e1, e2) {
  ( e2 * (-1)) + e1
}
setMethod("-", signature(e1 = "numeric", e2 = "MSplineBasis"), Dif.n.MS)



######################################################################
# TPSplineBasis
Prod.TPS.n <- function(e1, e2) {
    d <- length(e1@coef)
    le <- length(e2)
    if(le == 1 ) { # matching dim
        e1@coef <- e1@coef * e2
        e1
    } else if (le >= d){
      e1@coef <- (e1@coef * e2[1:d])
      e1
    }
    else stop ("length of args does not")
}
setMethod("*", signature(e1 = "TPSplineBasis", e2 = "numeric"), Prod.TPS.n)

Prod.n.TPS <- function(e1, e2) {
  e2 * e1
}
setMethod("*", signature(e1 ="numeric" , e2 = "TPSplineBasis"), Prod.n.TPS)


Sum.TPS.TPS <- function(e1, e2) {
  if(e1@knots== e2@knots && e1@degree == e2@degree && e1@min == e2@min && e1@max == e2@max && e1@degrees == e2@degrees && e1@type == e2@type){
    e1@coef <- e1@coef + e1@coef 
    e1
  }
  else stop("args cannot be added")
}
setMethod("+", signature(e1 = "TPSplineBasis", e2 = "TPSplineBasis"), Sum.TPS.TPS)


Sum.TPS.n <- function(e1, e2) {
    d <- length(e1@coef)
    le <- length(e2)
    if(le == 1 ) { # matching dim
        e1@coef[1] <- e1@coef[1]  +  e2
        e1
    } else if (le >= d){
      e1@coef <- e1@coef + e2[1:d]
      e1
    }
    else stop ("length of args does not match")
}
setMethod("+", signature(e1 = "TPSplineBasis", e2 = "numeric"), Sum.TPS.n)


Sum.n.TPS <- function(e1, e2) {
  e2 + e1
}
setMethod("+", signature(e1 ="numeric" , e2 = "TPSplineBasis"), Sum.n.TPS)


Dif.S.S <- function(e1, e2) {
  if(e1@knots== e2@knots && e1@degree == e2@degree && e1@min == e2@min && e1@max == e2@max && e1@degrees == e2@degrees && e1@type == e2@type){
    e1@coef <- e1@coef - e1@coef 
    e1
  }
  else stop("args cannot be added")
}
setMethod("-", signature(e1 = "SplineBasis", e2 = "SplineBasis"), Dif.S.S)

Dif.TPS.n <- function(e1, e2) {
  e1 + ((-1) * e2)
}
setMethod("-", signature(e1 = "TPSplineBasis", e2 = "numeric"), Dif.TPS.n)

Dif.n.TPS <- function(e1, e2) {
  ( e2 * (-1)) + e1
}
setMethod("-", signature(e1 = "numeric", e2 = "TPSplineBasis"), Dif.n.TPS)






