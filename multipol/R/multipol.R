if(FALSE){
"as.character" <- function(x, ...) {UseMethod("as.character")}
"as.character.default" <- base::as.character
}

"as.character.multipol" <- function(x, ..., xyz = getOption("xyz"), varname = getOption("varname")){

  if(is.null(xyz)){
    xyz <- TRUE
  }
  if(is.null(varname)){
    varname <- "x"
  }
  
  jj <- x != 0
  if(!any(jj)){return("0")}
  index <- which(jj,arr.ind=TRUE)-1
  val <- x[jj]
  negative <- val<0

  dx <- dim(x)
  sad <- seq_along(dx)
  
  if(xyz &  (length(dx) < 4)){
    u <- letters[23+sad]
  } else {
    u <- paste(varname, sad, sep="")
  }

  f <- function(ee){
    paste(paste(u,ee,sep="^"),collapse="*")
  }
  out <- paste(val,apply(index,1,f) ,sep="*")
  if(xyz & (length(dim(x)) < 4)){
    out <- gsub("\\*x\\^0","", out)
    out <- gsub("\\*y\\^0","", out)
    out <- gsub("\\*z\\^0","", out)
  } else {
    rgx <- paste("\\*",varname,"[0-9]\\^0",sep="")
    out <- gsub(rgx,"",out)
  }
  return(paste(out,collapse=" + "))
}

"polyprod" <- function(m1,m2,overlap=0){
  m1 <- unclass(m1)
  m2 <- unclass(m2)

  d1 <- dim(m1)
  d2 <- dim(m2)
  
  dim(m1) <- c(d1,rep(1,length(d2)-overlap)   )
  dim(m2) <- c(   rep(1,length(d1)-overlap),d2)

  as.multipol(m1)*as.multipol(m2)
}

"product" <- function(x){
  out <- array(0,x+1)
  out[length(out)] <- 1
  as.multipol(out)
}

"lone" <- function(d, x){
  stopifnot(all(d >= x))
  index <- rep(1,d)
  index[x] <- 2
  out <- array(0,index)
  out[length(out)] <- 1
  return(as.multipol(out))
}
  
"ooom" <- function(n, a, maxorder = NULL){
  nze <- which(a>0, arr.ind=TRUE)-1
  if(nrow(nze)==1){
    dims <- 1+(dim(a)-1)*n
    if(!is.null(maxorder)){
      dims <- pmin(dims,maxorder+1)
    } 
    out <- array(0,dims)
    index <- 1+kronecker(nze,0:n)
    if(!is.null(maxorder)){
      wanted <- apply(index,1,function(u){ !any(u>(maxorder+1))})
      out[index[wanted,]] <- 1
    } else {
      out[index] <- 1
    }
      
    return(as.multipol(out))
  }
  out <- 1
  temp <- 1
  for(i in seq_len(n)){
    temp <- taylor(temp*a,maxorder=maxorder)
    out <- out + temp
  }
  return(out)
}

"uni" <- function(d){product(rep(1,d))}

"single" <- function(d, e=d, power=1){
  stopifnot(e <= d)
  jj <- rep(1,d)
  jj[e] <- power+1
  out <- array(0,jj)
  out[power+1] <- 1
  as.multipol(out)
}

"linear" <- function(x,power=1){
  lx <- length(x)
  out <- as.multipol(array(0,rep(power+1,lx)))
  out[diag(lx)*power] <- x
  out
}

"ones" <- function(d,power=1){
  linear(rep(1,d),power=power)
}

"constant" <- function(d){
  as.multipol(array(1,rep(1,d)))
}

"zero" <- function(d){0*constant(d)}

"homog" <- function(d,n=d,value=1){
  jj <- apply(do.call("abind",c(lapply(seq_len(d),function(i){slice.index(array(0L,rep(n+1,d)),i)}) ,along=d+1)),seq_len(d),sum)-d
  out <- jj*0
  out[jj==n] <- value
  as.multipol(out)
}
  
"trim" <- function(a){
  a <- unclass(a)
  jj <- a != 0
  if(any(jj)){
    jj <- which(jj,arr.ind=TRUE)
    return(as.multipol(do.call("[",c(list(a), sapply(apply(jj,2,max),seq_len,simplify=FALSE),drop=FALSE))))
  } else {
    return(as.multipol(array(0,rep(1,length(dim(a))))))
  }
}

"is.zero" <- function(a, allow.untrimmed = TRUE){
  is.constant(a,allow.untrimmed=allow.untrimmed) & as.array(a)[1]==0
}

"is.constant" <- function(a,allow.untrimmed = TRUE){
  if(allow.untrimmed){
    return(Recall(trim(a), allow.untrimmed = FALSE))
  } else {
    return(all(dim(a)==1))
  }
}

".handleTheOffset" <- function(mc, dim, offset, dn)
{
  for (i in seq(along=dim)) {
    ii <- mc[[2+i]]

    if (missing(ii)) next

    if (is.symbol(ii) || is.call(ii))
      ii <- eval.parent(ii, 3)

    if (is.numeric(ii)) {

      if (!dn || all(ii>=0))
        ii <- ifelse(ii>=offset[i], ii - offset[i] + 1, dim[i]+1)
      else {
        if (all(ii <= -offset[i]))
          ii <- ii + offset[i] - 1
        else stop("subscript out of bounds")
      }

      mc[[2+i]] <- ii
    }
  }
  mc
}

"[.multipol" <- function(x, ...)
{
  mc <- match.call()
  k <- length(mc)
  offset <- rep(0,length(dim(x)))
  dn <- FALSE
  dim <- dim(x)

  if( k==3 ){
    if(mc[[3]] == ""){
      return(as.array(x))
    } 
    args <- list(...)
    index <- args[[1]]
    if(is.logical(index)){
      return(as.array(x)[index])
    }
    if(is.matrix(index)){
      return(as.array(x)[1+sweep(index,2,offset)])
    }
  }

  if (k < 2+length(dim))
    stop("incorrect number of dimensions")

  mc <- .handleTheOffset(mc, dim, offset, dn)
  mc[[1]] <- as.name("[")
  mc[[2]] <- as.name("x")
  x <- as.array(x)
  eval(mc)
}

"[<-.multipol" <- function(x, ..., value)
{
  mc <- match.call()
  k <- length(mc)
  offset <- rep(0,length(dim(x)))
  dn <- FALSE
  dim <- dim(x)

  if (k==4){
    if (mc[[3]] == ""){
      return(as.multipol(array(value, dim)))
    }
    args <- list(...)
    index <- args[[1]]
    att <- attributes(x)
    if(is.logical(index)){
      x <- as.array(x)
      x[index] <- value
      attributes(x) <- att
      return(x)
    }
    if(is.matrix(index)){
      x <- as.array(x)
      x[1+sweep(index,2,offset)] <- value
      attributes(x) <- att
      return(x)
    }
  }
  
  if (k < 3+length(dim))
    stop("incorrect number of dimensions")

  mc <- .handleTheOffset(mc, dim, offset, dn)
  mc[[1]] <- as.name("[<-")
  mc[[2]] <- as.name("x")
  x <- as.array(x)
  robj <- eval(mc)
  as.multipol(array(robj, dim))
}

".append2" <- 
function (x, value, after = length(x)) 
{
  out <- c(x,seq_along(after))
  out[+1+after] <- value
  out[-1-after] <- x
out
}

"multipol" <- function(x){
  stopifnot(is.array(x))
  class(x) <- c("multipol")
  return(x)
}

"print.multipol" <- function(x, ...){
  if(isTRUE(getOption("showchars"))){
    x <- noquote(as.character(x))
  } else {
    x <- do_dimnames(as.array(x))
  }
  return(invisible(print(x)))
}

"is.multipol" <- function(x){
  inherits(x,"multipol")
}

"as.multipol" <- function(x){
  if(is.character(x)) {
    out <- unlist(strsplit(x,"\\+|\\-"))
    stop("coercion from character strings not yet implemented")
  }  else if(is.vector(x)){
    return(multipol(as.array(x)))
  } else if(is.array(x)) {
    return(multipol(x))
  } else {
    stop("supply an array or a character string")
  }
}

if(getRversion() < "2.8.0"){
  "as.array" <- function(x, ...){UseMethod("as.array")}
  "as.array.default" <- function(x, ...){base::as.array(x)}
}

"as.array.multipol" <- function(x, ...){
  x <- unclass(x)
  NextMethod(x, ...)
}

".multipol.prod.multipol" <- function(... , trim=TRUE, maxorder=NULL){
  args <- list(...)
  if (length(args) == 1) {
    return(args[[1]])
  }
  if (length(args) > 2) {
    jj <-  do.call("Recall", c(args[-1]       ,          trim=trim,maxorder=maxorder))
    return(do.call("Recall", c(list(args[[1]]), list(jj),trim=trim,maxorder=maxorder)))
  }
  
  a <- args[[1]]
  b <- args[[2]]

  stopifnot(length(dim(a)) == length(dim(b)))

  if(trim){
    a <- trim(a)
    b <- trim(b)
  }
  
  a <- as.array(a)
  b <- as.array(b)

  if(is.null(maxorder)){  
    outDims <- dim(a)+dim(b)-1
  } else {
    stopifnot(length(maxorder) <= length(dim(a)))
    outDims <- pmin(dim(a)+dim(b)-1 , maxorder+1)
  }

  index <- as.matrix(expand.grid(lapply(outDims,seq_len)))
  
  f <- function(u){
    jj.a <- as.matrix(expand.grid(lapply(u,function(i){0:(i-1)})))
    jj.b <- -sweep(jj.a,2,u)-1L  # jj.a + jj.b "==" u
    
    wanted <-
    apply(jj.a,1,function(x){all(x < dim(a))}) &
    apply(jj.b,1,function(x){all(x < dim(b))}) &
    apply(jj.a,1,function(x){all(x >= 0)})        &
    apply(jj.b,1,function(x){all(x >= 0)})
    
    sum(a[1+jj.a[wanted,,drop=FALSE]] * b[1+jj.b[wanted,,drop=FALSE]])
  }
  out <- apply(index,1,f)
  dim(out) <- outDims
  out <- as.multipol(out)
  if(trim){
    out <- trim(out)
  }
  return(out)
} 

".multipol.prod.scalar" <- function(a,b) {
  as.multipol(as.array(a)*b)
}

".multipol.plus.multipol" <- function (..., trim=TRUE, maxorder=NULL) 
{
    args <- list(...)
  
    if (length(args) == 1) {
        return(args[[1]])
    }
    if (length(args) > 2) { 
        jj <- do.call("Recall", c(args[-1],trim=trim,maxorder=maxorder))
        return(do.call("Recall", c(list(args[[1]]), list(jj),trim=trim,maxorder=maxorder)))
    }
    a <- args[[1]]
    b <- args[[2]]
    if(trim){
      a <- trim(a)
      b <- trim(b)
    }
    dima <- dim(a)
    dimb <- dim(b)
    
    stopifnot(length(dima) == length(dimb))
    dd <- pmax(dima,dimb)
    out <- array(0, dd)
                   
    out <- do.call("[<-", c(list(out), lapply(dima, seq_len), list(a))) +
           do.call("[<-", c(list(out), lapply(dimb, seq_len), list(b)))

    out <- taylor(out,maxorder=maxorder)
    if(trim){
      out <- trim(out)
    }
    return(out)
  }

".multipol.neg" <- function(a, trim=TRUE, maxorder=NULL){
  if(trim){
    a <- trim(a)
  }
  return(taylor(-as.array(a),maxorder=maxorder))
}

".multipol.plus.scalar" <- function(a, b, trim=TRUE, maxorder=NULL){
  a <- as.array(a)
  a[1] <- a[1] + b
  a <- as.multipol(a)
  if(trim){
    a <- trim(a)
  }
  return(taylor(a,maxorder=maxorder))
}

".multipol.power.scalar"   <- function(a, n, trim=TRUE, maxorder=NULL){
  stopifnot(n==round(n))
  if(n==0){
    return(as.multipol(array(1,rep(1,length(dim(a))))))
  } else if (n==1) {
    if(trim){
      return(trim(a))
    } else {
      return(a)
    }
  } else {
    return(mprod(a,Recall(a,n-1,trim=trim,maxorder=maxorder),trim=trim,maxorder=maxorder))
  }
}

"mplus" <- .multipol.plus.multipol
"mprod" <- .multipol.prod.multipol
"mneg"  <- .multipol.neg
"mps"   <- .multipol.plus.scalar
"mpow"  <- .multipol.power.scalar

"taylor" <- function(a, maxorder=NULL){
  if(!is.null(maxorder)){
    a <- as.array(a)
    jj <- pmin(dim(a),maxorder+1)
    a <- do.call("[",c(list(a), lapply(jj,seq_len)))
  }
  return(as.multipol(a))
}
 
".unimplemented" <- function(e1=NULL,e2=NULL){
  stop("Operation reasonable but not yet implemented.  Try emailing me")  
}

".multipol.inv"            <- .unimplemented
".multipol.power.multipol" <- .unimplemented
".scalar.power.multipol"   <- .unimplemented

"do_dimnames" <- function(a, include.square.brackets=getOption("isb"), varname=getOption("varname"),xyz=getOption("xyz")){
  if(is.null(include.square.brackets)){
    include.square.brackets <- FALSE
  }
  if(is.null(varname)){
    varname <- "x"
  }
  if(is.null(xyz)){
    xyz <- TRUE
  }
  d <- dim(a)
  if ( (length(d) < 4) & xyz){
    f <- function(i){
      paste(letters[23+i],"^",seq_len(d[i])-1,sep="")}
  } else if(include.square.brackets){
    f <- function(i){
      paste("[",varname,i,"]^",seq_len(d[i])-1,sep="")}
  } else {
    f <- function(i){
      paste("",varname,i,"^",seq_len(d[i])-1,sep="")}
  }
  dimnames(a) <- sapply(seq_along(d),f,simplify=FALSE)
  class(a) <- "array"
  return(a)
}

"Ops.multipol" <-
  function (e1, e2 = NULL) 
{
  f <- function(...){stop("odd---neither argument has class multipol?")}
  unary <- nargs() == 1
  lclass <- inherits(e1,"multipol")
  rclass <- !unary && inherits(e2,"multipol")
  
  if(unary){
    if (.Generic == "+") {
      return(trim(e1))
    } else if (.Generic == "-") {
      return(.multipol.neg(e1))
    } else {
      stop("Unary operator '", .Generic, "' is not implemented for multipols")
    }
  }

  if (!is.element(.Generic,
                  c("+", "-", "*", "/", "^","==" , ">" , "<" , ">=" , "<=" , "!=" )
                  )){
    stop("operator '", .Generic, "' is not implemented for multipols")
  }
  
  if (.Generic == "*") {
    if (lclass && rclass) {
      return(.multipol.prod.multipol(e1, e2))
    } else if (lclass) {
      return(.multipol.prod.scalar  (e1, e2))
    } else if (rclass) {
      return(.multipol.prod.scalar  (e2, e1))
    } else {
      f()
    }
  } else if (.Generic == "/") {
    if (lclass && rclass) {
      return(.multipol.prod.multipol(e1, .multipol.inv(e2)))
    } else if (lclass) {
      return(.multipol.prod.scalar(e1, 1/e2))
    } else if (rclass) {
      return(.multipol.prod.scalar(.multipol.inv(e2), e1))
    } else {
      f()
    }
  } else if (.Generic == "^") {
    if(lclass && rclass){
      return(.multipol.power.multipol(e1,e2))
    } else if (lclass) {
      return(.multipol.power.scalar(e1,e2))
    } else if (rclass){
      return(.scalar.power.multipol(e1,e2))
    } else {
      f()
    }
  } else if (.Generic == "+") {
    if(lclass && rclass){
      return(.multipol.plus.multipol(e1,e2))
    } else if (lclass){
      return(.multipol.plus.scalar(e1,e2))
    } else if (rclass){
      return(.multipol.plus.scalar(e2,e1))
    } else {
      f()
    }
  } else if (.Generic == "-") {
    if(lclass && rclass){
      return(.multipol.plus.multipol(e1, .multipol.neg(e2)))
    } else if (lclass){
      return(.multipol.plus.scalar(e1, -e2))
    } else if (rclass){
      return(.multipol.plus.scalar(.multipol.neg(e2),e1))
    } else {
      f()
    }
  } else {
    return(NextMethod(.Generic))
  }
}

"deriv.multipol" <- function(expr, i, derivative=1, ...){
  a <- as.array(expr)
  if(i>length(dim(a))){
    stop("second argument exceeds arity of polynomial")
  }
  if(dim(a)[i]==1){return(as.multipol(array(0,rep(1,length(dim(a))))))}
  stopifnot(is.numeric(derivative))
  if(derivative != round(derivative)){
    stop("derivative not an integer")
  }
  if(derivative==0){
    return(a)
  } else if (derivative > 1){
    return(Recall(Recall(a,i,1),i,derivative-1))
  }
  
  jj <- rep(list(TRUE),length(dim(a)))
  jj[[i]] <- -1
  a <- do.call("[",c(list(a),jj,drop=FALSE))
  n <- seq_len(dim(a)[i])
  return(as.multipol(sweep(a,i,n,FUN="*")))
}

"put" <- function(a, i, value, keep=TRUE){
  if(i>length(dim(a))){
    stop("substituting to a variable that does not exist")
  }
  jj <- seq_along(dim(a))[-i]
  f <- function(x) {sum(x*value^(seq_along(x)-1))}
  out <- apply(a,jj,f)
  if(length(dim(a))>2){
    if(keep){
      dim(out) <- .append2(dim(out),1,after=i-1)
    }
  } else {
    if(i==1){
      out <- t(out)
    } else if (i==2) {
      out <- as.matrix(out)
    } else {
      stop("this cannot happen")
    }
  }
  return(as.multipol(out))
}

"as_function_multipol_vector" <- function(x,vec){
  stopifnot(is.vector(vec))
  x <- as.array(x)
  d <- length(dim(x))
  stopifnot(length(vec) == d)
  
  f <- function(i){vec[i]^(slice.index(x,i)-1)}
  jj <- do.call("abind", c(lapply(seq_len(d),f),along=d+1))
  sum(apply(jj,seq_len(d),prod)*x) 
}

"as_function_multipol" <- function(x,uu){
  if(is.vector(uu)){return(as_function_multipol_vector(x,uu))}
  f <- function(z){as_function_multipol_vector(x,z)}
  out <- apply(uu,1,f)
  return(out)
}

"as.function.multipol" <- function(x, ...){
  function(a){
    as_function_multipol(x,a)
  }
}
