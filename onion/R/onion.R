"biggest" <- function(...){
  a <-  unlist(lapply(list(...),type))
  if("octonion" %in% a){
    return("octonion")
  } else if("quaternion" %in% a)
    {return("quaternion")
   } else {
     return("scalar")
   }
}

"str.onion" <- function(object, vec.len=4, ...){
  string <- class(object)[2]
  if(!is.null(names(object))){
    string <- paste("Named",string,sep=" ")
  }
  cat(paste(string," [1:",length(object),"]\n",sep=""))
  l <- min(length(object),vec.len)
  if(l>0){
    cat(paste(condense(object[1:l]),collapse=", ",sep=""))
  }
  
  if(length(object) > l){
    if(l>0){
      cat(", ")
    }
    cat("...")
    
  }
  cat("\n")
}

"p3d" <- 
function (x, y, z, xlim = NULL, ylim = NULL, zlim = NULL, d0=0.2, h=1, ...) 
{
    if (is.matrix(x)) {
        z <- x[, 3]
        y <- x[, 2]
        x <- x[, 1]
    }
    if (missing(zlim)) {
        z.grid <- matrix(range(z), 2, 2)
    }
    else {
        z.grid <- matrix(zlim, 2, 2)
    }
    if (missing(xlim)) {
        xlim <- range(x)
    }
    if (missing(ylim)) {
        ylim <- range(y)
    }
    res <- persp(xlim, ylim, z.grid, col = NA, border = NA, ...)
    trans3d <- function(x, y, z, pmat) {
        tr <- cbind(x, y, z, 1) %*% pmat
        list(x = tr[, 1]/tr[, 4], y = tr[, 2]/tr[, 4])
    }
    depth3d <- function(x, y, z, pmat) {
        tr <- cbind(x, y, z, 1) %*% pmat
        return(tr[, 3]/tr[, 4])
    }
    rationalize <- function(x){(x-min(x))/(max(x)-min(x))}

    out <- trans3d(x, y, z, pm = res)
    depth <- rationalize(depth3d(x,y,z,pm=res))
    out$x <- out$x[order(depth,decreasing=FALSE)]
    out$y <- out$y[order(depth,decreasing=FALSE)]
    jj <- exp(-sort(depth,decreasing=TRUE)/d0)
    if(is.null(h)){ 
      colours <- hsv(h=1,s=0, v=1-jj)
    } else {
      colours <- hsv(h=h,s=jj)
    }
    points(out,col=colours, ...)
    return(invisible(out))
}

"rotate" <- function(x,H){
t(as.matrix(H*as.quaternion(t(cbind(0,x)))/H))[,-1]
}
  
"type" <- function(x){
  if(is.octonion(x)){
    out <- "octonion"
  } else if (is.quaternion(x)){
    out <- "quaternion"
  } else {
    out <- "scalar"
  }
  return(out)
}

"exp.onion" <- function(x){
  t <- Re(x)
  V <- Im(x)
  mV <- Mod(V)
  out <- exp(t)*cos(mV) + V*exp(t)*sin(mV)/mV
  i <- mV==0 | is.na(mV)
  Re(out[i]) <- exp(t[i])
  Im(out[i]) <- 0
  return(out)
}

"log.onion" <- function(x,base=exp(1)){
  if(!missing(base)){
    return(Recall(x)/Recall(log(base)))
  }
  t <- Re(x)
  V <- Im(x)
  mV <- Mod(V)
  mX <- Mod(x)
  out <- log(Norm(x))/2 + V*atan2(mV,t)/mV
  i <- mV==0  | is.na(mV)
  Re(out[i]) <- log(Norm(x[i]))/2
  Im(out[i]) <- 0
  return(out)
}

"sin.onion" <- function(x){
  t <- Re(x)
  V <- Im(x)
  mV <- Mod(V)
  out <- sin(t)*cosh(mV) + V*cos(t)*sinh(mV)/mV
  i <- mV==0 | is.na(mV)
  Re(out[i]) <- sin(t[i])
  Im(out[i]) <- 0
  return(out)
}

"cos.onion" <- function(x){
  t <- Re(x)
  V <- Im(x)
  mV <- Mod(V)
  out <- cos(t)*cosh(mV) - V*sin(t)*sinh(mV)/mV
  i <- mV==0 | is.na(mV)
  Re(out[i]) <- cos(t[i])
  Im(out[i]) <- 0
  return(out)
}

"tan.onion" <- function(x){
  return(sin(x)/cos(x))
}

"sinh.onion" <- function(x){
  t <- Re(x)
  V <- Im(x)
  mV <- Mod(V)
  out <- sinh(t)*cos(mV) + V*cosh(t)*sin(mV)/mV
  i <- mV==0 | is.na(mV)
  Re(out[i]) <- sinh(t[i])
  Im(out[i]) <- 0
  return(out)
}

"cosh.onion" <- function(x){
  t <- Re(x)
  V <- Im(x)
  mV <- Mod(V)
  out <- cosh(t)*cos(mV) + V*sinh(t)*sin(mV)/mV
  i <- mV==0 | is.na(mV)
  Re(out[i]) <- cosh(t[i])
  Im(out[i]) <- 0
  return(out)
}

"tanh.onion" <- function(x){
  return(sinh(x)/cosh(x))
}

"asinh.onion" <- function(x){
  log(x+sqrt(x*x+1))
}

"acosh.onion" <- function(x){
  log(x+sqrt(x*x-1))
}

"atanh.onion" <- function(x){
  log((1+x)/(1-x))/2
}

"asin.onion" <- function(x){
   V <- Im(x)
   mV <- Mod(V)
   v1 <- V/mV
   out <- x
   i <- mV==0 | is.na(mV)
   out[!i] <- -v1[!i]*asinh(x[!i]*v1[!i])
   Re(out[i]) <- Re(asin(0i+Re(x[i])))
   Im(out[i]) <- 0
   return(out)
 }

"acos.onion" <- function(x){
  V <- Im(x)
  mV <- Mod(V)
  v1 <- V/mV
  out <- x
  i <- mV==0 | is.na(mV)
  out[!i] <- -v1[!i]*acosh(x[!i])
  Re(out[i]) <- Re(acos(0i+Re(x[i])))
  Im(out[i]) <- 0
  return(out)
}
   
"atan.onion" <- function(x){
  V <- Im(x)
  mV <- Mod(V)
  v1 <- V/Mod(V)
  return(-v1*atanh(x*v1))
}

"sqrt.onion" <- function(x){exp(log(x)/2)}

"as.matrix.onion" <- function(x, ...){
  class(x) <- "matrix"
  NextMethod("as.matrix")
}

"roct" <- function(n, x=1:8, replace=TRUE, rand="sample", ...){
  return(switch(rand,
                sample = as.octonion(matrix(sample(x=x, size=8*n, replace=replace, ...),nrow=8)),
                norm = as.octonion(matrix(rnorm(n=8*n, ...),nrow=8)),
                unif = as.octonion(matrix(runif(n=8*n, ...),nrow=8)),
                binom = as.octonion(matrix(rbinom(n=8*n, ...),nrow=8)),
                pois = as.octonion(matrix(rpois(n=8*n, ...),nrow=8))
                )
         )
}

"rquat" <- function(n, x=1:4, replace=TRUE, rand="sample", ...){
  return(switch(rand,
                sample = as.quaternion(matrix(sample(x=x, size=4*n, replace=replace, ...),nrow=4)),
                norm = as.quaternion(matrix(rnorm(n=4*n, ...),nrow=4)),
                unif = as.quaternion(matrix(runif(n=4*n, ...),nrow=4)),
                binom = as.quaternion(matrix(rbinom(n=4*n, ...),nrow=4)),
                pois = as.quaternion(matrix(rpois(n=4*n, ...),nrow=4))
                )
         )
}

"associator" <- function(x1,x2,x3){
  return((x1*x2)*x3 - x1*(x2*x3))
}

"commutator" <- function(x1,x2){
  x1*x2 - x2*x1
}

"threeform" <- function(x1,x2,x3){
  Re(x1*(Conj(x2)*x3) - x3*(Conj(x2)*x1))/2
}

"harmonize" <- function(A,B=A){
  same <- all(class(A)==class(B))
  jj <- Amassage(A,B)
  A <- jj$u1
  B <- jj$u2
  if(is.octonion(A) | is.octonion(B)){
    A <- as.octonion(A)
    B <- as.octonion(B)
  } else if (is.quaternion(A) | is.quaternion(B)){
    A <- as.quaternion(A)
    B <- as.quaternion(B)
  } 
  return(
         list(
              A=A,
              B=B,
              names=jj$names,
              type=type(A),
              same=same
              )
         )
}

"Amassage" <- function(A,B){
  lA <- length(A)
  lB <- length(B)
  if( (lA >= lB) & (!is.null(names(A)))){
    names.out <- names(A)
  } else {
    names.out <- names(B)
  }
  jj <- rbind(seq(length.out=lA),seq(length.out=lB))
  return(list(u1=A[jj[1,]],u2=B[jj[2,]],names=names.out))
}

"AequalsA" <- function(A,B){
  jj <- harmonize(A,B)
  A <- as.matrix(jj$A)
  B <- as.matrix(jj$B)
  f <- function(i){all(A[,i] == B[,i])}
  out <- sapply(1:ncol(A),f)
  names(out) <- jj$names
  return(out)
}

"Ainv" <- function(A){
  out <- sweep(as.matrix(Conj(A)),2,Norm(A),FUN = "/")
  return(as.onion(out,type(A)))
}

"AprodA" <- function(A, B, ur=getOption("use.R")){
  jj <- harmonize(A,B)
  A <- jj$A
  B <- jj$B
  type <- jj$type
  if(isTRUE(ur)){
    out <- 
      switch(type,
             octonion=R_OprodO(A,B),
             quaternion=R_HprodH(A,B),
             A*B )
  } else {
    out <- 
      switch(type,
             octonion=OprodO(A,B),
             quaternion=HprodH(A,B),
             A*B )
  }
  names(out) <- jj$names
  return(out)
}
             

"OprodO" <- function(oct1, oct2){
  if(!is.octonion(oct1) | !is.octonion(oct2)){
    stop("Octonions may be multiplied by scalars or other octonions only")
  }
  jj <- Amassage(oct1,oct2)
  oct1 <- jj$u1
  oct2 <- jj$u2
  names.out <- jj$names
  n <- length(oct1)*8
  jj <- .C("octonion_prod",
           as.double(oct1),
           as.double(oct2),
           as.integer(n),
           z=as.double(oct1),
           PACKAGE="onion"
           )
  out <- jj$z
  dim(out) <- c(8,length(out)/8)
  colnames(out) <- names.out
  return(as.octonion(out))
}

"R_OprodO" <- function(oct1,oct2){
  x <- as.matrix(oct1)
  y <- as.matrix(oct2)
  out <- x
  out[1,] = +x[1,]*y[1,] -x[2,]*y[2,] -x[3,]*y[3,] -x[4,]*y[4,] -x[5,]*y[5,] -x[6,]*y[6,] -x[7,]*y[7,] -x[8,]*y[8,]
  out[2,] = +x[2,]*y[1,] +x[1,]*y[2,] -x[4,]*y[3,] +x[3,]*y[4,] -x[6,]*y[5,] +x[5,]*y[6,] +x[8,]*y[7,] -x[7,]*y[8,]
  out[3,] = +x[3,]*y[1,] +x[4,]*y[2,] +x[1,]*y[3,] -x[2,]*y[4,] -x[7,]*y[5,] -x[8,]*y[6,] +x[5,]*y[7,] +x[6,]*y[8,]
  out[4,] = +x[4,]*y[1,] -x[3,]*y[2,] +x[2,]*y[3,] +x[1,]*y[4,] -x[8,]*y[5,] +x[7,]*y[6,] -x[6,]*y[7,] +x[5,]*y[8,]
  out[5,] = +x[5,]*y[1,] +x[6,]*y[2,] +x[7,]*y[3,] +x[8,]*y[4,] +x[1,]*y[5,] -x[2,]*y[6,] -x[3,]*y[7,] -x[4,]*y[8,]
  out[6,] = +x[6,]*y[1,] -x[5,]*y[2,] +x[8,]*y[3,] -x[7,]*y[4,] +x[2,]*y[5,] +x[1,]*y[6,] +x[4,]*y[7,] -x[3,]*y[8,]
  out[7,] = +x[7,]*y[1,] -x[8,]*y[2,] -x[5,]*y[3,] +x[6,]*y[4,] +x[3,]*y[5,] -x[4,]*y[6,] +x[1,]*y[7,] +x[2,]*y[8,]
  out[8,] = +x[8,]*y[1,] +x[7,]*y[2,] -x[6,]*y[3,] -x[5,]*y[4,] +x[4,]*y[5,] +x[3,]*y[6,] -x[2,]*y[7,] +x[1,]*y[8,]
  return(as.octonion(out))
}

"HprodH" <- function(quat1, quat2){
  if(!is.quaternion(quat1) | !is.quaternion(quat2)){
    stop("Quaternions may be multiplied by scalars or other quaternions only in HprodH()")
  }
  jj <- Amassage(quat1,quat2)
  quat1 <- jj$u1
  quat2 <- jj$u2
  names.out <- jj$names
  n <- length(quat1)*4
  jj <- .C("quaternion_prod",
           as.double(quat1),
           as.double(quat2),
           as.integer(n),
           z=as.double(quat1),
           PACKAGE="onion"
           )
  out <- jj$z
  dim(out) <- c(4,length(out)/4)
  colnames(out) <- names.out
  return(as.quaternion(out))
}

"R_HprodH" <- function(quat1,quat2){
  x <- as.matrix(quat1)
  y <- as.matrix(quat2)
  out <- x
  out[1,] = +x[1,]*y[1,] -x[2,]*y[2,] -x[3,]*y[3,] -x[4,]*y[4,]
  out[2,] = +x[2,]*y[1,] +x[1,]*y[2,] -x[4,]*y[3,] +x[3,]*y[4,]
  out[3,] = +x[3,]*y[1,] +x[4,]*y[2,] +x[1,]*y[3,] -x[2,]*y[4,]
  out[4,] = +x[4,]*y[1,] -x[3,]*y[2,] +x[2,]*y[3,] +x[1,]*y[4,]
  return(as.quaternion(out))
}

"Apower" <- function(A,B){
  exp(log(A)*B)
}

"Aneg" <- function(A){
  out <- -as.matrix(A)
  out <- as.onion(out,type=type(A))
  return(out)
}

"AsumA" <- function(A,B){
  jj <- harmonize(A,B)
  A <- jj$A
  B <- jj$B
  class.out <- jj$class
  if(length(A)<length(B)){
    return(Recall(B,A))
  }
  out <- as.matrix(A) + as.matrix(B)
  return(as.onion(out,type=jj$type,names=jj$names))
}

"AprodS" <- function(A,scalar){
  x <- as.matrix(A)
  if(ncol(x)==1){
    out <- kronecker(x,t(scalar))
  } else {
    out <- sweep(x,2,scalar,"*")
  }
  class(out) <- class(A)
  out <- as.onion(out,type=type(A))
  return(out)
}

"[.onion" <- function(x, ...){
  out <- as.matrix(x)
  out <- out[, ..., drop=FALSE]
  out <- as.onion(out,type=type(x),names=colnames(out))
  return(out)
}

"[<-.onion" <- function(x,index,value){
  out <- as.matrix(x)
  if(is.vector(value)){
    value <- kronecker(t(value),c(1,rep(0,nrow(x)-1)))
  }
  out[,index] <- value
  out <- as.onion(out,type=type(x),names=names(x))
  return(out)
}

"Ops.onion" <-
  function (e1, e2 = NULL) 
{
  f <- function(...){stop("odd---neither argument has class octonion?")}
  unary <- nargs() == 1
  lclass <- nchar(.Method[1]) > 0
  rclass <- !unary && (nchar(.Method[2]) > 0)
  
  if(unary){
    if (.Generic == "+") {
      return(e1)
    } else if (.Generic == "-") {
      return(Aneg(e1))
    } else {
      stop("Unary operator '", .Generic, "' is not implemented for onions")
    }
  }
  if (!is.element(.Generic, c("+", "-", "*", "/", "^", "==", "!=")))
    stop("operator '", .Generic, "' is not implemented for onions")
  
  if (.Generic == "*") {
    if (lclass && rclass) {
      return(AprodA(e1, e2))
    } else if (lclass) {
      return(AprodS(e1, e2))
    } else if (rclass) {
      return(AprodS(e2, e1))
    } else {
      f()
    }
  } else if (.Generic == "/") {
    if (lclass && rclass) {
      return(AprodA(e1, Ainv(e2)))
    } else if (lclass) {
      return(AprodS(e1, 1/e2))
    } else if (rclass) {
      return(AprodS(Ainv(e2), e1))
    } else {
      f()
    }
  } else if (.Generic == "^") {
    return(Apower(e1, e2))
  } else if (.Generic == "+") { 
    return(AsumA(e1, e2)) 
  } else if (.Generic == "-") { 
    return(AsumA(e1, Aneg(e2)))
  } else if (.Generic == "==") {
    return("AequalsA"(e1,e2))
  } else if (.Generic == "!=") {
    return(!"AequalsA"(e1,e2))
  } else {
    stop("should not reach here")
  }
}

"octonion" <- function(length.out=NULL, names=NULL, Re=0, i=0, j=0, k=0, l=0, il=0, jl=0, kl=0){
  if (
      (missing(Re) | length(Re)==0) &
      (missing(i)  | length( i)==0) &
      (missing(j)  | length( j)==0) &
      (missing(k)  | length( k)==0) &
      (missing(l)  | length( l)==0) &
      (missing(il) | length(il)==0) &
      (missing(jl) | length(jl)==0) &
      (missing(kl) | length(kl)==0)
      )
    {
      if(missing(length.out)){
        return(Recall(Re=0)[0])
      } else {
        return(Recall(Re=rep(0,length.out)))
      }
    }
  
  out <- as.octonion(rbind(Re,i,j,k,l,il,jl,kl))
  if(!is.null(length.out)){
    if(ncol(out)==1){
      out <- as.octonion(kronecker(out,t(rep(1,length.out))))
    } else {
      length(out) <- length.out
    }
  }
  if(!is.null(names)){
    colnames(out) <- names
  }
  return(out)
}

"quaternion" <- function(length.out=NULL, names=NULL, Re=0, i=0, j=0, k=0){
  if (
      (missing(Re) | length(Re)==0) &
      (missing(i)  | length( i)==0) &
      (missing(j)  | length( j)==0) &
      (missing(k)  | length( k)==0) 
      )
    {
      if(missing(length.out)){
        return(Recall(Re=0)[0])
      } else {
        return(Recall(Re=rep(0,length.out)))
      }
    }
  
  out <- as.quaternion(rbind(Re,i,j,k))
  if(!is.null(length.out)){
    if(ncol(out)==1){
      out <- as.quaternion(kronecker(out,t(rep(1,length.out))))
    } else {
      length(out) <- length.out
    }
  }
  if(!is.null(names)){
    colnames(out) <- names
  }
  return(out)
}

"is.octonion" <- function(x){
  inherits(x,"octonion")
}

"is.quaternion" <- function(x){
  inherits(x,"quaternion")
}

"is.onion" <- function(x){
  inherits(x,"onion")
}

"as.onion" <- function(x, type, names=NULL, single=FALSE){
  if(single){
    return(Recall(as.matrix(as.vector(x)),type=type,names=names))
  }
  if(is.null(names)){names <- colnames(x)}
  switch(type,
         quaternion=as.quaternion(x, names=names),
         octonion=as.octonion(x, names=names),
         x
         )
}

"as.octonion" <- function(x, single=FALSE, names=NULL){
  if(is.octonion(x)){return(x)}
  if(is.quaternion(x)){return(as.octonion(rbind(x,x*0),names=names))}
  if(is.complex(x)){
    if(single){
      stop("single cannot be TRUE with complex x")
    } else {
      return(Recall(rbind(Re(x),Im(x),matrix(0,6,length(x))),names=names))
    }
  }
  if(is.matrix(x)){
    if(nrow(x) == 8){
      out <- x
    } else {
      stop("If matrix supplied, it must have eight rows")
    }
  } else {
    if(single){
      if(is.vector(x)){
        out <- as.matrix(x[1:8])
        if(length(x) != 8){
          warning("single set to TRUE, but a vector of length !=8 supplied.  Setting to length 8, Procrustes-style")
        }
      } else {
        stop("single set to TRUE, but nonvector supplied.")
      }
    } else {
      x <- as.vector(x)
      out <- kronecker(t(x),c(1,rep(0,7)))
    }
  }
    colnames(out) <- names
  class(out) <- c("onion","octonion")
  return(out)
}

"as.quaternion" <- function(x,single=FALSE,names=NULL){
  if(is.quaternion(x)){return(x)}
  if(is.complex(x)){
    if(single){
      stop("single cannot be TRUE with complex x")
    } else {
      return(Recall(rbind(Re(x),Im(x),matrix(0,2,length(x)))))
    }
  }
  if(is.matrix(x)){
    if(nrow(x) == 4){
      out <- x
    } else {
      stop("If matrix supplied, it must have four rows")
    }
  } else {
    if(single){
      if(is.vector(x)){
        out <- as.matrix(x[1:4])
        if(length(x) != 4){
          warning("single set to TRUE, but a vector of length !=4 supplied. Set to length 4, Procrustes-style")
        }
      } else {
        stop("single set to TRUE, but nonvector supplied.")
      }
    } else {
      x <- as.vector(x)
      out <- kronecker(t(x),c(1,rep(0,3)))
    }
  }
  class(out) <- c("onion","quaternion")
  rownames(out) <- NULL
  names(out) <- names
  return(out)
}

"names.onion" <- function(x){colnames(x)}

"names<-.onion" <- function(x,value){
  colnames(x) <- value
  return(x)
}

"condense" <- function(x){UseMethod("condense")}
"condense.onion" <- function(x){
  x <- as.matrix(x)
  out <- x
  out[x==0] <- "0"
  out[x> 0] <- "+"
  out[x< 0] <- "-"
  return(noquote(apply(out,2,paste,collapse="")))
  }

"c.onion" <- function(...){
  b <- biggest(...)
  as.onion(do.call("cbind",lapply(list(...) , as.onion,type=b)),type=b)
}

"t.onion" <- function(x){
  x <- as.matrix(x)
  NextMethod("t")
}

"sum.onion" <- function(... , na.rm=FALSE){
  a <- list(...)
  b <- biggest(...)
  a <- lapply(a,FUN=function(u){as.matrix(as.onion(u,type=b))})
  out <- do.call("cbind",a)
  out <- apply(out,1,base::sum)
  return(as.onion(out, type=b, single=TRUE))
}

"cumsum.onion" <- function(x){
  out <- t(apply(as.matrix(x),1,cumsum))
return(as.onion(out,type(x)))
}

"print.octonion" <- function(x, h=getOption("horiz"), ...){
  x <- as.matrix(x)
  rownames(x) <- c("Re","i","j","k","l","il","jl","kl")
  if(ncol(x)>0){
    if(is.null(colnames(x))){colnames(x) <- paste("[",1:ncol(x),"]",sep="")}
  }
  if(isTRUE(h)){
    return(invisible(print(t(x))))
  } else {
    return(invisible(print(x)))
  }
}

"print.quaternion" <- function(x, h=getOption("horiz"), ...){
  x <- as.matrix(x)
  rownames(x) <- c("Re","i","j","k")
  if(ncol(x)>0){
    if(is.null(colnames(x))){colnames(x) <- paste("[",1:ncol(x),"]",sep="")}
  }
  if(isTRUE(h)){
    return(invisible(print(t(x))))
  } else {
    return(invisible(print(x)))
  }
}

"set.comp<-" <- function(x, i, value){UseMethod("set.comp<-")}

"set.comp<-.onion" <- function(x,i,value){
  tx <- type(x)
  nx <- names(x)
  x <- as.matrix(x)
  x[i,] <- value
  return(as.onion(x, type=tx,names=nx))
}

"get.comp" <- function(x,i){UseMethod("get.comp")}
"get.comp.onion" <- function(x,i){
  as.matrix(x)[i,]
}

"seq.onion" <-
  function (from = 1, to = 1, by = ((to - from)/(length.out - 1)), length.out = NULL, slerp = FALSE, ...) 
{
  b <- biggest(from, to, by)
  if (identical(length.out,0)) {return(as.onion(0,type=b)[0])}

  from <- as.onion(from,type=b)
  to <- as.onion(to,type=b)
  by <- as.onion(by,type=b)

  if (!missing(length.out)){ 
    length.out <- ceiling(length.out)
  }

  if(missing(to)){
    to <- from + by*(length.out-1)
  }
  if(missing(from)){
    from <- to - by*(length.out-1)
  }
  del <- to - from
  if(missing(by)){
    by <- del/length.out
  }

  h <- seq(from=0,to=1,len=length.out)
  if(slerp){
    return(from*(to/from)^h)
  } else {
    return(from*(1-h)+ to*h)
  }
}

"length.onion" <- function(x){ncol(x)}

"length<-.onion" <- function(x,value){
  if(value <= length(x)){
    return(x[1:value])
  } else {
    out <- as.matrix(x)
    out <- cbind(out,matrix(NA,nrow(out),value-ncol(out)))
    return(as.onion(out,type=type(x)))
  }
}

"plot.onion" <- function(x, ...){plot(Re(x),Mod(Im(x)), ...)}

"rep.onion" <- function(x,  ...){
  u <- seq(length.out=length(x))
  return(x[rep(u, ...)])
}
    
"Re<-" <- function(x, value){UseMethod("Re<-")}
"Im<-" <- function(x, value){UseMethod("Im<-")}
 "i<-" <- function(x, value){UseMethod("i<-")}
 "j<-" <- function(x, value){UseMethod("j<-")}
 "k<-" <- function(x, value){UseMethod("k<-")}
 "l<-" <- function(x, value){UseMethod("l<-")}
"il<-" <- function(x, value){UseMethod("il<-")}
"jl<-" <- function(x, value){UseMethod("jl<-")}
"kl<-" <- function(x, value){UseMethod("kl<-")}

"Re<-.octonion" <- function(x,value){"set.comp<-"(x,1,value)}
 "i<-.octonion" <- function(x,value){"set.comp<-"(x,2,value)}
 "j<-.octonion" <- function(x,value){"set.comp<-"(x,3,value)}
 "k<-.octonion" <- function(x,value){"set.comp<-"(x,4,value)}
 "l<-.octonion" <- function(x,value){"set.comp<-"(x,5,value)}
"il<-.octonion" <- function(x,value){"set.comp<-"(x,6,value)}
"jl<-.octonion" <- function(x,value){"set.comp<-"(x,7,value)}
"kl<-.octonion" <- function(x,value){"set.comp<-"(x,8,value)}
"Im<-.octonion" <- function(x,value){
  nx <- names(x)
  x <- as.matrix(x)
  if(isTRUE(all.equal(nrow(value), 7)) || isTRUE(all.equal(length(value),7))) {
    x[-1,] <- value
  } else if (is.octonion(value)){
    if(max(Mod(Re(value)))>0){
      warning("Discarding nonzero real part of value")
    }
    x[-1,] <- as.matrix(value)[-1,]
  } else if (identical(value,0)){
    x[-1,] <- 0
  } else {
    stop("Error: when specifying imaginary part, give a seven-element vector, or a seven-row matrix, or an octonion with zero real part")
  }
  out <- as.octonion(x)
  names(out) <- nx
  return(out)
}

"Re<-.quaternion" <- function(x,value){
  "set.comp<-"(x,1,value)
}

 "i<-.quaternion" <- function(x,value){"set.comp<-"(x,2,value)}
 "j<-.quaternion" <- function(x,value){"set.comp<-"(x,3,value)}
 "k<-.quaternion" <- function(x,value){"set.comp<-"(x,4,value)}
"Im<-.quaternion" <- function(x,value){
  nx <- names(x)
  x <- as.matrix(x)
  if(isTRUE(all.equal(nrow(value), 3)) || isTRUE(all.equal(length(value),3))) {
    x[-1,] <- value
  } else if (is.quaternion(value)){
    if(max(Mod(Re(value)))>0){
      warning("Discarding nonzero real part of value")
    }
    x[-1,] <- as.matrix(value)[-1,]
  } else if (identical(value,0)){
    x[-1,] <- 0
  } else {
    stop("Error: when specifying imaginary part, give a three-element vector, or a three-row matrix, or an quaternion with zero real part")
  }
  out <- as.quaternion(x)
  names(out) <- nx
  return(out)
}

if(FALSE){
"Re.default" <- get("Re",pos="package:base",mode="function")
"Im.default" <- get("Im",pos="package:base",mode="function")
}
"i"  <- function(x){UseMethod("i" )}
"j"  <- function(x){UseMethod("j" )}
"k"  <- function(x){UseMethod("k" )}
"l"  <- function(x){UseMethod("l" )}
"il" <- function(x){UseMethod("il")}
"jl" <- function(x){UseMethod("jl")}
"kl" <- function(x){UseMethod("kl")}


"Re.octonion" <- function(z){get.comp(z,1)}
 "i.octonion" <- function(x){get.comp(x,2)}
 "j.octonion" <- function(x){get.comp(x,3)}
 "k.octonion" <- function(x){get.comp(x,4)}
 "l.octonion" <- function(x){get.comp(x,5)}
"il.octonion" <- function(x){get.comp(x,6)}
"jl.octonion" <- function(x){get.comp(x,7)}
"kl.octonion" <- function(x){get.comp(x,8)}
"Im.octonion" <- function(z){
  Re(z) <- 0
  return(z)
}

"Re.quaternion" <- function(z){get.comp(z,1)}
 "i.quaternion" <- function(x){get.comp(x,2)}
 "j.quaternion" <- function(x){get.comp(x,3)}
 "k.quaternion" <- function(x){get.comp(x,4)}
"Im.quaternion" <- function(z){
  Re(z) <- 0
  return(z)
}

"Conj.onion" <- function(z){
  Im(z) <- -Im(z)
  return(z)
}

"Mod.onion" <- function(z){
  return(sqrt(apply(as.matrix(z),2,function(u){sum(u^2)})))
}

"Mod<-" <- function(z, value){UseMethod("Mod<-")}

"Mod<-.onion" <- function(z,value){
  if(is.octonion(value)|is.complex(value)){
    stop("RHS must be real")
    } else {
      return(value*sign(z))
    }
}

"Norm" <- function(z){UseMethod("Norm")}

"Norm.onion" <- function(z){
  return(apply(z,2,function(x){sum(x^2)}))
}

"sign.onion" <- function(x){x/Mod(x)}

"g.even" <- function(x,y){UseMethod("g.even")}
 "g.odd" <- function(x,y){UseMethod("g.odd")}
"e.even" <- function(x,y){UseMethod("e.even")}
 "e.odd" <- function(x,y){UseMethod("e.odd")}

"g.even.onion" <- function(x,y){(x*y + y*x)/2}
 "g.odd.onion" <- function(x,y){(x*y - y*x)/2}
"e.even.onion" <- function(x,y){(Conj(x)*y + Conj(y)*x)/2}
 "e.odd.onion" <- function(x,y){(Conj(x)*y - Conj(y)*x)/2}

"%<*>%" <- function(x,y){UseMethod("%<*>%")}  #g.even
"%>*<%" <- function(x,y){UseMethod("%>*<%")}  #g.odd
"%<.>%" <- function(x,y){UseMethod("%<.>%")}  #e.even
"%>.<%" <- function(x,y){UseMethod("%>.<%")}  #e.odd

"%<*>%.onion" <- function(x,y){g.even.onion(x,y)}
"%>*<%.onion" <- function(x,y){ g.odd.onion(x,y)}
"%<.>%.onion" <- function(x,y){e.even.onion(x,y)}
"%>.<%.onion" <- function(x,y){ e.odd.onion(x,y)}

"dotprod" <- function(x,y){
  x <- as.matrix(x)
  y <- as.matrix(y)
  if(ncol(x)==1){
    x <- as.vector(x)
    return(apply(y,2,function(u){sum(u*x)}))
  } else if (ncol(y)==1){
    y <- as.vector(y)
    return(apply(x,2,function(u){sum(u*y)}))
  } else {
    return(apply(x*y,2,sum))
  }
}

"%.%" <- function(x,y){dotprod(x,y)}
