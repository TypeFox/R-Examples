## kernel functions
## Functions for computing a kernel value, matrix, matrix-vector
## product and quadratic form
##
## author : alexandros karatzoglou


## Define the kernel objects,
## functions with an additional slot for the kernel parameter list.
## kernel functions take two vector arguments and return a scalar (dot product)


rbfdot<- function(sigma=1)
  {

    rval <- function(x,y=NULL)
    {
       if(!is(x,"vector")) stop("x must be a vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must a vector")
      if (is(x,"vector") && is.null(y)){
        return(1)
      }
      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
        return(exp(sigma*(2*crossprod(x,y) - crossprod(x) - crossprod(y))))  
        # sigma/2 or sigma ??
      }
    }
     return(new("rbfkernel",.Data=rval,kpar=list(sigma=sigma)))
  }
setClass("rbfkernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

laplacedot<- function(sigma=1)
  {
    rval <- function(x,y=NULL)
    {
       if(!is(x,"vector")) stop("x must be a vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must a vector")
      if (is(x,"vector") && is.null(y)){
        return(1)
      }
      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
        return(exp(-sigma*sqrt(-(round(2*crossprod(x,y) - crossprod(x) - crossprod(y),9)))))
      }
    }
     return(new("laplacekernel",.Data=rval,kpar=list(sigma=sigma)))
  }

setClass("laplacekernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

besseldot<- function(sigma = 1, order = 1, degree = 1)
  {
    rval <- function(x,y=NULL)
    {
      if(!is(x,"vector")) stop("x must be a vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must a vector")
      if (is(x,"vector") && is.null(y)){
        return(1)
      }
      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
        lim <- 1/(gamma(order+1)*2^(order))
        bkt <- sigma*sqrt(-(2*crossprod(x,y) - crossprod(x) - crossprod(y)))
        if(bkt < 10e-5)
          res <- lim
        else
          res <- besselJ(bkt,order)*(bkt^(-order))
        return((res/lim)^degree)
      }
    }
     return(new("besselkernel",.Data=rval,kpar=list(sigma=sigma ,order = order ,degree = degree)))
  }

setClass("besselkernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

anovadot<- function(sigma = 1, degree = 1)
  {
    rval <- function(x,y=NULL)
    {
      if(!is(x,"vector")) stop("x must be a vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must a vector")
      if (is(x,"vector") && is.null(y)){
        return(1)
      }
      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")

          res <- sum(exp(- sigma * (x - y)^2))
        return((res)^degree)
      }
    }
     return(new("anovakernel",.Data=rval,kpar=list(sigma=sigma ,degree = degree)))
  }

setClass("anovakernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))


splinedot<- function()
  {
    rval <- function(x,y=NULL)
    {
      if(!is(x,"vector")) stop("x must be a vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must a vector")
      if (is(x,"vector") && is.null(y)){
        return(1)
      }
      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
        minv <- pmin(x,y)
        res <- 1 + x*y*(1+minv) - ((x+y)/2)*minv^2 + (minv^3)/3
          fres <- prod(res)
        return(fres)
      }
    }
     return(new("splinekernel",.Data=rval,kpar=list()))
  }

setClass("splinekernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))



fourierdot <- function(sigma = 1)
  {
    rval <- function(x,y=NULL)
    {
      if(!is(x,"vector")) stop("x must be a vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must a vector")
      if (is(x,"vector") && is.null(y)){
        return(1)
      }
      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
	res <- 	(1 - sigma^2)/2*(1 - 2*sigma*cos(x - y) + sigma^2)
	fres <- prod(res)
        return(fres)
      }
    }
     return(new("fourierkernel",.Data=rval,kpar=list()))
  }

setClass("fourierkernel",prototype=structure(.Data=function(){},kpar=list(sigma = 1)),contains=c("kernel"))





tanhdot <- function(scale = 1, offset = 1)
{
  rval<- function(x, y = NULL)
    {
      if(!is(x,"vector")) stop("x must be a  vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must be a vector")
      if (is(x,"vector") && is.null(y)){
        tanh(scale*crossprod(x)+offset)
      }
      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
        tanh(scale*crossprod(x,y)+offset)
      }
    }
  return(new("tanhkernel",.Data=rval,kpar=list(scale=scale,offset=offset)))
}
setClass("tanhkernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

setClass("polykernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

polydot <- function(degree = 1, scale = 1, offset = 1)
{
  rval<- function(x, y = NULL)
    {
      if(!is(x,"vector")) stop("x must be a vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must be a vector")
      if (is(x,"vector") && is.null(y)){
        (scale*crossprod(x)+offset)^degree
        }

      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
            stop("number of dimension must be the same on both data points")
        (scale*crossprod(x,y)+offset)^degree
          }

      }
  return(new("polykernel",.Data=rval,kpar=list(degree=degree,scale=scale,offset=offset)))
}

setClass("vanillakernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

vanilladot <- function( )
{
  rval<- function(x, y = NULL)
    {
      if(!is(x,"vector")) stop("x must be a vector")
      if(!is(y,"vector")&&!is.null(y)) stop("y must be a vector")
      if (is(x,"vector") && is.null(y)){
        crossprod(x)
        }

      if (is(x,"vector") && is(y,"vector")){
        if (!length(x)==length(y))
            stop("number of dimension must be the same on both data points")
        crossprod(x,y)
          }

      }
  return(new("vanillakernel",.Data=rval,kpar=list()))
}

setClass("stringkernel",prototype=structure(.Data=function(){},kpar=list(length = 4, lambda = 1.1, type = "spectrum", normalized = TRUE)),contains=c("kernel"))

stringdot <- function(length = 4, lambda = 1.1, type = "spectrum", normalized = TRUE)
{
  type <- match.arg(type,c("sequence","string","fullstring","exponential","constant","spectrum", "boundrange"))

  ## need to do this to set the length parameters
  if(type == "spectrum" | type == "boundrange")
    lambda <- length
  
   switch(type,
          "sequence" = {

            rval<- function(x, y = NULL)
              {
                if(!is(x,"vector")) stop("x must be a vector")
                if(!is(y,"vector")&&!is.null(y)) stop("y must be a vector")
                
                if (is(x,"vector") && is.null(y) && normalized == FALSE)
                  return(.Call("subsequencek",as.character(x), as.character(x), as.integer(nchar(x)), as.integer(nchar(x)), as.integer(length), as.double(lambda),PACKAGE = "kernlab"))
                if (is(x,"vector") && is(y,"vector") && normalized == FALSE)
                  return(.Call("subsequencek",as.character(x), as.character(y), as.integer(nchar(x)), as.integer(nchar(y)), as.integer(length), as.double(lambda),PACKAGE = "kernlab"))
                if (is(x,"vector") && is.null(y) && normalized == TRUE)
                  return(1)
                if (is(x,"vector") && is(y,"vector") && normalized == TRUE)
                  return(.Call("subsequencek",as.character(x), as.character(y), as.integer(nchar(x)), as.integer(nchar(y)), as.integer(length), as.double(lambda),PACKAGE = "kernlab")/sqrt(.Call("subsequencek",as.character(x), as.character(x), as.integer(nchar(x)), as.integer(nchar(x)), as.integer(length), as.double(lambda),PACKAGE = "kernlab")*.Call("subsequencek",as.character(y), as.character(y), as.integer(nchar(y)), as.integer(nchar(y)), as.integer(length), as.double(lambda),PACKAGE = "kernlab")))
              }
          },
          
          "exponential" = {
            rval <- function(x,y=NULL)
              {
                if(!is(x,"vector")) stop("x must be a vector")
                if(!is(y,"vector")&&!is.null(y)) stop("y must be a vector")

                x <- paste(x,"\n",sep="")
                if(!is.null(y))
                  y <- paste(y,"\n",sep="")

                if (normalized == FALSE){
                  if(is.null(y))
                    y <- x
                  return(.Call("stringtv",as.character(x),as.character(y),as.integer(1),as.integer(nchar(x)),as.integer(nchar(y)),as.integer(2),as.double(lambda)))}
                if (is(x,"vector") && is.null(y) && normalized == TRUE)
                  return(1)
                if (is(x,"vector") && is(y,"vector") && normalized == TRUE)
                  return(.Call("stringtv",as.character(x),as.character(y),as.integer(1),as.integer(nchar(x)),as.integer(nchar(y)),as.integer(2),as.double(lambda))/sqrt(.Call("stringtv",as.character(x),as.character(x),as.integer(1),as.integer(nchar(x)),as.integer(nchar(x)),as.integer(2),as.double(lambda))*.Call("stringtv",as.character(y),as.character(y),as.integer(1),as.integer(nchar(y)),as.integer(nchar(y)),as.integer(2),as.double(lambda))))
                
              }
          },

          "constant" = {
            rval <- function(x,y=NULL)
              {
                if(!is(x,"vector")) stop("x must be a vector")
                if(!is(y,"vector")&&!is.null(y)) stop("y must be a vector")

                x <- paste(x,"\n",sep="")
                if(!is.null(y))
                  y <- paste(y,"\n",sep="")
                
                if (normalized == FALSE){
                  if(is.null(y))
                    y <- x
                  return(.Call("stringtv",as.character(x),as.character(y),as.integer(1),as.integer(nchar(x)),as.integer(nchar(y)),as.integer(1),as.double(lambda)))}
                if (is(x,"vector") && is.null(y) && normalized == TRUE)
                  return(1)
                if (is(x,"vector") && is(y,"vector") && normalized == TRUE)
                  return(.Call("stringtv",as.character(x),as.character(y),as.integer(1),as.integer(nchar(x)),as.integer(nchar(y)),as.integer(1),as.double(lambda))/sqrt(.Call("stringtv",as.character(x),as.character(x),as.integer(1),as.integer(nchar(x)),as.integer(nchar(x)),as.integer(1),as.double(lambda))*.Call("stringtv",as.character(y),as.character(y),as.integer(1),as.integer(nchar(y)),as.integer(nchar(y)),as.integer(1),as.double(lambda))))
              }
          },
          
          "spectrum" = {
            rval <- function(x,y=NULL)
              {
                if(!is(x,"vector")) stop("x must be a vector")
                if(!is(y,"vector")&&!is.null(y)) stop("y must be a vector")

                x <- paste(x,"\n",sep="")
                if(!is.null(y))
                   y <- paste(y,"\n",sep="")
                
                n <- nchar(x)
                m <- nchar(y)
                if(n < length | m < length){
                  warning("String length smaller than length parameter value")
                  return(0)}

                if (normalized == FALSE){
                  if(is.null(y))
                    y <- x
                  return(.Call("stringtv",as.character(x),as.character(y),as.integer(1),as.integer(n),as.integer(m),as.integer(3),as.double(length)))}
                if (is(x,"vector") && is.null(y) && normalized == TRUE)
                  return(1)
                if (is(x,"vector") && is(y,"vector") && normalized == TRUE)
                  return(.Call("stringtv",as.character(x),as.character(y),as.integer(1),as.integer(n),as.integer(m),as.integer(3),as.double(length))/sqrt(.Call("stringtv",as.character(x),as.character(x),as.integer(1),as.integer(n),as.integer(n),as.integer(3),as.double(lambda))*.Call("stringtv",as.character(y),as.character(y),as.integer(1),as.integer(m),as.integer(m),as.integer(3),as.double(length))))
                
              }
          },
          "boundrange" = {
            rval <- function(x,y=NULL)
              {
                if(!is(x,"vector")) stop("x must be a vector")
                if(!is(y,"vector")&&!is.null(y)) stop("y must be a vector")

                x <- paste(x,"\n",sep="")
                if(!is.null(y))
                   y <- paste(y,"\n",sep="")
                
                if (normalized == FALSE){
                  if(is.null(y))
                    y <- x
                  return(.Call("stringtv",as.character(x),as.character(y),as.integer(1),as.integer(nchar(x)),as.integer(nchar(y)),as.integer(4),as.double(lambda)))}
                if (is(x,"vector") && is.null(y) && normalized == TRUE)
                  return(1)
                if (is(x,"vector") && is(y,"vector") && normalized == TRUE)
                  return(.Call("stringtv",as.character(x),as.character(y),as.integer(1),as.integer(nchar(x)),as.integer(nchar(y)),as.integer(4),as.double(lambda))/sqrt(.Call("stringtv",as.character(x),as.character(x),as.integer(1),as.integer(nchar(x)),as.integer(nchar(x)),as.integer(4),as.double(lambda))*.Call("stringtv",as.character(y),as.character(y),as.integer(1),as.integer(nchar(y)),as.integer(nchar(y)),as.integer(4),as.double(lambda))))
                
              }
          },

          "string" =
          {
            rval<- function(x, y = NULL)
              {
                if(!is(x,"vector")) stop("x must be a vector")
                if(!is(y,"vector")&&!is.null(y)) stop("y must be a vector")
                
                if (is(x,"vector") && is.null(y) && normalized == FALSE)
                  return(.Call("substringk",as.character(x), as.character(x), as.integer(nchar(x)), as.integer(nchar(x)), as.integer(length), as.double(lambda),PACKAGE = "kernlab"))
                
                if (is(x,"vector") && is(y,"vector") && normalized == FALSE)
                  return(.Call("substringk",as.character(x), as.character(y), as.integer(nchar(x)), as.integer(nchar(y)), as.integer(length), as.double(lambda),PACKAGE = "kernlab"))
                if (is(x,"vector") && is.null(y) && normalized == TRUE)
                  return(1)
                if (is(x,"vector") && is(y,"vector") && normalized == TRUE)
                  return(.Call("substringk",as.character(x), as.character(y), as.integer(nchar(x)), as.integer(nchar(y)), as.integer(length), as.double(lambda),PACKAGE = "kernlab")/sqrt(.Call("substringk",as.character(x), as.character(x), as.integer(nchar(x)), as.integer(nchar(x)), as.integer(length), as.double(lambda),PACKAGE = "kernlab")*.Call("substringk",as.character(y), as.character(y), as.integer(nchar(y)), as.integer(nchar(y)), as.integer(length), as.double(lambda),PACKAGE = "kernlab")))
              }
          },
          
          "fullstring" =
          {
            rval<- function(x, y = NULL)
              {
                if(!is(x,"vector")) stop("x must be a vector")
                if(!is(y,"vector")&&!is.null(y)) stop("y must be a vector")
                
                if (is(x,"vector") && is.null(y) && normalized == FALSE)
                  return(.Call("fullsubstringk",as.character(x), as.character(x), as.integer(nchar(x)), as.integer(nchar(x)), as.integer(length), as.double(lambda),PACKAGE = "kernlab"))
                
                if (is(x,"vector") && is(y,"vector") && normalized == FALSE)
                  return(.Call("fullsubstringk",as.character(x), as.character(y), as.integer(nchar(x)), as.integer(nchar(y)), as.integer(length), as.double(lambda),PACKAGE = "kernlab"))
                if (is(x,"vector") && is.null(y) && normalized == TRUE)
                  return(1)
                if (is(x,"vector") && is(y,"vector") && normalized == TRUE)
                  return(.Call("fullsubstringk",as.character(x), as.character(y), as.integer(nchar(x)), as.integer(nchar(y)), as.integer(length), as.double(lambda),PACKAGE = "kernlab")/sqrt(.Call("fullsubstringk",as.character(x), as.character(x), as.integer(nchar(x)), as.integer(nchar(x)), as.integer(length), as.double(lambda),PACKAGE = "kernlab")*.Call("fullsubstringk",as.character(y), as.character(y), as.integer(nchar(y)), as.integer(nchar(y)), as.integer(length), as.double(lambda),PACKAGE = "kernlab")))
              }
          })
  
  return(new("stringkernel",.Data=rval,kpar=list(length=length, lambda =lambda, type = type, normalized = normalized)))
}

## show method for kernel functions

setMethod("show",signature(object="kernel"),
          function(object)
          {
            switch(class(object),
                   "rbfkernel" = cat(paste("Gaussian Radial Basis kernel function.", "\n","Hyperparameter :" ,"sigma = ", kpar(object)$sigma,"\n")),
		   "laplacekernel" = cat(paste("Laplace kernel function.", "\n","Hyperparameter :" ,"sigma = ", kpar(object)$sigma,"\n")),
                   "besselkernel" = cat(paste("Bessel kernel function.", "\n","Hyperparameter :" ,"sigma = ", kpar(object)$sigma,"order = ",kpar(object)$order, "degree = ", kpar(object)$degree,"\n")),
                    "anovakernel" = cat(paste("Anova RBF kernel function.", "\n","Hyperparameter :" ,"sigma = ", kpar(object)$sigma, "degree = ", kpar(object)$degree,"\n")),
                   "tanhkernel" = cat(paste("Hyperbolic Tangent kernel function.", "\n","Hyperparameters :","scale = ", kpar(object)$scale," offset = ", kpar(object)$offset,"\n")),
                   "polykernel" = cat(paste("Polynomial kernel function.", "\n","Hyperparameters :","degree = ",kpar(object)$degree," scale = ", kpar(object)$scale," offset = ", kpar(object)$offset,"\n")),
                   "vanillakernel" = cat(paste("Linear (vanilla) kernel function.", "\n")),
                   "splinekernel" = cat(paste("Spline kernel function.", "\n")),

                   "stringkernel" = {
                     if(kpar(object)$type =="spectrum" | kpar(object)$type =="boundrange")
                       cat(paste("String kernel function.", " Type = ", kpar(object)$type, "\n","Hyperparameters :","sub-sequence/string length = ",kpar(object)$length, "\n"))
                     else
                       if(kpar(object)$type =="exponential" | kpar(object)$type =="constant")
                         cat(paste("String kernel function.", " Type = ", kpar(object)$type, "\n","Hyperparameters :"," lambda = ", kpar(object)$lambda, "\n"))
                       else
                         cat(paste("String kernel function.", " Type = ", kpar(object)$type, "\n","Hyperparameters :","sub-sequence/string length = ",kpar(object)$length," lambda = ", kpar(object)$lambda, "\n"))
                   if(kpar(object)$normalized == TRUE) cat(" Normalized","\n")
                   if(kpar(object)$normalized == FALSE) cat(" Not Normalized","\n")}
                   )
          })

## create accesor function as in "S4 Classses in 15 pages more or less", well..

if (!isGeneric("kpar")){
  if (is.function(kpar))
    fun <- kpar
  else fun <- function(object) standardGeneric("kpar")
  setGeneric("kpar",fun)
}

setMethod("kpar","kernel", function(object) object@kpar)




## Functions that return usefull kernel calculations (kernel matrix etc.)

## kernelMatrix function takes two or three arguments

kernelMatrix <- function(kernel, x, y=NULL)
{
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if(!is(x,"matrix")) stop("x must be a matrix")
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix")
  n <- nrow(x)
  res1 <- matrix(rep(0,n*n), ncol = n)
  if(is.null(y)){
    for(i in 1:n) {
      for(j in i:n) {
        res1[i,j]  <- kernel(x[i,],x[j,])
      }
    }
    res1 <- res1 + t(res1)
    diag(res1) <- diag(res1)/2
    
  
 }
  if (is(y,"matrix")){
    m<-dim(y)[1]
    res1 <- matrix(0,dim(x)[1],dim(y)[1])
    for(i in 1:n) {
      for(j in 1:m) {
        res1[i,j] <- kernel(x[i,],y[j,])
      }
    }
  }

  return(as.kernelMatrix(res1))
}

setGeneric("kernelMatrix",function(kernel, x, y = NULL) standardGeneric("kernelMatrix"))



kernelMatrix.rbfkernel <- function(kernel, x, y = NULL)
{
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  n <- dim(x)[1]
  dota <- rowSums(x*x)/2
  if (is(x,"matrix") && is.null(y)){
    res <- crossprod(t(x))
    for (i in 1:n)
      res[i,]<- exp(2*sigma*(res[i,] - dota - rep(dota[i],n)))
    return(as.kernelMatrix(res))
  }
  if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    m <- dim(y)[1]
    dotb <- rowSums(y*y)/2
    res <- x%*%t(y)
    for( i in 1:m)
      res[,i]<- exp(2*sigma*(res[,i] - dota - rep(dotb[i],n)))
    return(as.kernelMatrix(res))
  }
}
setMethod("kernelMatrix",signature(kernel="rbfkernel"),kernelMatrix.rbfkernel)

kernelMatrix.laplacekernel <- function(kernel, x, y = NULL)
{
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  n <- dim(x)[1]
  dota <- rowSums(x*x)/2
  if (is(x,"matrix") && is.null(y)){
    res <- crossprod(t(x))
    for (i in 1:n)
      res[i,]<- exp(-sigma*sqrt(round(-2*(res[i,] - dota - rep(dota[i],n)),9)))
    return(as.kernelMatrix(res))
  }
  if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    m <- dim(y)[1]
    dotb <- rowSums(y*y)/2
    res <- x%*%t(y)
    for( i in 1:m)
      res[,i]<- exp(-sigma*sqrt(round(-2*(res[,i] - dota - rep(dotb[i],n)),9)))
    return(as.kernelMatrix(res))
  }
}
setMethod("kernelMatrix",signature(kernel="laplacekernel"),kernelMatrix.laplacekernel)

kernelMatrix.besselkernel <- function(kernel, x, y = NULL)
{
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  nu = kpar(kernel)$order
  ni = kpar(kernel)$degree
  n <- dim(x)[1]
  lim <- 1/(gamma(nu+1)*2^(nu))
  dota <- rowSums(x*x)/2
  if (is(x,"matrix") && is.null(y)){
    res <- crossprod(t(x))
    for (i in 1:n){
      xx <- sigma*sqrt(round(-2*(res[i,] - dota - rep(dota[i],n)),9))
      res[i,] <- besselJ(xx,nu)*(xx^(-nu))
      res[i,which(xx<10e-5)] <- lim
    }
    return(as.kernelMatrix((res/lim)^ni))
  }
  if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    m <- dim(y)[1]
    dotb <- rowSums(y*y)/2
    res <- x%*%t(y)
    for( i in 1:m){
      xx <- sigma*sqrt(round(-2*(res[,i] - dota - rep(dotb[i],n)),9))
      res[,i] <- besselJ(xx,nu)*(xx^(-nu))
      res[which(xx<10e-5),i] <- lim
    }     
    return(as.kernelMatrix((res/lim)^ni))
  }
}
setMethod("kernelMatrix",signature(kernel="besselkernel"),kernelMatrix.besselkernel)


kernelMatrix.anovakernel <- function(kernel, x, y = NULL)
{
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  degree = kpar(kernel)$degree
  n <- dim(x)[1]
  if (is(x,"matrix") && is.null(y)){
    a <- matrix(0,  dim(x)[2], n)
    res <- matrix(0, n ,n)
    for (i in 1:n)
      {
        a[rep(TRUE,dim(x)[2]), rep(TRUE,n)] <- x[i,]
        res[i,]<- colSums(exp( - sigma*(a - t(x))^2))^degree
      }
    return(as.kernelMatrix(res))
  }
  if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    
    m <- dim(y)[1]
    b <- matrix(0, dim(x)[2],m)
    res <- matrix(0, dim(x)[1],m)
    for( i in 1:n)
      {
        b[rep(TRUE,dim(x)[2]), rep(TRUE,m)] <- x[i,]
        res[i,]<- colSums(exp( - sigma*(b - t(y))^2))^degree
      }
    return(as.kernelMatrix(res))
  }
}
setMethod("kernelMatrix",signature(kernel="anovakernel"),kernelMatrix.anovakernel)


kernelMatrix.polykernel <- function(kernel, x, y = NULL)
{
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix")
  scale = kpar(kernel)$scale
  offset = kpar(kernel)$offset
  degree = kpar(kernel)$degree
  if (is(x,"matrix") && is.null(y))
    {
      res <- (scale*crossprod(t(x))+offset)^degree
      return(as.kernelMatrix(res))
    }
      if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    res <- (scale*crossprod(t(x),t(y)) + offset)^degree
    return(as.kernelMatrix(res))
  }
}
setMethod("kernelMatrix",signature(kernel="polykernel"),kernelMatrix.polykernel)

kernelMatrix.vanilla <- function(kernel, x, y = NULL)
{
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix")
  if (is(x,"matrix") && is.null(y)){
    res <- crossprod(t(x))
    return(as.kernelMatrix(res))
  }
  if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    res <- crossprod(t(x),t(y))
    return(as.kernelMatrix(res))
  }
}
setMethod("kernelMatrix",signature(kernel="vanillakernel"),kernelMatrix.vanilla)

kernelMatrix.tanhkernel <- function(kernel, x, y = NULL)
{
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix")
  if (is(x,"matrix") && is.null(y)){
    scale = kpar(kernel)$scale
    offset = kpar(kernel)$offset
    res <- tanh(scale*crossprod(t(x)) + offset)
    return(as.kernelMatrix(res))
  }
  if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    res <- tanh(scale*crossprod(t(x),t(y)) + offset)
    return(as.kernelMatrix(res))
  }
}
setMethod("kernelMatrix",signature(kernel="tanhkernel"),kernelMatrix.tanhkernel)


kernelMatrix.splinekernel <- function(kernel, x, y = NULL)
{
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  degree = kpar(kernel)$degree
  n <- dim(x)[1]
  if (is(x,"matrix") && is.null(y)){
    a <- matrix(0,  dim(x)[2], n)
    res <- matrix(0, n ,n)
	x <- t(x)
    for (i in 1:n)
      {
	dr <- 	x + x[,i]	
	dp <-   x * x[,i]
	dm <-   pmin(x,x[,i])
	res[i,] <- apply((1 + dp + dp*dm - (dr/2)*dm^2 + (dm^3)/3),2, prod) 
      }
        return(as.kernelMatrix(res))
  }
  if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    
    m <- dim(y)[1]
    b <- matrix(0, dim(x)[2],m)
    res <- matrix(0, dim(x)[1],m)
	x <- t(x)
	y <- t(y)
    for( i in 1:n)
      {
	dr <- 	y + x[,i]	
	dp <-   y * x[,i]
	dm <-   pmin(y,x[,i])
	res[i,] <- apply((1 + dp + dp*dm - (dr/2)*dm^2 + (dm^3)/3),2, prod) 
   	}
    return(as.kernelMatrix(res))
  }
}
setMethod("kernelMatrix",signature(kernel="splinekernel"),kernelMatrix.splinekernel)

kernelMatrix.stringkernel <- function(kernel, x, y=NULL)
{ 
  n <- length(x)
  res1 <- matrix(rep(0,n*n), ncol = n)
  normalized = kpar(kernel)$normalized

  if(is(x,"list"))
    x <- sapply(x,paste,collapse="")
  
  if(is(y,"list"))
    y <- sapply(y,paste,collapse="")

  if (kpar(kernel)$type == "sequence" |kpar(kernel)$type == "string"|kpar(kernel)$type == "fullstring")
    {
      resdiag <- rep(0,n)
      if(normalized == TRUE)
        kernel <- stringdot(length = kpar(kernel)$length, type = kpar(kernel)$type, lambda = kpar(kernel)$lambda, normalized = FALSE)
      ## y is null   
      
      if(is.null(y)){
        if(normalized == TRUE){
          ## calculate diagonal elements first, and use them to normalize
          for (i in 1:n)
            resdiag[i] <- kernel(x[[i]],x[[i]])
      
          for(i in 1:n) {
            for(j in (i:n)[-1]) {
              res1[i,j]  <- kernel(x[[i]],x[[j]])/sqrt(resdiag[i]*resdiag[j])
            }
          }
          res1 <- res1 + t(res1)
          diag(res1) <- rep(1,n)
        }
        else{
          for (i in 1:n)
            resdiag[i] <- kernel(x[[i]],x[[i]])
          
          for(i in 1:n) {
            for(j in (i:n)[-1]) {
              res1[i,j]  <- kernel(x[[i]],x[[j]])
            }
          }
          res1 <- res1 + t(res1)
          diag(res1) <- resdiag
        }
      }
      
      if (!is.null(y)){
        m <- length(y)
        res1 <- matrix(0,n,m)
        resdiag1 <- rep(0,m)
        if(normalized == TRUE){     
          for(i in 1:n)
            resdiag[i] <- kernel(x[[i]],x[[i]])
          
          for(i in 1:m)
            resdiag1[i] <- kernel(y[[i]],y[[i]])
          
          for(i in 1:n) {
            for(j in 1:m) {
              res1[i,j] <- kernel(x[[i]],y[[j]])/sqrt(resdiag[i]*resdiag1[j])
            }
          }
        }
        else{
          for(i in 1:n) {
            for(j in 1:m) {
              res1[i,j] <- kernel(x[[i]],y[[j]])
            }
          }
        }
      }
      return(as.kernelMatrix(res1))
    }
  else {
  
      switch(kpar(kernel)$type,
             "exponential" = 
             sktype <- 2,
             "constant" =
             sktype <- 1,
             "spectrum" =
             sktype <- 3,
             "boundrange" =
             sktype <- 4)

        if(sktype==3 &(any(nchar(x) < kpar(kernel)$length)|any(nchar(x) < kpar(kernel)$length)))
          stop("spectral kernel does not accept strings shorter than the length parameter")
     
        if(is(x,"list"))
          x <- unlist(x)
        if(is(y,"list"))
          y <- unlist(y)

      
      x <- paste(x,"\n",sep="")
      if(!is.null(y))
        y <- paste(y,"\n",sep="")
        
      if(is.null(y))
          ret <- matrix(0, length(x),length(x))
        else
          ret <- matrix(0,length(x),length(y))

          if(is.null(y)){
            for(i in 1:length(x))
              ret[i,i:length(x)] <- .Call("stringtv",as.character(x[i]),as.character(x[i:length(x)]),as.integer(length(x) - i + 1),as.integer(nchar(x[i])),as.integer(nchar(x[i:length(x)])),as.integer(sktype),as.double(kpar(kernel)$lambda))
            ret <- ret + t(ret)
            diag(ret) <- diag(ret)/2
          }
          else
            for(i in 1:length(x))
              ret[i,] <- .Call("stringtv",as.character(x[i]),as.character(y),as.integer(length(y)),as.integer(nchar(x[i])),as.integer(nchar(y)),as.integer(sktype),as.double(kpar(kernel)$lambda))

        if(normalized == TRUE){
          if(is.null(y))
            ret <- t((1/sqrt(diag(ret)))*t(ret*(1/sqrt(diag(ret)))))
           else{
             norm1 <- rep(0,length(x))
             norm2 <- rep(0,length(y))
             for( i in 1:length(x))
               norm1[i] <- .Call("stringtv",as.character(x[i]),as.character(x[i]),as.integer(1),as.integer(nchar(x[i])),as.integer(nchar(x[i])),as.integer(sktype),as.double(kpar(kernel)$lambda))
             for( i in 1:length(y))
               norm2[i] <- .Call("stringtv",as.character(y[i]),as.character(y[i]),as.integer(1),as.integer(nchar(y[i])),as.integer(nchar(y[i])),as.integer(sktype),as.double(kpar(kernel)$lambda))
             ret <- t((1/sqrt(norm2))*t(ret*(1/sqrt(norm1))))
           }
        }
    }
  return(as.kernelMatrix(ret))
}

setMethod("kernelMatrix",signature(kernel="stringkernel"),kernelMatrix.stringkernel)



## kernelMult computes kernel matrix  - vector product
## function computing <x,x'> * z (<x,x'> %*% z)


kernelMult <- function(kernel, x, y=NULL, z, blocksize = 128)
{
#  if(is.function(kernel)) ker <- deparse(substitute(kernel))
#      kernel <- do.call(kernel, kpar)

  if(!is(x,"matrix")) stop("x must be a matrix")
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must ba a matrix or a vector")
  n <- nrow(x)

  if(is.null(y))
    {
      ## check if z,x match
      z <- as.matrix(z)
      if(is.null(y)&&!dim(z)[1]==n)
        stop("z columns/length do not match x columns")
      
      res1 <- matrix(rep(0,n*n), ncol = n)     
      
      for(i in 1:n)
        {
          for(j in i:n)
            {
              res1[j,i] <- kernel(x[i,],x[j,])
            }
        }
      res1 <- res1 + t(res1)
      diag(res1) <- diag(res1)/2
    }
  if (is(y,"matrix"))
    {
      
      m <- dim(y)[1]
      z <- as.matrix(z)

      if(!dim(z)[1] == m) stop("z has wrong dimension")
      res1 <- matrix(rep.int(0,m*n),ncol=m)
      for(i in 1:n)
        {
          for(j in 1:m)
            {
              res1[i,j] <- kernel(x[i,],y[j,])
            }
        }
    }
  return(res1%*%z)
}

setGeneric("kernelMult", function(kernel, x, y=NULL, z, blocksize = 256) standardGeneric("kernelMult"))

kernelMult.character <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  return(x%*%z)
}
setMethod("kernelMult",signature(kernel="character", x="kernelMatrix"),kernelMult.character)
  


kernelMult.rbfkernel <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  if(!is(y,"matrix")&&!is.null(y)&&!is(y,"vector")) stop("y must be a matrix or a vector")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
  if(is(x,"vector"))
    x <- as.matrix(x)
   if(is(y,"vector"))
    y <- as.matrix(y)
  sigma <- kpar(kernel)$sigma
  n <- dim(x)[1]
  m <- dim(x)[2]
  nblocks <- floor(n/blocksize)
  lowerl <- 1
  upperl <- 0
  dota <- as.matrix(rowSums(x^2))

  if (is.null(y))
    {
      z <- as.matrix(z)
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")

      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      if(nblocks > 0)
        {
          dotab <- rep(1,blocksize)%*%t(dota)
          for(i in 1:nblocks)
            {
              upperl = upperl + blocksize
              res[lowerl:upperl,] <- exp(sigma*(2*x[lowerl:upperl,]%*%t(x) - dotab - dota[lowerl:upperl]%*%t(rep.int(1,n))))%*%z
              lowerl <- upperl + 1
            }
        }
      if(lowerl <= n)
        res[lowerl:n,] <- exp(sigma*(2*x[lowerl:n,]%*%t(x) - rep.int(1,n+1-lowerl)%*%t(dota) - dota[lowerl:n]%*%t(rep.int(1,n))))%*%z

    }
  if(is(y,"matrix"))
    {
      n2 <- dim(y)[1]
      z <- as.matrix(z)
      
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      dotb <- as.matrix(rowSums(y*y))
      
       if(nblocks > 0)
         {
           dotbb <-  rep(1,blocksize)%*%t(dotb)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            res[lowerl:upperl,] <- exp(sigma*(2*x[lowerl:upperl,]%*%t(y) - dotbb - dota[lowerl:upperl]%*%t(rep.int(1,n2))))%*%z
            lowerl <- upperl + 1
          }
      }
      if(lowerl <= n)
        res[lowerl:n,] <- exp(sigma*(2*x[lowerl:n,]%*%t(y) - rep.int(1,n+1-lowerl)%*%t(dotb) - dota[lowerl:n]%*%t(rep.int(1,n2))))%*%z
    }
  return(res)
}  
setMethod("kernelMult",signature(kernel="rbfkernel"),kernelMult.rbfkernel)


kernelMult.laplacekernel <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  if(!is(y,"matrix")&&!is.null(y)&&!is(y,"vector")) stop("y must be a matrix or a vector")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
   if(is(x,"vector"))
    x <- as.matrix(x)
   if(is(y,"vector"))
    y <- as.matrix(y)
  sigma <- kpar(kernel)$sigma
  n <- dim(x)[1]
  m <- dim(x)[2]
  nblocks <- floor(n/blocksize)
  lowerl <- 1
  upperl <- 0
  dota <- as.matrix(rowSums(x^2))

  if (is.null(y))
    {
      z <- as.matrix(z)
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      if(nblocks > 0)
        {
          dotab <- rep(1,blocksize)%*%t(dota)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            res[lowerl:upperl,] <- exp(-sigma*sqrt(-round(2*x[lowerl:upperl,]%*%t(x) - dotab - dota[lowerl:upperl]%*%t(rep.int(1,n)),9)))%*%z
            lowerl <- upperl + 1
          
          }
      }
      if(lowerl <= n)
        res[lowerl:n,] <- exp(-sigma*sqrt(-round(2*x[lowerl:n,]%*%t(x) - rep.int(1,n+1-lowerl)%*%t(dota) - dota[lowerl:n]%*%t(rep.int(1,n)),9)))%*%z

    }
  if(is(y,"matrix"))
    {
      n2 <- dim(y)[1]
      z <- as.matrix(z)
      
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      dotb <- as.matrix(rowSums(y*y))
      
       if(nblocks > 0)
         {
           dotbb <-  rep(1,blocksize)%*%t(dotb)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            res[lowerl:upperl,] <- exp(-sigma*sqrt(-round(2*x[lowerl:upperl,]%*%t(y) - dotbb - dota[lowerl:upperl]%*%t(rep.int(1,n2)),9)))%*%z
            lowerl <- upperl + 1
          }
         }
      if(lowerl <= n)
        res[lowerl:n,] <- exp(-sigma*sqrt(-round(2*x[lowerl:n,]%*%t(y) - rep.int(1,n+1-lowerl)%*%t(dotb) - dota[lowerl:n]%*%t(rep.int(1,n2)),9)))%*%z
    }
  return(res)
}  
setMethod("kernelMult",signature(kernel="laplacekernel"),kernelMult.laplacekernel)



kernelMult.besselkernel <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  if(!is(y,"matrix")&&!is.null(y)&&!is(y,"vector")) stop("y must be a matrix")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
   if(is(x,"vector"))
    x <- as.matrix(x)
   if(is(y,"vector"))
    y <- as.matrix(y)
  sigma <- kpar(kernel)$sigma
  nu <- kpar(kernel)$order
  ni <- kpar(kernel)$degree
  n <- dim(x)[1]
  m <- dim(x)[2]
  nblocks <- floor(n/blocksize)
  lowerl <- 1
  upperl <- 0
  lim <- 1/(gamma(nu+1)*2^(nu))
  dota <- as.matrix(rowSums(x^2))

  if (is.null(y))
    {
      z <- as.matrix(z)
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      if(nblocks > 0)
        {
          dotab <- rep(1,blocksize)%*%t(dota)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            
            xx <- sigma*sqrt(-round(2*x[lowerl:upperl,]%*%t(x) - dotab - dota[lowerl:upperl]%*%t(rep.int(1,n)),9))
            res1 <- besselJ(xx,nu)*(xx^(-nu))
            res1[which(xx<10e-5)] <- lim
           
            res[lowerl:upperl,] <- ((res1/lim)^ni)%*%z
            lowerl <- upperl + 1
          }
      }
      if(lowerl <= n)
        {
          xx <- sigma*sqrt(-round(2*x[lowerl:n,]%*%t(x) - rep.int(1,n+1-lowerl)%*%t(dota) - dota[lowerl:n]%*%t(rep.int(1,n)),9))
          res1 <- besselJ(xx,nu)*(xx^(-nu))
          res1[which(xx<10e-5)] <- lim
          res[lowerl:n,] <- ((res1/lim)^ni)%*%z
      }
    }
  if(is(y,"matrix"))
    {
      n2 <- dim(y)[1]
      z <- as.matrix(z)
      
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      dotb <- as.matrix(rowSums(y*y))
      
       if(nblocks > 0)
         {
           dotbb <-  rep(1,blocksize)%*%t(dotb)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            xx <- sigma*sqrt(-round(2*x[lowerl:upperl,]%*%t(y) - dotbb - dota[lowerl:upperl]%*%t(rep.int(1,n2)),9))
            res1 <- besselJ(xx,nu)*(xx^(-nu))
            res1[which(xx < 10e-5)] <- lim

            res[lowerl:upperl,] <- ((res1/lim)^ni)%*%z
            lowerl <- upperl + 1
          }
      }
      if(lowerl <= n)
      {
        xx <- sigma*sqrt(-round(2*x[lowerl:n,]%*%t(y) - rep.int(1,n+1-lowerl)%*%t(dotb) - dota[lowerl:n]%*%t(rep.int(1,n2)),9))
        res1 <- besselJ(xx,nu)*(xx^(-nu))
        res1[which(xx < 10e-5)] <- lim
        res[lowerl:n,] <- ((res1/lim)^ni)%*%z 
      }
    }
  return(res)
}  
setMethod("kernelMult",signature(kernel="besselkernel"),kernelMult.besselkernel)

kernelMult.anovakernel <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  if(!is(y,"matrix")&&!is.null(y)&&!is(y,"vector")) stop("y must be a matrix or a vector")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  sigma <- kpar(kernel)$sigma
  degree <- kpar(kernel)$degree
  n <- dim(x)[1]
  m <- dim(x)[2]
  nblocks <- floor(n/blocksize)
  lowerl <- 1
  upperl <- 0
 

  if (is.null(y))
    {
      z <- as.matrix(z)
      
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      
      if(nblocks > 0)
        {
          a <- matrix(0,m,blocksize)
          re <- matrix(0, n, blocksize)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            for(j in 1:n)
              {
                a[rep(TRUE,m),rep(TRUE,blocksize)] <- x[j,]
                re[j,] <-  colSums(exp( - sigma*(a - t(x[lowerl:upperl,]))^2))^degree
              }
            res[lowerl:upperl,] <- t(re)%*%z
            lowerl <- upperl + 1
          
          }
        }
      if(lowerl <= n){
        a <- matrix(0,m,n-lowerl+1)
        re <- matrix(0,n,n-lowerl+1)
        for(j in 1:n)
          {
            a[rep(TRUE,m),rep(TRUE,n-lowerl+1)] <- x[j,]
            re[j,] <-  colSums(exp( - sigma*(a - t(x[lowerl:n,,drop=FALSE]))^2))^degree
          }
        res[lowerl:n,] <- t(re)%*%z
      }
    }
  if(is(y,"matrix"))
    {
      n2 <- dim(y)[1]
      nblocks <- floor(n2/blocksize)

      z <- as.matrix(z)
      
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])

      if(nblocks > 0)
        {
          b <- matrix(0, m, blocksize)
          re <- matrix(0, n, blocksize)
          for(i in 1:nblocks)
            {
              upperl = upperl + blocksize
              for(j in 1:n)
                {
                  b[rep(TRUE,dim(x)[2]), rep(TRUE,blocksize)] <- x[j,]
                  re[j,]<- colSums(exp( - sigma*(b - t(y[lowerl:upperl,]))^2)^degree)
                }
              res[,1] <- res[,1] + re %*%z[lowerl:upperl,]
              lowerl <- upperl + 1
            }
        }
      if(lowerl <= n)
        {
          b <- matrix(0, dim(x)[2], n2-lowerl+1)
          re <- matrix(0, n, n2-lowerl+1)
          for( i in 1:n)
            {
              b[rep(TRUE,dim(x)[2]),rep(TRUE,n2-lowerl+1)] <- x[i,]
              re[i,]<- colSums(exp( - sigma*(b - t(y[lowerl:n2,,drop=FALSE]))^2)^degree)
            }
       
          res[,1] <- res[,1] + re%*%z[lowerl:n2]
        }
    }
  return(res)
}  
setMethod("kernelMult",signature(kernel="anovakernel"),kernelMult.anovakernel)



kernelMult.splinekernel <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  if(!is(y,"matrix")&&!is.null(y)&&!is(y,"vector")) stop("y must be a matrix or a vector")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
  n <- dim(x)[1]
  m <- dim(x)[2]
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  nblocks <- floor(n/blocksize)
  lowerl <- 1
  upperl <- 0

  if (is.null(y))
    {
      z <- as.matrix(z)
      
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      	x <- t(x)
      if(nblocks > 0)
        {
          re <- matrix(0, dim(z)[1], blocksize)
          for(i in 1:nblocks)
            {
              upperl = upperl + blocksize
              
              for (j in lowerl:upperl)
                {	
                  dr <-  x + x[ , j]	
                  dp <-  x * x[ , j]
                  dm <-  pmin(x,x[,j])
                  re[,j-(i-1)*blocksize] <- apply((1 + dp + dp*dm - (dr/2)*dm^2 + (dm^3)/3),2, prod) 
                }
              res[lowerl:upperl,] <- crossprod(re,z)
              lowerl <- upperl + 1
            }
        }
      if(lowerl <= n){
        a <- matrix(0,m,n-lowerl+1)
        re <- matrix(0,dim(z)[1],n-lowerl+1)
        for(j in lowerl:(n-lowerl+1))
          {
            dr <- x + x[ , j]	
            dp <- x * x[ , j]
            dm <-  pmin(x,x[,j])
            re[,j-nblocks*blocksize] <- apply((1 + dp + dp*dm - (dr/2)*dm^2 + (dm^3)/3),2, prod) 
          }
        res[lowerl:n,] <- crossprod(re,z)
      }
    }
  if(is(y,"matrix"))
    {
      n2 <- dim(y)[1]
      nblocks <- floor(n2/blocksize)
      z <- as.matrix(z)
      
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
	x <- t(x)
	y <- t(y)
      if(nblocks > 0)
        {
          re <- matrix(0, dim(z)[1], blocksize)
          for(i in 1:nblocks)
            {
              upperl = upperl + blocksize
             
                for(j in lowerl:upperl)
                {
		dr <- y + x[ , j]	
		dp <- y * x[ , j]
		dm <- pmin(y,x[,j])
		re[,j-(i-1)*blocksize] <- apply((1 + dp + dp*dm - (dr/2)*dm^2 + (dm^3)/3),2, prod) 
           	}
              res[lowerl:upperl] <- crossprod(re, z)
              lowerl <- upperl + 1
            }
        }
      if(lowerl <= n)
        {
          b <- matrix(0, dim(x)[2], n-lowerl+1)
          re <- matrix(0, dim(z)[1], n-lowerl+1)
          for(j in lowerl:(n-lowerl+1))
            {
        	dr <- y + x[, j]	
		dp <- y * x[, j]
		dm <-  pmin(y,x[,j])
		re[,j-nblocks*blocksize] <- apply((1 + dp + dp*dm - (dr/2)*dm^2 + (dm^3)/3),2, prod) 
              }
       	  res[lowerl:n] <-  crossprod(re, z)
        }
    }
  return(res)
}  
setMethod("kernelMult",signature(kernel="splinekernel"),kernelMult.splinekernel)


kernelMult.polykernel <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  if(!is(y,"matrix")&&!is.null(y)&&!is(y,"vector")) stop("y must be a matrix")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  degree <- kpar(kernel)$degree
  scale <- kpar(kernel)$scale
  offset <- kpar(kernel)$offset
  n <- dim(x)[1]
  m <- dim(x)[2]
  nblocks <- floor(n/blocksize)
  lowerl <- 1
  upperl <- 0
  if (is.null(y))
    {
      z <- as.matrix(z)
      
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      if(nblocks > 0)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            res[lowerl:upperl,] <- ((scale*x[lowerl:upperl,]%*%t(x) + offset)^degree) %*% z
            lowerl <- upperl + 1
          }
      if(lowerl <= n)
        res[lowerl:n,] <- ((scale*x[lowerl:n,]%*%t(x) +offset)^degree)%*%z
    }
  if(is(y,"matrix"))
    {
      n2 <- dim(y)[1]
      z <- as.matrix(z)
      
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
 
       if(nblocks > 0)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            res[lowerl:upperl,] <- ((scale*x[lowerl:upperl,]%*%t(y) + offset)^degree)%*%z
            lowerl <- upperl + 1
          }
      if(lowerl <= n)
        res[lowerl:n,] <- ((scale*x[lowerl:n,]%*%t(y) + offset)^degree)%*%z
    }
  return(res)
} 
setMethod("kernelMult",signature(kernel="polykernel"),kernelMult.polykernel)


kernelMult.tanhkernel <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  if(!is(y,"matrix")&&!is.null(y)&&!is(y,"vector")) stop("y must be a matrix or a vector")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  scale <- kpar(kernel)$scale
  offset <- kpar(kernel)$offset
  n <- dim(x)[1]
  m <- dim(x)[2]
  nblocks <- floor(n/blocksize)
  lowerl <- 1
  upperl <- 0
 

  if (is.null(y))
    {
      z <- as.matrix(z)
      
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      if(nblocks > 0)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            res[lowerl:upperl,] <- tanh(scale*x[lowerl:upperl,]%*%t(x) + offset) %*% z
            lowerl <- upperl + 1
          }
      if(lowerl <= n)
        res[lowerl:n,] <- tanh(scale*x[lowerl:n,]%*%t(x) +offset)%*%z
    }
  if(is(y,"matrix"))
    {
      n2 <- dim(y)[1]
      z <- as.matrix(z)
      
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
 
       if(nblocks > 0)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            res[lowerl:upperl,] <- tanh(scale*x[lowerl:upperl,]%*%t(y) + offset)%*%z
            lowerl <- upperl + 1
          }
      if(lowerl <= n)
        res[lowerl:n,] <- tanh(scale*x[lowerl:n,]%*%t(y) + offset)%*%z
    }
  return(res)
} 
setMethod("kernelMult",signature(kernel="tanhkernel"),kernelMult.tanhkernel)


kernelMult.vanillakernel <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  if(!is(y,"matrix")&&!is.null(y)&&!is(y,"vector")) stop("y must be a matrix or vector")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
  n <- dim(x)[1]
  m <- dim(x)[2]
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if (is.null(y))
    {
      z <- as.matrix(z)
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")
      res <- t(crossprod(crossprod(x,z),t(x)))
    }
  if(is(y,"matrix"))
    {
      n2 <- dim(y)[1]
      z <- as.matrix(z)
      
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- t(crossprod(crossprod(y,z),t(x)))
    }
  return(res)
} 
setMethod("kernelMult",signature(kernel="vanillakernel"),kernelMult.vanillakernel)


kernelMult.stringkernel <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
  normalized = kpar(kernel)$normalized
      
  n <- length(x)
  res1 <- matrix(rep(0,n*n), ncol = n)
  resdiag <- rep(0,n)
  
  if(is(x,"list"))
    x <- sapply(x,paste,collapse="")
  
  if(is(y,"list"))
    y <- sapply(y,paste,collapse="")

  
  if (kpar(kernel)$type == "sequence" |kpar(kernel)$type == "string"|kpar(kernel)$type == "fullstring")
    {
      if(normalized == TRUE)
        kernel <- stringdot(length = kpar(kernel)$length, type = kpar(kernel)$type, lambda = kpar(kernel)$lambda, normalized = FALSE)
      
      ## y is null 
      if(is.null(y)){
        if(normalized == TRUE){
          z <- as.matrix(z)
          if(dim(z)[1]!= n) stop("z rows must be equal to x length")
          dz <- dim(z)[2]
          vres <- matrix(0,n,dz)
          ## calculate diagonal elements first, and use them to normalize
          for (i in 1:n)
            resdiag[i] <- kernel(x[[i]],x[[i]])
          
          for(i in 1:n) {
            for(j in (i:n)[-1]) {
              res1[i,j]  <- kernel(x[[i]],x[[j]])/sqrt(resdiag[i]*resdiag[j])
            }
          }
          res1 <- res1 + t(res1)
          diag(res1) <- rep(1,n)
          vres  <- res1 %*% z
        }
        else{
          z <- as.matrix(z)
          
          if(dim(z)[1]!= n) stop("z rows must be equal to x length")
          dz <- dim(z)[2]
          vres <- matrix(0,n,dz)
          ## calculate diagonal elements first, and use them to normalize
          for (i in 1:n)
            resdiag[i] <- kernel(x[[i]],x[[i]])
          
          for(i in 1:n) {
            for(j in (i:n)[-1]) {
              res1[i,j]  <- kernel(x[[i]],x[[j]])
            }
          }
          res1 <- res1 + t(res1)
          diag(res1) <- resdiag
          vres  <- res1 %*% z
        }
      }
      
      if (!is.null(y)){
        if(normalized == TRUE){
          nblocks <- floor(n/blocksize)
          lowerl <- 1
          upperl <- 0
          m <- length(y)
          z <- as.matrix(z)
          if(dim(z)[1]!= m) stop("z rows must be equal to y length")
          resdiag1 <- rep(0,m)
          dz <- dim(z)[2]
          vres <- matrix(0,n,dz)
          
          for(i in 1:n)
            resdiag[i] <- kernel(x[[i]],x[[i]])
          
          for(i in 1:m)
            resdiag1[i] <- kernel(y[[i]],y[[i]])
          
          if (nblocks > 0){
            res1 <- matrix(0,blocksize,m)
            for(k in 1:nblocks){
              upperl <- upperl + blocksize
              for(i in lowerl:(upperl)) {
                for(j in 1:m) {
                  res1[i - (k-1)*blocksize,j] <- kernel(x[[i]],y[[j]])/sqrt(resdiag[i]*resdiag1[j])
                }
              }
              vres[lowerl:upperl,] <- res1 %*% z
              lowerl <- upperl +1
            }
          }
          if(lowerl <= n)
            {
              res1 <- matrix(0,n-lowerl+1,m)
              for(i in lowerl:n) {
                for(j in 1:m) {
                  res1[i - nblocks*blocksize,j] <- kernel(x[[i]],y[[j]])/sqrt(resdiag[i]*resdiag1[j])
                }
              }
              vres[lowerl:n,] <- res1  %*% z
            }
        }
        else
          {
            nblocks <- floor(n/blocksize)
            lowerl <- 1
            upperl <- 0
            m <- length(y)
            z <- as.matrix(z)
            if(dim(z)[1]!= m) stop("z rows must be equal to y length")
            dz <- dim(z)[2]
            vres <- matrix(0,n,dz)
            
            if (nblocks > 0){
              res1 <- matrix(0,blocksize,m)
              
              for(k in 1:nblocks){
                
                upperl <- upperl + blocksize
                
                for(i in lowerl:(upperl)) {
                  for(j in 1:m) {
                    res1[i - (k-1)*blocksize, j] <- kernel(x[[i]],y[[j]])
                  }
                }
                vres[lowerl:upperl,] <- res1 %*% z
                lowerl <- upperl +1
              }
            }
            if(lowerl <= n)
              {
                res1 <- matrix(0,n-lowerl+1,m)
                for(i in lowerl:n) {
                  for(j in 1:m) {
                    res1[i - nblocks*blocksize,j] <- kernel(x[[i]],y[[j]])
                  }
                }
                vres[lowerl:n,] <- res1  %*% z
              }
          }
      }
    }
  else
    {
      switch(kpar(kernel)$type,
             "exponential" = 
             sktype <- 2,
             "constant" =
             sktype <- 1,
             "spectrum" =
             sktype <- 3,
             "boundrange" =
             sktype <- 4)

      if(sktype==3 &(any(nchar(x) < kpar(kernel)$length)|any(nchar(x) < kpar(kernel)$length)))
        stop("spectral kernel does not accept strings shorter than the length parameter")

      
      x <- paste(x,"\n",sep="")
      if(!is.null(y))
        y <- paste(y,"\n",sep="")
    
      
      ## y is null 
      if(is.null(y)){
        if(normalized == TRUE){
          nblocks <- floor(n/blocksize)
          lowerl <- 1
          upperl <- 0
          z <- as.matrix(z)
          if(dim(z)[1]!= n) stop("z rows must be equal to y length")
          dz <- dim(z)[2]
          vres <- matrix(0,n,dz)

          for (i in 1:n)
            resdiag[i] <- .Call("stringtv",as.character(x[i]),as.character(x[i]),as.integer(1),as.integer(nchar(x[i])),as.integer(nchar(x[i])),as.integer(sktype),as.double(kpar(kernel)$lambda))

          if (nblocks > 0){
            res1 <- matrix(0,blocksize,n)
            for(k in 1:nblocks){
              upperl <- upperl + blocksize
              for(i in lowerl:(upperl)) {
                res1[i - (k-1)*blocksize, ] <-  .Call("stringtv",as.character(x[i]),as.character(x),as.integer(length(x)),as.integer(nchar(x[i])),as.integer(nchar(x)),as.integer(sktype),as.double(kpar(kernel)$lambda))/sqrt(resdiag[i]*resdiag)
              }
              vres[lowerl:upperl,] <- res1 %*% z
              lowerl <- upperl +1
            }
          }
          if(lowerl <= n)
            {
              res1 <- matrix(0,n-lowerl+1,n)
              for(i in lowerl:n) {
                res1[i - nblocks*blocksize,] <- .Call("stringtv",as.character(x[i]),as.character(x),as.integer(length(x)),as.integer(nchar(x[i])),as.integer(nchar(x)),as.integer(sktype),as.double(kpar(kernel)$lambda))/sqrt(resdiag[i]*resdiag)
              }
              vres[lowerl:n,] <- res1  %*% z
            }
        }
        else
          {
            nblocks <- floor(n/blocksize)
            lowerl <- 1
            upperl <- 0
            z <- as.matrix(z)
            if(dim(z)[1]!= n) stop("z rows must be equal to y length")
            dz <- dim(z)[2]
            vres <- matrix(0,n,dz)
            
            if (nblocks > 0){
              res1 <- matrix(0,blocksize,n)
              for(k in 1:nblocks){
                upperl <- upperl + blocksize
                for(i in lowerl:(upperl)) {
                  res1[i - (k-1)*blocksize, ] <- .Call("stringtv",as.character(x[i]),as.character(x),as.integer(length(x)),as.integer(nchar(x[i])),as.integer(nchar(x)),as.integer(sktype),as.double(kpar(kernel)$lambda))
                }
                vres[lowerl:upperl,] <- res1 %*% z
                lowerl <- upperl +1
              }
            }
            if(lowerl <= n)
              {
                res1 <- matrix(0,n-lowerl+1,n)
                for(i in lowerl:n) {
                  res1[i - nblocks*blocksize,] <- .Call("stringtv",as.character(x[i]),as.character(x),as.integer(length(x)),as.integer(nchar(x[i])),as.integer(nchar(x)),as.integer(sktype),as.double(kpar(kernel)$lambda))
                }
                vres[lowerl:n,] <- res1  %*% z
              }
          }
      }
      
      if (!is.null(y)){
        if(normalized == TRUE){
          nblocks <- floor(n/blocksize)
          lowerl <- 1
          upperl <- 0
          m <- length(y)
          z <- as.matrix(z)
          if(dim(z)[1]!= m) stop("z rows must be equal to y length")
          resdiag1 <- rep(0,m)
          dz <- dim(z)[2]
          vres <- matrix(0,n,dz)
          
          for(i in 1:n)
            resdiag[i] <- .Call("stringtv",as.character(x[i]),as.character(x[i]),as.integer(1),as.integer(nchar(x[i])),as.integer(nchar(x[i])),as.integer(sktype),as.double(kpar(kernel)$lambda))

          
          for(i in 1:m)
            resdiag1[i] <- .Call("stringtv",as.character(y[i]),as.character(y[i]),as.integer(1),as.integer(nchar(y[i])),as.integer(nchar(y[i])),as.integer(sktype),as.double(kpar(kernel)$lambda))
          
          if (nblocks > 0){
            res1 <- matrix(0,blocksize,m)
            for(k in 1:nblocks){
              upperl <- upperl + blocksize
              for(i in lowerl:(upperl)) {
                res1[i - (k-1)*blocksize, ] <-  .Call("stringtv",as.character(x[i]),as.character(y),as.integer(length(y)),as.integer(nchar(x[i])),as.integer(nchar(y)),as.integer(sktype),as.double(kpar(kernel)$lambda))/sqrt(resdiag[i]*resdiag1)
              }
              vres[lowerl:upperl,] <- res1 %*% z
              lowerl <- upperl +1
            }
          }
          if(lowerl <= n)
            {
              res1 <- matrix(0,n-lowerl+1,m)
              for(i in lowerl:n) {
                res1[i - nblocks*blocksize,] <- .Call("stringtv",as.character(x[i]),as.character(y),as.integer(length(y)),as.integer(nchar(x[i])),as.integer(nchar(y)),as.integer(sktype),as.double(kpar(kernel)$lambda))/sqrt(resdiag[i]*resdiag1)
              }
              vres[lowerl:n,] <- res1  %*% z
            }
        }
        else
          {
            nblocks <- floor(n/blocksize)
            lowerl <- 1
            upperl <- 0
            m <- length(y)
            z <- as.matrix(z)
            if(dim(z)[1]!= m) stop("z rows must be equal to y length")
            dz <- dim(z)[2]
            vres <- matrix(0,n,dz)
            
            if (nblocks > 0){
              res1 <- matrix(0,blocksize,m)
              
              for(k in 1:nblocks){
                
                upperl <- upperl + blocksize
                
                for(i in lowerl:(upperl)) {
                  res1[i - (k-1)*blocksize, ] <- .Call("stringtv",as.character(x[i]),as.character(y),as.integer(length(y)),as.integer(nchar(x[i])),as.integer(nchar(y)),as.integer(sktype),as.double(kpar(kernel)$lambda))
               
                }
                vres[lowerl:upperl,] <- res1 %*% z
                lowerl <- upperl +1
              }
            }
            if(lowerl <= n)
              {
                res1 <- matrix(0,n-lowerl+1,m)
                for(i in lowerl:n) {
                  res1[i - nblocks*blocksize,] <- .Call("stringtv",as.character(x[i]),as.character(y),as.integer(length(y)),as.integer(nchar(x[i])),as.integer(nchar(y)),as.integer(sktype),as.double(kpar(kernel)$lambda))
                  }
                vres[lowerl:n,] <- res1  %*% z
              }
          }
      }
    }
  return(vres)
  
}
setMethod("kernelMult",signature(kernel="stringkernel"),kernelMult.stringkernel)


## kernelPol return the quadratic form of a kernel matrix
## kernelPol returns the scalar product of x y componentwise with polarities
## of z and k 

kernelPol <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is(x,"matrix")) stop("x must be a matrix")
  if(!is(y,"matrix")&&!is.null(y)) stop("y must be a matrix")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must ba a matrix or a vector")
  n <- nrow(x)
  z <- as.matrix(z)


  if(!dim(z)[1]==n)
    stop("z must have the length equal to x colums")
  res1 <- matrix(rep(0,n*n), ncol = n)
  if (is.null(y))
    {
      for(i in 1:n)
        {
          for(j in i:n)
            {
              res1[i,j] <- kernel(x[i,],x[j,])*z[j]*z[i]
      }
    }
      res1 <- res1 + t(res1)
      diag(res1) <- diag(res1)/2
  }
  if (is(x,"matrix") && is(y,"matrix")){
    m <- dim(y)[1]
    if(is.null(k)) stop("k not specified!")
    k <- as.matrix(k)
    if(!dim(x)[2]==dim(y)[2])
      stop("matrixes must have the same number of columns")
    if(!dim(z)[2]==dim(k)[2])
      stop("z and k vectors must have the same number of columns")
    if(!dim(x)[1]==dim(z)[1])
      stop("z and x must have the same number of rows")
    if(!dim(y)[1]==dim(k)[1])
      stop("y and k must have the same number of rows")
    res1 <- matrix(0,dim(x)[1],dim(y)[1])
    for(i in 1:n)
      {
        for(j in 1:m)
          {
            res1[i,j] <- kernel(x[i,],y[j,])*z[i]*k[j]
          }
      }
  }
  return(res1)
}

setGeneric("kernelPol", function(kernel, x, y=NULL, z, k = NULL) standardGeneric("kernelPol"))


kernelPol.rbfkernel <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is(y,"matrix")&&!is.null(y)&&!is(y,"vector")) stop("y must be a matrix a vector or NULL")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
  if(!is(k,"matrix")&&!is(k,"vector")&&!is.null(k)) stop("k must be a matrix or a vector")
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  sigma <- kpar(kernel)$sigma
  n <- dim(x)[1]
  dota <- rowSums(x*x)/2
  z <- as.matrix(z)
  if(!dim(z)[1]==n)
    stop("z must have the length equal to x colums")
  if (is.null(y))
    {
      if(is(z,"matrix")&&!dim(z)[1]==n)
       stop("z must have size equal to x colums")
      res <- crossprod(t(x))
      for (i in 1:n)
        res[i,] <- z[i,]*(exp(2*sigma*(res[i,] - dota - rep(dota[i],n)))*z)
      return(res)
    }
  if (is(y,"matrix"))
    {
      if(is.null(k)) stop("k not specified!")
      m <- dim(y)[1]
      k <- as.matrix(k)
      if(!dim(k)[1]==m)
        stop("k must have equal rows to y")
      if(!dim(x)[2]==dim(y)[2])
        stop("matrixes must have the same number of columns")
      dotb <- rowSums(y*y)/2
      res <- x%*%t(y)
      for( i in 1:m)#2*sigma or sigma
        res[,i]<- k[i,]*(exp(2*sigma*(res[,i] - dota - rep(dotb[i],n)))*z)
      return(res)
    }
}
setMethod("kernelPol",signature(kernel="rbfkernel"),kernelPol.rbfkernel)

kernelPol.laplacekernel <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is(y,"matrix")&&!is.null(y)&&!is(y,"vector")) stop("y must be a matrix, vector or NULL")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
  if(!is(k,"matrix")&&!is(k,"vector")&&!is.null(k)) stop("k must be a matrix or a vector")
  sigma <- kpar(kernel)$sigma
 if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  n <- dim(x)[1]
  dota <- rowSums(x*x)/2
  z <- as.matrix(z)
  if(!dim(z)[1]==n)
    stop("z must have the length equal to x colums")
  if (is.null(y))
    {
      if(is(z,"matrix")&&!dim(z)[1]==n)
       stop("z must have size equal to x colums")
      res <- crossprod(t(x))
      for (i in 1:n)
        res[i,] <- z[i,]*(exp(-sigma*sqrt(-round(2*(res[i,] - dota - rep(dota[i],n)),9)))*z)
      return(res)
    }
  if (is(y,"matrix"))
    {
      if(is.null(k)) stop("k not specified!")
      m <- dim(y)[1]
      k <- as.matrix(k)
      if(!dim(k)[1]==m)
        stop("k must have equal rows to y")
      if(!dim(x)[2]==dim(y)[2])
        stop("matrixes must have the same number of columns")
      dotb <- rowSums(y*y)/2
      res <- x%*%t(y)
      for( i in 1:m)#2*sigma or sigma
        res[,i]<- k[i,]*(exp(-sigma*sqrt(-round(2*(res[,i] - dota - rep(dotb[i],n)),9)))*z)
      return(res)
    }
}
setMethod("kernelPol",signature(kernel="laplacekernel"),kernelPol.laplacekernel)


kernelPol.besselkernel <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is(y,"matrix")&&!is.null(y)&&!is(y,"vector")) stop("y must be a matrix or NULL")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
  if(!is(k,"matrix")&&!is(k,"vector")&&!is.null(k)) stop("k must be a matrix or a vector")
  sigma <- kpar(kernel)$sigma
   if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  nu <- kpar(kernel)$order
  ni <- kpar(kernel)$degree
  n <- dim(x)[1]
  lim <- 1/(gamma(nu + 1)*2^nu)
  dota <- rowSums(x*x)/2
  z <- as.matrix(z)

  if(!dim(z)[1]==n)
    stop("z must have the length equal to x colums")
  if (is.null(y))
    {
      if(is(z,"matrix")&&!dim(z)[1]==n)
       stop("z must have size equal to x colums")
      res <- crossprod(t(x))
      for (i in 1:n)
        {
        xx <- sigma*sqrt(-round(2*(res[i,] - dota - rep(dota[i],n)),9))
        res[i,] <- besselJ(xx,nu)*(xx^(-nu))
        res[i,which(xx < 10e-5)] <- lim
        res[i,] <- z[i,]*(((res[i,]/lim)^ni)*z)
      }
      return(res)
    }
  if (is(y,"matrix"))
    {
      if(is.null(k)) stop("k not specified!")
      m <- dim(y)[1]
      if(!dim(k)[1]==m)
        stop("k must have equal rows to y")
      k <- as.matrix(k)
      if(!dim(x)[2]==dim(y)[2])
        stop("matrixes must have the same number of columns")
      dotb <- rowSums(y*y)/2
      res <- x%*%t(y)
      for( i in 1:m){#2*sigma or sigma
        xx <- sigma*sqrt(-round(2*(res[,i] - dota - rep(dotb[i],n)),9))
        res[,i] <- besselJ(xx,nu)*(xx^(-nu))
        res[which(xx<10e-5),i] <- lim
        res[,i]<- k[i,]*(((res[,i]/lim)^ni)*z)
      }
      return(res)
    }
}
setMethod("kernelPol",signature(kernel="besselkernel"),kernelPol.besselkernel)


kernelPol.anovakernel <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is(y,"matrix")&&!is.null(y)&&!is(y,"vector")) stop("y must be a matrix or NULL")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
  if(!is(k,"matrix")&&!is(k,"vector")&&!is.null(k)) stop("k must be a matrix or a vector")
  sigma <- kpar(kernel)$sigma
  degree <- kpar(kernel)$degree
 if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  n <- dim(x)[1]
  z <- as.matrix(z)
  if(!dim(z)[1]==n)
    stop("z must have the length equal to x colums")
  if (is.null(y))
    {
      if(is(z,"matrix")&&!dim(z)[1]==n)
       stop("z must have size equal to x colums")
      a <- matrix(0, dim(x)[2], n)
      res <- matrix(0,n,n)
      for (i in 1:n)
        {
          a[rep(TRUE,dim(x)[2]), rep(TRUE,n)] <- x[i,]
          res[i,]<- z[i,]*((colSums(exp( - sigma*(a - t(x))^2))^degree)*z)
        }
      return(res)
    }
  if (is(y,"matrix"))
    {
      if(is.null(k)) stop("k not specified!")
      m <- dim(y)[1]
      k <- as.matrix(k)
      if(!dim(k)[1]==m)
        stop("k must have equal rows to y")
      if(!dim(x)[2]==dim(y)[2])
        stop("matrixes must have the same number of columns")

      b <- matrix(0, dim(x)[2],m)
      res <- matrix(0, dim(x)[1],m)
      for( i in 1:n)
        {
          b[rep(TRUE,dim(x)[2]), rep(TRUE,m)] <- x[i,]
          res[i,] <- z[i,]*((colSums(exp( - sigma*(b - t(y))^2))^degree)*k)
        }
      return(res)
    }
}
setMethod("kernelPol",signature(kernel="anovakernel"),kernelPol.anovakernel)


kernelPol.splinekernel <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is(y,"matrix")&&!is.null(y)&&!is(y,"vector")) stop("y must be a matrix or NULL")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
  if(!is(k,"matrix")&&!is(k,"vector")&&!is.null(k)) stop("k must be a matrix or a vector")
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  sigma <- kpar(kernel)$sigma
  degree <- kpar(kernel)$degree
  n <- dim(x)[1]
  z <- as.vector(z)
  if(!(length(z)==n))
    stop("z must have the length equal to x colums")
  if (is.null(y))
    {
      
      res <- kernelMatrix(kernel,x)
      return(unclass(z*t(res*z)))
    }
  if (is(y,"matrix"))
    {
      if(is.null(k)) stop("k not specified!")
      m <- dim(y)[1]
      k <- as.vector(k)
      if(!(length(k)==m))
        stop("k must have length equal to rows of y")

      res <- kernelMatrix(kernel,x,y)
      return(unclass(k*t(res*z)))
    }
}
setMethod("kernelPol",signature(kernel="splinekernel"),kernelPol.splinekernel)


kernelPol.polykernel <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is(y,"matrix")&&!is.null(y)&&!is(y,"vector")) stop("y must be a matrix or NULL")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
  if(!is(k,"matrix")&&!is(k,"vector")&&!is.null(k)) stop("k must be a matrix or a vector")
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  degree <- kpar(kernel)$degree
  scale <- kpar(kernel)$scale
  offset <- kpar(kernel)$offset
  n <- dim(x)[1]

  if(is(z,"matrix"))
    {
      z <- as.vector(z)
    }
  m <- length(z)
  
  if(!(m==n))
    stop("z must have the length equal to x colums")
  if (is.null(y))
    {
      res <- z*t(((scale*crossprod(t(x))+offset)^degree)*z)
      return(res)
    }
  if (is(y,"matrix"))
    {
      if(is.null(k)) stop("k not specified!")
      m <- dim(y)[1]
      k <- as.vector(k)
      if(!(length(k)==m))
        stop("k must have length equal to rows of y")
      if(!dim(x)[2]==dim(y)[2])
        stop("matrixes must have the same number of columns")
      res<- k*t(((scale*x%*%t(y) + offset)^degree)*z)
      return(res)
    }
}
setMethod("kernelPol",signature(kernel="polykernel"),kernelPol.polykernel)


kernelPol.tanhkernel <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is(y,"matrix")&&!is.null(y)&&!is(y,"vector")) stop("y must be a matrix, vector or NULL")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
  if(!is(k,"matrix")&&!is(k,"vector")&&!is.null(k)) stop("k must be a matrix or a vector")
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  scale <- kpar(kernel)$scale
  offset <- kpar(kernel)$offset
  n <- dim(x)[1]
  if(is(z,"matrix"))
    {
      z <- as.vector(z)
    }
  m <- length(z)

  if(!(m==n))
    stop("z must have the length equal to x colums")
  if (is.null(y))
    {
      res <- z*t(tanh(scale*crossprod(t(x))+offset)*z)
      return(res)
    }
  if (is(y,"matrix"))
    {
      if(is.null(k)) stop("k not specified!")
      m <- dim(y)[1]
      k <- as.vector(k)
      if(!(length(k)==m))
        stop("k must have length equal rows to y")
      if(!dim(x)[2]==dim(y)[2])
        stop("matrixes x, y must have the same number of columns")
      res<- k*t(tanh(scale*x%*%t(y) + offset)*z)
      return(res)
    }
}
setMethod("kernelPol",signature(kernel="tanhkernel"),kernelPol.tanhkernel)


kernelPol.vanillakernel <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is(y,"matrix")&&!is.null(y)&&!is(y,"vector")) stop("y must be a matrix, vector or NULL")
  if(!is(z,"matrix")&&!is(z,"vector")) stop("z must be a matrix or a vector")
  if(!is(k,"matrix")&&!is(k,"vector")&&!is.null(k)) stop("k must be a matrix or a vector")
  n <- dim(x)[1]
  if(is(z,"matrix"))
    {
      z <- as.vector(z)
    }
      m <- length(z)
      
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if(!(m==n))
    stop("z must have the length equal to x colums")
  if (is.null(y))
    {
        res <- z*t(crossprod(t(x))*z)
      return(res)
    }
  if (is(y,"matrix"))
    {
      if(is.null(k)) stop("k not specified!")
      m <- dim(y)[1]
      k <- as.vector(k)
      if(!length(k)==m)
        stop("k must have length equal rows to y")
      if(!dim(x)[2]==dim(y)[2])
        stop("matrixes x, y must have the same number of columns")
      for( i in 1:m)
        res<- k*t(x%*%t(y)*z)
      return(res)
    }
}
setMethod("kernelPol",signature(kernel="vanillakernel"),kernelPol.vanillakernel)

kernelPol.stringkernel <- function(kernel, x, y=NULL ,z ,k=NULL)
{
  n <- length(x)
  res1 <- matrix(rep(0,n*n), ncol = n)
  resdiag <- rep(0,n)

  if(is(x,"list"))
    x <- sapply(x,paste,collapse="")
  
  if(is(y,"list"))
    y <- sapply(y,paste,collapse="")

  
  normalized = kpar(kernel)$normalized
  if(normalized == TRUE)
    kernel <- stringdot(length = kpar(kernel)$length, type = kpar(kernel)$type, lambda = kpar(kernel)$lambda, normalized = FALSE)
  
  z <- as.matrix(z)
  ## y is null 
  if (kpar(kernel)$type == "sequence" |kpar(kernel)$type == "string"|kpar(kernel)$type == "fullstring")
    {

      if(is.null(y)){
        if(normalized == TRUE){
          ## calculate diagonal elements first, and use them to normalize
          for (i in 1:n)
            resdiag[i] <- kernel(x[[i]],x[[i]])
      
          for(i in 1:n) {
            for(j in (i:n)[-1]) {
              res1[i,j]  <- (z[i,]*kernel(x[[i]],x[[j]])*z[j,])/sqrt(resdiag[i]*resdiag[j])
            }
          }
          res1 <- res1 + t(res1)
          diag(res1) <- z^2
        }
        else
          {
            for (i in 1:n)
              resdiag[i] <- kernel(x[[i]],x[[i]])
            
            for(i in 1:n) {
              for(j in (i:n)[-1]) {
                res1[i,j]  <- (z[i,]*kernel(x[[i]],x[[j]])*z[j,])
              }
            }
            res1 <- res1 + t(res1)
            diag(res1) <- resdiag * z^2
          }
      }
      
      if (!is.null(y)){
        if(normalized == TRUE){
          m <- length(y)
          res1 <- matrix(0,n,m)
          resdiag1 <- rep(0,m)
          
          k <- as.matrix(k)
          for(i in 1:n)
            resdiag[i] <- kernel(x[[i]],x[[i]])
          
          for(i in 1:m)
            resdiag1[i] <- kernel(y[[i]],y[[i]])
          
          for(i in 1:n) {
            for(j in 1:m) {
              res1[i,j] <- (z[i,]*kernel(x[[i]],y[[j]])*k[j,])/sqrt(resdiag[i]*resdiag1[j])
            }
          }
        }
      }
      else{
        m <- length(y)
        res1 <- matrix(0,n,m)
        k <- as.matrix(k)
        
        for(i in 1:n) {
          for(j in 1:m) {
            res1[i,j] <- (z[i,]*kernel(x[[i]],y[[j]])*k[j,])
          }
        }
      }
    }
  else {
    
    switch(kpar(kernel)$type,
           "exponential" = 
           sktype <- 2,
           "constant" =
           sktype <- 1,
           "spectrum" =
           sktype <- 3,
           "boundrange" =
           sktype <- 4)
    
    
    if(is(x,"list"))
      x <- unlist(x)
    if(is(y,"list"))
      y <- unlist(y)

    
    x <- paste(x,"\n",seq="")
    if(!is.null(y))
      y <- paste(y,"\n",seq="")  


    if(is.null(y))
      ret <- matrix(0, length(x),length(x))
    else
      ret <- matrix(0,length(x),length(y))
    
    if(is.null(y)){
      for( i in 1:length(x))
        ret[i,] <- .Call("stringtv",as.character(x[i]),as.character(x),as.integer(length(x)),as.integer(nchar(x[i])),as.integer(nchar(x)),as.integer(sktype),as.double(kpar(kernel)$lambda))
      res1 <- k*ret*k
    }
    else{
      for( i in 1:length(x))
        ret[i,] <- .Call("stringtv",as.character(x[i]),as.character(y),as.integer(length(x)),as.integer(nchar(x[i])),as.integer(nchar(y)),as.integer(sktype),as.double(kpar(kernel)$lambda))
      res1 <- k*ret*z
    }
    if(normalized == TRUE){
      if(is.null(y)){
        ret <- t((1/sqrt(diag(ret)))*t(ret*(1/sqrt(diag(ret)))))
        res1 <- k*ret*k
      }
      else{
        norm1 <- rep(0,length(x))
        norm2 <- rep(0,length(y))
        for( i in 1:length(x))
          norm1[i] <- .Call("stringtv",as.character(x[i]),as.character(x[i]),as.integer(1),as.integer(nchar(x[i])),as.integer(nchar(x[i])),as.integer(sktype),as.double(kpar(kernel)$lambda))
        for( i in 1:length(y))
          norm2[i] <- .Call("stringtv",as.character(y[i]),as.character(y[i]),as.integer(1),as.integer(nchar(y[i])),as.integer(nchar(y[i])),as.integer(sktype),as.double(kpar(kernel)$lambda))
        ret <- t((1/sqrt(norm2))*t(ret*(1/sqrt(norm1))))
        res1 <- k*ret*z 
      }
    }
  }

  return(res1)
}
setMethod("kernelPol",signature(kernel="stringkernel"),kernelPol.stringkernel)

## kernelFast returns the kernel matrix, its usefull in algorithms
## which require iterative kernel matrix computations

kernelFast <- function(kernel, x, y, a)
{
  return(kernelMatrix(kernel,x,y))
}
setGeneric("kernelFast",function(kernel, x, y, a) standardGeneric("kernelFast"))



kernelFast.rbfkernel <- function(kernel, x, y, a)
{
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if(!is(y,"matrix")) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  n <- dim(x)[1]
  dota <- a/2
   if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    m <- dim(y)[1]
    dotb <- rowSums(y*y)/2
    res <- x%*%t(y)
    for( i in 1:m)
      res[,i]<- exp(2*sigma*(res[,i] - dota - rep(dotb[i],n)))
    return(res)
  }
}
setMethod("kernelFast",signature(kernel="rbfkernel"),kernelFast.rbfkernel)

kernelFast.laplacekernel <- function(kernel, x, y, a)
{
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if(!is(y,"matrix")) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  n <- dim(x)[1]
  dota <- a/2
   if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    m <- dim(y)[1]
    dotb <- rowSums(y*y)/2
    res <- x%*%t(y)
    for( i in 1:m)
      res[,i]<- exp(-sigma*sqrt(round(-2*(res[,i] - dota - rep(dotb[i],n)),9)))
    return(res)
  }
}
setMethod("kernelFast",signature(kernel="laplacekernel"),kernelFast.laplacekernel)

kernelFast.besselkernel <- function(kernel, x, y, a)
{
  if(is(x,"vector"))
    x <- as.matrix(x)
  if(is(y,"vector"))
    y <- as.matrix(y)
  if(!is(y,"matrix")) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  nu = kpar(kernel)$order
  ni = kpar(kernel)$degree
  n <- dim(x)[1]
  lim <- 1/(gamma(nu+1)*2^(nu))
  dota <- a/2
  if (is(x,"matrix") && is(y,"matrix")){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    m <- dim(y)[1]
    dotb <- rowSums(y*y)/2
    res <- x%*%t(y)
    for( i in 1:m){
      xx <- sigma*sqrt(round(-2*(res[,i] - dota - rep(dotb[i],n)),9))
      res[,i] <- besselJ(xx,nu)*(xx^(-nu))
      res[which(xx<10e-5),i] <- lim
    }     
    return((res/lim)^ni)
  }
}
setMethod("kernelFast",signature(kernel="besselkernel"),kernelFast.besselkernel)


kernelFast.anovakernel <- function(kernel, x, y, a)
{
  return(kernelMatrix(kernel,x,y))
}
setMethod("kernelFast",signature(kernel="anovakernel"),kernelFast.anovakernel)


kernelFast.polykernel <- function(kernel, x, y, a)
{
  return(kernelMatrix(kernel,x,y))
}
setMethod("kernelFast",signature(kernel="polykernel"),kernelFast.polykernel)

kernelFast.vanilla <- function(kernel, x, y, a)
{
  return(kernelMatrix(kernel,x,y))
}
setMethod("kernelFast",signature(kernel="vanillakernel"),kernelFast.vanilla)

kernelFast.tanhkernel <- function(kernel, x, y, a)
{
  return(kernelMatrix(kernel,x,y))
}
setMethod("kernelFast",signature(kernel="tanhkernel"),kernelFast.tanhkernel)

kernelFast.stringkernel <- function(kernel, x, y, a)
{
  return(kernelMatrix(kernel,x,y))
}
setMethod("kernelFast",signature(kernel="stringkernel"),kernelFast.stringkernel)

kernelFast.splinekernel <- function(kernel, x, y, a)
{
  return(kernelMatrix(kernel,x,y))
}
setMethod("kernelFast",signature(kernel="splinekernel"),kernelFast.splinekernel)




