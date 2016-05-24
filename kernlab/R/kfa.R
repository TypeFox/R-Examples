
## This code takes the set x of vectors from the input space
## and does projection pursuit to find a good basis for x.
##
## The algorithm is described in Section 14.5 of
## Learning with Kernels by B. Schoelkopf and A. Smola, entitled 
## Kernel Feature Analysis.
## 
## created : 17.09.04 alexandros
## updated :

setGeneric("kfa",function(x, ...) standardGeneric("kfa"))
setMethod("kfa", signature(x = "formula"),
function(x, data = NULL, na.action = na.omit, ...)
{
    mt <- terms(x, data = data)
    if(attr(mt, "response") > 0) stop("response not allowed in formula")
    attr(mt, "intercept") <- 0
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf$formula <- mf$x
    mf$... <- NULL
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    Terms <- attr(mf, "terms")
    na.act <- attr(mf, "na.action")
    x <- model.matrix(mt, mf)
    res <- kfa(x, ...)
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("kfa")
    kcall(res) <- cl
    attr(Terms,"intercept") <- 0
    terms(res) <- Terms
    if(!is.null(na.act))
        n.action(res) <- na.act
  
   return(res)
  })

setMethod("kfa",signature(x="matrix"),
function(x, kernel="rbfdot", kpar=list(sigma=0.1), features = 0, subset = 59, normalize = TRUE, na.action = na.omit)
{  
  if(!is.matrix(x))
    stop("x must be a matrix")
  
  x <- na.action(x)  
  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }
  
  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

  ## initialize variables
  m <- dim(x)[1]

  if(subset > m)
    subset <- m

  if (features==0)
      features <-  subset 
  
  alpha <- matrix(0,subset,features)
  alphazero <- rep(1,subset)
  alphafeat <- matrix(0,features,features)
  idx <- -(1:subset)
  randomindex <- sample(1:m, subset)
  K <- kernelMatrix(kernel,x[randomindex,,drop=FALSE],x)
  
  ## main loop
 for (i in 1:features)
    {
      K.cols <- K[-idx, , drop = FALSE] 

      if(i > 1)
        projections <- K.cols *  (alphazero[-idx]%*%t(rep(1,m))) + crossprod(t(alpha[-idx,1:(i-1),drop=FALSE]),K[idx, ,drop = FALSE])
      else
        projections <- K.cols *  (alphazero%*%t(rep(1,m))) 
      
      Q <- apply(projections, 1, sd) 
      Q.tmp <- rep(0,subset)
      Q.tmp[-idx] <- Q
      Qidx <- which.max(Q.tmp)
      Qmax <- Q.tmp[Qidx]

      if(i > 1)
        alphafeat[i,1:(i-1)] <- alpha[Qidx,1:(i-1)]

      alphafeat[i,i] <- alphazero[Qidx]

      if (i > 1)
        idx <- c(idx,Qidx)
      else
        idx <- Qidx

      if (i > 1)
        Qfeat <- c(Qfeat, Qmax)
      else
        Qfeat <- Qmax

      Ksub <- K[idx, idx, drop = FALSE]
      alphasub <- alphafeat[i,1:i]
      phisquare <- alphasub %*% Ksub %*% t(t(alphasub))
      dotprod <- (alphazero * (K[,idx, drop = FALSE] %*% t(t(alphasub))) + alpha[,1:i]%*%(Ksub%*%t(t(alphasub))))/drop(phisquare)
      alpha[,1:i] <- alpha[,1:i] - dotprod %*%alphasub

      if(normalize){
        sumalpha <- alphazero + rowSums(abs(alpha))
        alphazero <- alphazero / sumalpha
        alpha <- alpha/ (sumalpha %*% t(rep(1,features)))
      }
    }

  obj <- new("kfa")
  alpha(obj) <- alphafeat
  alphaindex(obj) <- randomindex[idx]
  xmatrix(obj) <- x[alphaindex(obj),]
  kernelf(obj) <- kernel
  kcall(obj) <- match.call()
  return(obj)
})


## project a new matrix into the feature space 

setMethod("predict",signature(object="kfa"),
function(object , x)
  {
    if (!is.null(terms(object)))
      {
        if(!is.matrix(x))
          x <- model.matrix(delete.response(terms(object)), as.data.frame(x), na.action = n.action(object))
      }
    else
      x  <- if (is.vector(x)) t(t(x)) else as.matrix(x)
    
    if (!is.matrix(x)) stop("x must be a matrix a vector or a data frame")
    tmpres <- kernelMult(kernelf(object), x, xmatrix(object), alpha(object))
    return(tmpres - matrix(colSums(tmpres)/dim(tmpres)[1],dim(tmpres)[1],dim(tmpres)[2],byrow=TRUE))

    
  })

setMethod("show",signature(object="kfa"),
function(object)
  {
    cat(paste("Number of features :",dim(alpha(object))[2],"\n"))
    show(kernelf(object))
  })
  


