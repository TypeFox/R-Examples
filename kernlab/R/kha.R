

#Kernel Hebbian Algorithm function

setGeneric("kha",function(x, ...) standardGeneric("kha"))
setMethod("kha", signature(x = "formula"),
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
    na.act <- attr(mf, "na.action")
    Terms <- attr(mf, "terms")
    x <- model.matrix(mt, mf)
    res <- kha(x, ...)
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("kha")
    kcall(res) <- cl
    attr(Terms,"intercept") <- 0
    terms(res) <- Terms
    if(!is.null(na.act)) 
        n.action(res) <- na.act
    return(res)
  })



setMethod("kha",signature(x="matrix"),
          function(x, kernel = "rbfdot", kpar = list(sigma = 0.1),
          features = 5, eta = 0.005, th = 1e-4, maxiter = 10000, verbose = FALSE, na.action = na.omit, ...)
{
  x <- na.action(x)
  x <- as.matrix(x)
  m <- nrow(x)
  ret <- new("kha")
  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }
  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

  ## Initialize A dual variables
  A <- matrix(runif(features*m),m,features)*2 - 1
  AOld <- A

  ## compute square norm of data
  a <- rowSums(x^2)

  ## initialize the empirical sum kernel map
  eskm <- rep(0,m)

  for (i in 1:m)
    eskm[i] <- sum(kernelFast(kernel,x,x[i,,drop=FALSE], a))

  eks <- sum(eskm)

  counter <- 0
  step <- th + 1
  Aold <- A
  
  while(step > th && counter < maxiter)
    {
      y <- rep(0, features)
      ot <- rep(0,m)
      
      ## Hebbian Iteration
      for (i in 1:m)
        { 
          ## compute y output
          etkm <- as.vector(kernelFast(kernel,x,x[i,,drop=FALSE], a))
          sum1 <- as.vector(etkm %*% A)
          sum2 <- as.vector(eskm%*%A)/m
          asum <- colSums(A)
          sum3 <- as.vector(eskm[i]*asum)/m
          sum4 <- as.vector(eks * asum)/m^2
          y <- sum1 - sum2 - sum3 + sum4

          ## update A
          yy <- y%*%t(y)
          yy[upper.tri(yy)] <- 0
          tA <- t(A)
          A <- t(tA - eta * yy%*%tA)
          A[i,] <- A[i,] + eta * y
        }

      if (counter %% 100 == 0 )
        {
         step = mean(abs(Aold - A))
         Aold <- A
         if(verbose)
           cat("Iteration :", counter, "Converged :", step,"\n")
       }
      counter <- counter + 1
    }

  ## Normalize in Feature space
  cA <- t(A) - colSums(A)
  Fnorm <- rep(0,features)
  for (j in 1:m)
    Fnorm <- Fnorm + colSums(t(cA[,j] * cA) * as.vector(kernelFast(kernel,x,x[j,,drop=FALSE],a)))
  
  
  if(any(Fnorm==0))
    {
      warning("Normalization vector contains zeros, replacing them with ones")
      Fnorm[which(Fnorm==0)] <- 1
    }

  A <- t(t(A)/sqrt(Fnorm))
  
  pcv(ret) <- A
  eig(ret) <- Fnorm
  names(eig(ret)) <- paste("Comp.", 1:features, sep = "")
  eskm(ret) <- eskm
  kcall(ret) <- match.call()
  kernelf(ret) <- kernel
  xmatrix(ret) <- x
  return(ret)
})


## Project a new matrix into the feature space 
setMethod("predict",signature(object="kha"),
function(object , x)
  {
    if (!is.null(terms(object)))
      {
        if(!is.matrix(x))
          x <- model.matrix(delete.response(terms(object)), as.data.frame(x), na.action = n.action(object))
      }
    else
      x  <- if (is.vector(x)) t(t(x)) else as.matrix(x)

    if (is.vector(x)||is.data.frame(x))
      x<-as.matrix(x)
    if (!is.matrix(x)) stop("x must be a matrix a vector or a data frame")
    n <- nrow(x)
    m <- nrow(xmatrix(object))
    A <- pcv(object)
    y <- matrix(0,n,dim(A)[2])
    eks <- sum(eskm(object))
    a <- rowSums(xmatrix(object)^2)
    
    ## Project data
    sum2 <- as.vector(eskm(object)%*%A)/m
    asum <- colSums(A)

    sum4 <- as.vector(eks * asum)/m^2
        
    for (i in 1:n)
      {
        ## compute y output
        etkm <- as.vector(kernelFast(kernelf(object),xmatrix(object),x[i,,drop=FALSE], a))
        sum1 <- as.vector(etkm %*% A)
        sum3 <- sum(etkm)*asum/m
        y[i,] <- sum1 - sum2 - sum3 + sum4
      }

    return(y)
  })


  
