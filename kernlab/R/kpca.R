## kpca function
## author : alexandros

setGeneric("kpca",function(x, ...) standardGeneric("kpca"))
setMethod("kpca", signature(x = "formula"),
function(x,  data = NULL, na.action = na.omit, ...)
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
    res <- kpca(x, ...)
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("kpca")
    kcall(res) <- cl
    attr(Terms,"intercept") <- 0
    terms(res) <- Terms
    if(!is.null(na.act)) 
        n.action(res) <- na.act
  
    return(res)
  })


## Matrix Interface
setMethod("kpca",signature(x="matrix"),
          function(x, kernel = "rbfdot", kpar = list(sigma = 0.1), features = 0, th = 1e-4, na.action = na.omit, ...)
{
  x <- na.action(x)
  x <- as.matrix(x)
  m <- nrow(x)
  ret <- new("kpca")
  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }
  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

  km <- kernelMatrix(kernel,x)

  ## center kernel matrix
  kc <- t(t(km - colSums(km)/m) -  rowSums(km)/m) + sum(km)/m^2

  ## compute eigenvectors
  res <- eigen(kc/m,symmetric=TRUE)
  
  if(features == 0)
    features <- sum(res$values > th)
  else 
    if(res$values[features] < th)
      warning(paste("eigenvalues of the kernel matrix are below threshold!"))
 
  pcv(ret) <- t(t(res$vectors[,1:features])/sqrt(res$values[1:features]))
  eig(ret) <- res$values[1:features]
  names(eig(ret)) <- paste("Comp.", 1:features, sep = "")
  rotated(ret) <- kc %*% pcv(ret)
  kcall(ret) <- match.call()
  kernelf(ret) <- kernel
  xmatrix(ret) <- x
  return(ret)
})

## List Interface
setMethod("kpca",signature(x="list"),
          function(x, kernel = "stringdot", kpar = list(length = 4, lambda = 0.5), features = 0, th = 1e-4, na.action = na.omit, ...)
{
  x <- na.action(x)
  m <- length(x)
  ret <- new("kpca")

  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }
  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

  km <- kernelMatrix(kernel,x)
  ## center kernel matrix
  kc <- t(t(km - colSums(km)/m) -  rowSums(km)/m) + sum(km)/m^2
 
  ## compute eigenvectors
  res <- eigen(kc/m,symmetric=TRUE)
  
  if(features == 0)
    features <- sum(res$values > th)
  else 
    if(res$values[features] < th)
      warning(paste("eigenvalues of the kernel matrix are below threshold!"))
 
  pcv(ret) <- t(t(res$vectors[,1:features])/sqrt(res$values[1:features]))
  eig(ret) <- res$values[1:features]
  names(eig(ret)) <- paste("Comp.", 1:features, sep = "")
  rotated(ret) <- kc %*% pcv(ret)
  kcall(ret) <- match.call()
  kernelf(ret) <- kernel
  xmatrix(ret) <- x
  return(ret)
})

## Kernel Matrix Interface
setMethod("kpca",signature(x= "kernelMatrix"),
          function(x, features = 0, th = 1e-4, ...)
{
  ret <- new("kpca")
  m <- dim(x)[1]
  if(m!= dim(x)[2])
    stop("Kernel matrix has to be symetric, and positive semidefinite")

  ## center kernel matrix
   kc <- t(t(x - colSums(x)/m) -  rowSums(x)/m) + sum(x)/m^2
  
  ## compute eigenvectors
  res <- eigen(kc/m,symmetric=TRUE)
  
  if(features == 0)
    features <- sum(res$values > th)
  else 
    if(res$values[features] < th)
      warning(paste("eigenvalues of the kernel matrix are below threshold!"))
 
  pcv(ret) <- t(t(res$vectors[,1:features])/sqrt(res$values[1:features]))
  eig(ret) <- res$values[1:features]
  names(eig(ret)) <- paste("Comp.", 1:features, sep = "")
  rotated(ret) <- kc %*% pcv(ret)
  kcall(ret) <- match.call()
  xmatrix(ret) <- x
  kernelf(ret) <- " Kernel matrix used."
  return(ret)
})


## project a new matrix into the feature space 
setMethod("predict",signature(object="kpca"),
function(object , x)
  {
    if (!is.null(terms(object)))
      {
        if(!is.matrix(x) || !is(x,"list"))
          x <- model.matrix(delete.response(terms(object)), as.data.frame(x), na.action = n.action(object))
      }
    else
      x  <- if (is.vector(x)) t(t(x)) else if (!is(x,"list")) x <- as.matrix(x)

    if (is.vector(x) || is.data.frame(x))
      x <- as.matrix(x)
    if (!is.matrix(x) && !is(x,"list")) stop("x must be a matrix a vector, a data frame, or a list")

    if(is(x,"matrix"))
      {
        n <- nrow(x)
        m <- nrow(xmatrix(object))}
    else
      {
        n <- length(x)
        m <- length(xmatrix(object))
      }

    if(is.character(kernelf(object)))
      {
        knc <- x
        ka <- xmatrix(object)
      }
    else
      {
        knc <- kernelMatrix(kernelf(object),x,xmatrix(object))
        ka <- kernelMatrix(kernelf(object),xmatrix(object))
      }
    ## center
    ret <- t(t(knc - rowSums(knc)/m) - rowSums(ka)/m) + sum(ka)/(m*n) 

    return(ret %*% pcv(object))
  })



  
