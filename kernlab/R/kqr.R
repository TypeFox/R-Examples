setGeneric("kqr", function(x, ...) standardGeneric("kqr"))
setMethod("kqr",signature(x="formula"),
function (x, data=NULL, ..., subset, na.action = na.omit, scaled = TRUE){
  cl <- match.call()
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m$formula <- m$x
  m$x <- NULL
  m[[1L]] <- quote(stats::model.frame)
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  attr(Terms, "intercept") <- 0
  x <- model.matrix(Terms, m)
  y <- model.extract(m, "response")

  if (length(scaled) == 1)
    scaled <- rep(scaled, ncol(x))
  if (any(scaled)) {
    remove <- unique(c(which(labels(Terms) %in% names(attr(x, "contrasts"))),
                       which(!scaled)
                       )
                     )
    scaled <- !attr(x, "assign") %in% remove
  }
  
  ret <- kqr(x, y, scaled = scaled, ...)
  kcall(ret) <- cl
  terms(ret) <- Terms
  if (!is.null(attr(m, "na.action")))
    n.action(ret) <- attr(m, "na.action")
  return (ret)
})

setMethod("kqr",signature(x="vector"),
function(x,...)
  {
    x <- t(t(x))
    ret <- kqr(x, ...)
    ret
  })
    
setMethod("kqr",signature(x="matrix"),
function (x, y, scaled = TRUE, tau = 0.5, C = 0.1, kernel = "rbfdot", kpar = "automatic", reduced = FALSE, rank = dim(x)[1]/6, fit = TRUE, cross = 0, na.action = na.omit)
  {
    if((tau > 1)||(tau < 0 )) stop("tau has to be strictly between 0 and 1")
    
    ret <- new("kqr")
    param(ret) <- list(C = C, tau = tau)
    if (is.null(y))
      x <- na.action(x)
    else {
      df <- na.action(data.frame(y, x))
      y <- df[,1]
      x <- as.matrix(df[,-1])
    }
    ncols <- ncol(x)
    m <- nrows <- nrow(x)
    tmpsc <- NULL
  x.scale <- y.scale <- NULL
 ## scaling
  if (length(scaled) == 1)
    scaled <- rep(scaled, ncol(x))
  if (any(scaled)) {
    co <- !apply(x[,scaled, drop = FALSE], 2, var)
    if (any(co)) {
      scaled <- rep(FALSE, ncol(x))
      warning(paste("Variable(s)",
                    paste("`",colnames(x[,scaled, drop = FALSE])[co],
                          "'", sep="", collapse=" and "),
                    "constant. Cannot scale data.")
              )
    } else {
      xtmp <- scale(x[,scaled])
      x[,scaled] <- xtmp
      x.scale <- attributes(xtmp)[c("scaled:center","scaled:scale")]
      y <- scale(y)
      y.scale <- attributes(y)[c("scaled:center","scaled:scale")]
      y <- as.vector(y)
      tmpsc <- list(scaled = scaled, x.scale = x.scale,y.scale = y.scale)
    }
  }

    ## Arrange all the kernel mambo jumpo
    
    if(is.character(kernel)){
      kernel <- match.arg(kernel,c("rbfdot","polydot","tanhdot","vanilladot","laplacedot","besseldot","anovadot","splinedot"))
      if(is.character(kpar))
        if((kernel == "tanhdot" || kernel == "vanilladot" || kernel == "polydot"|| kernel == "besseldot" || kernel== "anovadot"|| kernel=="splinedot") &&  kpar=="automatic" )
          {
            cat (" Setting default kernel parameters ","\n")
            kpar <- list()
          }
    }
    
    if (!is.function(kernel))
      if (!is.list(kpar)&&is.character(kpar)&&(class(kernel)=="rbfkernel" || class(kernel) =="laplacedot" || kernel == "laplacedot"|| kernel=="rbfdot")){
        kp <- match.arg(kpar,"automatic")
        if(kp=="automatic")
          kpar <- list(sigma=mean(sigest(x,scaled=FALSE,frac=1)[c(1,3)]))
        cat("Using automatic sigma estimation (sigest) for RBF or laplace kernel","\n")
        
      }
    
    if(!is(kernel,"kernel"))
      {
        if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
      }
    if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")
    
    ## Setup QP problem and call ipop
    if(!reduced)
      H = kernelMatrix(kernel,x)
    else
      H = csi(x, kernel = kernel, rank = rank)
    c = -y
    A = rep(1,m)
    b = 0
    r = 0
    l = matrix(C * (tau-1),m,1)
    u = matrix(C * tau ,m,1)    
                       
    qpsol = ipop(c, H, A, b, l, u, r)
    alpha(ret)= coef(ret) = primal(qpsol)
    b(ret) = dual(qpsol)[1]
    
    ## Compute training error/loss
    xmatrix(ret) <- x
    ymatrix(ret) <- y
    kernelf(ret) <- kernel
    kpar(ret) <- kpar
    type(ret) <- ("Quantile Regresion")
    
    if (fit){
      fitted(ret) <- predict(ret, x)
      if (!is.null(scaling(ret)$y.scale))
        fitted(ret) <- fitted(ret) * tmpsc$y.scale$"scaled:scale" + tmpsc$y.scale$"scaled:center"
      error(ret) <- c(pinloss(y, fitted(ret), tau), ramploss(y,fitted(ret),tau))

    }
    else fitted(ret) <- NULL

    if(any(scaled))
      scaling(ret) <- tmpsc
    
    ## Crossvalidation
    cross(ret) <- -1
  if(cross == 1)
    cat("\n","cross should be >1 no cross-validation done!","\n","\n")
  else if (cross > 1)
    {
      pinloss <- 0
      ramloss <- 0
      crescs <- NULL
      suppressWarnings(vgr<-split(sample(1:m,m),1:cross))
      for(i in 1:cross)
        {
          cind <- unsplit(vgr[-i],factor(rep((1:cross)[-i],unlist(lapply(vgr[-i],length)))))
          cret <- kqr(x[cind,],y[cind], tau = tau, C = C, scale = FALSE, kernel = kernel, cross = 0, fit = FALSE)
          cres <- predict(cret, x[vgr[[i]],])
          crescs <- c(crescs,cres)
        }
      if (!is.null(scaling(ret)$y.scale)){
        crescs <- crescs * tmpsc$y.scale$"scaled:scale" + tmpsc$y.scale$"scaled:center"
        ysvgr <- y[unlist(vgr)] * tmpsc$y.scale$"scaled:scale" + tmpsc$y.scale$"scaled:center"
      }
      else
        ysvgr <-  y[unlist(vgr)]
        
      pinloss <- drop(pinloss(ysvgr, crescs, tau))
      ramloss <- drop(ramloss(ysvgr, crescs, tau))
      cross(ret) <- c(pinloss, ramloss)
    }
  
    return(ret)
})


setMethod("kqr",signature(x="list"),
function (x, y, tau = 0.5, C = 0.1, kernel = "strigdot", kpar = list(length=4, C=0.5), fit = TRUE, cross = 0)
  {
    if((tau > 1)||(tau < 0 )) stop("tau has to be strictly between 0 and 1")
    
  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }
  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

  K <- kernelMatrix(kernel,x)

    ret <- kqr(K,y = y,tau = tau, C = C, fit = fit, cross = cross)

    kernelf(ret) <- kernel
    kpar(ret) <- kpar
    
  return(ret)

})



setMethod("kqr",signature(x="kernelMatrix"),
function (x, y, tau = 0.5, C = 0.1, fit = TRUE, cross = 0)
  {
    if((tau > 1)||(tau < 0 )) stop("tau has to be strictly between 0 and 1")
    ret <- new("kqr")
    param(ret) <- list(C = C, tau = tau)
    ncols <- ncol(x)
    m <- nrows <- nrow(x)

    y <- as.vector(y)
    
    ## Setup QP problem and call ipop
    
    H = x
    c = -y
    A = rep(1,m)
    b = 0
    r = 0
    l = matrix(C * (tau-1),m,1)
    u = matrix(C * tau ,m,1)    
                       
    qpsol = ipop(c, H, A, b, l, u, r)
    alpha(ret)= coef(ret) = primal(qpsol)
    b(ret) = dual(qpsol)[1]
    
    ## Compute training error/loss
    ymatrix(ret) <- y
    kernelf(ret) <- "Kernel Matrix used."
    type(ret) <- ("Quantile Regresion")
    
    if (fit){
      fitted(ret) <- predict(ret, x)
      error(ret) <- c(pinloss(y, fitted(ret), tau), ramploss(y,fitted(ret),tau))

    }
    else NA 

    ## Crossvalidation
    cross(ret) <- -1
  if(cross == 1)
    cat("\n","cross should be >1 no cross-validation done!","\n","\n")
  else if (cross > 1)
    {
      pinloss <- 0
      ramloss <- 0
      crescs <- NULL
      suppressWarnings(vgr<-split(sample(1:m,m),1:cross))
      for(i in 1:cross)
        {
          cind <- unsplit(vgr[-i],factor(rep((1:cross)[-i],unlist(lapply(vgr[-i],length)))))
          cret <- kqr(x[cind,cind],y[cind], tau = tau, C = C, scale = FALSE, cross = 0, fit = FALSE)
          cres <- predict(cret, x[vgr[[i]],vgr[[i]]])
          crescs <- c(crescs,cres)
        }
      ysvgr <-  y[unlist(vgr)]
      
      pinloss <- drop(pinloss(ysvgr, crescs, tau))
      ramloss <- drop(ramloss(ysvgr, crescs, tau))
      cross(ret) <- c(pinloss, ramloss)
    }
  
    return(ret)
  })


pinloss <- function(y,f,tau)
  {

    if(is.vector(y)) m <- length(y)
    else m <-  dim(y)[1]
    tmp <- y - f
    return((tau *sum(tmp*(tmp>=0)) +  (tau-1) * sum(tmp * (tmp<0)))/m)

  }

ramploss <- function(y,f,tau)
  {
    if(is.vector(y)) m <- length(y)
    else m <-  dim(y)[1]
    
    return(sum(y<=f)/m)
  }


setMethod("predict", signature(object = "kqr"),
function (object, newdata)
{
  sc <- 0
  if (missing(newdata))
    if(!is.null(fitted(object)))
      return(fitted(object))
    else
      stop("newdata is missing and no fitted values found.")

  if(!is(newdata,"kernelMatrix")){
    ncols <- ncol(xmatrix(object))
    nrows <- nrow(xmatrix(object))
    oldco <- ncols
    
    if (!is.null(terms(object)))
      {  
        newdata <- model.matrix(delete.response(terms(object)), as.data.frame(newdata), na.action = na.action)
      }
    else
      newdata  <- if (is.vector (newdata)) t(t(newdata)) else as.matrix(newdata)
    
    newcols  <- 0
    newnrows <- nrow(newdata)
    newncols <- ncol(newdata)
    newco    <- newncols
    
    if (oldco != newco) stop ("test vector does not match model !")
    
    if (is.list(scaling(object)) && sc != 1)
      newdata[,scaling(object)$scaled] <-
        scale(newdata[,scaling(object)$scaled, drop = FALSE],
              center = scaling(object)$x.scale$"scaled:center",
              scale  = scaling(object)$x.scale$"scaled:scale"
              )
    
    predres <- kernelMult(kernelf(object),newdata,xmatrix(object),as.matrix(alpha(object))) - b(object)
    
    if (!is.null(scaling(object)$y.scale))
      return(predres * scaling(object)$y.scale$"scaled:scale" + scaling(object)$y.scale$"scaled:center")
    else
      return(predres)
  }
     else
     {
       return(newdata%*%alpha(object) - b(object))
     }

})


setMethod("show","kqr",
function(object){
  cat("Kernel Quantile Regression object of class \"kqr\"","\n")
  cat("\n")
  show(kernelf(object))
  cat("\n")
  cat("Regularization Cost Parameter C: ",round(param(object)[[1]],9))
  cat(paste("\nNumber of training instances learned :", dim(xmatrix(object))[1],"\n"))
  if(!is.null(fitted(object)))
    cat(paste("Train error :"," pinball loss : ", round(error(object)[1],9)," rambloss :", round(error(object)[2],9),"\n"))
  ##train error & loss
  if(cross(object)!=-1)
    cat("Cross validation error :", " pinballoss : ", round(cross(object)[1],9)," rambloss :", round(cross(object)[2],9),"\n")
})
