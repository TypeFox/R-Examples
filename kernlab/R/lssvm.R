## reduced least squares support vector machines
## author : alexandros

setGeneric("lssvm", function(x, ...) standardGeneric("lssvm"))
setMethod("lssvm",signature(x="formula"),
function (x, data=NULL, ..., subset, na.action = na.omit, scaled = TRUE){
    
  cl <- match.call()
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m$formula <- m$x
  m$x <- NULL
  m$scaled <- NULL
  m[[1L]] <- quote(stats::model.frame)    
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
   attr(Terms, "intercept") <- 0    ## no intercept
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
   ret <- lssvm(x, y, scaled = scaled, ...)
  kcall(ret) <- cl
  attr(Terms,"intercept") <- 0 ## no intercept
  terms(ret) <- Terms
  if (!is.null(attr(m, "na.action")))
    n.action(ret) <- attr(m, "na.action")
  return (ret)
})

setMethod("lssvm",signature(x="vector"),
function(x,...)
          { x <- t(t(x))
            ret <- lssvm(x, ...)
            return(ret)
          })
    
setMethod("lssvm",signature(x="matrix"),
function (x,
          y,        
          scaled    = TRUE,
          kernel    = "rbfdot",
          kpar      = "automatic",
          type      = NULL,
          tau       = 0.01,
          reduced   = TRUE, 
          tol       = 0.0001,
          rank      = floor(dim(x)[1]/3),
          delta      = 40,
          ##          prob.model = FALSE, 
          cross     = 0,
          fit       = TRUE,
          ...,
          subset,
          na.action = na.omit)
{ 
  ## subsetting and na-handling for matrices
  ret <- new("lssvm")
  if (!missing(subset)) x <- x[subset,]
 
  df <- unique(na.action(data.frame(y, x)))
  y <- df[,1]
  x <- as.matrix(df[,-1])
 
  n.action(ret) <- na.action

  if(!is.null(type))
    type(ret) <- match.arg(type,c("classification","regression"))
  
  if (is.null(type)) type(ret) <- if (is.factor(y)) "classification" else "regression"
  else type(ret) <- type
  
  ## scaling, subsetting, and NA handling
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
    }
  }
  ncols <- ncol(x)
  m <- nrows <- nrow(x)


    if(is.character(kernel)){
    kernel <- match.arg(kernel,c("rbfdot","polydot","tanhdot","vanilladot","laplacedot","besseldot","anovadot","splinedot","matrix"))

    if(kernel == "matrix")
      if(dim(x)[1]==dim(x)[2])
        return(lssvm(as.kernelMatrix(x), y = y,type = NULL,
          tau       = 0.01,
          tol       = 0.0001,
          rank      = floor(dim(x)[1]/3),
          delta      = 40,
          cross     = 0,
          fit       = TRUE,
          ...))
      else
        stop(" kernel matrix not square!")
    
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
      kpar <- list(sigma=mean(sigest(x,scaled=FALSE)[c(1,3)]))
   cat("Using automatic sigma estimation (sigest) for RBF or laplace kernel","\n")
   
  }
  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }

  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")


  
  if(type(ret)=="classification")
    {
      if (!is.vector(y) && !is.factor (y)) stop("y must be a vector or a factor.")
      if(is(y,"vector"))
        {
          y <- as.matrix(y)
          if (nrows != nrow(y)) stop("x and y don't match.")
        }
      
      if (is.factor(y)) {
        lev(ret) <- levels (y)
        y <- as.integer (y)
        if (nrows != length(y)) stop("x and y don't match.")
      }
      else if (is.numeric(y))
        {
          y <- as.integer(y)
          lev(ret) <- unique (y)  
        }
      else 
        stop ("dependent variable has to be of factor or integer type for classification mode.")

      ## initialize    
      nclass(ret) <- length (unique(y))
      p <- 0

      svindex <- NULL
      
      ## create multidimensional y matrix
      yind <- t(matrix(1:nclass(ret),nclass(ret),m))
      ymat <- matrix(0, m, nclass(ret))
      ymat[yind==y] <- 1

      if(reduced == FALSE)
        {
          K <- kernelMatrix(kernel,x)

          KP <- K - (1/m)*colSums(K) 
          beta <- solve((KP%*%K + m * tau * K), KP%*%ymat)
          b <- colMeans(ymat) - colMeans(K%*%beta)
          alphaindex(ret) <- 1:m
        }
      else
        {
          G <- csi(x, ymat, rank = rank ,kernel= kernel, delta = delta , tol = tol)
          rep <- sort(pivots(G),index.return=TRUE)$ix
          G <- G[rep,]
          GtP <- t(G) - matrix(rowSums(t(G))/dim(G)[1],dim(G)[2],dim(G)[1])
          Gtalpha <- (GtP)%*%G
          diag(Gtalpha) <-  diag(Gtalpha) + tau
          Gtalpha <- solve(Gtalpha) %*% GtP %*% ymat[rep,,drop=FALSE]
          beta <- solve(t(G[1:dim(G)[2],]), Gtalpha)
          b <- colMeans(ymat) - colMeans(G%*%Gtalpha) 
          alphaindex(ret) <- rep[1:dim(G)[2]]
        }
      
      alpha(ret) <- beta
      ## nonzero alpha*y
      coef(ret) <- alpha(ret)
      ## store SV indexes from current problem for later use in predict
  
      ## save the indexes from all the SV in a vector (use unique?)
      svindex <- alphaindex(ret)
      ## store betas in a vector 
      b(ret) <- b
      ##store C  in return object
      param(ret)$tau <- tau

      ## calculate class prob.
  ##    if (prob.model& reduced== TRUE)
  #      warning("Class Probapilities not supported for reduced model.)
      
 ##     if(prob.model & reduced == FALSE)
 ##       {
 ##         pos <- as.vector(ymat)==1
 ##         neg <- as.vector(ymat)==-1
 ##         ones <- rep(1,dim(x)[1])
 ##         onesneg <- ones[pos] <- 0
  ##        ones <- rep(1,dim(x)[1])
  ##        onespos <- ones[neg] <- 0
          ##Kpos <- kernelMult(kernel,x,x[pos,],rep(1,sum(pos)))
          ##Kneg <- kernelMult(kernel,x,x[neg,],rep(1,sum(neg)))
  ##        Kpos <- K[,pos]%*%rep(1,sum(pos))
  ##        Kneg <- K[,neg]%*%rep(1,sum(neg))
  ##        classmeans <- c(sum( Kpos * coef(ret)[pos] * as.vector(ymat)[pos]),sum( Kneg * coef(ret)[pos] * as.vector(ymat)[pos]))
  ##        kneg <- K%*%onesneg
  ##        kpos <- K%*%onespos
  ##        M <- (diag(dim(x)[1])- (1/dim(x)[1])*rep(1,dim(x)[1])%*%t(rep(1,dim(x)[1])))
  ##        kcentered <- M%*%solve(diag(dim(x)[1]) - tau*M%*%K%*%M)%*%M
          
  ##        prob.model(ret) <- list(Kpos=Kpos, Kneg=Kneg, kcentered=kcentered, classmeans=classmeans)
   ##     }
    }

  if(type(ret)=="regression")
    {
      if (nrows != nrow(x)) stop("x and y don't match.")
      
      ## initialize    
      p <- 0
      svindex <- NULL

      ymat <- y
            
      G <- csi(x, ymat, rank = rank ,kernel= kernel, delta = delta , tol = tol)
      
      GtP <- t(G) - matrix(rowSums(t(G))/dim(G)[1],dim(G)[2],dim(G)[1])
      Gtalpha <- (GtP)%*%G
      diag(Gtalpha) <-  diag(Gtalpha) + tau
      Gtalpha <- solve(Gtalpha) %*% GtP %*% ymat
      beta <- solve(t(G[1:dim(G)[2],]), Gtalpha)
      b <- colMeans(ymat) - colMeans(G%*%Gtalpha) 
          
      alpha(ret) <- beta
      ## nonzero alpha*y
      coef(ret) <- alpha(ret)
      ## store SV indexes from current problem for later use in predict
      alphaindex(ret) <- pivots(G)[1:dim(G)[2]]
      ## save the indexes from all the SV in a vector (use unique?)
      svindex <- alphaindex(ret)
      ## store betas in a vector 
      b(ret) <- b
      ##store C  in return object
      param(ret)$tau <- tau
    }
      
  kcall(ret) <- match.call()
  kernelf(ret) <- kernel
  ## param(ret) <- list(C=C, nu = nu, epsilon = epsilon)
  xmatrix(ret) <- x[alphaindex(ret),,drop = FALSE]
  ymatrix(ret) <- y
  nSV(ret)  <- length(svindex)
  if(nSV(ret)==0)
    stop("No Support Vectors found. You may want to change your parameters")
  fitted(ret)  <- if (fit)
    predict(ret, x) else NA

  scaling(ret) <- list(scaled = scaled, x.scale = x.scale)

  if (fit){
    if(type(ret)=="classification")
      error(ret) <- 1 - .classAgreement(table(y,as.integer(fitted(ret))))
    if(type(ret)=="regression")
      error(ret) <- drop(crossprod(fitted(ret) - y)/m)
  }
  
  cross(ret) <- -1
  if(cross == 1)
    cat("\n","cross should be >1 no cross-validation done!","\n","\n")
  else if (cross > 1)
    {
     
      cerror <- 0
      suppressWarnings(vgr<-split(sample(1:m,m),1:cross))
      for(i in 1:cross)
        {
          cind <-  unsplit(vgr[-i],factor(rep((1:cross)[-i],unlist(lapply(vgr[-i],length)))))
          cret <- lssvm(x[cind,],y[cind],type = type(ret),kernel=kernel,kpar = NULL,reduced = reduced,
      tau=tau,  tol=tol, rank = floor(rank/cross), delta = floor(delta/cross), scaled=FALSE, cross = 0, fit = FALSE)
          cres <- predict(cret, x[vgr[[i]],])
          cerror <- (1 - .classAgreement(table(y[vgr[[i]]],as.integer(cres))))/cross + cerror
        }
      cross(ret) <- cerror
    }
 
  return(ret)
})



## kernelMatrix interface

setMethod("lssvm",signature(x="kernelMatrix"),
function (x,
          y,
          type      = NULL,
          tau       = 0.01,
          tol       = 0.0001,
          rank      = floor(dim(x)[1]/3),
          delta      = 40,
          cross     = 0,
          fit       = TRUE,
          ...)
{ 
  ## subsetting and na-handling for matrices
  ret <- new("lssvm")
  
  if(!is.null(type))
    type(ret) <- match.arg(type,c("classification","regression"))
  
  if (is.null(type)) type(ret) <- if (is.factor(y)) "classification" else "regression"
  else type(ret) <- type

  ncols <- ncol(x)
  m <- nrows <- nrow(x)

  if(type(ret)=="classification")
    {
      if (!is.vector(y) && !is.factor (y)) stop("y must be a vector or a factor.")
      if (is(y,"vector")) {
        y <- as.matrix(y)
      if (nrows != nrow(y)) stop("x and y don't match.")}
      
      if (is.factor(y)) {
        lev(ret) <- levels (y)
        y <- as.integer (y)
        if (nrows != length(y)) stop("x and y don't match.")
      }
      else if (is.numeric(y))
        {
          y <- as.integer(y)
          lev(ret) <- unique (y)  
        }
      else 
        stop ("dependent variable has to be of factor or integer type for classification mode.")

      ## initialize    
      nclass(ret) <- length (unique(y))
      p <- 0

      svindex <- NULL
      
      ## create multidimensional y matrix
      yind <- t(matrix(1:nclass(ret),nclass(ret),m))
      ymat <- matrix(0, m, nclass(ret))
      ymat[yind==y] <- 1

         
      KP <- x - (1/m)*colSums(x) 
      beta <- solve((KP%*%x + m * tau * x), KP%*%ymat)
      b <- colMeans(ymat) - colMeans(x%*%beta)
      alphaindex(ret) <- 1:m
   
      
      alpha(ret) <- beta
      ## nonzero alpha*y
      coef(ret) <- alpha(ret)
      ## store SV indexes from current problem for later use in predict
  
      ## save the indexes from all the SV in a vector (use unique?)
      svindex <- alphaindex(ret)
      ## store betas in a vector 
      b(ret) <- b
      ##store C  in return object
      param(ret)$tau <- tau
    }

  if(type(ret)=="regression")
    {
      if (nrows != nrow(x)) stop("x and y don't match.")
      
      ## initialize    
      p <- 0

      svindex <- NULL

      ymat <- y
            
      G <- csi(x, ymat, rank = rank , delta = delta , tol = tol)
      
      GtP <- t(G) - matrix(rowSums(t(G))/dim(G)[1],dim(G)[2],dim(G)[1])
      Gtalpha <- (GtP)%*%G
      diag(Gtalpha) <-  diag(Gtalpha) + tau
      Gtalpha <- solve(Gtalpha) %*% GtP %*% ymat[pivots(G),,drop=FALSE]
      beta <- solve(t(G[1:dim(G)[2],]), Gtalpha)
      b <- colMeans(ymat) - colMeans(G%*%Gtalpha) 
          
      alpha(ret) <- beta
      ## nonzero alpha*y
      coef(ret) <- alpha(ret)
      ## store SV indexes from current problem for later use in predict
      alphaindex(ret) <- pivots(G)[1:dim(G)[2]]
      ## save the indexes from all the SV in a vector (use unique?)
      svindex <- alphaindex(ret)
      ## store betas in a vector 
      b(ret) <- b
      ##store C  in return object
      param(ret)$tau <- tau
    }
      
  kcall(ret) <- match.call()
  ## param(ret) <- list(C=C, nu = nu, epsilon = epsilon)
  xmatrix(ret) <- x
  ymatrix(ret) <- y
  kernelf(ret) <- "Kernel matrix used for training."
  nSV(ret)  <- length(svindex)
  if(nSV(ret)==0)
    stop("No Support Vectors found. You may want to change your parameters")
  fitted(ret)  <- if (fit)
    predict(ret, x) else NA


  if (fit){
    if(type(ret)=="classification")
      error(ret) <- 1 - .classAgreement(table(y,as.integer(fitted(ret))))
    if(type(ret)=="regression")
      error(ret) <- drop(crossprod(fitted(ret) - y)/m)
  }
  
  cross(ret) <- -1
  if(cross == 1)
    cat("\n","cross should be >1 no cross-validation done!","\n","\n")
  else if (cross > 1)
    {
     
      cerror <- 0
      suppressWarnings(vgr<-split(sample(1:m,m),1:cross))
      for(i in 1:cross)
        {
          cind <-  unsplit(vgr[-i],factor(rep((1:cross)[-i],unlist(lapply(vgr[-i],length)))))
          cret <- lssvm(x[cind,cind],y[cind],type = type(ret), tau=tau, rank = floor(rank/cross), delta = floor(delta/cross),  cross = 0, fit = FALSE)
          cres <- predict(cret, as.kernelMatrix(x[vgr[[i]], cind,drop = FALSE][,svindex,drop=FALSE]))
          cerror <- (1 - .classAgreement(table(y[vgr[[i]]],as.integer(cres))))/cross + cerror
        }
      cross(ret) <- cerror
    }
 
  return(ret)
})



## list interface
setMethod("lssvm",signature(x="list"),
function (x,
          y,
          scaled    = TRUE,
          kernel    = "stringdot",
          kpar      = list(length=4, lambda = 0.5),
          type      = NULL,
          tau       = 0.01,
          reduced   = TRUE, 
          tol       = 0.0001,
          rank      = floor(dim(x)[1]/3),
          delta      = 40,
          cross     = 0,
          fit       = TRUE,
          ...,
          subset)
{ 
  ## subsetting and na-handling for matrices
  ret <- new("lssvm")
  if (!missing(subset)) x <- x[subset]
  
  if(!is.null(type))
    type(ret) <- match.arg(type,c("classification","regression"))
  
  if (is.null(type)) type(ret) <- if (is.factor(y)) "classification" else "regression"
  else type(ret) <- type

  
  m <- nrows <- length(x)

 if(is.character(kernel)){
    kernel <- match.arg(kernel,c("rbfdot","polydot","tanhdot","vanilladot","laplacedot","besseldot","anovadot","splinedot","stringdot"))

    if(is.character(kpar))
       if(kernel == "tanhdot" || kernel == "vanilladot" || kernel == "polydot"|| kernel == "besseldot" || kernel== "anovadot"|| kernel=="splinedot" || kernel == "rbfdot" || kernel == "laplacedot" )
       {
         stop("List interface supports only the stringdot kernel.")
       }
     }
  
    if(is(kernel,"kernel")) 
    
    if(!is(kernel,"kernel"))
      {
        if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
        kernel <- do.call(kernel, kpar)
      }

    if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

  if(type(ret)=="classification")
    {
      if (!is.vector(y) && !is.factor (y)) stop("y must be a vector or a factor.")
      if (nrows != nrow(x)) stop("x and y don't match.")
      
      if (is.factor(y)) {
        lev(ret) <- levels (y)
        y <- as.integer (y)
      }
      else if (is.numeric(y))
        {
          y <- as.integer(y)
          lev(ret) <- unique (y)  
        }
      else 
        stop ("dependent variable has to be of factor or integer type for classification mode.")

      ## initialize    
      nclass(ret) <- length (unique(y))
      p <- 0

      svindex <- NULL
      
      ## create multidimensional y matrix
      yind <- t(matrix(1:nclass(ret),nclass(ret),m))
      ymat <- matrix(0, m, nclass(ret))
      ymat[yind==y] <- 1

      if(reduced == FALSE)
        {
          K <- kernelMatrix(kernel,x)

          KP <- K - (1/m)*colSums(K) 
          beta <- solve((KP%*%K + m * tau * K), KP%*%ymat)
          b <- colMeans(ymat) - colMeans(K%*%beta)
          alphaindex(ret) <- 1:m
        }
      else
        {
          G <- csi(x, ymat, rank = rank ,kernel= kernel, delta = delta , tol = tol)
          
          GtP <- t(G) - matrix(rowSums(t(G))/dim(G)[1],dim(G)[2],dim(G)[1])
          Gtalpha <- (GtP)%*%G
          diag(Gtalpha) <-  diag(Gtalpha) + tau
          Gtalpha <- solve(Gtalpha) %*% GtP %*% ymat[pivots(G),,drop=FALSE]
          beta <- solve(t(G[1:dim(G)[2],]), Gtalpha)
          b <- colMeans(ymat) - colMeans(G%*%Gtalpha) 
          alphaindex(ret) <- pivots(G)[1:dim(G)[2]]
        }
      
      alpha(ret) <- beta
      ## nonzero alpha*y
      coef(ret) <- alpha(ret)
      ## store SV indexes from current problem for later use in predict
  
      ## save the indexes from all the SV in a vector (use unique?)
      svindex <- alphaindex(ret)
      ## store betas in a vector 
      b(ret) <- b
      ##store C  in return object
      param(ret)$tau <- tau
    }

  if(type(ret)=="regression")
    {
      if (nrows != nrow(x)) stop("x and y don't match.")
      
      ## initialize    
      p <- 0

      svindex <- NULL

      ymat <- y
            
      G <- csi(x, ymat, rank = rank ,kernel= kernel, delta = delta , tol = tol)
      
      GtP <- t(G) - matrix(rowSums(t(G))/dim(G)[1],dim(G)[2],dim(G)[1])
      Gtalpha <- (GtP)%*%G
      diag(Gtalpha) <-  diag(Gtalpha) + tau
      Gtalpha <- solve(Gtalpha) %*% GtP %*% ymat[pivots(G),,drop=FALSE]
      beta <- solve(t(G[1:dim(G)[2],]), Gtalpha)
      b <- colMeans(ymat) - colMeans(G%*%Gtalpha) 
          
      alpha(ret) <- beta
      ## nonzero alpha*y
      coef(ret) <- alpha(ret)
      ## store SV indexes from current problem for later use in predict
      alphaindex(ret) <- pivots(G)[1:dim(G)[2]]
      ## save the indexes from all the SV in a vector (use unique?)
      svindex <- alphaindex(ret)
      ## store betas in a vector 
      b(ret) <- b
      ##store C  in return object
      param(ret)$tau <- tau
    }
      
  kcall(ret) <- match.call()
  kernelf(ret) <- kernel
  ## param(ret) <- list(C=C, nu = nu, epsilon = epsilon)
  xmatrix(ret) <- x[alphaindex(ret)]
  ymatrix(ret) <- y
  SVindex(ret) <- svindex
  nSV(ret)  <- length(svindex)
  if(nSV(ret)==0)
    stop("No Support Vectors found. You may want to change your parameters")
  fitted(ret)  <- if (fit)
    predict(ret, x) else NA

  

  if (fit){
    if(type(ret)=="classification")
      error(ret) <- 1 - .classAgreement(table(y,as.integer(fitted(ret))))
    if(type(ret)=="regression")
      error(ret) <- drop(crossprod(fitted(ret) - y)/m)
  }
  
  cross(ret) <- -1
  if(cross == 1)
    cat("\n","cross should be >1 no cross-validation done!","\n","\n")
  else if (cross > 1)
    {
     
      cerror <- 0
      suppressWarnings(vgr<-split(sample(1:m,m),1:cross))
      for(i in 1:cross)
        {
          cind <-  unsplit(vgr[-i],factor(rep((1:cross)[-i],unlist(lapply(vgr[-i],length)))))
          cret <- lssvm(x[cind,],y[cind],type = type(ret),kernel=kernel,kpar = NULL,reduced = reduced,
      tau=tau,  tol=tol, rank = floor(rank/cross), delta = floor(delta/cross), scaled=FALSE, cross = 0, fit = FALSE )
          cres <- predict(cret, x[vgr[[i]],])
          cerror <- (1 - .classAgreement(table(y[vgr[[i]]],as.integer(cres))))/cross + cerror
        }
      cross(ret) <- cerror
    }
 

  return(ret)
})


#**************************************************************#

setMethod("predict", signature(object = "lssvm"),
function (object, newdata, type = "response", coupler = "minpair")
{
  sc <- 0
  type <- match.arg(type,c("response","probabilities","decision"))
  if (missing(newdata) && type!="response")
    return(fitted(object))
  else if(missing(newdata))
    {
      newdata <- xmatrix(object)
      sc <- 1
    }
  
  ncols <- ncol(xmatrix(object))
  nrows <- nrow(xmatrix(object))
  oldco <- ncols

  if (!is.null(terms(object)))
    {
      if(!is.matrix(newdata))
        newdata <- model.matrix(delete.response(terms(object)), as.data.frame(newdata), na.action = n.action(object))
    }
  else
    newdata  <- if (is.vector(newdata)) t(t(newdata)) else as.matrix(newdata)

  newcols  <- 0
  newnrows <- nrow(newdata)
  newncols <- ncol(newdata)
  newco    <- newncols
    
  if (oldco != newco) stop ("test vector does not match model !")
  p<-0

  if (!is.null(scaling(object)$x.scale) && sc != 1)
    newdata[,scaling(object)$scaled] <-
      scale(newdata[,scaling(object)$scaled, drop = FALSE],
            center = scaling(object)$x.scale$"scaled:center",
            scale  = scaling(object)$x.scale$"scaled:scale"
            )

  if(is(newdata,"kernelMatrix"))
        res <- newdata %*% coef(object) - b(object)
      else
        res <- t(t(kernelMult(kernelf(object), newdata,xmatrix(object), alpha(object))) + b(object))

  if(type == "response" && type(object)=="classification"){
    predres <- max.col(res)
    return(factor (lev(object)[predres], levels = lev(object)))
  }
  
  if (type == "decision" || type(object)=="regression")
    return(res)

  if (type =="probabilities" && type(object)=="classification")
    {
      res - prob.model(object)$classmeans

      
    return(res)
  }
})

#****************************************************************************************#

setMethod("show","lssvm",
function(object){
  cat("Least Squares Support Vector Machine object of class \"lssvm\"","\n")
  cat("\n")
  cat(paste("problem type :",type(object), "\n"))
  cat(paste(" parameter : tau =",param(object)$tau, "\n"))
         
  cat("\n")
 show(kernelf(object))
  cat(paste("\nNumber of data points used for training :", nSV(object),"\n"))

  if(!is.null(fitted(object)))
    cat(paste("Training error :", round(error(object),6),"\n"))
  if(cross(object)!= -1)
    cat("Cross validation error :",round(cross(object),6),"\n")
})

##.partopro <- function(z,s,m){
##return(2*pi*(1/sqrt((1/z)+s^2))*exp(-(m^2)/(2*((1/z)+s^2))))
##}



