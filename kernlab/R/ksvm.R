## Support Vector Machines
## author : alexandros karatzoglou
## updated : 08.02.06

setGeneric("ksvm", function(x, ...) standardGeneric("ksvm"))
setMethod("ksvm",signature(x="formula"),
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
   ret <- ksvm(x, y, scaled = scaled, ...)
  kcall(ret) <- cl
  attr(Terms,"intercept") <- 0 ## no intercept
  terms(ret) <- Terms
  if (!is.null(attr(m, "na.action")))
    n.action(ret) <- attr(m, "na.action")
  return (ret)
})

setMethod("ksvm",signature(x="vector"),
function(x, ...)
          { x <- t(t(x))
            ret <- ksvm(x, ...)
            return(ret)
          })
    
setMethod("ksvm",signature(x="matrix"),
function (x,
          y         = NULL,
          scaled    = TRUE,
          type      = NULL,
          kernel    = "rbfdot",
          kpar      = "automatic",
          C         = 1,
          nu        = 0.2,
          epsilon   = 0.1,
          prob.model = FALSE,
          class.weights = NULL,
          cross     = 0,
          fit       = TRUE,
          cache     = 40,
          tol       = 0.001,
          shrinking = TRUE,
          ...
          ,subset 
         ,na.action = na.omit)
{ 
  ## Comment out sparse code, future impl. will be based on "Matrix"
  ##  sparse  <- inherits(x, "matrix.csr")
  ##  if (sparse) {
  ##    if (!require(SparseM))
  ##      stop("Need SparseM package for handling of sparse structures!")
  ##  }
  sparse <- FALSE
  
  if(is.character(kernel)){
    kernel <- match.arg(kernel,c("rbfdot","polydot","tanhdot","vanilladot","laplacedot","besseldot","anovadot","splinedot","matrix"))

    if(kernel == "matrix")
      if(dim(x)[1]==dim(x)[2])
        return(ksvm(as.kernelMatrix(x), y = y, type = type, C = C, nu = nu, epsilon  = epsilon, prob.model = prob.model, class.weights = class.weights, cross = cross, fit = fit, cache = cache, tol = tol, shrinking = shrinking, ...))
      else
        stop(" kernel matrix not square!")
    
    if(is.character(kpar))
      if((kernel == "tanhdot" || kernel == "vanilladot" || kernel == "polydot"|| kernel == "besseldot" || kernel== "anovadot"|| kernel=="splinedot") &&  kpar=="automatic" )
        {
          cat (" Setting default kernel parameters ","\n")
          kpar <- list()
        }
  }

  ## subsetting and na-handling for matrices
  ret <- new("ksvm")
  if (!missing(subset)) x <- x[subset,]
  if (is.null(y))
    x <- na.action(x)
  else {
    df <- na.action(data.frame(y, x))
    y <- df[,1]
    x <- as.matrix(df[,-1])
  }
  n.action(ret) <- na.action
  
 if (is.null(type)) type(ret) <- if (is.null(y)) "one-svc" else if (is.factor(y)) "C-svc" else "eps-svr"
  
  if(!is.null(type))
  type(ret) <- match.arg(type,c("C-svc",
                                "nu-svc",
                                "kbb-svc",
                                "spoc-svc",
                                "C-bsvc",
                                "one-svc",
                                "eps-svr",
                                "eps-bsvr",
                                "nu-svr"))

  ## ## scaling, subsetting, and NA handling
  ##  if (sparse) {
  ##    scale <- rep(FALSE, ncol(x))
  ##    if(!is.null(y)) na.fail(y)
  ##    x <- t(t(x)) ## make shure that col-indices are sorted
  ##  }

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
      if (is.numeric(y)&&(type(ret)!="C-svc"&&type(ret)!="nu-svc"&&type(ret)!="C-bsvc"&&type(ret)!="spoc-svc"&&type(ret)!="kbb-svc")) {
        y <- scale(y)
        y.scale <- attributes(y)[c("scaled:center","scaled:scale")]
        y <- as.vector(y)
      }
    }
  }
  ncols <- ncol(x)
  m <- nrows <- nrow(x)
  
  if (!is.function(kernel))
  if (!is.list(kpar)&&is.character(kpar)&&(class(kernel)=="rbfkernel" || class(kernel) =="laplacedot" || kernel == "laplacedot"|| kernel=="rbfdot")){
    kp <- match.arg(kpar,"automatic")
    if(kp=="automatic")
      kpar <- list(sigma=mean(sigest(x,scaled=FALSE)[c(1,3)]))
   #cat("Using automatic sigma estimation (sigest) for RBF or laplace kernel","\n") 
  }
  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }

  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

  if (!is(y,"vector") && !is.factor (y) & is(y,"matrix") & !(type(ret)=="one-svc")) stop("y must be a vector or a factor.")

  if(!(type(ret)=="one-svc"))
    if(is(y,"vector") | is(y,"factor") ) ym <- length(y) else if(is(y,"matrix")) ym <-  dim(y)[1] else stop("y must be a matrix or a vector")
    
    if ((type(ret) != "one-svc") && ym != m) stop("x and y don't match.")

  if(nu > 1|| nu <0) stop("nu must be between 0 an 1.")
   
  weightlabels <- NULL
  nweights <- 0
  weight <- 0
  wl <- 0
  ## in case of classification: transform factors into integers
  if (type(ret) == "one-svc") # one class classification --> set dummy
    y <- 1
  else
    if (is.factor(y)) {
      lev(ret) <- levels (y)
      y <- as.integer (y)
      if (!is.null(class.weights)) {
        weightlabels <- match (names(class.weights),lev(ret))
        if (any(is.na(weightlabels)))
          stop ("At least one level name is missing or misspelled.")
      }
    }
    else {
      if ((type(ret) =="C-svc" || type(ret) == "nu-svc" ||type(ret) == "C-bsvc" || type(ret) == "spoc-svc" || type(ret) == "kbb-svc") && any(as.integer (y) != y))
        stop ("dependent variable has to be of factor or integer type for classification mode.")

      if (type(ret) != "eps-svr" || type(ret) != "nu-svr"|| type(ret)!="eps-bsvr")
        lev(ret) <- sort(unique (y))
    }
 ## initialize    
  nclass(ret) <- length (unique(y))
  p <- 0
  K <- 0 
  svindex <- problem <- NULL
  sigma <- 0.1
  degree <- offset <- scale <- 1

  switch(is(kernel)[1],
         "rbfkernel" =
         {
           sigma <- kpar(kernel)$sigma
           ktype <- 2
         },
         "tanhkernel" =
         {
           sigma <- kpar(kernel)$scale
           offset <- kpar(kernel)$offset
           ktype <- 3
         },
         "polykernel" =
         {
           degree <- kpar(kernel)$degree
           sigma <- kpar(kernel)$scale
           offset <- kpar(kernel)$offset
           ktype <- 1
         },
         "vanillakernel" =
         {
           ktype <- 0
         },
	 "laplacekernel" =
	 {
	 ktype <- 5
	 sigma <- kpar(kernel)$sigma
	 },
         "besselkernel" =
         {
           ktype <- 6
           sigma <- kpar(kernel)$sigma
           degree <- kpar(kernel)$order
           offset <- kpar(kernel)$degree
         },
         "anovakernel" =
         {
           ktype <- 7
           sigma <- kpar(kernel)$sigma
           degree <-  kpar(kernel)$degree
         },
         "splinekernel" =
         {
           ktype <- 8
         },
         {
           ktype <- 4
         }
         )
  prior(ret) <- list(NULL)

## C classification
  if(type(ret) == "C-svc"){

    indexes <- lapply(sort(unique(y)), function(kk) which(y == kk))
    for (i in 1:(nclass(ret)-1)) {
      jj <- i+1
      for(j in jj:nclass(ret)) {
        p <- p+1
        ## prepare the data
        li <- length(indexes[[i]])
        lj <- length(indexes[[j]])
      
        if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
          {
            yd <- c(rep(-1,li),rep(1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(j,i)]]
            wl <- c(1,0)
            nweights <- 2
          }
          }
        else
          {
            yd <- c(rep(1,li),rep(-1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(i,j)]]
            wl <- c(0,1)
            nweigths <- 2
          }
          }
       
        boolabel <- yd >= 0
        prior1 <- sum(boolabel)

        md <- length(yd)
        prior0 <- md - prior1
        prior(ret)[[p]] <- list(prior1 = prior1, prior0 = prior0) 

        if(ktype==4)
          K <- kernelMatrix(kernel,x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE])
               
        resv <- .Call("smo_optim",
                      as.double(t(x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE])),
                      as.integer(li+lj),
                      as.integer(ncol(x)),
                      as.double(yd),
                      as.double(K),
                      
                      as.integer(if (sparse) x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE]@ia else 0),
                      as.integer(if (sparse) x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE]@ja else 0),
                      as.integer(sparse),
                      
                      as.double(matrix(rep(-1,m))), ##linear term
                      as.integer(ktype),
                      as.integer(0), 
                      as.double(C),
                      as.double(nu),
                      as.double(epsilon),
                      as.double(sigma),
                      as.integer(degree),
                      as.double(offset),
                      as.integer(wl), ##weightlabel
                      as.double(weight),
                      as.integer(nweights),
                      as.double(cache), 
                      as.double(tol),
                      as.integer(shrinking),
                      PACKAGE="kernlab")

        reind <- sort(c(indexes[[i]],indexes[[j]]),method="quick",index.return=TRUE)$ix
        tmpres <- resv[c(-(li+lj+1),-(li+lj+2))][reind]
        ## alpha
        svind <- tmpres > 0
        alpha(ret)[p] <- list(tmpres[svind])
        ## coefficients alpha*y
        coef(ret)[p] <- list(alpha(ret)[[p]]*yd[reind][svind])
        ## store SV indexes from current problem for later use in predict
        alphaindex(ret)[p] <- list(c(indexes[[i]],indexes[[j]])[reind][svind])
        ## store Support Vectors
        xmatrix(ret)[p] <- list(x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE][reind,,drop=FALSE][svind, ,drop=FALSE])
        ## save the indexes from all the SV in a vector (use unique?)
        svindex <- c(svindex,alphaindex(ret)[[p]])
        ## store betas in a vector 
        b(ret) <- c(b(ret), resv[li+lj+1])
	## store objective function values in a vector
	obj(ret) <- c(obj(ret), resv[li+lj+2])
        ## used to reconstruct indexes for the patterns matrix x from "indexes" (really usefull ?)
        problem[p] <- list(c(i,j))
        ##store C  in return object
        param(ret)$C <- C
        ## margin(ret)[p] <- (min(kernelMult(kernel,xd[1:li,],,alpha(ret)[[p]][1:li])) - max(kernelMult(kernel,xd[li:(li+lj),],,alpha(ret)[[p]][li:(li+lj)])))/2
      }
    }
  } 

## nu classification
if(type(ret) == "nu-svc"){
  indexes <- lapply(sort(unique(y)), function(kk) which(y == kk))
    for (i in 1:(nclass(ret)-1)) {
      jj <- i+1
      for(j in jj:nclass(ret)) {
        p <- p+1
       ##prepare data
        li <- length(indexes[[i]])
        lj <- length(indexes[[j]])

        if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
          {
            yd <- c(rep(-1,li),rep(1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(j,i)]]
            wl <- c(1,0)
            nweights <- 2
          }
          }
        else
          {
            yd <- c(rep(1,li),rep(-1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(i,j)]]
            wl <- c(0,1)
            nweigths <- 2
          }
          }
       
        boolabel <- yd >= 0
        prior1 <- sum(boolabel)
        md <- length(yd)
        prior0 <- md - prior1
        prior(ret)[[p]] <- list(prior1 = prior1, prior0 = prior0)

        if(ktype==4)
             K <- kernelMatrix(kernel,x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE])
        
        resv <- .Call("smo_optim",
                      as.double(t(x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE])),
                      as.integer(li+lj),
                      as.integer(ncol(x)),
                      as.double(yd),
                      as.double(K),
                      as.integer(if (sparse) x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE]@ia else 0),
                      as.integer(if (sparse) x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE]@ja else 0),
                      as.integer(sparse),
                      
                      as.double(matrix(rep(-1,m))), #linear term
                      as.integer(ktype),
                      as.integer(1),
                      as.double(C),
                      as.double(nu),
                      as.double(epsilon),
                      as.double(sigma),
                      as.integer(degree),
                      as.double(offset),
                      as.integer(wl), #weightlabl.
                      as.double(weight),
                      as.integer(nweights),
                      as.double(cache),
                      as.double(tol), 
                      as.integer(shrinking),
                      PACKAGE="kernlab")
        
        reind <- sort(c(indexes[[i]],indexes[[j]]),method="quick",index.return=TRUE)$ix
        tmpres <- resv[c(-(li+lj+1),-(li+lj+2))][reind]
        svind <- tmpres != 0
        alpha(ret)[p] <- coef(ret)[p] <- list(tmpres[svind])
        ##store SV indexes from current problem for later use in predict
       	alphaindex(ret)[p] <- list(c(indexes[[i]],indexes[[j]])[reind][svind])
        ## store Support Vectors
        xmatrix(ret)[p] <- list(x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE][reind,,drop=FALSE][svind,,drop=FALSE])
        ##save the indexes from all the SV in a vector (use unique!)
        svindex <- c(svindex,alphaindex(ret)[[p]])
        ## store betas in a vector 
        b(ret) <- c(b(ret), resv[li+lj+1])
	## store objective function values in a vector
	obj(ret) <- c(obj(ret), resv[li+lj+2])
        ## used to reconstruct indexes for the patterns matrix x from "indexes"
        problem[p] <- list(c(i,j))
        param(ret)$nu <- nu
        ## margin(ret)[p] <- (min(kernelMult(kernel,xd[1:li,],,alpha(ret)[[p]][1:li])) - max(kernelMult(kernel,xd[li:(li+lj),],,alpha(ret)[[p]][li:(li+lj)])))/2
      }
    }
}  

## Bound constraint C classification
  if(type(ret) == "C-bsvc"){
     if(!is.null(class.weights))
     weightedC <- class.weights[weightlabels] * rep(C,nclass(ret))
    else
      weightedC <- rep(C,nclass(ret)) 
    indexes <- lapply(sort(unique(y)), function(kk) which(y == kk))
    for (i in 1:(nclass(ret)-1)) {
      jj <- i+1
      for(j in jj:nclass(ret)) {
        p <- p+1
        ##prepare data
        li <- length(indexes[[i]])
        lj <- length(indexes[[j]])

        if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
          {
            yd <- c(rep(-1,li),rep(1,lj))
            if(!is.null(class.weights)){
              weight <- class.weights[weightlabels[c(j,i)]]
              wl <- c(1,0)
              nweights <- 2
            }
          }
        else
          {
            yd <- c(rep(1,li),rep(-1,lj))
            if(!is.null(class.weights)){
              weight <- class.weights[weightlabels[c(i,j)]]
              wl <- c(0,1)
              nweigths <- 2
            }
          }
       
        boolabel <- yd >= 0
        prior1 <- sum(boolabel)
        md <- length(yd)
        prior0 <- md - prior1
        prior(ret)[[p]] <- list(prior1 = prior1, prior0 = prior0) 

           if(ktype==4)
             K <- kernelMatrix(kernel,x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE])
        
        resv <- .Call("tron_optim",
                      as.double(t(x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE])),
                      as.integer(li+lj),
                      as.integer(ncol(x)),
                      as.double(yd),
                      as.double(K),
                      as.integer(if (sparse) x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE]@ia else 0),
                      as.integer(if (sparse) x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE]@ja else 0),
                      as.integer(sparse),
                      as.integer(2),
                      as.double(0), ##countc
                      as.integer(ktype),
                      as.integer(5), 
                      as.double(C),
                      as.double(epsilon),
                      as.double(sigma),
                      as.integer(degree),
                      as.double(offset),
                      as.double(1),  ##  cost value of alpha seeding
                      as.double(2),  ## step value of alpha seeding
                      as.integer(wl), ##weightlabel
                      as.double(weight),
                      as.integer(nweights),
                      as.double(weightedC),
                      as.double(cache), 
                      as.double(tol),
                      as.integer(10), ##qpsize
                      as.integer(shrinking),
                      PACKAGE="kernlab")
        
        reind <- sort(c(indexes[[i]],indexes[[j]]),method="quick",index.return=TRUE)$ix
        svind <- resv[-(li+lj+1)][reind] > 0
        alpha(ret)[p] <- list(resv[-(li+lj+1)][reind][svind])
        ## nonzero alpha*y
        coef(ret)[p] <- list(alpha(ret)[[p]] * yd[reind][svind])
        ## store SV indexes from current problem for later use in predict
       	alphaindex(ret)[p] <- list(c(indexes[[i]],indexes[[j]])[reind][svind])
        ## store Support Vectors
        xmatrix(ret)[p] <- list(x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE][reind,,drop = FALSE][svind,,drop = FALSE])
        ## save the indexes from all the SV in a vector (use unique?)
        svindex <- c(svindex,alphaindex(ret)[[p]])
        ## store betas in a vector 
        b(ret) <- - sapply(coef(ret),sum) 
	## store obj. values in vector 
	obj(ret) <- c(obj(ret), resv[(li+lj+1)])
        ## used to reconstruct indexes for the patterns matrix x from "indexes" (really usefull ?)
        problem[p] <- list(c(i,j))
        ##store C  in return object
        param(ret)$C <- C
##        margin(ret)[p] <- (min(kernelMult(kernel,xd[1:li,],,alpha(ret)[[p]][1:li])) - max(kernelMult(kernel,xd[li:(li+lj),],,alpha(ret)[[p]][li:(li+lj)])))/2
      }
    }
  } 

## SPOC multiclass classification 
if(type(ret) =="spoc-svc")
  {
    if(!is.null(class.weights))
     weightedC <- class.weights[weightlabels] * rep(C,nclass(ret))
    else
      weightedC <- rep(C,nclass(ret)) 
    yd <- sort(y,method="quick", index.return = TRUE)
    xd <- matrix(x[yd$ix,],nrow=dim(x)[1])
    count <- 0

       if(ktype==4)
          K <- kernelMatrix(kernel,x)
    
    resv <- .Call("tron_optim",
                  as.double(t(xd)),
                  as.integer(nrow(xd)),
                  as.integer(ncol(xd)),
                  as.double(rep(yd$x-1,2)),
                  as.double(K),
                  as.integer(if (sparse) xd@ia else 0),
                  as.integer(if (sparse) xd@ja else 0),
                  as.integer(sparse),
                  as.integer(nclass(ret)),
                  as.integer(count),
                  as.integer(ktype),
                  as.integer(7), 
                  as.double(C),
                  as.double(epsilon),
                  as.double(sigma),
                  as.integer(degree),
                  as.double(offset),
                  as.double(C), 
                  as.double(2), #Cstep
                  as.integer(0), #weightlabel
                  as.double(0),
                  as.integer(0),
                  as.double(weightedC),
                  as.double(cache), 
                  as.double(tol),
                  as.integer(10), #qpsize
                  as.integer(shrinking),
                  PACKAGE="kernlab")
    
    reind <- sort(yd$ix,method="quick",index.return=TRUE)$ix
    alpha(ret) <- t(matrix(resv[-(nclass(ret)*nrow(xd) + 1)],nclass(ret)))[reind,,drop=FALSE]
    coef(ret) <- lapply(1:nclass(ret), function(x) alpha(ret)[,x][alpha(ret)[,x]!=0])
    names(coef(ret)) <- lev(ret)
    alphaindex(ret) <-  lapply(sort(unique(y)), function(x) which(alpha(ret)[,x]!=0))
    xmatrix(ret) <- x
    obj(ret) <- resv[(nclass(ret)*nrow(xd) + 1)]
    names(alphaindex(ret)) <- lev(ret)
    svindex <- which(rowSums(alpha(ret)!=0)!=0)
    b(ret) <- 0
    param(ret)$C <- C
  }

## KBB multiclass classification  
if(type(ret) =="kbb-svc")
  {
    if(!is.null(class.weights))
      weightedC <- weightlabels * rep(C,nclass(ret))
    else
      weightedC <- rep(C,nclass(ret)) 
    yd <- sort(y,method="quick", index.return = TRUE)
    x <- x[yd$ix,,drop=FALSE]
    count <-  sapply(unique(yd$x), function(c) length(yd$x[yd$x==c]))
    if(ktype==4)
      K <- kernelMatrix(kernel,x)
    resv <- .Call("tron_optim",
                  as.double(t(x)),
                  as.integer(nrow(x)),
                  as.integer(ncol(x)),
                  as.double(yd$x-1),
                  as.double(K),
                  as.integer(if (sparse) x@ia else 0),
                  as.integer(if (sparse) x@ja else 0),
                  as.integer(sparse),
                  as.integer(nclass(ret)),
                  as.integer(count),
                  as.integer(ktype),
                  as.integer(8),
                  as.double(C),
                  as.double(epsilon),
                  as.double(sigma),
                  as.integer(degree),
                  as.double(offset),
                  as.double(C), #Cbegin
                  as.double(2), #Cstep
                  as.integer(0), #weightlabl.
                  as.double(0),
                  as.integer(0),
                  as.double(weightedC),
                  as.double(cache),
                  as.double(tol),
                  as.integer(10), #qpsize
                  as.integer(shrinking),
                  PACKAGE="kernlab")

    reind <- sort(yd$ix,method="quick",index.return=TRUE)$ix
    alpha(ret) <- matrix(resv[-(nrow(x)*(nclass(ret)-1)+1)],nrow(x))[reind,,drop=FALSE]
    xmatrix(ret) <- x<- x[reind,,drop=FALSE]
    coef(ret) <-  lapply(1:(nclass(ret)-1), function(x) alpha(ret)[,x][alpha(ret)[,x]!=0])
    alphaindex(ret) <-  lapply(sort(unique(y)), function(x) which((y == x) & (rowSums(alpha(ret))!=0)))
    svindex <- which(rowSums(alpha(ret)!=0)!=0)
    b(ret) <- - sapply(coef(ret),sum)
    obj(ret) <- resv[(nrow(x)*(nclass(ret)-1)+1)]
    param(ret)$C <- C
  }

  ## Novelty detection
  if(type(ret) =="one-svc")
  {
    if(ktype==4)
      K <- kernelMatrix(kernel,x)
       
    resv <- .Call("smo_optim",
                  as.double(t(x)),
                  as.integer(nrow(x)),
                  as.integer(ncol(x)),
                  as.double(matrix(rep(1,m))),
                  as.double(K),
                  as.integer(if (sparse) x@ia else 0),
                  as.integer(if (sparse) x@ja else 0),
                  as.integer(sparse),
                  as.double(matrix(rep(-1,m))),
                  as.integer(ktype),
                  as.integer(2),
                  as.double(C),
                  as.double(nu),
                  as.double(epsilon),
                  as.double(sigma),
                  as.integer(degree),
                  as.double(offset),
                  as.integer(0), #weightlabl.
                  as.double(0),
                  as.integer(0),
                  as.double(cache),
                  as.double(tol),
                  as.integer(shrinking),
                  PACKAGE="kernlab")

       tmpres <- resv[c(-(m+1),-(m+2))]
       alpha(ret) <- coef(ret) <- tmpres[tmpres != 0]
       svindex <-  alphaindex(ret) <- which(tmpres != 0) 
       xmatrix(ret) <- x[svindex,,drop=FALSE]
       b(ret) <- resv[(m+1)]
       obj(ret) <- resv[(m+2)]
       param(ret)$nu <- nu
  }

  ## epsilon regression
  if(type(ret) =="eps-svr")
    {
      if(ktype==4)
        K <- kernelMatrix(kernel,x)
      
      resv <- .Call("smo_optim",
                    as.double(t(x)),
                    as.integer(nrow(x)),
                    as.integer(ncol(x)),
                    as.double(y),
                    as.double(K),
                    as.integer(if (sparse) x@ia else 0),
                    as.integer(if (sparse) x@ja else 0),
                    as.integer(sparse),
                    as.double(matrix(rep(-1,m))),
                    as.integer(ktype),
                    as.integer(3),
                    as.double(C),
                    as.double(nu),
                    as.double(epsilon),
                    as.double(sigma),
                    as.integer(degree),
                    as.double(offset),
                    as.integer(0), #weightlabl.
                    as.double(0),
                    as.integer(0),
                    as.double(cache), 
                    as.double(tol), 
                    as.integer(shrinking), 
                    PACKAGE="kernlab")
      tmpres <- resv[c(-(m+1),-(m+2))]
      alpha(ret) <- coef(ret) <- tmpres[tmpres != 0]
      svindex <-  alphaindex(ret) <- which(tmpres != 0) 
      xmatrix(ret) <- x[svindex, ,drop=FALSE]
      b(ret) <- resv[(m+1)]
      obj(ret) <- resv[(m+2)]
      param(ret)$epsilon <- epsilon
      param(ret)$C <- C
    }

  ## nu regression
  if(type(ret) =="nu-svr")
    {
      if(ktype==4)
        K <- kernelMatrix(kernel,x)
      
      resv <- .Call("smo_optim",
                    as.double(t(x)),
                    as.integer(nrow(x)),
                    as.integer(ncol(x)),
                    as.double(y),
                    as.double(K),
                    as.integer(if (sparse) x@ia else 0),
                    as.integer(if (sparse) x@ja else 0),
                    as.integer(sparse),
                    as.double(matrix(rep(-1,m))),
                    as.integer(ktype),
                    as.integer(4),
                    as.double(C),
                    as.double(nu),
                    as.double(epsilon),
                    as.double(sigma),
                    as.integer(degree),
                    as.double(offset),
                    as.integer(0),
                    as.double(0),
                    as.integer(0),
                    as.double(cache), 
                    as.double(tol), 
                    as.integer(shrinking), 
                    PACKAGE="kernlab")
      tmpres <- resv[c(-(m+1),-(m+2))]
      alpha(ret) <- coef(ret) <- tmpres[tmpres!=0]
      svindex <-  alphaindex(ret) <- which(tmpres != 0)
      xmatrix(ret) <- x[svindex,,drop=FALSE]
      b(ret) <- resv[(m+1)]
      obj(ret) <- resv[(m+2)]
      param(ret)$epsilon <- epsilon
      param(ret)$nu <- nu
    }

  ## bound constraint eps regression
  if(type(ret) =="eps-bsvr")
    {
      if(ktype==4)
        K <- kernelMatrix(kernel,x)
      
      resv <- .Call("tron_optim",
                    as.double(t(x)),
                    as.integer(nrow(x)),
                    as.integer(ncol(x)),
                    as.double(y),
                    as.double(K),
                    as.integer(if (sparse) x@ia else 0),
                    as.integer(if (sparse) x@ja else 0),
                    as.integer(sparse),
                    as.integer(2),
                    as.integer(0),
                    as.integer(ktype),
                    as.integer(6),
                    as.double(C),
                    as.double(epsilon),
                    as.double(sigma),
                    as.integer(degree),
                    as.double(offset),
                    as.double(1),  #Cbegin
                    as.double(2), #Cstep
                    as.integer(0), #weightlabl.
                    as.double(0),
                    as.integer(0),
                    as.double(0),
                    as.double(cache), 
                    as.double(tol),
                    as.integer(10), #qpsize
                    as.integer(shrinking), 
                   PACKAGE="kernlab")
      tmpres <- resv[-(m + 1)]
      alpha(ret) <- coef(ret) <- tmpres[tmpres!=0]
      svindex <-  alphaindex(ret) <- which(tmpres != 0)
      xmatrix(ret) <- x[svindex,,drop=FALSE]
      b(ret) <- -sum(alpha(ret))
      obj(ret) <-  resv[(m + 1)]
      param(ret)$epsilon <- epsilon
      param(ret)$C <- C
  }

  
  kcall(ret) <- match.call()
  kernelf(ret) <- kernel
  ymatrix(ret) <- y
  SVindex(ret) <- sort(unique(svindex),method="quick")
  nSV(ret)  <- length(unique(svindex))
  if(nSV(ret)==0)
    stop("No Support Vectors found. You may want to change your parameters")
  
  fitted(ret)  <- if (fit)
    predict(ret, x) else NULL


  if(any(scaled))
    scaling(ret) <- list(scaled = scaled, x.scale = x.scale, y.scale = y.scale)

  
  if (fit){
    if(type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="spoc-svc"||type(ret)=="kbb-svc"||type(ret)=="C-bsvc")
      error(ret) <- 1 - .classAgreement(table(y,as.integer(fitted(ret))))
    if(type(ret)=="one-svc")
      error(ret) <- sum(!fitted(ret))/m
    if(type(ret)=="eps-svr"||type(ret)=="nu-svr"||type(ret)=="eps-bsvr"){
      if (!is.null(scaling(ret)$y.scale)){
        scal <- scaling(ret)$y.scale$"scaled:scale"
        fitted(ret) <- fitted(ret) # / scaling(ret)$y.scale$"scaled:scale" + scaling(ret)$y.scale$"scaled:center"
      }
      else
        scal <- 1
     
      error(ret) <- drop(crossprod(fitted(ret) - y)/m)
    }
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
          if(type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="spoc-svc"||type(ret)=="kbb-svc"||type(ret)=="C-bsvc")
            {
              if(is.null(class.weights))
                cret <- ksvm(x[cind,],y[cind],type = type(ret),kernel=kernel,kpar = NULL, C=C, nu=nu, tol=tol, scaled=FALSE, cross = 0, fit = FALSE ,cache = cache)
              else
                cret <- ksvm(x[cind,],as.factor(lev(ret)[y[cind]]),type = type(ret),kernel=kernel,kpar = NULL, C=C, nu=nu, tol=tol, scaled=FALSE, cross = 0, fit = FALSE, class.weights = class.weights,cache = cache)
               cres <- predict(cret, x[vgr[[i]],,drop=FALSE])
            cerror <- (1 - .classAgreement(table(y[vgr[[i]]],as.integer(cres))))/cross + cerror
            }
          if(type(ret)=="one-svc")
            {
              cret <- ksvm(x[cind,],type=type(ret),kernel=kernel,kpar = NULL,C=C,nu=nu,epsilon=epsilon,tol=tol,scaled=FALSE, cross = 0, fit = FALSE, cache = cache, prob.model = FALSE)
              cres <- predict(cret, x[vgr[[i]],, drop=FALSE])
              cerror <- (1 - sum(cres)/length(cres))/cross + cerror
            }
           
          if(type(ret)=="eps-svr"||type(ret)=="nu-svr"||type(ret)=="eps-bsvr")
            {
              cret <- ksvm(x[cind,],y[cind],type=type(ret),kernel=kernel,kpar = NULL,C=C,nu=nu,epsilon=epsilon,tol=tol,scaled=FALSE, cross = 0, fit = FALSE, cache = cache, prob.model = FALSE)
              cres <- predict(cret, x[vgr[[i]],,drop=FALSE])
              if (!is.null(scaling(ret)$y.scale))
                scal <- scaling(ret)$y.scale$"scaled:scale"
              else
                scal <- 1
              cerror <- drop((scal^2)*crossprod(cres - y[vgr[[i]]])/m) + cerror
            }
        }
      cross(ret) <- cerror
    }

  prob.model(ret) <- list(NULL)
  
  if(prob.model)
    {
      if(type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="C-bsvc")
        {
          p <- 0
          for (i in 1:(nclass(ret)-1)) {
            jj <- i+1
            for(j in jj:nclass(ret)) {
              p <- p+1
              ##prepare data
              li <- length(indexes[[i]])
              lj <- length(indexes[[j]])

              if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
                {
                  yd <- c(rep(-1,li),rep(1,lj))
                  if(!is.null(class.weights)){
                    weight <- weightlabels[c(j,i)]
                    wl <- c(1,0)
                    nweights <- 2
                  }
                }
              else
                {
                  yd <- c(rep(1,li),rep(-1,lj))
                  if(!is.null(class.weights)){
                    weight <- weightlabels[c(i,j)]
                    wl <- c(0,1)
                    nweigths <- 2
                  }
                }
              m <- li+lj
              suppressWarnings(vgr <- split(c(sample(1:li,li),sample((li+1):(li+lj),lj)),1:3)) 
              pres <- yres <- NULL
              for(k in 1:3)
                {
                  cind <- unsplit(vgr[-k],factor(rep((1:3)[-k],unlist(lapply(vgr[-k],length)))))
                  if(is.null(class.weights))
                    cret <- ksvm(x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE][cind,],yd[cind],type = type(ret),kernel=kernel,kpar = NULL, C=C, nu=nu, tol=tol, scaled=FALSE, cross = 0, fit = FALSE ,cache = cache, prob.model = FALSE)
                  else
                    cret <- ksvm(x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE][cind,],as.factor(lev(ret)[y[c(indexes[[i]],indexes[[j]])][cind]]),type = type(ret),kernel=kernel,kpar = NULL, C=C, nu=nu, tol=tol, scaled=FALSE, cross = 0, fit = FALSE, class.weights = class.weights,cache = cache, prob.model = FALSE)
                  
                                    
                  yres <- c(yres, yd[vgr[[k]]])
                  pres <- rbind(pres, predict(cret, x[c(indexes[[i]],indexes[[j]]), ,drop=FALSE][vgr[[k]],],type="decision"))
                }
              prob.model(ret)[[p]] <- .probPlatt(pres,yres)
            }
          }
        }
      if(type(ret) == "eps-svr"||type(ret) == "nu-svr"||type(ret)=="eps-bsvr"){
        suppressWarnings(vgr<-split(sample(1:m,m),1:3))
        pres <- NULL
        for(i in 1:3)
          {
            cind <- unsplit(vgr[-i],factor(rep((1:3)[-i],unlist(lapply(vgr[-i],length)))))

            cret <- ksvm(x[cind,],y[cind],type=type(ret),kernel=kernel,kpar = NULL,C=C,nu=nu,epsilon=epsilon,tol=tol,scaled=FALSE, cross = 0, fit = FALSE, cache = cache, prob.model = FALSE)
            cres <- predict(cret, x[vgr[[i]],])
            if (!is.null(scaling(ret)$y.scale))
              cres <- cres * scaling(ret)$y.scale$"scaled:scale" + scaling(ret)$y.scale$"scaled:center"
            pres <- rbind(pres, cres)
          }
        pres[abs(pres) > (5*sd(pres))] <- 0
        prob.model(ret) <- list(sum(abs(pres))/dim(pres)[1])
      }
    }

  return(ret)
})




## kernelmatrix interface

setMethod("ksvm",signature(x="kernelMatrix"),
function (x,
          y         = NULL,
          type      = NULL,
          C         = 1,
          nu        = 0.2,
          epsilon   = 0.1,
          prob.model = FALSE,
          class.weights = NULL,
          cross     = 0,
          fit       = TRUE,
          cache     = 40,
          tol       = 0.001,
          shrinking = TRUE,
          ...)
{ 
  sparse <- FALSE
  ## subsetting and na-handling for matrices
  ret <- new("ksvm")

 if (is.null(type)) type(ret) <- if (is.null(y)) "one-svc" else if (is.factor(y)) "C-svc" else "eps-svr"
  
  if(!is.null(type))
  type(ret) <- match.arg(type,c("C-svc",
                                "nu-svc",
                                "kbb-svc",
                                "spoc-svc",
                                "C-bsvc",
                                "one-svc",
                                "eps-svr",
                                "eps-bsvr",
                                "nu-svr"))


  ncols <- ncol(x)
  m <- nrows <- nrow(x)
  

  if (!is(y,"vector") && !is.factor (y) & !is(y,"matrix") & !(type(ret)=="one-svc")) stop("y must be a vector or a factor.")

  if(!(type(ret)=="one-svc"))
    if(is(y,"vector") | is(y,"factor")) ym <- length(y) else if(is(y,"matrix")) ym <-  dim(y)[1] else stop("y must be a matrix or a vector")
    
    if ((type(ret) != "one-svc") && ym != m) stop("x and y don't match.")

  if(nu > 1|| nu <0) stop("nu must be between 0 an 1.")
  
  weightlabels <- NULL
  nweights <- 0
  weight <- 0
  wl <- 0
  ## in case of classification: transform factors into integers
  if (type(ret) == "one-svc") # one class classification --> set dummy
    y <- 1
  else
    if (is.factor(y)) {
      lev(ret) <- levels (y)
      y <- as.integer (y)
      if (!is.null(class.weights)) {
        if (is.null(names (class.weights)))
          stop ("Weights have to be specified along with their according level names !")
        weightlabels <- match (names(class.weights),lev(ret))
        if (any(is.na(weightlabels)))
          stop ("At least one level name is missing or misspelled.")
      }
    }
    else {
      if ((type(ret) =="C-svc" || type(ret) == "nu-svc" ||type(ret) == "C-bsvc" || type(ret) == "spoc-svc" || type(ret) == "kbb-svc") && any(as.integer (y) != y))
        stop ("dependent variable has to be of factor or integer type for classification mode.")

      if (type(ret) != "eps-svr" || type(ret) != "nu-svr"|| type(ret)!="eps-bsvr")
        lev(ret) <- sort(unique (y))
    }
 ## initialize    
  nclass(ret) <- length (unique(y))
  p <- 0
  svindex <- problem <- NULL
  sigma <- 0.1
  degree <- offset <- scale <- 1

  ktype <- 4
  
  prior(ret) <- list(NULL)

## C classification
  if(type(ret) == "C-svc"){
    indexes <- lapply(sort(unique(y)), function(kk) which(y == kk))
    for (i in 1:(nclass(ret)-1)) {
      jj <- i+1
      for(j in jj:nclass(ret)) {
        p <- p+1
        ##prepare data
        li <- length(indexes[[i]])
        lj <- length(indexes[[j]])
        
        if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
          {
            yd <- c(rep(-1,li),rep(1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(j,i)]]
            wl <- c(1,0)
            nweights <- 2
          }
          }
        else
          {
            yd <- c(rep(1,li),rep(-1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(i,j)]]
            wl <- c(0,1)
            nweigths <- 2
          }
          }
       
        boolabel <- yd >= 0
        prior1 <- sum(boolabel)

        md <- length(yd)
        prior0 <- md - prior1
        prior(ret)[[p]] <- list(prior1 = prior1, prior0 = prior0) 

        xdd <- matrix(1,li+lj,1)
        
        resv <- .Call("smo_optim",
                      as.double(t(xdd)),
                      as.integer(nrow(xdd)),
                      as.integer(ncol(xdd)),
                      as.double(yd),
                      as.double(as.vector(x[c(indexes[[i]],indexes[[j]]),c(indexes[[i]],indexes[[j]]),drop=FALSE])),
                      
                      as.integer(if (sparse) x[c(indexes[[i]],indexes[[j]]),c(indexes[[i]],indexes[[j]]),drop=FALSE]@ia else 0),
                      as.integer(if (sparse) x[c(indexes[[i]],indexes[[j]]),c(indexes[[i]],indexes[[j]]),drop=FALSE]@ja else 0),
                      as.integer(sparse),
                      
                      as.double(matrix(rep(-1,m))), ##linear term
                      as.integer(ktype),
                      as.integer(0), 
                      as.double(C),
                      as.double(nu),
                      as.double(epsilon),
                      as.double(sigma),
                      as.integer(degree),
                      as.double(offset),
                      as.integer(wl), ##weightlabel
                      as.double(weight),
                      as.integer(nweights),
                      as.double(cache), 
                      as.double(tol),
                      as.integer(shrinking),
                      PACKAGE="kernlab")
        
        reind <- sort(c(indexes[[i]],indexes[[j]]),method="quick",index.return=TRUE)$ix        
        tmpres <- resv[c(-(li+lj+1),-(li+lj+2))][reind]
        ## alpha
        svind <- tmpres > 0
        alpha(ret)[p] <- list(tmpres[svind])
        ## coefficients alpha*y
        coef(ret)[p] <- list(alpha(ret)[[p]]*yd[reind][svind])
        ## store SV indexes from current problem for later use in predict
        alphaindex(ret)[p] <- list(c(indexes[[i]],indexes[[j]])[reind][svind])
        ## store Support Vectors
        ## xmatrix(ret)[p] <- list(xd[svind, svind,drop=FALSE])
        ## save the indexes from all the SV in a vector (use unique?)
        svindex <- c(svindex,alphaindex(ret)[[p]])
        ## store betas in a vector 
        b(ret) <- c(b(ret), resv[li+lj+1])
	## store objective function values in vector 
	obj(ret) <- c(obj(ret), resv[li+lj+2])
        ## used to reconstruct indexes for the patterns matrix x from "indexes" (really usefull ?)
        problem[p] <- list(c(i,j))
        ##store C  in return object
        param(ret)$C <- C
      }
    }
  } 

## nu classification
if(type(ret) == "nu-svc"){
  indexes <- lapply(sort(unique(y)), function(kk) which(y == kk))
    for (i in 1:(nclass(ret)-1)) {
      jj <- i+1
      for(j in jj:nclass(ret)) {
        p <- p+1
       ##prepare data
        li <- length(indexes[[i]])
        lj <- length(indexes[[j]])

        ##xd <- matrix(0,(li+lj),(li+lj))
        ##xdi <- 1:(li+lj) <= li
        ##xd[xdi,rep(TRUE,li+lj)] <- x[indexes[[i]],c(indexes[[i]],indexes[[j]])]
        ##xd[xdi == FALSE,rep(TRUE,li+lj)] <- x[indexes[[j]],c(indexes[[i]],indexes[[j]])]
        
        if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
          {
            yd <- c(rep(-1,li),rep(1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(j,i)]]
            wl <- c(1,0)
            nweights <- 2
          }
          }
        else
          {
            yd <- c(rep(1,li),rep(-1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(i,j)]]
            wl <- c(0,1)
            nweigths <- 2
          }
          }
       
        boolabel <- yd >= 0
        prior1 <- sum(boolabel)
        md <- length(yd)
        prior0 <- md - prior1
        prior(ret)[[p]] <- list(prior1 = prior1, prior0 = prior0)

        xdd <- matrix(1,li+lj,1)
        
        resv <- .Call("smo_optim",
                      as.double(t(xdd)),
                      as.integer(nrow(xdd)),
                      as.integer(ncol(xdd)),
                      as.double(yd),
                      as.double(x[c(indexes[[i]],indexes[[j]]),c(indexes[[i]],indexes[[j]]),drop=FALSE]),
                      as.integer(if (sparse) x[c(indexes[[i]],indexes[[j]]),c(indexes[[i]],indexes[[j]]),drop=FALSE]@ia else 0),
                      as.integer(if (sparse) x[c(indexes[[i]],indexes[[j]]),c(indexes[[i]],indexes[[j]]),drop=FALSE]@ja else 0),
                      as.integer(sparse),
                      
                      as.double(matrix(rep(-1,m))), #linear term
                      as.integer(ktype),
                      as.integer(1),
                      as.double(C),
                      as.double(nu),
                      as.double(epsilon),
                      as.double(sigma),
                      as.integer(degree),
                      as.double(offset),
                      as.integer(wl), #weightlabl.
                      as.double(weight),
                      as.integer(nweights),
                      as.double(cache),
                      as.double(tol), 
                      as.integer(shrinking),
                      PACKAGE="kernlab")

        reind <- sort(c(indexes[[i]],indexes[[j]]),method="quick",index.return=TRUE)$ix
        tmpres <- resv[c(-(li+lj+1),-(li+lj+2))][reind]
        alpha(ret)[p] <- coef(ret)[p] <- list(tmpres[tmpres != 0])
        ##store SV indexes from current problem for later use in predict
        alphaindex(ret)[p] <- list(c(indexes[[i]],indexes[[j]])[reind][tmpres != 0])
	## store Support Vectors
        ## xmatrix(ret)[p] <- list(xd[tmpres != 0,tmpres != 0,drop=FALSE])
        ##save the indexes from all the SV in a vector (use unique!)
        svindex <- c(svindex,alphaindex(ret)[[p]])
        ## store betas in a vector 
        b(ret) <- c(b(ret), resv[li+lj+1])
	## store objective function values in vector
	obj(ret) <- c(obj(ret), resv[li+lj+2])
        ## used to reconstruct indexes for the patterns matrix x from "indexes"
        problem[p] <- list(c(i,j))
        param(ret)$nu <- nu
        ## margin(ret)[p] <- (min(kernelMult(kernel,xd[1:li,],,alpha(ret)[[p]][1:li])) - max(kernelMult(kernel,xd[li:(li+lj),],,alpha(ret)[[p]][li:(li+lj)])))/2
      }
    }
}  

## Bound constraint C classification
  if(type(ret) == "C-bsvc"){
     if(!is.null(class.weights))
     weightedC <- class.weights[weightlabels] * rep(C,nclass(ret))
    else
      weightedC <- rep(C,nclass(ret)) 
    indexes <- lapply(sort(unique(y)), function(kk) which(y == kk))
    for (i in 1:(nclass(ret)-1)) {
      jj <- i+1
      for(j in jj:nclass(ret)) {
        p <- p+1
        ##prepare data
        li <- length(indexes[[i]])
        lj <- length(indexes[[j]])
        
        if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
          {
            yd <- c(rep(-1,li),rep(1,lj))
            if(!is.null(class.weights)){
              weight <- class.weights[weightlabels[c(j,i)]]
              wl <- c(1,0)
              nweights <- 2
            }
          }
        else
          {
            yd <- c(rep(1,li),rep(-1,lj))
            if(!is.null(class.weights)){
              weight <- class.weights[weightlabels[c(i,j)]]
              wl <- c(0,1)
              nweigths <- 2
            }
          }
       
        boolabel <- yd >= 0
        prior1 <- sum(boolabel)
        md <- length(yd)
        prior0 <- md - prior1
        prior(ret)[[p]] <- list(prior1 = prior1, prior0 = prior0) 

        xdd <- matrix(rnorm(li+lj),li+lj,1)
        
        resv <- .Call("tron_optim",
                      as.double(t(xdd)),
                      as.integer(nrow(xdd)),
                      as.integer(ncol(xdd)),
                      as.double(yd),
                      as.double(x[c(indexes[[i]],indexes[[j]]),c(indexes[[i]],indexes[[j]]),drop=FALSE]),
                      as.integer(if (sparse) x[c(indexes[[i]],indexes[[j]]),c(indexes[[i]],indexes[[j]]),drop=FALSE]@ia else 0),
                      as.integer(if (sparse) x[c(indexes[[i]],indexes[[j]]),c(indexes[[i]],indexes[[j]]),drop=FALSE]@ja else 0),
                      as.integer(sparse),
                      as.integer(2),
                      as.double(0), ##countc
                      as.integer(ktype),
                      as.integer(5), 
                      as.double(C),
                      as.double(epsilon),
                      as.double(sigma),
                      as.integer(degree),
                      as.double(offset),
                      as.double(1),  ##  cost value of alpha seeding
                      as.double(2),  ## step value of alpha seeding
                      as.integer(wl), ##weightlabel
                      as.double(weight),
                      as.integer(nweights),
                      as.double(weightedC),
                      as.double(cache), 
                      as.double(tol),
                      as.integer(10), ##qpsize
                      as.integer(shrinking),
                      PACKAGE="kernlab")

        reind <- sort(c(indexes[[i]],indexes[[j]]),method="quick",index.return=TRUE)$ix
        alpha(ret)[p] <- list(resv[-(li+lj+1)][reind][resv[-(li+lj+1)][reind] > 0])
        ## nonzero alpha*y
        coef(ret)[p] <- list(alpha(ret)[[p]] * yd[reind][resv[-(li+lj+1)][reind] > 0])
        ## store SV indexes from current problem for later use in predict
	alphaindex(ret)[p] <- list(c(indexes[[i]],indexes[[j]])[reind][resv[-(li+lj+1)][reind] > 0])	
        ## store Support Vectors
        ## xmatrix(ret)[p] <- list(xd[resv > 0 ,resv > 0,drop = FALSE])
        ## save the indexes from all the SV in a vector (use unique?)
        svindex <- c(svindex,alphaindex(ret)[[p]])
        ## store betas in a vector 
        b(ret) <- - sapply(coef(ret),sum)
        ## store objective function values vector
        obj(ret) <- c(obj(ret), resv[(li+lj+1)])
        ## used to reconstruct indexes for the patterns matrix x from "indexes" (really usefull ?)
        problem[p] <- list(c(i,j))
        ##store C  in return object
        param(ret)$C <- C
      }
    }
  } 

## SPOC multiclass classification 
if(type(ret) =="spoc-svc")
  {
    if(!is.null(class.weights))
     weightedC <- class.weights[weightlabels] * rep(C,nclass(ret))
    else
      weightedC <- rep(C,nclass(ret)) 
    yd <- sort(y,method="quick", index.return = TRUE)
    x <- matrix(x[yd$ix,yd$ix],nrow=dim(x)[1])
    count <- 0
    
    xdd <- matrix(1,m,1)
    
    resv <- .Call("tron_optim",
                  as.double(t(xdd)),
                  as.integer(nrow(xdd)),
                  as.integer(ncol(xdd)),
                  as.double(rep(yd$x-1,2)),
                  as.double(x),
                  as.integer(if (sparse) x@ia else 0),
                  as.integer(if (sparse) x@ja else 0),
                  as.integer(sparse),
                  as.integer(nclass(ret)),
                  as.integer(count),
                  as.integer(ktype),
                  as.integer(7), 
                  as.double(C),
                  as.double(epsilon),
                  as.double(sigma),
                  as.integer(degree),
                  as.double(offset),
                  as.double(C), 
                  as.double(2), #Cstep
                  as.integer(0), #weightlabel
                  as.double(0),
                  as.integer(0),
                  as.double(weightedC),
                  as.double(cache), 
                  as.double(tol),
                  as.integer(10), #qpsize
                  as.integer(shrinking),
                  PACKAGE="kernlab")
    reind <- sort(yd$ix,method="quick",index.return=TRUE)$ix
    alpha(ret) <- t(matrix(resv[-(nclass(ret)*nrow(xdd)+1)],nclass(ret)))[reind,,drop=FALSE]
    coef(ret) <- lapply(1:nclass(ret), function(x) alpha(ret)[,x][alpha(ret)[,x]!=0])
    names(coef(ret)) <- lev(ret)
    alphaindex(ret) <-  lapply(sort(unique(y)), function(x) which(alpha(ret)[,x]!=0))
    ## xmatrix(ret) <- x
    names(alphaindex(ret)) <- lev(ret)
    svindex <- which(rowSums(alpha(ret)!=0)!=0)
    b(ret) <- 0
    obj(ret) <- resv[(nclass(ret)*nrow(xdd)+1)]
    param(ret)$C <- C
  }

## KBB multiclass classification  
if(type(ret) =="kbb-svc")
  {
     if(!is.null(class.weights))
     weightedC <- class.weights[weightlabels] * rep(C,nclass(ret))
    else
      weightedC <- rep(C,nclass(ret)) 
     yd <- sort(y,method="quick", index.return = TRUE)
     x <- matrix(x[yd$ix,yd$ix],nrow=dim(x)[1])
     count <-  sapply(unique(yd$x), function(c) length(yd$x[yd$x==c]))
     
     xdd <- matrix(1,m,1)

    resv <- .Call("tron_optim",
                  as.double(t(xdd)),
                  as.integer(nrow(xdd)),
                  as.integer(ncol(xdd)),
                  as.double(yd$x-1),
                  as.double(x),
                  as.integer(if (sparse) x@ia else 0),
                  as.integer(if (sparse) x@ja else 0),
                  as.integer(sparse),
                  as.integer(nclass(ret)),
                  as.integer(count),
                  as.integer(ktype),
                  as.integer(8),
                  as.double(C),
                  as.double(epsilon),
                  as.double(sigma),
                  as.integer(degree),
                  as.double(offset),
                  as.double(1), #Cbegin
                  as.double(2), #Cstep
                  as.integer(0), #weightlabl.
                  as.double(0),
                  as.integer(0),
                  as.double(weightedC),
                  as.double(cache),
                  as.double(tol),
                  as.integer(10), #qpsize
                  as.integer(shrinking),
                  PACKAGE="kernlab")
     
     reind <- sort(yd$ix,method="quick",index.return=TRUE)$ix
     alpha(ret) <- matrix(resv[-(nrow(x)*(nclass(ret)-1) + 1)],nrow(x))[reind,,drop=FALSE]
     coef(ret) <-  lapply(1:(nclass(ret)-1), function(x) alpha(ret)[,x][alpha(ret)[,x]!=0])
     alphaindex(ret) <-  lapply(sort(unique(y)), function(x) which((y == x) & (rowSums(alpha(ret))!=0)))
     svindex <- which(rowSums(alpha(ret)!=0)!=0)
     b(ret) <- - sapply(coef(ret),sum)
     obj(ret) <- resv[(nrow(x)*(nclass(ret)-1) + 1)]
     param(ret)$C <- C
   }

  ## Novelty detection
  if(type(ret) =="one-svc")
  {
    xdd <- matrix(1,m,1)
       
    resv <- .Call("smo_optim",
                  as.double(t(xdd)),
                  as.integer(nrow(xdd)),
                  as.integer(ncol(xdd)),
                  as.double(matrix(rep(1,m))),
                  as.double(x),
                  as.integer(if (sparse) x@ia else 0),
                  as.integer(if (sparse) x@ja else 0),
                  as.integer(sparse),
                  as.double(matrix(rep(-1,m))),
                  as.integer(ktype),
                  as.integer(2),
                  as.double(C),
                  as.double(nu),
                  as.double(epsilon),
                  as.double(sigma),
                  as.integer(degree),
                  as.double(offset),
                  as.integer(0), #weightlabl.
                  as.double(0),
                  as.integer(0),
                  as.double(cache),
                  as.double(tol),
                  as.integer(shrinking),
                  PACKAGE="kernlab")

       tmpres <- resv[c(-(m+1),-(m+2))]
       alpha(ret) <- coef(ret) <- tmpres[tmpres != 0]
       svindex <-  alphaindex(ret) <- which(tmpres != 0) 
       ## xmatrix(ret) <- x[svindex,svindex,drop=FALSE]
       b(ret) <- resv[(m+1)]
       obj(ret) <- resv[(m+2)]
       param(ret)$nu <- nu
  }

  ## epsilon regression
  if(type(ret) =="eps-svr")
    {
      xdd <- matrix(1,m,1)
      resv <- .Call("smo_optim",
                    as.double(t(xdd)),
                    as.integer(nrow(xdd)),
                    as.integer(ncol(xdd)),
                    as.double(y),
                    as.double(x),
                    as.integer(if (sparse) x@ia else 0),
                    as.integer(if (sparse) x@ja else 0),
                    as.integer(sparse),
                    as.double(matrix(rep(-1,m))),
                    as.integer(ktype),
                    as.integer(3),
                    as.double(C),
                    as.double(nu),
                    as.double(epsilon),
                    as.double(sigma),
                    as.integer(degree),
                    as.double(offset),
                    as.integer(0), #weightlabl.
                    as.double(0),
                    as.integer(0),
                    as.double(cache), 
                    as.double(tol), 
                    as.integer(shrinking), 
                    PACKAGE="kernlab")

      tmpres <- resv[c(-(m+1),-(m+2))]
      alpha(ret) <- coef(ret) <- tmpres[tmpres != 0]
      svindex <-  alphaindex(ret) <- which(tmpres != 0) 
      ## xmatrix(ret) <- x[svindex,svindex ,drop=FALSE]
      b(ret) <- resv[(m+1)]
      obj(ret) <- resv[(m+2)]
      param(ret)$epsilon <- epsilon
      param(ret)$C <- C
    }

  ## nu regression
  if(type(ret) =="nu-svr")
    {
      xdd <- matrix(1,m,1)
      resv <- .Call("smo_optim",
                    as.double(t(xdd)),
                    as.integer(nrow(xdd)),
                    as.integer(ncol(xdd)),
                    as.double(y),
                    as.double(x),
                    as.integer(if (sparse) x@ia else 0),
                    as.integer(if (sparse) x@ja else 0),
                    as.integer(sparse),
                    as.double(matrix(rep(-1,m))),
                    as.integer(ktype),
                    as.integer(4),
                    as.double(C),
                    as.double(nu),
                    as.double(epsilon),
                    as.double(sigma),
                    as.integer(degree),
                    as.double(offset),
                    as.integer(0),
                    as.double(0),
                    as.integer(0),
                    as.double(cache), 
                    as.double(tol), 
                    as.integer(shrinking), 
                    PACKAGE="kernlab")
      tmpres <- resv[c(-(m+1),-(m+2))]
      alpha(ret) <- coef(ret) <- tmpres[tmpres!=0]
      svindex <-  alphaindex(ret) <- which(tmpres != 0)
      ## xmatrix(ret) <- x[svindex,svindex,drop=FALSE]
      b(ret) <- resv[(m+1)]
      obj(ret) <- resv[(m+2)]
      param(ret)$epsilon <- epsilon
      param(ret)$nu <- nu
    }

  ## bound constraint eps regression
  if(type(ret) =="eps-bsvr")
    {
      xdd <- matrix(1,m,1)
      resv <- .Call("tron_optim",
                    as.double(t(xdd)),
                    as.integer(nrow(xdd)),
                    as.integer(ncol(xdd)),
                    as.double(y),
                    as.double(x),
                    as.integer(if (sparse) x@ia else 0),
                    as.integer(if (sparse) x@ja else 0),
                    as.integer(sparse),
                    as.integer(2),
                    as.integer(0),
                    as.integer(ktype),
                    as.integer(6),
                    as.double(C),
                    as.double(epsilon),
                    as.double(sigma),
                    as.integer(degree),
                    as.double(offset),
                    as.double(1),  #Cbegin
                    as.double(2), #Cstep
                    as.integer(0), #weightlabl.
                    as.double(0),
                    as.integer(0),
                    as.double(0),
                    as.double(cache), 
                    as.double(tol),
                    as.integer(10), #qpsize
                    as.integer(shrinking), 
                   PACKAGE="kernlab")
      tmpres <- resv[-(m+1)]
      alpha(ret) <- coef(ret) <- tmpres[tmpres!=0]
      svindex <-  alphaindex(ret) <- which(tmpres != 0)
      ## xmatrix(ret) <- x[svindex,,drop=FALSE]
      b(ret) <- -sum(alpha(ret))
      obj(ret) <- resv[(m+1)]
      param(ret)$epsilon <- epsilon
      param(ret)$C <- C
  }

  
  kcall(ret) <- match.call()
  kernelf(ret) <- " Kernel matrix used as input."
  ymatrix(ret) <- y
  SVindex(ret) <- unique(sort(svindex,method="quick"))
  nSV(ret)  <- length(unique(svindex))
  if(nSV(ret)==0)
    stop("No Support Vectors found. You may want to change your parameters")
  fitted(ret)  <- if (fit)
    predict(ret, as.kernelMatrix(x[,SVindex(ret),drop = FALSE])) else NULL

  if (fit){
    if(type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="spoc-svc"||type(ret)=="kbb-svc"||type(ret)=="C-bsvc")
      error(ret) <- 1 - .classAgreement(table(y,as.integer(fitted(ret))))
    if(type(ret)=="one-svc")
      error(ret) <- sum(!fitted(ret))/m
    if(type(ret)=="eps-svr"||type(ret)=="nu-svr"||type(ret)=="eps-bsvr")
      error(ret) <- drop(crossprod(fitted(ret) - y)/m)
  }

  cross(ret) <- -1
  if(cross == 1)
    cat("\n","cross should be >1 no cross-validation done!","\n","\n")
  else if (cross > 1)
    {
     
      cerror <- 0
      suppressWarnings(vgr <- split(sample(1:m,m),1:cross))

      for(i in 1:cross)
        {
         
          cind <-  unsplit(vgr[-i],factor(rep((1:cross)[-i],unlist(lapply(vgr[-i],length)))))
          if(type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="spoc-svc"||type(ret)=="kbb-svc"||type(ret)=="C-bsvc")
            {
              if(is.null(class.weights))
                cret <- ksvm(as.kernelMatrix(x[cind,cind]),y[cind],type = type(ret), C=C, nu=nu, tol=tol, cross = 0, fit = FALSE ,cache = cache)
              else
                cret <- ksvm(as.kernelMatrix(x[cind,cind]), as.factor(lev(ret)[y[cind]]),type = type(ret), C=C, nu=nu, tol=tol, cross = 0, fit = FALSE, class.weights = class.weights,cache = cache)
              cres <- predict(cret, as.kernelMatrix(x[vgr[[i]], cind,drop = FALSE][,SVindex(cret),drop=FALSE]))
              cerror <- (1 - .classAgreement(table(y[vgr[[i]]],as.integer(cres))))/cross + cerror
            }
          if(type(ret)=="one-svc")
            {
              cret <- ksvm(as.kernelMatrix(x[cind,cind]),type = type(ret), C=C, nu=nu, tol=tol, cross = 0, fit = FALSE ,cache = cache)
              cres <- predict(cret, as.kernelMatrix(x[vgr[[i]], cind,drop = FALSE][,SVindex(cret),drop=FALSE]))
              cerror <- (1 - sum(cres)/length(cres))/cross + cerror
            }
          if(type(ret)=="eps-svr"||type(ret)=="nu-svr"||type(ret)=="eps-bsvr")
            {
              cret <- ksvm(as.kernelMatrix(x[cind,cind]),y[cind],type=type(ret), C=C,nu=nu,epsilon=epsilon,tol=tol, cross = 0, fit = FALSE, cache = cache, prob.model = FALSE)
              cres <- predict(cret, as.kernelMatrix(x[vgr[[i]], cind,drop = FALSE][,SVindex(cret),drop=FALSE]))
              cerror <- drop(crossprod(cres - y[vgr[[i]]])/m) + cerror
            }
         }
      cross(ret) <- cerror
    }

  prob.model(ret) <- list(NULL)
  
  if(prob.model)
    {
      if(type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="C-bsvc")
        {
          p <- 0
          for (i in 1:(nclass(ret)-1)) {
            jj <- i+1
            for(j in jj:nclass(ret)) {
              p <- p+1
              ##prepare data
              li <- length(indexes[[i]])
              lj <- length(indexes[[j]])
              
              if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
                {
                  yd <- c(rep(-1,li),rep(1,lj))
                  if(!is.null(class.weights)){
                    weight <- weightlabels[c(j,i)]
                    wl <- c(1,0)
                    nweights <- 2
                  }
                }
              else
                {
                  yd <- c(rep(1,li),rep(-1,lj))
                  if(!is.null(class.weights)){
                    weight <- weightlabels[c(i,j)]
                    wl <- c(0,1)
                    nweigths <- 2
                  }
                }
              m <- li+lj
             suppressWarnings(vgr <- split(c(sample(1:li,li),sample((li+1):(li+lj),lj)),1:3)) 

              pres <- yres <- NULL
              for(k in 1:3)
                {
                  cind <- unsplit(vgr[-k],factor(rep((1:3)[-k],unlist(lapply(vgr[-k],length)))))
                  if(is.null(class.weights))
                    cret <- ksvm(as.kernelMatrix(x[c(indexes[[i]],indexes[[j]]),c(indexes[[i]],indexes[[j]]),drop=FALSE][cind,cind]),yd[cind],type = type(ret), C=C, nu=nu, tol=tol, cross = 0, fit = FALSE ,cache = cache, prob.model=FALSE)
                  else
                    cret <- ksvm(as.kernelMatrix(x[c(indexes[[i]],indexes[[j]]),c(indexes[[i]],indexes[[j]]),drop=FALSE][cind,cind]), as.factor(lev(ret)[y[c(indexes[[i]],indexes[[j]])][cind]]),type = type(ret), C=C, nu=nu, tol=tol, cross = 0, fit = FALSE, class.weights = class.weights,cache = cache, prob.model=FALSE)
                  yres <- c(yres,yd[vgr[[k]]])
                  pres <- rbind(pres,predict(cret, as.kernelMatrix(x[c(indexes[[i]],indexes[[j]]),c(indexes[[i]],indexes[[j]]),drop=FALSE][vgr[[k]], cind,drop = FALSE][,SVindex(cret),drop = FALSE]),type="decision"))
                }
              prob.model(ret)[[p]] <- .probPlatt(pres,yres)
            }
          }
        }
      if(type(ret) == "eps-svr"||type(ret) == "nu-svr"||type(ret)=="eps-bsvr"){
        suppressWarnings(vgr<-split(sample(1:m,m),1:3))
        pres <-  NULL
        for(i in 1:3)
          {
            cind <- unsplit(vgr[-i],factor(rep((1:3)[-i],unlist(lapply(vgr[-i],length)))))
            cret <- ksvm(as.kernelMatrix(x[cind,cind]),y[cind],type=type(ret), C=C, nu=nu, epsilon=epsilon, tol=tol, cross = 0, fit = FALSE, cache = cache, prob.model = FALSE)
            cres <- predict(cret, as.kernelMatrix(x[vgr[[i]], cind, drop = FALSE][,SVindex(cret), drop = FALSE]))
            pres <- rbind(pres,predict(cret, as.kernelMatrix(x[vgr[[i]],cind , drop = FALSE][,SVindex(cret) ,drop = FALSE]),type="decision"))
          }
        pres[abs(pres) > (5*sd(pres))] <- 0
        prob.model(ret) <- list(sum(abs(pres))/dim(pres)[1])
      }
    }

  return(ret)
})



.classAgreement <- function (tab) {
  n <- sum(tab)
  if (!is.null(dimnames(tab))) {
    lev <- intersect(colnames(tab), rownames(tab))
    p0 <- sum(diag(tab[lev, lev])) / n
  } else {
    m <- min(dim(tab))
    p0 <- sum(diag(tab[1:m, 1:m])) / n
  }
  return(p0)
}

## List Interface


setMethod("ksvm",signature(x="list"),
function (x,
          y         = NULL,
          type      = NULL,
          kernel    = "stringdot",
          kpar      = list(length = 4, lambda = 0.5),
          C         = 1,
          nu        = 0.2,
          epsilon   = 0.1,
          prob.model = FALSE,
          class.weights = NULL,
          cross     = 0,
          fit       = TRUE,
          cache     = 40,
          tol       = 0.001,
          shrinking = TRUE,
          ...
         ,na.action = na.omit)
{ 
  ret <- new("ksvm")

  if (is.null(y))
    x <- na.action(x)

  n.action(ret) <- na.action
  sparse <- FALSE 
  if (is.null(type)) type(ret) <- if (is.null(y)) "one-svc" else if (is.factor(y)) "C-svc" else "eps-svr"
  
  if(!is.null(type))
  type(ret) <- match.arg(type,c("C-svc",
                                "nu-svc",
                                "kbb-svc",
                                "spoc-svc",
                                "C-bsvc",
                                "one-svc",
                                "eps-svr",
                                "eps-bsvr",
                                "nu-svr"))

  m <- length(x)
  
  if(is.character(kernel)){
    kernel <- match.arg(kernel,c("rbfdot","polydot","tanhdot","vanilladot","laplacedot","besseldot","anovadot","splinedot","stringdot"))

    if(is.character(kpar))
       if(kernel == "tanhdot" || kernel == "vanilladot" || kernel == "polydot"|| kernel == "besseldot" || kernel== "anovadot"|| kernel=="splinedot" || kernel == "rbfdot" || kernel == "laplacedot" )
       {
         stop("List interface supports only the stringdot kernel.")
       }
     }
  
    if(is(kernel,"kernel") & !is(kernel,"stringkernel")) stop("List interface supports only the stringdot kernel.")
    
    if(!is(kernel,"kernel"))
      {
        if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
        kernel <- do.call(kernel, kpar)
      }

    if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

    if (!is(y,"vector") && !is.factor(y) & !is(y,"matrix") & !(type(ret)=="one-svc")) stop("y must be a vector or a factor.")

    if(!(type(ret)=="one-svc"))
      if(is(y,"vector") | is(y,"factor")) ym <- length(y) else if(is(y,"matrix")) ym <-  dim(y)[1] else stop("y must be a matrix or a vector")
    
    if ((type(ret) != "one-svc") && ym != m) stop("x and y don't match.")

    if(nu > 1|| nu <0) stop("nu must be between 0 an 1.")
  
  weightlabels <- NULL
  nweights <- 0
  weight <- 0
  wl <- 0
  ## in case of classification: transform factors into integers
  if (type(ret) == "one-svc") # one class classification --> set dummy
    y <- 1
  else
    if (is.factor(y)) {
      lev(ret) <- levels (y)
      y <- as.integer (y)
      if (!is.null(class.weights)) {
        if (is.null(names (class.weights)))
          stop ("Weights have to be specified along with their according level names !")
        weightlabels <- match (names(class.weights),lev(ret))
        if (any(is.na(weightlabels)))
          stop ("At least one level name is missing or misspelled.")
      }
    }
    else {
      if ((type(ret) =="C-svc" || type(ret) == "nu-svc" ||type(ret) == "C-bsvc" || type(ret) == "spoc-svc" || type(ret) == "kbb-svc") && any(as.integer (y) != y))
        stop ("dependent variable has to be of factor or integer type for classification mode.")
      
      if (type(ret) != "eps-svr" || type(ret) != "nu-svr"|| type(ret)!="eps-bsvr")
        lev(ret) <- sort(unique (y))
    }
 ## initialize
    if (type(ret) =="C-svc" || type(ret) == "nu-svc" ||type(ret) == "C-bsvc" || type(ret) == "spoc-svc" || type(ret) == "kbb-svc")
      nclass(ret) <- length (unique(y))

    p <- 0
    K <- 0 
    svindex <- problem <- NULL
    ktype <- 4
    prior(ret) <- list(NULL)
    sigma <- 0.1
    degree <- offset <- scale <- 1

## C classification
  if(type(ret) == "C-svc"){
    indexes <- lapply(sort(unique(y)), function(kk) which(y == kk))
    for (i in 1:(nclass(ret)-1)) {
      jj <- i+1
      for(j in jj:nclass(ret)) {
        p <- p+1
        ##prepare data
        li <- length(indexes[[i]])
        lj <- length(indexes[[j]])
    
        if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
          {
            yd <- c(rep(-1,li),rep(1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(j,i)]]
            wl <- c(1,0)
            nweights <- 2
          }
          }
        else
          {
            yd <- c(rep(1,li),rep(-1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(i,j)]]
            wl <- c(0,1)
            nweigths <- 2
          }
          }
       
        boolabel <- yd >= 0
        prior1 <- sum(boolabel)

        md <- length(yd)
        prior0 <- md - prior1
        prior(ret)[[p]] <- list(prior1 = prior1, prior0 = prior0) 
        
        K <- kernelMatrix(kernel,x[c(indexes[[i]],indexes[[j]])])
        xdd <- matrix(1,li+lj,1) 
        resv <- .Call("smo_optim",
                      as.double(t(xdd)),
                      as.integer(nrow(xdd)),
                      as.integer(ncol(xdd)),
                      as.double(yd),
                      as.double(K),
                      
                      as.integer(if (sparse) x@ia else 0),
                      as.integer(if (sparse) x@ja else 0),
                      as.integer(sparse),
                      
                      as.double(matrix(rep(-1,m))), ##linear term
                      as.integer(ktype),
                      as.integer(0), 
                      as.double(C),
                      as.double(nu),
                      as.double(epsilon),
                      as.double(sigma),
                      as.integer(degree),
                      as.double(offset),
                      as.integer(wl), ##weightlabel
                      as.double(weight),
                      as.integer(nweights),
                      as.double(cache), 
                      as.double(tol),
                      as.integer(shrinking),
                      PACKAGE="kernlab")

        reind <- sort(c(indexes[[i]],indexes[[j]]),method="quick",index.return=TRUE)$ix
        tmpres <- resv[c(-(li+lj+1),-(li+lj+2))][reind]
        ## alpha
        alpha(ret)[p] <- list(tmpres[tmpres > 0])
        ## coefficients alpha*y
        coef(ret)[p] <- list(alpha(ret)[[p]]*yd[reind][tmpres > 0])
        ## store SV indexes from current problem for later use in predict
        alphaindex(ret)[p] <- list(c(indexes[[i]],indexes[[j]])[reind][tmpres>0])
        ## store Support Vectors
        xmatrix(ret)[p] <- list(x[c(indexes[[i]],indexes[[j]])][reind][tmpres > 0])
        ## save the indexes from all the SV in a vector (use unique?)
        svindex <- c(svindex,alphaindex(ret)[[p]])
        ## store betas in a vector 
        b(ret) <- c(b(ret), resv[li+lj+1])
        obj(ret) <- c(obj(ret),resv[li+lj+2])
        ## used to reconstruct indexes for the patterns matrix x from "indexes" (really usefull ?)
        problem[p] <- list(c(i,j))
        ##store C  in return object
        param(ret)$C <- C
        ## margin(ret)[p] <- (min(kernelMult(kernel,xd[1:li,],,alpha(ret)[[p]][1:li])) - max(kernelMult(kernel,xd[li:(li+lj),],,alpha(ret)[[p]][li:(li+lj)])))/2
      }
    }
  } 

## nu classification
if(type(ret) == "nu-svc"){
  indexes <- lapply(sort(unique(y)), function(kk) which(y == kk))
    for (i in 1:(nclass(ret)-1)) {
      jj <- i+1
      for(j in jj:nclass(ret)) {
        p <- p+1
       ##prepare data
        li <- length(indexes[[i]])
        lj <- length(indexes[[j]])

        if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
          {
            yd <- c(rep(-1,li),rep(1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(j,i)]]
            wl <- c(1,0)
            nweights <- 2
          }
          }
        else
          {
            yd <- c(rep(1,li),rep(-1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(i,j)]]
            wl <- c(0,1)
            nweigths <- 2
          }
          }
       
        boolabel <- yd >= 0
        prior1 <- sum(boolabel)
        md <- length(yd)
        prior0 <- md - prior1
        prior(ret)[[p]] <- list(prior1 = prior1, prior0 = prior0)

        K <- kernelMatrix(kernel,x[c(indexes[[i]],indexes[[j]])])
        xdd <- matrix(1,li+lj,1)
        resv <- .Call("smo_optim",
                      as.double(t(xdd)),
                      as.integer(nrow(xdd)),
                      as.integer(ncol(xdd)),
                      as.double(yd),
                      as.double(K),
                      as.integer(if (sparse) x@ia else 0),
                      as.integer(if (sparse) x@ja else 0),
                      as.integer(sparse),
                      
                      as.double(matrix(rep(-1,m))), #linear term
                      as.integer(ktype),
                      as.integer(1),
                      as.double(C),
                      as.double(nu),
                      as.double(epsilon),
                      as.double(sigma),
                      as.integer(degree),
                      as.double(offset),
                      as.integer(wl), #weightlabl.
                      as.double(weight),
                      as.integer(nweights),
                      as.double(cache),
                      as.double(tol), 
                      as.integer(shrinking),
                      PACKAGE="kernlab")
        reind <- sort(c(indexes[[i]],indexes[[j]]),method="quick",index.return=TRUE)$ix
        tmpres <- resv[c(-(li+lj+1),-(li+lj+2))][reind]
        alpha(ret)[p] <- coef(ret)[p] <- list(tmpres[tmpres != 0])
        ##store SV indexes from current problem for later use in predict
        alphaindex(ret)[p] <- list(c(indexes[[i]],indexes[[j]])[reind][tmpres!=0])
        ## store Support Vectors
        xmatrix(ret)[p] <- list(x[c(indexes[[i]],indexes[[j]])][reind][tmpres != 0])
        ##save the indexes from all the SV in a vector (use unique!)
        svindex <- c(svindex,alphaindex(ret)[[p]])
        ## store betas in a vector 
        b(ret) <- c(b(ret), resv[li+lj+1])
        obj(ret) <- c(obj(ret), resv[li+lj+2])
        ## used to reconstruct indexes for the patterns matrix x from "indexes"
        problem[p] <- list(c(i,j))
        param(ret)$nu <- nu
        ## margin(ret)[p] <- (min(kernelMult(kernel,xd[1:li,],,alpha(ret)[[p]][1:li])) - max(kernelMult(kernel,xd[li:(li+lj),],,alpha(ret)[[p]][li:(li+lj)])))/2
      }
    }
}  

## Bound constraint C classification
  if(type(ret) == "C-bsvc"){
     if(!is.null(class.weights))
     weightedC <- class.weights[weightlabels] * rep(C,nclass(ret))
    else
      weightedC <- rep(C,nclass(ret)) 
    indexes <- lapply(sort(unique(y)), function(kk) which(y == kk))
    for (i in 1:(nclass(ret)-1)) {
      jj <- i+1
      for(j in jj:nclass(ret)) {
        p <- p+1
        ##prepare data
        li <- length(indexes[[i]])
        lj <- length(indexes[[j]])

        if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
          {
            yd <- c(rep(-1,li),rep(1,lj))
            if(!is.null(class.weights)){
              weight <- class.weights[weightlabels[c(j,i)]]
              wl <- c(1,0)
              nweights <- 2
            }
          }
        else
          {
            yd <- c(rep(1,li),rep(-1,lj))
            if(!is.null(class.weights)){
              weight <- class.weights[weightlabels[c(i,j)]]
              wl <- c(0,1)
              nweigths <- 2
            }
          }
       
        boolabel <- yd >= 0
        prior1 <- sum(boolabel)
        md <- length(yd)
        prior0 <- md - prior1
        prior(ret)[[p]] <- list(prior1 = prior1, prior0 = prior0) 

        K <- kernelMatrix(kernel,x[c(indexes[[i]],indexes[[j]])])
        xdd <- matrix(1,li+lj,1) 

        resv <- .Call("tron_optim",
                      as.double(t(xdd)),
                      as.integer(nrow(xdd)),
                      as.integer(ncol(xdd)),
                      as.double(yd),
                      as.double(K),
                      as.integer(if (sparse) x@ia else 0),
                      as.integer(if (sparse) x@ja else 0),
                      as.integer(sparse),
                      as.integer(2),
                      as.double(0), ##countc
                      as.integer(ktype),
                      as.integer(5), 
                      as.double(C),
                      as.double(epsilon),
                      as.double(sigma),
                      as.integer(degree),
                      as.double(offset),
                      as.double(1),  ##  cost value of alpha seeding
                      as.double(2),  ## step value of alpha seeding
                      as.integer(wl), ##weightlabel
                      as.double(weight),
                      as.integer(nweights),
                      as.double(weightedC),
                      as.double(cache), 
                      as.double(tol),
                      as.integer(10), ##qpsize
                      as.integer(shrinking),
                      PACKAGE="kernlab")
                
        reind <- sort(c(indexes[[i]],indexes[[j]]),method="quick",index.return=TRUE)$ix
        alpha(ret)[p] <- list(resv[-(li+lj+1)][reind][resv[-(li+lj+1)][reind] > 0])
        ## nonzero alpha*y
        coef(ret)[p] <- list(alpha(ret)[[p]] * yd[reind][resv[-(li+lj+1)][reind] > 0])
        ## store SV indexes from current problem for later use in predict
        alphaindex(ret)[p] <- list(c(indexes[[i]],indexes[[j]])[reind][resv[-(li+lj+1)][reind] > 0])
        ## store Support Vectors
        xmatrix(ret)[p] <- list(x[c(indexes[[i]],indexes[[j]])][reind][resv[-(li+lj+1)][reind] > 0])
        ## save the indexes from all the SV in a vector (use unique?)
        svindex <- c(svindex,alphaindex(ret)[[p]])
        ## store betas in a vector 
        b(ret) <- - sapply(coef(ret),sum)
	obj(ret) <- c(obj(ret),resv[(li+lj+1)])
        ## used to reconstruct indexes for the patterns matrix x from "indexes" (really usefull ?)
        problem[p] <- list(c(i,j))
        ##store C  in return object
        param(ret)$C <- C
##        margin(ret)[p] <- (min(kernelMult(kernel,xd[1:li,],,alpha(ret)[[p]][1:li])) - max(kernelMult(kernel,xd[li:(li+lj),],,alpha(ret)[[p]][li:(li+lj)])))/2
      }
    }
  } 

## SPOC multiclass classification 
if(type(ret) =="spoc-svc")
  {
    if(!is.null(class.weights))
     weightedC <- class.weights[weightlabels] * rep(C,nclass(ret))
    else
      weightedC <- rep(C,nclass(ret)) 
    yd <- sort(y,method="quick", index.return = TRUE)
    x <- x[yd$ix]
    count <- 0

    K <- kernelMatrix(kernel,x)
    xdd <- matrix(1,length(x),1) 
    resv <- .Call("tron_optim",
                  as.double(t(xdd)),
                  as.integer(nrow(xdd)),
                  as.integer(ncol(xdd)),
                  as.double(rep(yd$x-1,2)),
                  as.double(K),
                  as.integer(if (sparse) x@ia else 0),
                  as.integer(if (sparse) x@ja else 0),
                  as.integer(sparse),
                  as.integer(nclass(ret)),
                  as.integer(count),
                  as.integer(ktype),
                  as.integer(7), 
                  as.double(C),
                  as.double(epsilon),
                  as.double(sigma),
                  as.integer(degree),
                  as.double(offset),
                  as.double(C), 
                  as.double(2), #Cstep
                  as.integer(0), #weightlabel
                  as.double(0),
                  as.integer(0),
                  as.double(weightedC),
                  as.double(cache), 
                  as.double(tol),
                  as.integer(10), #qpsize
                  as.integer(shrinking),
                  PACKAGE="kernlab")

    reind <- sort(yd$ix,method="quick",index.return=TRUE)$ix
    alpha(ret) <- t(matrix(resv[-(nclass(ret)*nrow(xdd) + 1)],nclass(ret)))[reind,,drop=FALSE]
    coef(ret) <- lapply(1:nclass(ret), function(x) alpha(ret)[,x][alpha(ret)[,x]!=0])
    names(coef(ret)) <- lev(ret)
    alphaindex(ret) <-  lapply(1:nclass(ret), function(x) which(alpha(ret)[,x]!=0))
    names(alphaindex(ret)) <- lev(ret)
    xmatrix(ret) <- x
    svindex <- which(rowSums(alpha(ret)!=0)!=0)
    b(ret) <- 0
    obj(ret) <- resv[(nclass(ret)*nrow(xdd) + 1)]
    param(ret)$C <- C
  }

## KBB multiclass classification  
if(type(ret) =="kbb-svc")
  {
    if(!is.null(class.weights))
      weightedC <- weightlabels * rep(C,nclass(ret))
    else
      weightedC <- rep(C,nclass(ret)) 
    yd <- sort(y,method="quick", index.return = TRUE)
    x <- x[yd$ix]
    count <-  sapply(unique(yd$x), function(c) length(yd$x[yd$x==c]))

    K <- kernelMatrix(kernel,x)
    xdd <- matrix(1,length(x),1)

    resv <- .Call("tron_optim",
                  as.double(t(xdd)),
                  as.integer(nrow(xdd)),
                  as.integer(ncol(xdd)),
                  as.double(yd$x-1),
                  as.double(K),
                  as.integer(if (sparse) x@ia else 0),
                  as.integer(if (sparse) x@ja else 0),
                  as.integer(sparse),
                  as.integer(nclass(ret)),
                  as.integer(count),
                  as.integer(ktype),
                  as.integer(8),
                  as.double(C),
                  as.double(epsilon),
                  as.double(sigma),
                  as.integer(degree),
                  as.double(offset),
                  as.double(1), #Cbegin
                  as.double(2), #Cstep
                  as.integer(0), #weightlabl.
                  as.double(0),
                  as.integer(0),
                  as.double(weightedC),
                  as.double(cache),
                  as.double(tol),
                  as.integer(10), #qpsize
                  as.integer(shrinking),
                  PACKAGE="kernlab")
    reind <- sort(yd$ix,method="quick",index.return=TRUE)$ix
    alpha(ret) <- matrix(resv[-((nclass(ret)-1)*length(x)+1)],length(x))[reind,,drop=FALSE]
    xmatrix(ret) <- x<- x[reind]
    coef(ret) <-  lapply(1:(nclass(ret)-1), function(x) alpha(ret)[,x][alpha(ret)[,x]!=0])
    alphaindex(ret) <-  lapply(sort(unique(y)), function(x) which((y == x) & (rowSums(alpha(ret))!=0)))
    svindex <- which(rowSums(alpha(ret)!=0)!=0)
    b(ret) <- - sapply(coef(ret),sum)
    obj(ret) <- resv[((nclass(ret)-1)*length(x)+1)]
    param(ret)$C <- C
  }

  ## Novelty detection
  if(type(ret) =="one-svc")
  {
    K <- kernelMatrix(kernel,x)
    xdd <- matrix(1,length(x),1) 
    resv <- .Call("smo_optim",
                  as.double(t(xdd)),
                  as.integer(nrow(xdd)),
                  as.integer(ncol(xdd)),
                  as.double(matrix(rep(1,m))),
                  as.double(K),
                  as.integer(if (sparse) x@ia else 0),
                  as.integer(if (sparse) x@ja else 0),
                  as.integer(sparse),
                  as.double(matrix(rep(-1,m))),
                  as.integer(ktype),
                  as.integer(2),
                  as.double(C),
                  as.double(nu),
                  as.double(epsilon),
                  as.double(sigma),
                  as.integer(degree),
                  as.double(offset),
                  as.integer(0), #weightlabl.
                  as.double(0),
                  as.integer(0),
                  as.double(cache),
                  as.double(tol),
                  as.integer(shrinking),
                  PACKAGE="kernlab")

       tmpres <- resv[c(-(m+1),-(m+2))]
       alpha(ret) <- coef(ret) <- tmpres[tmpres != 0]
       svindex <- alphaindex(ret) <- which(tmpres !=0)
       xmatrix(ret) <- x[svindex]
       b(ret) <- resv[(m+1)]
       obj(ret) <- resv[(m+2)]
       param(ret)$nu <- nu
  }

  ## epsilon regression
  if(type(ret) =="eps-svr")
    {
      K <- kernelMatrix(kernel,x)
      xdd <- matrix(1,length(x),1)  
      resv <- .Call("smo_optim",
                    as.double(t(xdd)),
                    as.integer(nrow(xdd)),
                    as.integer(ncol(xdd)),
                    as.double(y),
                    as.double(K),
                    as.integer(if (sparse) x@ia else 0),
                    as.integer(if (sparse) x@ja else 0),
                    as.integer(sparse),
                    as.double(matrix(rep(-1,m))),
                    as.integer(ktype),
                    as.integer(3),
                    as.double(C),
                    as.double(nu),
                    as.double(epsilon),
                    as.double(sigma),
                    as.integer(degree),
                    as.double(offset),
                    as.integer(0), #weightlabl.
                    as.double(0),
                    as.integer(0),
                    as.double(cache), 
                    as.double(tol), 
                    as.integer(shrinking), 
                    PACKAGE="kernlab")
      tmpres <- resv[c(-(m+1),-(m+2))]
      alpha(ret) <- coef(ret) <- tmpres[tmpres != 0]
      svindex <-  alphaindex(ret) <- which(tmpres != 0)
      xmatrix(ret) <- x[svindex]
      b(ret) <- resv[(m+1)]
      obj(ret) <- resv[(m+2)]
      param(ret)$epsilon <- epsilon
      param(ret)$C <- C
    }

  ## nu regression
  if(type(ret) =="nu-svr")
    {
      K <- kernelMatrix(kernel,x)
      xdd <- matrix(1,length(x),1)
      resv <- .Call("smo_optim",
                    as.double(t(xdd)),
                    as.integer(nrow(xdd)),
                    as.integer(ncol(xdd)),
                    as.double(y),
                    as.double(K),
                    as.integer(if (sparse) x@ia else 0),
                    as.integer(if (sparse) x@ja else 0),
                    as.integer(sparse),
                    as.double(matrix(rep(-1,m))),
                    as.integer(ktype),
                    as.integer(4),
                    as.double(C),
                    as.double(nu),
                    as.double(epsilon),
                    as.double(sigma),
                    as.integer(degree),
                    as.double(offset),
                    as.integer(0),
                    as.double(0),
                    as.integer(0),
                    as.double(cache), 
                    as.double(tol), 
                    as.integer(shrinking), 
                    PACKAGE="kernlab")
      tmpres <- resv[c(-(m+1),-(m+2))]
      alpha(ret) <- coef(ret) <- tmpres[tmpres!=0]
      svindex <-  alphaindex(ret) <- which(tmpres != 0)
      xmatrix(ret) <- x[svindex]
      b(ret) <- resv[(m+1)]
      obj(ret) <- resv[(m+2)]
      param(ret)$epsilon <- epsilon
      param(ret)$nu <- nu
    }

  ## bound constraint eps regression
  if(type(ret) =="eps-bsvr")
    {
      K <- kernelMatrix(kernel,x)
      xdd <- matrix(1,length(x),1)	 
      resv <- .Call("tron_optim",
                    as.double(t(xdd)),
                    as.integer(nrow(xdd)),
                    as.integer(ncol(xdd)),
                    as.double(y),
                    as.double(K),
                    as.integer(if (sparse) x@ia else 0),
                    as.integer(if (sparse) x@ja else 0),
                    as.integer(sparse),
                    as.integer(2),
                    as.integer(0),
                    as.integer(ktype),
                    as.integer(6),
                    as.double(C),
                    as.double(epsilon),
                    as.double(sigma),
                    as.integer(degree),
                    as.double(offset),
                    as.double(1),  #Cbegin
                    as.double(2), #Cstep
                    as.integer(0), #weightlabl.
                    as.double(0),
                    as.integer(0),
                    as.double(0),
                    as.double(cache), 
                    as.double(tol),
                    as.integer(10), #qpsize
                    as.integer(shrinking), 
                   PACKAGE="kernlab")
      tmpres <- resv[-(m+1)]
      alpha(ret) <- coef(ret) <- tmpres[tmpres!=0]
      svindex <-  alphaindex(ret) <- which(tmpres != 0)
      xmatrix(ret) <- x[svindex]
      b(ret) <- -sum(alpha(ret))
      obj(ret) <- resv[(m+1)]
      param(ret)$epsilon <- epsilon
      param(ret)$C <- C
  }

  kcall(ret) <- match.call()
  kernelf(ret) <- kernel
  ymatrix(ret) <- y
  SVindex(ret) <- unique(svindex)
  nSV(ret)  <- length(unique(svindex))

  if(type(ret)=="eps-svr"||type(ret)=="nu-svr"||type(ret)=="eps-bsvr")
    nclass(ret) <- m

  if(type(ret)=="one-svc")
    nclass(ret) <- 1
  
  if(nSV(ret)==0)
    stop("No Support Vectors found. You may want to change your parameters")
  fitted(ret) <- if (fit) {
    if((type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="C-bsvc") & nclass(ret) > 2)
      predict(ret, x)
    else
      if((type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="C-bsvc"||type(ret)=="spoc-bsvc"||type(ret)=="kbb-bsvc"))
        predict(ret,as.kernelMatrix(K[reind,reind][,SVindex(ret), drop=FALSE]))
      else
        predict(ret,as.kernelMatrix(K[,SVindex(ret), drop=FALSE]))
  }
  else NULL

  if (fit){
    if(type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="spoc-svc"||type(ret)=="kbb-svc"||type(ret)=="C-bsvc")
      error(ret) <- 1 - .classAgreement(table(y,as.integer(fitted(ret))))
    if(type(ret)=="one-svc")
      error(ret) <- sum(!fitted(ret))/m
    if(type(ret)=="eps-svr"||type(ret)=="nu-svr"||type(ret)=="eps-bsvr")
      error(ret) <- drop(crossprod(fitted(ret) - y)/m)
  }

  cross(ret) <- -1
  if(!((type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="C-bsvc") & nclass(ret) > 2))
    {
      if((type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="C-bsvc"||type(ret)=="spoc-bsvc"||type(ret)=="kbb-bsvc"))
        K <- as.kernelMatrix(K[reind,reind])
      if(cross == 1)
        cat("\n","cross should be >1 no cross-validation done!","\n","\n")
      else if (cross > 1)
        {
          cerror <- 0
          suppressWarnings(vgr <- split(sample(1:dim(K)[1],dim(K)[1]),1:cross))
          for(i in 1:cross)
            {
              cind <-  unsplit(vgr[-i],factor(rep((1:cross)[-i],unlist(lapply(vgr[-i],length)))))
              if(type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="spoc-svc"||type(ret)=="kbb-svc"||type(ret)=="C-bsvc")
                {
                  if(is.null(class.weights))
                    cret <- ksvm(as.kernelMatrix(K[cind,cind]),y[cind],type = type(ret), C=C, nu=nu, tol=tol, cross = 0, fit = FALSE ,cache = cache)
                  else
                    cret <- ksvm(as.kernelMatrix(K[cind,cind]),as.factor(lev(ret)[y[cind]]),type = type(ret), C=C, nu=nu, tol=tol, cross = 0, fit = FALSE, class.weights = class.weights,cache = cache)
                  cres <- predict(cret, as.kernelMatrix(K[vgr[[i]], cind,drop = FALSE][,SVindex(cret),drop=FALSE]))
                  cerror <- (1 - .classAgreement(table(y[vgr[[i]]],as.integer(cres))))/cross + cerror
                }
              if(type(ret)=="one-svc")
                {
                  cret <- ksvm(as.kernelMatrix(K[cind,cind]), type = type(ret), C=C, nu=nu, tol=tol, cross = 0, fit = FALSE ,cache = cache)
                  cres <- predict(cret, as.kernelMatrix(K[vgr[[i]], cind,drop = FALSE][,SVindex(cret),drop=FALSE]))
                  cerror <- (1 - sum(cres)/length(cres))/cross + cerror
            }

              if(type(ret)=="eps-svr"||type(ret)=="nu-svr"||type(ret)=="eps-bsvr")
                {
                  cret <- ksvm(as.kernelMatrix(K[cind,cind]),y[cind],type=type(ret), C=C,nu=nu,epsilon=epsilon,tol=tol, cross = 0, fit = FALSE, cache = cache, prob.model = FALSE)
                  cres <- predict(cret, as.kernelMatrix(K[vgr[[i]], cind,drop = FALSE][,SVindex(cret),drop=FALSE]))
                  cerror <- drop(crossprod(cres - y[vgr[[i]]])/m) + cerror
                }
            }
          cross(ret) <- cerror
        }
      prob.model(ret) <- list(NULL)
      if(prob.model)
        {
          if(type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="C-bsvc")
            {
              p <- 0
              for (i in 1:(nclass(ret)-1)) {
                jj <- i+1
                for(j in jj:nclass(ret)) {
                  p <- p+1
                  ##prepare data
                  li <- length(indexes[[i]])
                  lj <- length(indexes[[j]])
                  
                  if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
                    {
                      yd <- c(rep(-1,li),rep(1,lj))
                      if(!is.null(class.weights)){
                        weight <- weightlabels[c(j,i)]
                        wl <- c(1,0)
                        nweights <- 2
                      }
                    }
                  else
                    {
                      yd <- c(rep(1,li),rep(-1,lj))
                      if(!is.null(class.weights)){
                        weight <- weightlabels[c(i,j)]
                        wl <- c(0,1)
                        nweigths <- 2
                      }
                    }
                  m <- li+lj
                 suppressWarnings(vgr <- split(c(sample(1:li,li),sample((li+1):(li+lj),lj)),1:3)) 
                  
                  pres <- yres <- NULL
                  for(k in 1:3)
                    {
                      cind <- unsplit(vgr[-k],factor(rep((1:3)[-k],unlist(lapply(vgr[-k],length)))))
                      cret <- ksvm(as.kernelMatrix(as.kernelMatrix(K[c(indexes[[i]],indexes[[j]]),c(indexes[[i]],indexes[[j]]),drop=FALSE][cind,cind])), yd[cind], type = type(ret),  C=C, nu=nu, tol=tol, cross = 0, fit = FALSE, cache = cache, prob.model=FALSE)
                      yres <- c(yres,yd[vgr[[k]]])
                      pres <- rbind(pres,predict(cret, as.kernelMatrix(K[c(indexes[[i]],indexes[[j]]),c(indexes[[i]],indexes[[j]]),drop=FALSE][vgr[[k]], cind,drop = FALSE][,SVindex(cret),drop = FALSE]),type="decision"))
                      
                    }
                  prob.model(ret)[[p]] <- .probPlatt(pres,yres)
                }
              }
            }
          if(type(ret) == "eps-svr"||type(ret) == "nu-svr"||type(ret)=="eps-bsvr"){
            suppressWarnings(vgr<-split(sample(1:m,m),1:3))
            pres <-  NULL
            for(i in 1:3)
              {
                cind <- unsplit(vgr[-i],factor(rep((1:3)[-i],unlist(lapply(vgr[-i],length)))))
                cret <- ksvm(as.kernelMatrix(K[cind,cind]),y[cind],type=type(ret), C=C, nu=nu, epsilon=epsilon, tol=tol, cross = 0, fit = FALSE, cache = cache, prob.model = FALSE)

               cres <- predict(cret, as.kernelMatrix(K[vgr[[i]], cind, drop = FALSE][,SVindex(cret), drop = FALSE]))
                pres <- rbind(pres,predict(cret, as.kernelMatrix(K[vgr[[i]],cind , drop = FALSE][,SVindex(cret) ,drop = FALSE]),type="decision"))
              }
            pres[abs(pres) > (5*sd(pres))] <- 0
            prob.model(ret) <- list(sum(abs(pres))/dim(pres)[1])
          }
        }
    }
  else{
    if(cross == 1)
      cat("\n","cross should be >1 no cross-validation done!","\n","\n")
    else if (cross > 1)
      {
        cerror <- 0
        suppressWarnings(vgr<-split(sample(1:m,m),1:cross))
        for(i in 1:cross)
          {
            cind <-  unsplit(vgr[-i],factor(rep((1:cross)[-i],unlist(lapply(vgr[-i],length)))))
            if(type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="spoc-svc"||type(ret)=="kbb-svc"||type(ret)=="C-bsvc")
              {
                if(is.null(class.weights))
                  cret <- ksvm(x[cind],y[cind],type = type(ret),kernel=kernel,kpar = NULL, C=C, nu=nu, tol=tol, cross = 0, fit = FALSE ,cache = cache)
                else
                  cret <- ksvm(x[cind],as.factor(lev(ret)[y[cind]]),type = type(ret),kernel=kernel,kpar = NULL, C=C, nu=nu, tol=tol, cross = 0, fit = FALSE, class.weights = class.weights,cache = cache)
                cres <- predict(cret, x[vgr[[i]]])
                cerror <- (1 - .classAgreement(table(y[vgr[[i]]],as.integer(cres))))/cross + cerror
              }
            if(type(ret)=="eps-svr"||type(ret)=="nu-svr"||type(ret)=="eps-bsvr")
              {
                cret <- ksvm(x[cind],y[cind],type=type(ret),kernel=kernel,kpar = NULL,C=C,nu=nu,epsilon=epsilon,tol=tol, cross = 0, fit = FALSE, cache = cache, prob.model = FALSE)
                cres <- predict(cret, x[vgr[[i]]])
                cerror <- drop(crossprod(cres - y[vgr[[i]]])/m)/cross + cerror
              }
          }
        cross(ret) <- cerror
      }
    
    prob.model(ret) <- list(NULL)
    
    if(prob.model)
      {
        if(type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="C-bsvc")
          {
            p <- 0
            for (i in 1:(nclass(ret)-1)) {
              jj <- i+1
              for(j in jj:nclass(ret)) {
                p <- p+1
                ##prepare data
                li <- length(indexes[[i]])
                lj <- length(indexes[[j]])
                                
                if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
                  {
                    yd <- c(rep(-1,li),rep(1,lj))
                    if(!is.null(class.weights)){
                      weight <- weightlabels[c(j,i)]
                      wl <- c(1,0)
                      nweights <- 2
                    }
                  }
                else
                  {
                    yd <- c(rep(1,li),rep(-1,lj))
                    if(!is.null(class.weights)){
                      weight <- weightlabels[c(i,j)]
                      wl <- c(0,1)
                      nweigths <- 2
                    }
                  }
                m <- li+lj
                suppressWarnings(vgr <- split(c(sample(1:li,li),sample((li+1):(li+lj),lj)),1:3)) 
                
                pres <- yres <- NULL
                for(k in 1:3)
                  { 
                    cind <- unsplit(vgr[-k],factor(rep((1:3)[-k],unlist(lapply(vgr[-k],length)))))


                if(is.null(class.weights))
                  cret <- ksvm(x[c(indexes[[i]], indexes[[j]])][cind],yd[cind],type = type(ret),kernel=kernel,kpar = NULL, C=C, nu=nu, tol=tol, cross = 0, fit = FALSE ,cache = cache, prob.model=FALSE)
                else
                  cret <- ksvm(x[c(indexes[[i]], indexes[[j]])][cind],as.factor(lev(ret)[y[cind]]),type = type(ret),kernel=kernel,kpar = NULL, C=C, nu=nu, tol=tol, cross = 0, fit = FALSE, class.weights = class.weights,cache = cache, prob.model=FALSE)
                    yres <- c(yres,yd[vgr[[k]]])
                    pres <- rbind(pres,predict(cret, x[c(indexes[[i]], indexes[[j]])][vgr[[k]]],type="decision"))
                  }
                prob.model(ret)[[p]] <- .probPlatt(pres,yres)
              }
            }
          }
        if(type(ret) == "eps-svr"||type(ret) == "nu-svr"||type(ret)=="eps-bsvr"){
          suppressWarnings(vgr<-split(sample(1:m,m),1:3))
          for(i in 1:3)
            {
              cind <- unsplit(vgr[-i],factor(rep((1:3)[-i],unlist(lapply(vgr[-i],length)))))
              
              cret <- ksvm(x[cind],y[cind],type=type(ret),kernel=kernel,kpar = NULL,C=C,nu=nu,epsilon=epsilon,tol=tol, cross = 0, fit = FALSE, cache = cache, prob.model = FALSE)
              cres <- predict(cret, x[vgr[[i]]])
              pres <- rbind(pres,predict(cret, x[vgr[[i]]],type="decision"))
            }
          pres[abs(pres) > (5*sd(pres))] <- 0
          prob.model(ret) <- list(sum(abs(pres))/dim(pres)[1])
        }
      }
  }

  return(ret)
})

##**************************************************************#
## predict for matrix, data.frame input

setMethod("predict", signature(object = "ksvm"),
function (object, newdata, type = "response", coupler = "minpair")
{
  type <- match.arg(type,c("response","probabilities","votes","decision"))
  if (missing(newdata) && type=="response" & !is.null(fitted(object)))
    return(fitted(object))
  else if(missing(newdata))
    stop("Missing data !")

 if(!is(newdata,"list")){ 
  if (!is.null(terms(object)) & !is(newdata,"kernelMatrix"))
    {
      if(!is.matrix(newdata))
        newdata <- model.matrix(delete.response(terms(object)), as.data.frame(newdata), na.action = n.action(object))
    }
  else
    newdata  <- if (is.vector(newdata)) t(t(newdata)) else as.matrix(newdata)

  
    newnrows <- nrow(newdata)
    newncols <- ncol(newdata)
    if(!is(newdata,"kernelMatrix") && !is.null(xmatrix(object))){
      if(is(xmatrix(object),"list") && is(xmatrix(object)[[1]],"matrix")) oldco <- ncol(xmatrix(object)[[1]])
      if(is(xmatrix(object),"matrix")) oldco <- ncol(xmatrix(object))
      if (oldco != newncols) stop ("test vector does not match model !")
    }
  }
  else
   newnrows <- length(newdata)
  
  p <- 0
  
  if (is.list(scaling(object)))
    newdata[,scaling(object)$scaled] <-
      scale(newdata[,scaling(object)$scaled, drop = FALSE],
            center = scaling(object)$x.scale$"scaled:center", scale  = scaling(object)$x.scale$"scaled:scale")

  if(type == "response" || type =="decision" || type=="votes")
    {
      if(type(object)=="C-svc"||type(object)=="nu-svc"||type(object)=="C-bsvc")
        {
          predres <- 1:newnrows
          if(type=="decision")
	   votematrix <- matrix(0,nclass(object)*(nclass(object)-1)/2,newnrows)
	  else
	   votematrix <- matrix(0,nclass(object),newnrows)
          
	  for(i in 1:(nclass(object)-1))
            {
              jj <- i+1
              for(j in jj:nclass(object))
                {
                  p <- p+1

                  if(is(newdata,"kernelMatrix"))
                    ret <- newdata[,which(SVindex(object)%in%alphaindex(object)[[p]]), drop=FALSE] %*% coef(object)[[p]] - b(object)[p]
                  else
                    ret <- kernelMult(kernelf(object),newdata,xmatrix(object)[[p]],coef(object)[[p]]) - b(object)[p]

                  if(type=="decision")
                    votematrix[p,] <- ret
                  else{
                    votematrix[i,ret<0] <- votematrix[i,ret<0] + 1
                    votematrix[j,ret>0] <- votematrix[j,ret>0] + 1
                  }
                }
            }
          if(type == "decision")
            predres <-  t(votematrix)
          else 
            predres <- sapply(predres, function(x) which.max(votematrix[,x]))
        }
      
  if(type(object) == "spoc-svc")
    {
      predres <- 1:newnrows
      votematrix <- matrix(0,nclass(object),newnrows)
      for(i in 1:nclass(object)){
        if(is(newdata,"kernelMatrix"))
          votematrix[i,] <- newdata[,which(SVindex(object)%in%alphaindex(object)[[i]]), drop=FALSE] %*% coef(object)[[i]] 
        else if (is(newdata,"list"))
          votematrix[i,] <- kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]]],coef(object)[[i]])
        else
          votematrix[i,] <- kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]],,drop=FALSE],coef(object)[[i]])
      }
      predres <- sapply(predres, function(x) which.max(votematrix[,x]))
    }

  if(type(object) == "kbb-svc")
    { 
      predres <- 1:newnrows
      votematrix <- matrix(0,nclass(object),newnrows)
      A <- rowSums(alpha(object))

      for(i in 1:nclass(object))
        {
          for(k in (1:i)[-i])
            if(is(newdata,"kernelMatrix"))
              votematrix[k,] <- votematrix[k,] - (newdata[,which(SVindex(object)%in%alphaindex(object)[[i]]), drop=FALSE] %*% alpha(object)[,k][alphaindex(object)[[i]]] + sum(alpha(object)[,k][alphaindex(object)[[i]]]))
            else if (is(newdata,"list"))
              votematrix[k,] <- votematrix[k,] - (kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]]],alpha(object)[,k][alphaindex(object)[[i]]]) + sum(alpha(object)[,k][alphaindex(object)[[i]]]))
            else
              votematrix[k,] <- votematrix[k,] - (kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]],,drop=FALSE],alpha(object)[,k][alphaindex(object)[[i]]]) + sum(alpha(object)[,k][alphaindex(object)[[i]]]))

          if(is(newdata,"kernelMatrix"))
            votematrix[i,] <- votematrix[i,] + (newdata[,which(SVindex(object)%in%alphaindex(object)[[i]]), drop=FALSE] %*% A[alphaindex(object)[[i]]] + sum(A[alphaindex(object)[[i]]]))
          else if (is(newdata,"list"))
            votematrix[i,] <- votematrix[i,] + (kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]]],A[alphaindex(object)[[i]]]) + sum(A[alphaindex(object)[[i]]]))
          else
            votematrix[i,] <- votematrix[i,] + (kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]],,drop=FALSE],A[alphaindex(object)[[i]]]) + sum(A[alphaindex(object)[[i]]]))

          if(i <= (nclass(object)-1))
            for(kk in i:(nclass(object)-1))
              if(is(newdata,"kernelMatrix"))
                votematrix[kk+1,] <- votematrix[kk+1,] - (newdata[,which(SVindex(object)%in%alphaindex(object)[[i]]), drop=FALSE] %*% alpha(object)[,kk][alphaindex(object)[[i]]] + sum(alpha(object)[,kk][alphaindex(object)[[i]]]))
              else if (is(newdata,"list"))
                votematrix[kk+1,] <- votematrix[kk+1,] - (kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]]],alpha(object)[,kk][alphaindex(object)[[i]]]) + sum(alpha(object)[,kk][alphaindex(object)[[i]]]))
              else
                votematrix[kk+1,] <- votematrix[kk+1,] - (kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]],,drop=FALSE],alpha(object)[,kk][alphaindex(object)[[i]]]) + sum(alpha(object)[,kk][alphaindex(object)[[i]]]))
        }
      predres <- sapply(predres, function(x) which.max(votematrix[,x]))
    }
}

  if(type == "probabilities")
    { 
      if(is.null(prob.model(object)[[1]]))
        stop("ksvm object contains no probability model. Make sure you set the paramater prob.model in ksvm during training.")
      
      if(type(object)=="C-svc"||type(object)=="nu-svc"||type(object)=="C-bsvc")
        {
          binprob <- matrix(0, newnrows, nclass(object)*(nclass(object) - 1)/2)
          for(i in 1:(nclass(object)-1))
            {
              jj <- i+1
              for(j in jj:nclass(object))
                {
                  p <- p+1
                  if(is(newdata,"kernelMatrix"))
                    binprob[,p] <- 1 - .SigmoidPredict(as.vector(newdata[,which(SVindex(object)%in%alphaindex(object)[[p]]), drop=FALSE] %*% coef(object)[[p]] - b(object)[p]), prob.model(object)[[p]]$A, prob.model(object)[[p]]$B)
                  else
                    binprob[,p] <- 1 - .SigmoidPredict(as.vector(kernelMult(kernelf(object),newdata,xmatrix(object)[[p]],coef(object)[[p]]) - b(object)[p]), prob.model(object)[[p]]$A, prob.model(object)[[p]]$B)
                }
            }
          multiprob <- couple(binprob, coupler = coupler)
        }
      else
        stop("probability estimates only supported for C-svc, C-bsvc and nu-svc")
    }
  
  if(type(object) == "one-svc")
    {
      if(is(newdata,"kernelMatrix"))
        ret <- newdata %*% coef(object) - b(object)
      else
        ret <- kernelMult(kernelf(object),newdata,xmatrix(object),coef(object)) - b(object)
      ##one-class-classification: return TRUE/FALSE (probabilities ?)
      if(type=="decision")
      	return(ret)
	else
	{
	ret[ret>0]<-1
      	return(ret == 1)
      }      
    }
  else {
    if(type(object)=="eps-svr"||type(object)=="nu-svr"||type(object)=="eps-bsvr")
      {
        if(is(newdata,"kernelMatrix"))
          predres <- newdata %*% coef(object) - b(object)
       else
         predres <- kernelMult(kernelf(object),newdata,xmatrix(object),coef(object)) - b(object)
      }
    else {
      ##classification & votes : return votematrix
      if(type == "votes")
        return(votematrix)
      
      ##classification & probabilities : return probability matrix
      if(type == "probabilities")
        {
          colnames(multiprob) <- lev(object)
          return(multiprob)
        }

      if(is.numeric(lev(object)) && type == "response")
         return(lev(object)[predres])
      
      if (is.character(lev(object)) && type!="decision")
        {
          ##classification & type response: return factors
          if(type == "response")
            return(factor (lev(object)[predres], levels = lev(object)))
        }
    }
  }
 
  if (!is.null(scaling(object)$y.scale) & !is(newdata,"kernelMatrix") & !is(newdata,"list"))
    ## return raw values, possibly scaled back
    return(predres * scaling(object)$y.scale$"scaled:scale" + scaling(object)$y.scale$"scaled:center")
  else
    ##else: return raw values
    return(predres)
})



#****************************************************************************************#

setMethod("show","ksvm",
function(object){
  cat("Support Vector Machine object of class \"ksvm\"","\n")
  cat("\n")
   cat(paste("SV type:", type(object)))

  
  switch(type(object),
         "C-svc" = cat(paste("  (classification)", "\n")),
         "nu-svc" = cat(paste("  (classification)", "\n")),
         "C-bsvc" = cat(paste("  (classification)", "\n")),
         "one-svc" = cat(paste("  (novelty detection)", "\n")),
         "spoc-svc" = cat(paste("  (classification)", "\n")),
         "kbb-svc" = cat(paste("  (classification)", "\n")),
         "eps-svr" = cat(paste("  (regression)","\n")),
         "nu-svr" = cat(paste("  (regression)","\n"))
         )
  
  switch(type(object),
         "C-svc" = cat(paste(" parameter : cost C =",param(object)$C, "\n")),
         "nu-svc" = cat(paste(" parameter : nu =", param(object)$nu, "\n")),
         "C-bsvc" = cat(paste(" parameter : cost C =",param(object)$C, "\n")),
         "one-svc" = cat(paste(" parameter : nu =", param(object)$nu, "\n")),
         "spoc-svc" = cat(paste(" parameter : cost C =",param(object)$C, "\n")),
         "kbb-svc" = cat(paste(" parameter : cost C =",param(object)$C, "\n")),
         "eps-svr" = cat(paste(" parameter : epsilon =",param(object)$epsilon, " cost C =", param(object)$C,"\n")),
         "nu-svr" = cat(paste(" parameter : epsilon =", param(object)$epsilon, " nu =", param(object)$nu,"\n"))
         )
  cat("\n")

  
  show(kernelf(object))
  cat(paste("\nNumber of Support Vectors :", nSV(object),"\n"))

  cat("\nObjective Function Value :", round(obj(object),4),"\n")
  

##    if(type(object)=="C-svc" || type(object) == "nu-svc")
##      cat(paste("Margin width :",margin(object),"\n"))
  if(!is.null(fitted(object)))
    cat(paste("Training error :", round(error(object),6),"\n"))
  if(cross(object)!= -1)
    cat("Cross validation error :",round(cross(object),6),"\n")
  if(!is.null(prob.model(object)[[1]])&&(type(object)=="eps-svr" ||type(object)=="nu-svr"||type(object)=="eps-bsvr"))
    cat("Laplace distr. width :",round(prob.model(object)[[1]],6),"\n")
  if(!is.null(prob.model(object)[[1]]) & (type(object) == "C-svc"| type(object) == "nu-svc"| type(object) == "C-bsvc"))
    cat("Probability model included.","\n")

  ##train error & loss
})


setMethod("plot", signature(x = "ksvm", y = "missing"),
function(x, data = NULL, grid = 50, slice = list(), ...) {

  if (type(x) =="C-svc" || type(x) == "nu-svc") {
    if(nclass(x) > 2)
      stop("plot function only supports binary classification")
  
    if (!is.null(terms(x))&&!is.null(data))
      {
        if(!is.matrix(data))
          sub <- model.matrix(delete.response(terms(x)), as.data.frame(data), na.action = n.action(x))
      }
    else if(!is.null(data))
      sub <-  as.matrix(data)
    else
      sub <- xmatrix(x)[[1]]

##    sub <- sub[,!colnames(xmatrix(x)[[1]])%in%names(slice)]
    xr <- seq(min(sub[,2]), max(sub[,2]), length = grid)
    yr <- seq(min(sub[,1]), max(sub[,1]), length = grid)
    sc <- 0

#    if(is.null(data))
#      {
#        sc  <- 1
#        data <- xmatrix(x)[[1]]
#      }

    if(is.data.frame(data) || !is.null(terms(x))){
      lis <- c(list(yr), list(xr), slice)
      names(lis)[1:2] <- setdiff(colnames(sub),names(slice))
      new <- expand.grid(lis)[,labels(terms(x))]
    }
    else
      new <- expand.grid(xr,yr)
    
    if(sc== 1) 
      scaling(x) <- NULL

    preds <- predict(x, new ,type = "decision")
    
    if(is.null(terms(x)))
      xylb <- colnames(sub)
    else
      xylb <- names(lis)
    lvl <- 37
    
    mymax <- max(abs(preds))
    mylevels <- pretty(c(0, mymax), 15)
    nl <- length(mylevels)-2
    
    mycols <- c(hcl(0, 100 * (nl:0/nl)^1.3, 90 - 40 *(nl:0/nl)^1.3),
    		rev(hcl(260, 100 * (nl:0/nl)^1.3, 90 - 40 *(nl:0/nl)^1.3)))

    mylevels <- c(-rev(mylevels[-1]), mylevels)

    index <- max(which(mylevels < min(preds))):min(which(mylevels > max(preds)))
    mycols <- mycols[index]
    mylevels <- mylevels[index]
    
    #FIXME# previously the plot code assumed that the y values are either
    #FIXME# -1 or 1, but this is not generally true. If generated from a
    #FIXME# factor, they are typically 1 and 2. Maybe ymatrix should be
    #FIXME# changed?
    ymat <- ymatrix(x)
    ymean <- mean(unique(ymat))

    filled.contour(xr, yr, matrix(as.numeric(preds), nrow = length(xr), byrow = TRUE),
                   col = mycols, levels = mylevels,
		   plot.axes = {
                     axis(1)
                     axis(2)
                     if(!is.null(data)){
                       points(sub[-SVindex(x),2], sub[-SVindex(x),1],
                             pch = ifelse(ymat[-SVindex(x)] < ymean, 2, 1))
                       points(sub[SVindex(x),2], sub[SVindex(x),1],
                              pch = ifelse(ymat[SVindex(x)] < ymean, 17, 16))}
                     else{
                       ## points(sub[-SVindex(x),], pch = ifelse(ymat[-SVindex(x)] < ymean, 2, 1))
                       points(sub,
		         pch = ifelse(ymat[SVindex(x)] < ymean, 17, 16))
                     }},
                   nlevels = lvl,
                   plot.title = title(main = "SVM classification plot", xlab = xylb[2], ylab = xylb[1]),
                   ...
                   )
  } else {
    stop("Only plots of classification ksvm objects supported")
  }
})


setGeneric(".probPlatt", function(deci, yres) standardGeneric(".probPlatt"))
setMethod(".probPlatt",signature(deci="ANY"),
function(deci,yres)
  { 
    if (is.matrix(deci))
      deci <- as.vector(deci)
    if (!is.vector(deci))
      stop("input should be matrix or vector")
    yres <- as.vector(yres)
    ## Create label and count priors
    boolabel <- yres >= 0
    prior1 <- sum(boolabel)
    m <- length(yres)
    prior0 <- m - prior1

    ## set parameters (should be on the interface I guess)
    maxiter <- 100
    minstep <- 1e-10
    sigma <- 1e-3
    eps <- 1e-5

    ## Construct target support
    hiTarget <- (prior1 + 1)/(prior1 + 2)
    loTarget <- 1/(prior0 + 2)
    length <- prior1 + prior0
    t <- rep(loTarget, m)
    t[boolabel] <- hiTarget

    ##Initial Point & Initial Fun Value
    A <- 0
    B <- log((prior0 + 1)/(prior1 + 1))
    fval <- 0

    fApB <- deci*A + B
    bindex <- fApB >= 0
    p <- q <- rep(0,m)
    
    fval <- sum(t[bindex]*fApB[bindex] + log(1 + exp(-fApB[bindex])))
    fval <- fval + sum((t[!bindex] - 1)*fApB[!bindex] + log(1+exp(fApB[!bindex])))

    for (it in 1:maxiter)
      {
        h11 <- h22 <- sigma
        h21 <- g1 <- g2 <- 0
        fApB <- deci*A + B
        
        bindex <- fApB >= 0
        p[bindex] <- exp(-fApB[bindex])/(1 + exp(-fApB[bindex]))
        q[bindex] <- 1/(1+exp(-fApB[bindex]))

        bindex <- fApB < 0
        p[bindex] <- 1/(1 + exp(fApB[bindex]))
        q[bindex] <- exp(fApB[bindex])/(1 + exp(fApB[bindex]))
                  
        d2 <- p*q
        h11 <- h11 + sum(d2*deci^2)
        h22 <- h22 + sum(d2)
        h21 <- h21 + sum(deci*d2)
        d1 <- t - p
        g1 <- g1 + sum(deci*d1)
        g2 <- g2 + sum(d1)
        
        ## Stopping Criteria
        if (abs(g1) < eps && abs(g2) < eps)
          break

        ## Finding Newton Direction -inv(t(H))%*%g
        det <- h11*h22 - h21^2
        dA <- -(h22*g1 - h21*g2) / det
        dB <- -(-h21*g1 + h11*g2) / det
        gd <- g1*dA + g2*dB

        ## Line Search
        stepsize <- 1

        while(stepsize >= minstep)
          {
            newA <- A + stepsize * dA
            newB <- B + stepsize * dB

            ## New function value
            newf <- 0
            fApB <- deci * newA + newB
            bindex <- fApB >= 0 
            newf <- sum(t[bindex] * fApB[bindex] + log(1 + exp(-fApB[bindex])))
            newf <- newf + sum((t[!bindex] - 1)*fApB[!bindex] + log(1 + exp(fApB[!bindex])))

           ## Check decrease
            if (newf < (fval + 0.0001 * stepsize * gd))
              {
                A <- newA
                B <- newB
                fval <- newf
                break
              }
            else
              stepsize <- stepsize/2
          }
        if (stepsize < minstep)
          {
            cat("line search fails", A, B, g1, g2, dA, dB, gd)
            ret <- .SigmoidPredict(deci, A, B)
            return(ret)
          }
      }
    if(it >= maxiter -1)
      cat("maximum number of iterations reached",g1,g2)

    ret <- list(A=A, B=B)
    return(ret)
  })

 ## Sigmoid predict function

.SigmoidPredict <- function(deci, A, B)
  {
fApB <- deci*A +B
k <- length(deci)
ret <- rep(0,k)
bindex <- fApB >= 0
ret[bindex] <- exp(-fApB[bindex])/(1 + exp(-fApB[bindex]))
ret[!bindex] <- 1/(1 + exp(fApB[!bindex]))
return(ret)
}

