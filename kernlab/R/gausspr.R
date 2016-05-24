## Gaussian Processes implementation. Laplace approximation for classification.
## author : alexandros karatzoglou

setGeneric("gausspr", function(x, ...) standardGeneric("gausspr"))
setMethod("gausspr",signature(x="formula"),
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
  
  ret <- gausspr(x, y, scaled = scaled, ...)
  kcall(ret) <- cl
  terms(ret) <- Terms
  if (!is.null(attr(m, "na.action")))
    n.action(ret) <- attr(m, "na.action")
  return (ret)
})

setMethod("gausspr",signature(x="vector"),
function(x,...)
  {
    x <- t(t(x))
    ret <- gausspr(x, ...)
    ret
  })
    
setMethod("gausspr",signature(x="matrix"),
function (x,
          y,
          scaled    = TRUE, 
          type      = NULL,
          kernel    = "rbfdot",
          kpar      = "automatic",
          var       = 1,
          variance.model = FALSE,
          tol       = 0.0005,  
          cross     = 0,
          fit       = TRUE,
          ...
          ,subset 
         ,na.action = na.omit)
{

## should become an option
  reduced <- FALSE
## subsetting and na-handling for matrices
  ret <- new("gausspr")
  if (!missing(subset)) x <- x[subset,]
  if (is.null(y))
    x <- na.action(x)
  else {
    df <- na.action(data.frame(y, x))
    y <- df[,1]
    x <- as.matrix(df[,-1])
  }
  ncols <- ncol(x)
  m <- nrows <- nrow(x)

 if (is.null (type)) type(ret) <-
   if (is.factor(y)) "classification"
    else "regression"
  else type(ret) <- type
  
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
      if (is.numeric(y)&&(type(ret)!="classification")) {
        y <- scale(y)
        y.scale <- attributes(y)[c("scaled:center","scaled:scale")]
        y <- as.vector(y)
      }
      tmpsc <- list(scaled = scaled, x.scale = x.scale, y.scale = y.scale)
    }
  }

  
  if (var < 10^-3)
    stop("Noise variance parameter var has to be greater than 10^-3")

 # in case of classification: transform factors into integers
  if (is.factor(y)) {
    lev(ret) <- levels (y)
    y <- as.integer (y)
  }
  else {
    if (type(ret) == "classification" && any(as.integer (y) != y))
        stop ("dependent variable has to be of factor or integer type for classification mode.")
    if(type(ret) == "classification")
      lev(ret) <- unique (y)
    }
 # initialize    
  nclass(ret) <- length (lev(ret))
  
  if(!is.null(type))
    type(ret) <- match.arg(type,c("classification", "regression"))

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
      kpar <- list(sigma=mean(sigest(x,scaled=FALSE)[c(1,3)]))
   cat("Using automatic sigma estimation (sigest) for RBF or laplace kernel","\n")
   
  }

  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }
  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

  p <- 0
  
  if (type(ret) == "classification")
    {
      indexes <- lapply(1:nclass(ret), function(kk) which(y == kk))
      for (i in 1:(nclass(ret)-1)) {
        jj <- i+1
        for(j in jj:nclass(ret)) {
          p <- p+1
          ##prepare data
          li <- length(indexes[[i]])
          lj <- length(indexes[[j]])
          xd <- matrix(0,(li+lj),dim(x)[2])
          xdi <- 1:(li+lj) <= li
          xd[xdi,rep(TRUE,dim(x)[2])] <- x[indexes[[i]],]
          xd[xdi == FALSE,rep(TRUE,dim(x)[2])] <- x[indexes[[j]],]
          if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
            yd <- c(rep(1,li),rep(-1,lj))
          else
            yd <- c(rep(-1,li),rep(1,lj))
          if(reduced == FALSE){
            K <- kernelMatrix(kernel,xd)
            gradnorm <- 1 
            alphag <- solut <- rep(0,li+lj)
            while (gradnorm > tol)
              {
                f <- crossprod(K,alphag)
                grad <- -yd/(1 + exp(yd*f))
                hess <- exp(yd*f)
                hess <- hess / ((1 + hess)^2)

                ## We use solveiter instead of solve to speed up things
                ## A <- t(t(K)*as.vector(hess))
                ## diag(A) <- diag(A) + 1
                ## alphag <- alphag - solve(A,(grad + alphag))

                solut <- solveiter(K, hess, (grad + alphag), solut)
                alphag <- alphag - solut
                gradnorm <- sqrt(sum((grad + alphag)^2))
              }
          }
          else if (reduced ==TRUE)
            {
             
              yind <- t(matrix(unique(yd),2,length(yd)))
              ymat <- matrix(0, length(yd), 2)
              ymat[yind==yd] <- 1
              ##Z <- csi(xd, ymat, kernel = kernel, rank = dim(yd)[1])
              ##Z <- Z[sort(pivots(Z),index.return = TRUE)$ix, ,drop=FALSE]
              Z <- inchol(xd, kernel = kernel)
              gradnorm <- 1 
              alphag <- rep(0,li+lj)
              m1 <- dim(Z)[1]
              n1 <- dim(Z)[2]
              Ksub <- diag(rep(1,n1))
           
           while (gradnorm > tol)
              {
                f <- drop(Z%*%crossprod(Z,alphag))
                f[which(f>20)] <- 20
                grad <- -yd/(1 + exp(yd*f))
                hess <- exp(yd*f)
                hess <- as.vector(hess / ((1 + hess)^2))
                
                alphag <- alphag - (- Z %*%solve(Ksub + (t(Z)*hess)%*%Z) %*% (t(Z)*hess))%*%(grad + alphag) + (grad + alphag) 
                
                gradnorm <- sqrt(sum((grad + alphag)^2))
              }
              
            }
              alpha(ret)[[p]] <- alphag
              alphaindex(ret)[[p]] <- c(indexes[[i]],indexes[[j]])
        }
      }
    }
  
  if (type(ret) == "regression")
    {
      K <- kernelMatrix(kernel,x)
      if(variance.model)
        {
          sol <- solve(K + diag(rep(var, length = m)))
          rm(K)
          alpha(ret) <- sol%*%y
        }
      else
        alpha(ret) <- solve(K + diag(rep(var, length = m))) %*% y
      
    }

  kcall(ret) <- match.call()
  kernelf(ret) <- kernel
  xmatrix(ret) <- x
  if(variance.model)
      sol(ret) <- sol
  
  fitted(ret)  <- if (fit)
    predict(ret, x) else NA

  if (fit){
    if(type(ret)=="classification")
      error(ret) <- 1 - .classAgreement(table(y,as.integer(fitted(ret))))
    if(type(ret)=="regression"){
      if (!is.null(scaling(ret)$y.scale))
        fitted(ret) <- fitted(ret) * tmpsc$y.scale$"scaled:scale" + tmpsc$y.scale$"scaled:center"
      
    error(ret) <- drop(crossprod(fitted(ret) - y)/m)
    }
  }
  if(any(scaled))
    scaling(ret) <- tmpsc
  
  cross(ret) <- -1
  if(cross == 1)
    cat("\n","cross should be >1 no cross-validation done!","\n","\n")
  else if (cross > 1)
    {
      cerror <- 0
      suppressWarnings(vgr<-split(sample(1:m,m),1:cross))
      for(i in 1:cross)
        {
          cind <- unsplit(vgr[-i],factor(rep((1:cross)[-i],unlist(lapply(vgr[-i],length)))))
          if(type(ret)=="classification")
            {
              cret <- gausspr(x[cind,], y[cind], scaled = FALSE, type=type(ret),kernel=kernel,var = var, cross = 0, fit = FALSE)
               cres <- predict(cret, x[vgr[[i]],])
            cerror <- (1 - .classAgreement(table(y[vgr[[i]]],as.integer(cres))))/cross + cerror
            }
          if(type(ret)=="regression")
            {
              cret <- gausspr(x[cind,],y[cind],type=type(ret),scaled = FALSE, kernel=kernel,var = var,tol=tol, cross = 0, fit = FALSE)
              cres <- predict(cret, x[vgr[[i]],])
              if (!is.null(scaling(ret)$y.scale))
                scal <- scaling(ret)$y.scale$"scaled:scale"
              cerror <- drop((scal^2)*crossprod(cres - y[vgr[[i]]])/m) + cerror
            }
        }
      cross(ret) <- cerror
    }


  
  return(ret)
})


setMethod("predict", signature(object = "gausspr"),
function (object, newdata, type = "response", coupler = "minpair")
{
  sc <- 0
  type <- match.arg(type,c("response","probabilities","votes", "variance", "sdeviation"))
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
  
  p <- 0
  if(type == "response")
    {
  if(type(object)=="classification")
    {
      predres <- 1:newnrows
      votematrix <- matrix(0,nclass(object),nrows)
      for(i in 1:(nclass(object)-1))
        {
        jj <- i+1
        for(j in jj:nclass(object))
          {
            p <- p+1
            ret <- kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[p]],],alpha(object)[[p]])
            votematrix[i,ret>0] <- votematrix[i,ret>0] + 1
            votematrix[j,ret<0] <- votematrix[j,ret<0] + 1
          }
      }
      predres <- sapply(predres, function(x) which.max(votematrix[,x]))
    }

}

  if(type == "probabilities")
    {
      if(type(object)=="classification")
        {
          binprob <- matrix(0, newnrows, nclass(object)*(nclass(object) - 1)/2)
          for(i in 1:(nclass(object)-1))
            {
              jj <- i+1
              for(j in jj:nclass(object))
                {
                  p <- p+1
                  binprob[,p] <-  1/(1+exp(-kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[p]],],alpha(object)[[p]])))

                }
            }
          ## multiprob <- sapply(1:newnrows, function(x) couple(binprob[x ,],coupler = coupler))
          multiprob <- couple(binprob, coupler = coupler)
        }
    }
    
  
  if(type(object) == "regression")
    {
      if (type == "variance"||type == "sdeviation")
        { 
          Ktest <- kernelMatrix(kernelf(object),xmatrix(object), newdata)
          predres <- diag(kernelMatrix(kernelf(object),newdata) - t(Ktest)  %*% sol(object) %*% Ktest)
          if (type== "sdeviation")
            predres <- sqrt(predres)
          if (!is.null(scaling(object)$y.scale))
           predres <- predres * scaling(object)$y.scale$"scaled:scale" + scaling(object)$y.scale$"scaled:center"
        }

      else
        {
      
      predres <- kernelMult(kernelf(object),newdata,xmatrix(object),as.matrix(alpha(object)))

      if (!is.null(scaling(object)$y.scale))
        predres <- predres * scaling(object)$y.scale$"scaled:scale" + scaling(object)$y.scale$"scaled:center"
    }
      
   }


 if (is.character(lev(object)))
    {
      ##classification & probabilities : return probabilitie matrix
      if(type == "probabilities")
        {
          colnames(multiprob) <- lev(object)
          return(multiprob)
        }
      ##classification & type response: return factors
      if(type == "response")
        return(factor (lev(object)[predres], levels = lev(object)))
      ##classification & votes : return votematrix
      if(type == "votes")
        return(votematrix)
    }
  else
    ##else: return raw values
    return(predres)

})


setMethod("show","gausspr",
function(object){
  cat("Gaussian Processes object of class \"gausspr\"","\n")
  cat(paste("Problem type:", type(object),"\n"))
  cat("\n")
  show(kernelf(object))
  cat(paste("\nNumber of training instances learned :", dim(xmatrix(object))[1],"\n"))
  if(!is.null(fitted(object)))
    cat(paste("Train error :", round(error(object),9),"\n"))
  ##train error & loss
  if(cross(object)!=-1)
    cat("Cross validation error :",round(cross(object),9),"\n")
})


solveiter <- function(B,noiseproc,b,x,itmax = 50,tol = 10e-4 ,verbose = FALSE) {

## ----------------------------
## Preconditioned Biconjugate Gradient method
## solves linear system Ax <- b for general A
## ------------------------------------------
## x : initial guess
## itmax : max # iterations
## iterates while mean(abs(Ax-b)) > tol
##
## Simplified form of Numerical Recipes: linbcg
## 
## The preconditioned matrix is set to inv(diag(A))
## A defined through A <- I + N*B

diagA <- matrix(1,dim(B)[1],1) + colSums(B)+ diag(B)*(noiseproc-1)
## diags of A

cont <- 0
iter <- 0
r <- .Amul2(x,B,noiseproc)
r <- b - r
rr <- r
znrm <- 1

bnrm <- sqrt(sum((b)^2))
z <- r/diagA

err <- sqrt(sum((.Amul2(x,B,noiseproc) - b)^2))/bnrm

while (iter <= itmax){
  iter <- iter + 1
  zm1nrm <- znrm
  zz <- rr/diagA
  bknum<- drop(crossprod(z,rr))
  if (iter == 1)
    {
      p <- z
      pp <- zz
    }
  else
    {
      bk <- bknum/bkden
      p <- bk*p + z
      pp <- bk*pp + zz
    } 
  
  bkden <- bknum
  z <- .Amul2(p,B,noiseproc)
  akden <- drop(crossprod(z,pp))
  ak <- bknum/akden
  zz <- .Amul2T(pp,B,noiseproc)
  
  x <- x + ak*p
  r <- r - ak*z
  rr <- rr - ak*zz
  z <- r/diagA
  znrm <- 1
  
  err <- mean(abs(r))
  
  if (err<tol)
    break
}    
return(x)
}

.Amul2 <- function(d, B, noiseproc){
ee <- B%*%d
return(d + noiseproc*ee)
}
.Amul2T <- function(d, B, noiseproc){
ee <- noiseproc*d
return(d + B%*%ee)
}
