## kernel based on-line learning algorithms for classification, novelty detection and regression. 
##
## created 15.09.04 alexandros
## updated

setGeneric("onlearn",function(obj, x, y = NULL, nu = 0.2, lambda = 1e-4) standardGeneric("onlearn"))
setMethod("onlearn", signature(obj = "onlearn"),
          function(obj , x, y = NULL, nu = 0.2, lambda = 1e-4)
          {
            if(onstart(obj) == 1 && onstop(obj) < buffer(obj))
              buffernotfull  <- TRUE
            else
              buffernotfull <- FALSE

            if(is.vector(x))
              x <- matrix(x,,length(x))  
            d <- dim(x)[2]
            
         
            for (i in 1:dim(x)[1])
              {
                xt <- x[i,,drop=FALSE]
                yt <- y[i]
            if(type(obj)=="novelty")
              {
                phi <- fit(obj)
                if(phi < 0)
                  {
                    alpha(obj) <- (1-lambda) * alpha(obj)
                    if(buffernotfull)
                      onstop(obj) <- onstop(obj) + 1
                    else{
                      onstop(obj) <- onstop(obj)%%buffer(obj) + 1
                      onstart(obj) <- onstart(obj)%%buffer(obj) +1
                    }
                    alpha(obj)[onstop(obj)] <- lambda
                    xmatrix(obj)[onstop(obj),] <- xt
                    rho(obj) <- rho(obj) + lambda*(nu-1)
                  }
                else
                  rho(obj) <- rho(obj) + lambda*nu
                
                rho(obj) <- max(rho(obj), 0)

                if(onstart(obj) == 1 && onstop(obj) < buffer(obj))
                  fit(obj) <- drop(kernelMult(kernelf(obj), xt, matrix(xmatrix(obj)[1:onstop(obj),],ncol=d), matrix(alpha(obj)[1:onstop(obj)],ncol=1)) - rho(obj)) 
                else
                  fit(obj) <- drop(kernelMult(kernelf(obj), xt, xmatrix(obj), matrix(alpha(obj),ncol=1)) - rho(obj))
              }
            if(type(obj)=="classification")
              { 
                if(is.null(pattern(obj)) && is.factor(y))
                  pattern(obj) <- yt
                if(!is.null(pattern(obj)))
                  if(pattern(obj) == yt)
                    yt <- 1
                  else yt <-  -1
                
                phi <- fit(obj)
                
                alpha(obj) <- (1-lambda) * alpha(obj)

                if(yt*phi < rho(obj))
                  {
                    if(buffernotfull)
                      onstop(obj) <- onstop(obj) + 1
                    else{
                      onstop(obj) <- onstop(obj)%%buffer(obj) + 1
                      onstart(obj) <- onstart(obj)%%buffer(obj) +1
                    }
                    alpha(obj)[onstop(obj)] <- lambda*yt
                    b(obj) <- b(obj) + lambda*yt
                    xmatrix(obj)[onstop(obj),] <- xt
                    rho(obj) <- rho(obj) + lambda*(nu-1) ## (1-nu) ??
                  }
                else
                  rho(obj) <- rho(obj) + lambda*nu
                
                rho(obj) <- max(rho(obj), 0)

                if(onstart(obj) == 1 && onstop(obj) < buffer(obj))
                  fit(obj) <- drop(kernelMult(kernelf(obj), xt, xmatrix(obj)[1:onstop(obj),,drop=FALSE], matrix(alpha(obj)[1:onstop(obj)],ncol=1)) + b(obj))
                else
                  fit(obj) <-drop(kernelMult(kernelf(obj), xt, xmatrix(obj), matrix(alpha(obj),ncol=1)) + b(obj))
          
              }

            if(type(obj)=="regression")
              {
                alpha(obj) <- (1-lambda) * alpha(obj)
                phi <- fit(obj)
                
                if(abs(-phi) < rho(obj))
                  {
                    if(buffernotfull)
                      onstop(obj) <- onstop(obj) + 1
                    else{
                      onstop(obj) <- onstop(obj)%%buffer(obj) + 1
                      onstart(obj) <- onstart(obj)%% buffer(obj) +1
                    }
                    alpha(obj)[onstop(obj)] <- sign(yt-phi)*lambda
                    xmatrix(obj)[onstop(obj),] <- xt
                    rho(obj) <- rho(obj) + lambda*(1-nu) ## (1-nu) ??
                  }
                else{
                  rho(obj) <- rho(obj) - lambda*nu
                  alpha(obj)[onstop(obj)] <- sign(yt-phi)/rho(obj)
                }
                if(onstart(obj) == 1 && onstop(obj) < buffer(obj))
                  fit(obj) <- drop(kernelMult(kernelf(obj), xt, matrix(xmatrix(obj)[1:onstop(obj),],ncol=d), matrix(alpha(obj)[1:onstop(obj)],ncol=1)) + b(obj)) 
                else
                  fit(obj) <- drop(kernelMult(kernelf(obj), xt, xmatrix(obj), matrix(alpha(obj),ncol=1)) + b(obj)) 
              }
          }
            return(obj)
          })


setGeneric("inlearn",function(d, kernel = "rbfdot", kpar = list(sigma=0.1), type = "novelty", buffersize = 1000) standardGeneric("inlearn"))
setMethod("inlearn", signature(d = "numeric"),
          function(d ,kernel = "rbfdot", kpar = list(sigma=0.1), type = "novelty", buffersize = 1000)
          {
            obj <- new("onlearn")

            if(!is(kernel,"kernel"))
              {
                if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
                kernel <- do.call(kernel, kpar)
              }
            if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

            type(obj) <- match.arg(type,c("novelty","classification","regression"))
            xmatrix(obj) <- matrix(0,buffersize,d)
            kernelf(obj) <- kernel
            onstart(obj) <- 1
            onstop(obj) <- 1
            fit(obj) <- 0 
            b(obj) <- 0
            alpha(obj) <- rep(0, buffersize)
            rho(obj) <- 0
            buffer(obj) <- buffersize
            return(obj)
          })


setMethod("show","onlearn",
function(object){
  cat("On-line learning object of class \"onlearn\"","\n")
  cat("\n")
  cat(paste("Learning problem :", type(object), "\n"))
  cat
  cat(paste("Data dimensions :", dim(xmatrix(object))[2], "\n"))
  cat(paste("Buffersize :", buffer(object), "\n"))
  cat("\n")
 show(kernelf(object))
})


setMethod("predict",signature(object="onlearn"),
function(object, x)
  {
    if(is.vector(x))
      x<- matrix(x,1)

    d <- dim(xmatrix(object))[2]
          
    if(type(object)=="novelty")
      {
        if(onstart(object) == 1 && onstop(object) < buffer(object))
          res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object)[1:onstop(object),],ncol= d), matrix(alpha(object)[1:onstop(object)],ncol=1)) - rho(object)) 
        else
          res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object),ncol=d), matrix(alpha(object),ncol=1)) - rho(object))
      }

    if(type(object)=="classification")
      {
        if(onstart(object) == 1 && onstop(object) < buffer(object))
          res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object)[1:onstop(object),],ncol=d), matrix(alpha(object)[1:onstop(object)],ncol=1)) + b(object))
        else
          res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object),ncol=d), matrix(alpha(object),ncol=1)) + b(object))
       
      }

    if(type(object)=="regression")
      {
        if(onstart(object) == 1 && onstop(object) < buffer(object))
          res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object)[1:onstop(object),],ncol=d), matrix(alpha(object)[1:onstop(object)],ncol=1)) + b(object)) 
        else
          res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object),ncol=d), matrix(alpha(object),ncol=1)) + b(object)) 
      }

    return(res)
    
  })

  
