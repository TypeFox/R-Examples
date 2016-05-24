predict.ordPen <- function(object, newx, newu = NULL, newz = NULL,
  offset = rep(0,nrow(as.matrix(newx))), type = c("link", "response"), ...)
  {
    type <- match.arg(type)
    TYPE <- switch(type, link="link", response="response")

    ## Check the newx matrix
    if(!(is.matrix(newx) | is.numeric(newx) | is.data.frame(newx)))
      stop("newx has to be a matrix, numeric vector or data.frame")

    if(any(is.na(newx)))
      stop("Missing values in newx are not allowed")
     
    newx <- as.matrix(newx)
    if(any(!apply(newx,2,is.numeric)))
      stop("Entries of newx have to be of type 'numeric'")
      
    tol <- .Machine$double.eps^0.5  
    if(any(abs(newx - round(newx)) > tol) | any(newx < 1))
      stop("newx has to contain positive integers")

    if(ncol(newx) != length(object$xlevels))
      stop("inadequate ncol(newx)")

    if(any(apply(newx,2,max) > object$xlevels))
      stop("too high class levels in newx")

    ## Check the other arguments
    if(!is.null(object$ulevels) & is.null(newu))
      stop("no newu given")
      
    if(object$zcovars > 0 & is.null(newz))
      stop("no newz given")


    if (!is.null(object$ulevels))
      {
        if(!(is.matrix(newu) | is.numeric(newu) | is.data.frame(newu)))
          stop("newu has to be a matrix, numeric vector or data.frame")
        if(any(is.na(newu)))
          stop("Missing values in newu are not allowed")
        
        newu <- as.matrix(newu)
        if(nrow(newu) != nrow(newx))
          stop("newu and newx do not have correct dimensions")
        if(any(!apply(newu,2,is.numeric)))
          stop("Entries of newu have to be of type 'numeric'")
        if(any(abs(newu - round(newu)) > tol) | any(newu < 1))
          stop("newu has to contain positive integers")
        if(ncol(newu) != length(object$ulevels))
          stop("inadequate ncol(newu)")
        if(any(apply(newu,2,max) > object$ulevels))
          stop("too high class levels in newu")
          
        newcat <- rbind(c(object$xlevels,object$ulevels),
        cbind(newx,newu))
      }
    else
      {
        if (!is.null(newu))
          warning("newu not used")
        newcat <- rbind(object$xlevels,newx)
      }
    xuz <- rbind(coding(newcat + 1, constant=TRUE, splitcod=FALSE)[-1,])
    if (object$zcovars > 0)
      {
        if(!(is.matrix(newz) | is.numeric(newz) | is.data.frame(newz)))
          stop("newz has to be a matrix, numeric vector or data.frame")
        if(any(is.na(newz)))
          stop("Missing values in newz are not allowed")
     
        newz <- rbind(as.matrix(newz))
        if(nrow(newz) != nrow(newx))
          stop("newz and newx do not have correct dimensions")
        if(any(!apply(newz,2,is.numeric)))
          stop("Entries of newz have to be of type 'numeric'")
        if(ncol(newz) != object$zcovars)
          stop("inadequate ncol(newz)")
   
        xuz <- cbind(xuz,newz)
      }
    else if (!is.null(newz))
      warning("newz not used")
    
    rownames(xuz) <- rownames(newx)
    eta <- xuz%*%object$coef + matrix(rep(offset,ncol(object$coef)),
    nrow(xuz),ncol(object$coef))
    if (object$model == "linear" | TYPE == "link")
      {
        fit <- eta
      }
    else
      {
        if (object$model == "logit")
          fit <- plogis(eta)
        else
          fit <- exp(eta)
      }
    return(fit)
  }