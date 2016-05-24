##' Class "blockSeg"
##'
##' Class of object returned by the \code{blockSeg} function.
##'
##' @slot Beta a Matrix object of type \code{dgCMatrix},
##' encoding the solution path of the underlying LARS algorithm. Ommited
##' if the blockSeg function was called with the option
##' \code{Beta=FALSE}.
##'
##' @slot Lambda a numeric vector with the successive values
##' of \code{Lambda}, that is, the value of the penalty parameter
##' corresponding to a new event in the path (either a variable
##' activation or deactivation).
##' 
##' @slot RowBreaks a list of vectors, one per step of the 
##' LARS algorithm. Each vector contains the breaks currently identified
##' along the ROWS of the 2-dimensional signalat the current step.
##'
##' @slot ColBreaks a list of vectors, one per step of the
##' LARS algorithm. Each vector contains the breaks currently identified
##' along the COLUMNS of the 2-dimensional signal at the current step.
##'
##' @slot Actions a list with the successive actions at each
##' step of the LARS algorithm.
##' 
##' @param Y the original data matrix
##'
##' @param object an object with class blockSeg
##'
##' @param x in the print method, a blockSeg object 
##' 
##' @param ... in the print method, additional parameters (ignored)
##' 
##' @seealso See also \code{\link{plot,blockSeg-method}}, \code{\link{predict,blockSeg-method}}
##' and \code{\link{blockSeg}}.
##'
##' @rdname blockSeg-class
##'
##' @exportClass blockSeg
##' @exportMethod getBreaks
##' @exportMethod getComplexity
##' @exportMethod print
##' @exportMethod show
##' @exportMethod deviance
##' @exportMethod residuals
##'
setClass("blockSeg",
  representation = representation(
    Beta      = "dgCMatrix",
    Lambda    = "numeric",
    RowBreaks = "list",
    ColBreaks = "list",
    Actions   = "list")
)

##' @rdname blockSeg-class
setMethod("print", "blockSeg", definition =
   function(x, ...) {
     cat("Two-dimensional change points detection (blockSeg)\n")
     cat("- number of steps :", length(slot(x,"Lambda")), "\n")
     cat("- final number of breaks detected (rows) :", length(x@RowBreaks[[length(slot(x,"Lambda"))]]), "\n")
     cat("- final number of breaks detected (columns) :", length(x@ColBreaks[[length(slot(x,"Lambda"))]]), "\n")
     invisible(x)
   }
)

##' @rdname blockSeg-class
setMethod("show", "blockSeg", definition =
   function(object) {print(object)}
)

##' @rdname blockSeg-class
setGeneric("getComplexity",function(object){standardGeneric("getComplexity")})
##' @rdname blockSeg-class
setMethod("getComplexity", "blockSeg", definition =
   function (object)  {
       return(cumsum(sapply(object@Actions, function(x) sum(sign(x)))))
   })

getExpandedYhat <- function(Y.hat) {
    return(Y.hat$mu.hat[rep(1:nrow(Y.hat$mu.hat), Y.hat$RowCounts), rep(1:ncol(Y.hat$mu.hat), Y.hat$ColCounts)])
}


##' @rdname blockSeg-class
setMethod("residuals", "blockSeg", definition =
   function (object, Y)  {
       return(lapply(predict(object, Y), function(Y.hat) {
           return(getExpandedYhat(Y.hat) - Y)
       }))
       
   })

##' @rdname blockSeg-class
setMethod("deviance", "blockSeg", definition =
   function (object, Y)  {
       return(sapply(residuals(object, Y), function(x) sum(x^2)))
   })

##' @rdname blockSeg-class
setGeneric("getBreaks",function(object){standardGeneric("getBreaks")})
##' @rdname blockSeg-class
setMethod("getBreaks", "blockSeg", definition =
   function(object) {
    return(lapply(seq(length(object@Actions)), function(k) {
        return(list(RowBreaks = sort(unique(c(1,object@RowBreaks[[k]]))),
                    ColBreaks = sort(unique(c(1,object@ColBreaks[[k]])))))        
    }))
})

##' @rdname blockSeg-class
setGeneric("getCompressYhat",function(object, Y){standardGeneric("getCompressYhat")})
##' @rdname blockSeg-class
setMethod("getCompressYhat", "blockSeg", definition =
   function(object,Y) {
    Breaks <- getBreaks(object)
    return(lapply(seq(length(object@Actions)), function(k) {
        RowCounts <- diff(c(Breaks[[k]]$RowBreaks,nrow(Y)+1))
        ColCounts <- diff(c(Breaks[[k]]$ColBreaks,ncol(Y)+1))
        RowGroups <- rep(Breaks[[k]]$RowBreaks, RowCounts)
        ColGroups <- rep(Breaks[[k]]$ColBreaks, ColCounts)
        GrpSizes  <- tcrossprod(RowCounts,ColCounts)
        return(list(
            mu.hat = t(rowsum( t(rowsum(Y, RowGroups)), ColGroups)) / GrpSizes ,
            RowBreaks = Breaks[[k]]$RowBreaks, ColBreaks = Breaks[[k]]$ColBreaks,
            RowCounts = RowCounts, ColCounts = ColCounts))
    }))
})

##' Predict method for a blockSeg object
##'
##' Produce a prediction for a vector of \code{lambda} parameter and an array of \code{class}.
##'
##' @param object an object of class \code{blockSeg}.
##' @param Y matrix of observations.
##' @param lambda a numeric vector giving the list of \eqn{\lambda}{lambda} for which to predict.
##' By default, \code{NULL}. If \code{NULL}, it is set to the \code{lambdalist} slot
##' of \code{object}. If this slot is empty, \code{lambda} is set to the fusion times detected in
##' the \code{blockSeg} function.
##' @seealso \code{\linkS4class{blockSeg}}.
##' @rdname predict
##'
##' @examples
##' require(blockseg)
##' n <- 100
##' K <- 5
##' mu <- suppressWarnings(matrix(rep(c(1,0),ceiling(K**2/2)), K,K))
##' Y <- rblockdata(n,mu,sigma=.5)$Y
##' res <- blockSeg(Y, 100)
##' predict(res, Y, lambda=slot(res, "Lambda")[1:3])
##'
##' @exportMethod predict
setMethod("predict", "blockSeg", definition =
   function (object, Y, lambda=NULL)  {
       
       ## recover all the mu.hat 
       mu.hat <- getCompressYhat(object, Y)
       lambda.path <- object@Lambda

       if (is.null(lambda)){ # no new grid
           return(setNames(mu.hat, round(lambda.path,3)))           
       } else { # linear interpolation
           
           ## adjust the request value of lambda
           lambda.min <- min(lambda.path)
           lambda.max <- max(lambda.path)
           lambda[lambda >= lambda.max] <- lambda.max
           lambda[lambda <= lambda.min] <- lambda.min
           lambda <- sort(unique(lambda), decreasing=TRUE)
           
           ## compute the linear approximation along the path
           frac.hat  <- (lambda - lambda.max)/(lambda.min - lambda.max)
           frac.path <- (lambda.path - lambda.max)/(lambda.min - lambda.max)
           useq <- which(!duplicated(frac.path))
           frac.path <- frac.path[useq]
           coord <- approx(frac.path, seq(frac.path), frac.hat)$y
           left  <- floor(coord) ; right <- ceiling(coord)

           mu.hat.left  <- mu.hat[useq[left ]]
           mu.hat.right <- mu.hat[useq[right]]

           delta1 <- (frac.path[right] - frac.hat)
           delta2 <- (frac.hat - frac.path[left])
           delta3 <- (frac.path[right] - frac.path[left])

           ## the tricky thing is to expand the "left" mu.hat to match the dimension of the "right" one
           ## This is beacause we use a compress version of the fit
           new.mu.hat <- lapply(1:length(lambda), function(k) {
               ## find the row (or column) to duplicate to have left and right mu.hat comparable
               row.left <- match(mu.hat.right[[k]]$RowBreaks, mu.hat.left[[k]]$RowBreaks)
               col.left <- match(mu.hat.right[[k]]$ColBreaks, mu.hat.left[[k]]$ColBreaks)
               while(anyNA(row.left)) {
                   row.left[which(is.na(row.left))[1]] <- row.left[which(is.na(row.left))[1]-1]
               }
               while(anyNA(col.left)) {
                   col.left[which(is.na(col.left))[1]] <- row.left[which(is.na(col.left))[1]-1]
               }
               new.mu <- mu.hat.right[[k]]
               if (delta3[k]!=0){
                 new.mu$mu.hat <- (delta1[k] * mu.hat.left[[k]]$mu.hat[row.left,col.left] + delta2[k] * mu.hat.right[[k]]$mu.hat)/delta3[k]
               }
               return(new.mu)
           })
           
           #new.mu.hat[left == right] <- new.mu.hat[left[left == right]]

           ## dirty handling of the boundaries of the path (lambda.min/lambda.max)
           if (!is.na(match(lambda.max,lambda))) {
               new.mu.hat[[match(lambda.max,lambda)]] <- mu.hat[[1]]
           }

           if (!is.na(match(lambda.min,lambda))) {
               new.mu.hat[[match(lambda.min,lambda)]] <- mu.hat[[length(mu.hat)]]
           }

           return(setNames(new.mu.hat, round(lambda)))
       }       
})

##' Plot method for a blockSeg object
##'
##' Produce a plot of two-dimensional segmentation of a \code{blockSeg} fit.
##'
##' @param x an object of class \code{blockSeg}.
##' @param y used  for S4 compatibility.
##' @param lambda parameter used in the LASSO.
##' @param ask If \code{TRUE}, to hit will be necessary to see next plot.
##' @param col for the colors of the representations
##' @param ... used for S4 compatibility.
##'
##' @return a \pkg{ggplot2} object which can be plotted via the \code{print} method.
##' 
##' @rdname plot-blockSeg
##' @seealso \code{\linkS4class{blockSeg}}.
##'
##' @exportMethod plot
setMethod(
  f="plot",
  signature="blockSeg",
  definition=function(x,y,lambda=NULL,ask=TRUE,col="GrayLevel",...){

    if (is.null(lambda)) {
        lambda <- x@Lambda[length(x@Lambda)]
        Y.hat <- predict(x,y,lambda=lambda)
    } else if (!is.numeric(lambda)){
      stop("lambda must be numeric")
    } else {
        Y.hat <- predict(x,y,lambda)
    }
    if (!is.character(col)){
      stop("col must be a character")
    }else{
      if (length(col)==1){
        if (col=="GrayLevel"){
          mypalette=gray(seq(0,1, length=256))
        }else{
          mypalette=c(rgb((0:201)/201*200,(0:201)/201*200,255,maxColorValue = 255),"white",
                    rgb(255,(0:201)/201*200,(0:201)/201*200,maxColorValue = 255)[201:0])
        }
      }else{
        mypalette=col
      }
    }
    brek=seq(min(y),max(y),length=length(mypalette)+1)
    
    if (!is.logical(ask)){
      stop("ask must be equal to TRUE or FALSE")
    }else{
      if (length(lambda)==1){
        ask=FALSE
      }
    }
    
    n=nrow(y)
    d=ncol(y)
    
    par(ask=ask,mfrow=c(1,2),oma=c(0,0,3,0))
    for (k in 1:length(lambda)){
        ## Diplay of the original matrix
      image(1:d,1:n,t(y)[,n:1],xlab="",ylab="",xaxt="n",yaxt="n",main="Original data",col=mypalette,breaks=brek)
        
        ## Diplay of the estimated matrix
      ty=sort(n+1-unique(c(1,Y.hat[[k]]$RowBreaks,nrow(y))))
      tx=unique(c(1,Y.hat[[k]]$ColBreaks,ncol(y)))
      if ((length(ty)==2)||(length(tx)==2)){
        image(y = ty,
              x = tx,
              z = t(as.matrix(Y.hat[[k]]$mu.hat[nrow(Y.hat[[k]]$mu.hat):1,])),xlab="",ylab="",xaxt="n",yaxt="n", main="Estimated matrix",
              col=mypalette)
    }else{
      image(y = ty,
            x = tx,
            z = as.matrix(t(Y.hat[[k]]$mu.hat)[,nrow(Y.hat[[k]]$mu.hat):1]),xlab="",ylab="",xaxt="n",yaxt="n", main="Estimated matrix",
            col=mypalette)
    }
        abline(v=tx,col="purple")
        abline(h=ty,col="purple")
        title(main=paste("Lambda = ",round(lambda[k],3)," with (",
              length(Y.hat[[k]]$RowBreaks)-1,",", length(Y.hat[[k]]$ColBreaks)-1,
              ") breaks along the (rows,columns).", sep=""),outer=TRUE)
    }
    par(ask=FALSE)
})

