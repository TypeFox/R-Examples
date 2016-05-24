methods::setMethod(f = "predict",
                   signature = c(object="ramsvm"),
                   definition = function(object, 
                                         newdata = NULL, 
                                         lambda = NULL, ...) {

  if( is(newdata, "NULL") ) newdata <- object@x

  if( is(lambda, "NULL") ) lambda <- object@lambda

  if( class(newdata) == "data.frame" ) newdata <- data.matrix(newdata)
  if( class(newdata) != "matrix" ) newdata <- as.matrix(newdata)

  if( ncol(newdata) != ncol(object@x) ) {
    stop(paste("The number of columns in new ",
               "covariates matrix/data.frame is incorrect.", sep=""),
         call. = FALSE)
  }

  if( !is(lambda, "numeric") ) {
    stop("All lambda should be numeric.", call. = FALSE)
  }

  if( any(lambda < 0.0) ) {
    stop("All lambda should be non-negative.", call. = FALSE)
  }

  pred.y <- numeric(0L)

  nol <- length(object@lambda)

  for( i in 1L:length(lambda) ) {
    temp <- lambda[i]
    index <- which(sapply(X = object@lambda, 
                          FUN = function(x){isTRUE(all.equal(x,temp))}))
    if( length(index) == 1L ) {
      temp.beta <- object@beta[[index]]
      temp.beta0 <- object@beta0[[index]]
    } else if( length(index) == 0L ) {
      if( temp > object@lambda[1L] ) {
        temp.beta <- object@beta[[1L]]
        temp.beta0 <- object@beta0[[1L]]
        cat(paste("Lambda",
                  temp,
                  "is bigger than the largest lambda in the",
                  "solution path.\nUsing the parameters",
                  "corresponding to the largest lambda in the",
                  "solution path.\n"))
      } else if( temp < object@lambda[nol] ) {
        temp.beta <- object@beta[[nol]]
        temp.beta0 <- object@beta0[[nol]]
        cat(paste("Lambda",
                  temp,
                  "is less than the smallest lambda in the",
                  "solution path.\nUsing the parameters",
                  "corresponding to the smallest lambda in the",
                  "solution path.\n"))
      } else if( temp > object@lambda[nol] && 
                 temp < object@lambda[1L] ) {
        idx <- max(which(object@lambda > temp))
        temp.beta <- object@beta[[idx]] + 
                     (temp - object@lambda[idx]) / 
                     (object@lambda[(idx+1L)] - object@lambda[idx]) * 
                     (object@beta[[(idx+1L)]] - object@beta[[idx]])
        temp.beta0 <- object@beta0[[idx]] + 
                      (temp - object@lambda[idx]) / 
                      (object@lambda[(idx+1L)] - object@lambda[idx]) * 
                      (object@beta0[[(idx+1L)]] - object@beta0[[idx]])
      }
    }

    temp.pred.y <- Pred(object = object,
                        x = newdata,
                        beta = temp.beta,
                        beta0 = temp.beta0)

    temp.pred.y2 <- character(nrow(newdata))
    for( i in 1L:object@k ) {
      temp.pred.y2[temp.pred.y == i] <- object@y.name[i]
    }

    if( is(object@y, "numeric") ) {
      temp.pred.y2 <- as.numeric(temp.pred.y2)
    }
  
    pred.y <- c(pred.y, list(temp.pred.y2))
  }

  names(pred.y) <- lambda

  return(pred.y)
})

