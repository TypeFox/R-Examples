setClass("FLXMRmultinom",
         contains = "FLXMRglm")

FLXMRmultinom <- function(formula=.~., ...)
{
    z <- new("FLXMRmultinom", weighted=TRUE, formula=formula,
             family = "multinom", name=paste("FLXMRglm", "multinom", sep=":"))
    z@preproc.y <- function(x){
      x <- as.integer(factor(x))
      if (min(x) < 1 | length(unique(x)) != max(x))
        stop("x needs to be coercible to an integer vector containing all numbers from 1 to max(x)")
      y <- matrix(0, nrow = length(x), ncol = max(x))
      y[cbind(seq_along(x), x)] <- 1
      y
    }
    z@defineComponent <- expression({
      predict <- function(x) {
        p <- tcrossprod(x, coef)
        eta <- cbind(1, exp(p))
        eta/rowSums(eta)
      }
      logLik <- function(x, y) {
        log(predict(x))[cbind(seq_len(nrow(y)), max.col(y))]
      }
      
      new("FLXcomponent",
          parameters=list(coef=coef),
          logLik=logLik, predict=predict,
          df=df)
    })
    
    z@fit <- function(x, y, w, component){
      r <- ncol(x)
      p <- ncol(y)
      if (p < 2) stop("Multinom requires at least two components.")
      mask <- c(rep(0, r + 1), rep(c(0, rep(1, r)), p - 1))
      fit <- nnet::nnet.default(x, y, w, mask = mask, size = 0, 
                                skip = TRUE, softmax = TRUE, censored = FALSE, 
                                rang = 0, trace=FALSE, ...)
      fit$coefnames <- colnames(x)
      fit$weights <- w
      fit$vcoefnames <- fit$coefnames[seq_len(ncol(x))]
      fit$lab <- seq_len(ncol(y))
      class(fit) <- c("multinom", "nnet")
      coef <- coef(fit)
      with(list(coef = coef, df = length(coef)),
           eval(z@defineComponent))
    }
    z
}

setMethod("existGradient", signature(object = "FLXMRmultinom"),
          function(object) FALSE)

