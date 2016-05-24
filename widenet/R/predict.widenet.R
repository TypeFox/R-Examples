################################################################################

predict.widenet <- function(object,
                            newx,
                            order = object$which.order.min,
                            alpha.val = object$which.alpha.min,
                            type = c("link", "response",
                              "coefficients",
                              "nonzero", "class"),
                            ...) {

  type <- match.arg(type)
  
  if(length(alpha.val) != 1)
    stop("alpha.val must be a single value")

  if(length(order) != 1)
    stop("order must be a single value")

  alpha.index <- which(object$alpha == alpha.val)

  order.index <- which(object$order == order)

  if(length(alpha.index) == 0)
    stop("The specified alpha.val is not part of this widenet object")

  if(length(order.index) == 0)
    stop("The specified order is not part of this widenet object")

  if(type %in% c("link", "response", "class")) {

    if(missing(newx)) stop("must specify newx")

    ## expand the basis if necessary

    if(object$screen.method == "none") {

      screened.in.index <- 1:ncol(newx)

    } else screened.in.index <- object$screened.in.index

    newx.to.use <- newx[, screened.in.index, drop=FALSE]

    if (order == 2 || order == 3) {

      colsBinary.screened <- object$colsBinary[screened.in.index]

      numBinary.screened <- sum(colsBinary.screened == 2)
      numNotBinary.screened <- sum(colsBinary.screened == 3)

      x2.screened <- expand.to.order.2(newx.to.use,
                                       colsBinary.screened,
                                       numBinary.screened,
                                       numNotBinary.screened)

      beta.rownames <-
        rownames(object$cv.relaxnet.results[[order.index]][[alpha.index]]$relaxnet.fit$main.glmnet.fit$beta)

      if(order == 2) {

        newx.to.use <- cbind(newx.to.use,
                             x2.screened)[, beta.rownames]

      } else  { ## order == 3

        x3.screened <- expand.to.order.3(newx.to.use,
                                         x2.screened,
                                         colsBinary.screened,
                                         numBinary.screened,
                                         numNotBinary.screened)
          
        newx.to.use <-
          cbind(newx.to.use,
                x2.screened,
                x3.screened)[, beta.rownames]
      }
    }
  }

  predict(object$cv.relaxnet.results[[order.index]][[alpha.index]],
          newx = newx.to.use,
          type = type,
          ...)
}
