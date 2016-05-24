
setClass("unmarkedLinComb",
         representation(parentEstimate = "unmarkedEstimate",
                        estimate = "numeric",
                        covMat = "matrix",
                        covMatBS = "optionalMatrix",
                        coefficients = "matrix"))

setClass("unmarkedBackTrans",
         representation(parentLinComb = "unmarkedLinComb",
                        estimate = "numeric",
                        covMat = "matrix",
                        covMatBS = "optionalMatrix"))

setClassUnion("linCombOrBackTrans",
              c("unmarkedLinComb", "unmarkedBackTrans"))

setMethod("show",
          signature(object = "unmarkedLinComb"),
          function(object) {

    cat("Linear combination(s) of", object@parentEstimate@name, "estimate(s)\n\n")

    df.coefs <- as.data.frame(object@coefficients)
    colnames(df.coefs) <- names(object@parentEstimate@estimates)
    lcTable <- data.frame(Estimate = object@estimate, SE = SE(object))
    lcTable <- cbind(lcTable, df.coefs)
    if(nrow(lcTable) > 1) {
        print(lcTable, digits = 3, row.names = 1:nrow(lcTable))
    } else {
        print(lcTable, digits = 3, row.names = FALSE)
    }
    cat("\n")

})


setMethod("show",
          signature(object = "unmarkedBackTrans"),
          function(object) {
    cat("Backtransformed linear combination(s) of",object@parentLinComb@parentEstimate@name,"estimate(s)\n\n")

    lcTable <- data.frame(LinComb = object@parentLinComb@estimate)

    df.coefs <- as.data.frame(object@parentLinComb@coefficients)
    colnames(df.coefs) <-
        names(object@parentLinComb@parentEstimate@estimates)

    lcTable <- cbind(lcTable, df.coefs)

    btTable <- data.frame(Estimate = object@estimate, SE = SE(object))

    btTable <- cbind(btTable, lcTable)

    if(nrow(btTable) > 1) {
        print(btTable, digits = 3, row.names = 1:nrow(btTable))
    } else {
        print(btTable, digits = 3, row.names = FALSE)
    }

    cat("\nTransformation:", object@parentLinComb@parentEstimate@invlink,"\n")
})


setMethod("backTransform",
          signature(obj = "unmarkedLinComb"),
          function(obj) {

    ## In general, MV delta method is Var=J*Sigma*J^T where J is Jacobian
    ## In this case, J is diagonal with elements = gradient
    ## Reduces to scaling the rows then columns of Sigma by the gradient
    e <- do.call(obj@parentEstimate@invlink,list(obj@estimate))
    grad <- do.call(obj@parentEstimate@invlinkGrad,list(obj@estimate))

    if(length(obj@estimate) > 1) {
        v <- diag(grad) %*% obj@covMat %*% diag(grad)
    } else {
        v <- grad^2 * obj@covMat
    }

    if (!is.null(obj@covMatBS)) {
        if(length(obj@estimate) > 1) {
            v.bs <- diag(grad) %*% obj@covMatBS %*% diag(grad)
        } else {
            v.bs <- grad^2 * obj@covMatBS
        }
    } else {
        v.bs <- NULL
    }

    umbt <- new("unmarkedBackTrans", parentLinComb = obj,
                estimate = e, covMat = v, covMatBS = v.bs)
    umbt
})


setMethod("SE", "linCombOrBackTrans",
          function(obj, ...) {
              v <- vcov(obj, ...)
              sqrt(diag(v))
          })

setMethod("coef", "linCombOrBackTrans",
          function(object) {
              object@estimate
          })

setMethod("vcov", "linCombOrBackTrans",
          function(object, method="hessian") {
              method <- match.arg(method, c("hessian", "nonparboot"))
              switch(method,
                     hessian = return(object@covMat),
                     nonparboot = return(object@covMatBS))
          })

setMethod("confint", "unmarkedLinComb",
          function(object, parm, level = 0.95) {
              if(missing(parm)) parm <- 1:length(object@estimate)
              ests <- object@estimate[parm]
              ses <- SE(object)[parm]
              z <- qnorm((1-level)/2, lower.tail = FALSE)
              lower.lim <- ests - z*ses
              upper.lim <- ests + z*ses
              ci <- as.matrix(cbind(lower.lim, upper.lim))
              colnames(ci) <- c((1-level)/2, 1- (1-level)/2)
              rownames(ci) <- 1:nrow(ci)
              if(nrow(ci) == 1) rownames(ci) <- ""
              ci
          })

setMethod("confint", "unmarkedBackTrans",
          function(object, parm, level = 0.95) {
              if(missing(parm)) parm <- 1:length(object@estimate)
              ci.orig <- callGeneric(object@parentLinComb, parm, level)
              invlink <- object@parentLinComb@parentEstimate@invlink
              lower.lim <- do.call(invlink, list(ci.orig[,1]))
              upper.lim <- do.call(invlink, list(ci.orig[,2]))
              ci <- as.matrix(cbind(lower.lim, upper.lim))
              colnames(ci) <- c((1-level)/2, 1- (1-level)/2)
              rownames(ci) <- 1:nrow(ci)
              if(nrow(ci) == 1) rownames(ci) <- ""
              ci
          })

setAs("unmarkedBackTrans", "data.frame",
      def=function(from) {
          lcTable <- data.frame(LinComb = from@parentLinComb@estimate)
          df.coefs <- as.data.frame(from@parentLinComb@coefficients)
          colnames(df.coefs) <-
              names(from@parentLinComb@parentEstimate@estimates)
          lcTable <- cbind(lcTable, df.coefs)
          btTable <- data.frame(Estimate = from@estimate, SE = SE(from))
          btTable <- cbind(btTable, lcTable)
          btTable
      })

setAs("unmarkedLinComb", "data.frame",
      def=function(from) {
          df.coefs <- as.data.frame(from@coefficients)
          colnames(df.coefs) <- names(from@parentEstimate@estimates)
          lcTable <- data.frame(Estimate = from@estimate, SE = SE(from))
          lcTable <- cbind(lcTable, df.coefs)
          lcTable
      })

