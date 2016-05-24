if(!isGeneric("summary"))
  setGeneric(name = "summary", def = function(object, ...) standardGeneric("summary"))

## summary method for glimML objects
setMethod("summary", signature = "glimML",
  function(object){
# function to insert NAs in se when there are NAs in b
    insNA <- function(b, se){
      if(any(is.na(b))){
        nb <- length(b)
        SE <- rep(NA, nb)
        j <- 1
        for(i in seq(nb))
          if(!is.na(b[i])){
            SE[i] <- se[j]
            j <- j + 1
            }
        se <- SE
        }
      se
      }

# check whether any fixed-effect coef was set to a fixed value
    nb <- length(coef(object))
    param <- object@param
    vpar <- object@varparam
# Checks whether any estimated variance was negative and replaces these elements with NA's
    diagv <- diag(vpar)
    if(!all(is.na(diagv))){
      if(any(diagv[!is.na(diagv)] < 0)){
        diagv[diagv < 0] <- NA
        warning("At least one variance was < 0 in the var-cov matrix. Any such element was replaced with NA.\n")
        }
      }

# position of fixed-effect coefficients
    pos1 <- seq(nb)

# position of parameters set to a fixed value, if any
    fp <- match("fixpar", table = names(object@CALL))
    pos2 <- NA
    if(!is.na(fp))
      pos2 <- eval(object@CALL$fixpar[[2]])

# remove fixed parameters, if any
    pos3 <- setdiff(pos1, pos2)
    Coef <- data.frame()

# compute new var-cov mat, coef vector and position of term(s) to be tested
    if(length(pos3) > 0){
      b3 <- param[pos3]
      v3 <- if(object@singular.hessian == 0 & !all(is.na(diagv))) diagv[pos3] else diagv

# coef, se, z and t test
      se3 <- sqrt(v3)
      se3 <- insNA(b3, se3)
      Coef <- data.frame(b = b3, se = se3, z = b3 / se3, P = 2 * (1 - pnorm(abs(b3) / se3)))
      nam <- names(b3)
      rownames(Coef) <- nam
      colnames(Coef) <- c("Estimate", "Std. Error", "z value", "Pr(> |z|)")
      }
# fixed-effect coefficients which were set to a fixed value, if any
    pos4 <- setdiff(pos1, pos3)
    FixedCoef <- data.frame()
    if(length(pos4) > 0){
      FixedCoef <- data.frame(Value = param[pos4])
      }

# position of overdispersion parameters
    pos1 <- (nb + 1):length(param)

# position of parameters set to a fixed value, if any
    fp <- match("fixpar", table = names(object@CALL))
    pos2 <- NA
    if(!is.na(fp))
      pos2 <- eval(object@CALL$fixpar[[2]])

# remove fixed parameters, if any
    pos3 <- setdiff(pos1, pos2)
    Phi <- data.frame()

# compute new var-cov mat, coef vector and position of term(s) to be tested
    if(length(pos3) > 0){
      b3 <- param[pos3]
      if(object@singular.hessian == 0 & !all(is.na(diagv))){
        va3 <- diagv[pos3]
        va3[va3 < 0] <- NA
        se3 <- sqrt(va3)
        se3 <- insNA(b3, se3)
        }
      else
        se3 <- rep(NA, length(b3))
# coef, se, z and t test for phi
# beware: unilateral test for phi because it cannot be negative
      if(any(b3 < 0))
        warning("Negative values for phi.")

## Modif R Lancelot 26/08/2008 suite à remarque de F Bonnot et proposition de M Lesnoff
      Phi <- data.frame(b  = b3,
                        se = se3,
                        z = ifelse(se3 <= 2e-13, 0, b3 / se3),
                        P = ifelse(se3 <= 2e-13, 1, 1 - pnorm(abs(b3) / se3)))
## fin modif

      nam <- names(b3)
      rownames(Phi) <- nam
      colnames(Phi) <- c("Estimate", "Std. Error", "z value", "Pr(> z)")
      }

# print random coefficients which were set to a fixed value, if any
    pos4 <- setdiff(pos1, pos3)
    FixedPhi <- data.frame()
    if(length(pos4) > 0){
      FixedPhi <- data.frame(Value = param[pos4])
      }
    res <- new(Class = "summary.glimML",
               object = object, Coef = Coef, FixedCoef = FixedCoef, Phi = Phi, FixedPhi = FixedPhi)
    res
    })


### show method for objects of class summary.glimML
setMethod("show", signature = "summary.glimML",
  function(object){
  Object <- object@object

# Title and call
    switch(Object@method,
      BB = cat("Beta-binomial model\n", "-------------------\n", sep = ""),
      NB = cat("Negative-binomial model\n", "-----------------------\n", sep = ""))
    print(Object@CALL)

# Checks whether convergence problems occurred
    n <- Object@code
    iter <- Object@iterations
    msg <- Object@msg
    if(n == 0)
      cat("\nConvergence was obtained after " , iter, " iterations.\n", sep = "")

# Print estimated fixed effects, if any
    Coef <- object@Coef
    if(nrow(Coef) > 0){
      nam <- rownames(Coef)
      List <- vector(mode = "list", length = 4)
      for(i in 1:4){
        x <- Coef[,i]
        List[[i]] <- format(x, digits = 4, scientific = TRUE)
        }
      Coeftext <- as.data.frame(t(do.call("rbind", List)))
      rownames(Coeftext) <- nam
      colnames(Coeftext) <- c("Estimate", "Std. Error", "z value", "Pr(> |z|)")
      cat("\nFixed-effect coefficients:\n")
      print(Coeftext)
      }

# print fixed-effect coefficients which were set to a fixed value, if any
    FixedCoef <- object@FixedCoef
    if(nrow(FixedCoef) > 0){
      cat("\nFixed-effect coefficients set to fixed values:\n")
      print(FixedCoef)
      }

# overdispersion coefficients phi
# compute new var-cov mat, coef vector and position of term(s) to be tested
    Phi <- object@Phi
    if(nrow(Phi) > 0){
      nam <- rownames(Phi)
      List <- vector(mode = "list", length = 4)
#      for(i in 1:4){
#        x <- Phi[,i]
#        List[[i]] <- if(i < 4)
#                       format(x, scientific = TRUE)
#                     else
#                       ifelse(x < 1e-4, "< 1e-4", format(x, scientific = TRUE))
      for(i in 1:4){
        x <- Phi[,i]
        List[[i]] <- format(x, digits = 4, scientific = TRUE)
        }
      Phitext <- as.data.frame(t(do.call("rbind", List)))
      rownames(Phitext) <- nam
      colnames(Phitext) <- c("Estimate", "Std. Error", "z value", "Pr(> z)")
      cat("\nOverdispersion coefficients:\n")
      print(Phitext)
      }

# print overdispersion coefficients which were set to a fixed value, if any
    FixedPhi <- object@FixedPhi
    if(nrow(FixedPhi) > 0){
      cat("\nOverdispersion coefficients set to fixed values:\n")
      print(FixedPhi)
      }
    akic <- AIC(Object)@istats; aic <- akic[,2]; aicc <- akic[,3]

    ll    <- format(Object@logL, digits = 4, scientific = TRUE)
    nbpar <- format(Object@nbpar)
    dfres <- format(df.residual(Object))
    dev   <- format(deviance(Object), digits = 4, scientific = TRUE)
    aic   <- format(akic[,2], digits = 4, scientific = TRUE)
    aicc  <- format(akic[,3], digits = 4, scientific = TRUE)
    res   <- c(ll, nbpar, dfres, dev, aic, aicc)
    names(res) <- c("Log-lik", "nbpar", "df res.", "Deviance", "AIC", "AICc")
    cat("\nLog-likelihood statistics\n")
    print(res, quote = FALSE)
    invisible(object)
    })


### show method for glimML objects
setMethod("show", signature = "glimML",  function(object) show(summary(object)))
