summary.aodml <- function(object, ...) {
	## Vector b
	# function to insert NAs in se when there are NAs
  insNA <- function(b, se) {
    if(any(is.na(b))){
      nb <- length(b)
      SE <- rep(NA, nb)
      j <- 1
      for(i in seq(nb))
        if(!is.na(b[i])) {
          SE[i] <- se[j]
          j <- j + 1
        }
      se <- SE
    }
    se
  }

	# check whether any coef was set to a fixed value
  nb <- length(coef(object))
  param <- object$param
  vpar <- object$varparam
	
	# check whether any estimated var was negative and replaces these elements with NA's
  diagv <- diag(vpar)
  if(!all(is.na(diagv))){
    if(any(diagv[!is.na(diagv)] < 0)) {
      diagv[diagv < 0] <- NA
      warning("At least one variance was < 0 in the var-cov matrix. Any such element was replaced with NA.\n")
    }
  }

	# position of coef
  pos1 <- seq(nb)

	# position of coef set to a fixed value, if any
  fp <- match("fixpar", table = names(object$call))
  pos2 <- NA
  if(!is.na(fp))
    pos2 <- eval(object$call$fixpar[[2]])

	# remove coef set to a fixed value, if any
  pos3 <- setdiff(pos1, pos2)
  BCoef <- data.frame()

	# compute new var-cov mat, coef vector and position of term(s) to be tested
  if(length(pos3) > 0) {
    b3 <- param[pos3]
    v3 <- if(object$singular.hessian == 0 & !all(is.na(diagv))) diagv[pos3] else diagv

	  # coef, se, z and t test
    se3 <- sqrt(v3)
    se3 <- insNA(b3, se3)
    BCoef <- data.frame(b = b3, se = se3, z = b3 / se3,
                        P = 2 * (1 - pnorm(abs(b3) / se3)))
    nam <- names(b3)
    rownames(BCoef) <- nam
    colnames(BCoef) <- c("Estimate", "Std. Error", "z value", "Pr(> |z|)")
  }
	
	# coef set to a fixed value, if any
  pos4 <- setdiff(pos1, pos3)
  FixedBCoef <- data.frame()
  if(length(pos4) > 0)
    FixedBCoef <- data.frame(Value = param[pos4])

	## Vector phi	
	
	# position of scoef
  pos1 <- (nb + 1):length(param)

	# position of coef set to a fixed value, if any
  fp <- match("fixpar", table = names(object$call))
  pos2 <- NA
  if(!is.na(fp)) pos2 <- eval(object$call$fixpar[[2]])

	# remove coef set to a fixed value, if any
  pos3 <- setdiff(pos1, pos2)
  Phi <- data.frame()

	# compute new var-cov mat, coef vector and position of term(s) to be tested
  if(length(pos3) > 0){
    b3 <- param[pos3]
    if(object$singular.hessian == 0 & !all(is.na(diagv))){
      va3 <- diagv[pos3]
      va3[va3 < 0] <- NA
      se3 <- sqrt(va3)
      se3 <- insNA(b3, se3)
      }
    else
      se3 <- rep(NA, length(b3))

    Phi <- data.frame(b  = b3, se = se3)
  	nam <- names(b3)
  	rownames(Phi) <- nam
    colnames(Phi) <- c("Estimate", "Std. Error") #, "z value", "Pr(> z)"
    }

	# print coef set to a fixed value, if any
  pos4 <- setdiff(pos1, pos3)
  FixedPhi <- data.frame()
  if(length(pos4) > 0)
    FixedPhi <- data.frame(Value = param[pos4])
	
	structure(
		list(object = object, BCoef = BCoef, FixedBCoef = FixedBCoef,
         Phi = Phi, FixedPhi = FixedPhi),
		class = "summary.aodml")
	}
