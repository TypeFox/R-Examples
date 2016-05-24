#
### Simple utility functions for building linear models.
#
lmforwardsequentialAICc <- function (y, x, object) {
  AuxTrait <- !missing(x)
  included <- numeric(0)
  candidates <- 1L:ncol(object$u)
  df1 <- if(AuxTrait) {x <- cbind(x) ; data.frame(x,object)} else object
  p1 <- if(AuxTrait) paste("y~",paste(colnames(x),collapse="+"),sep="") else "y~"
  while (TRUE) {
    p2 <- paste(if(length(included)) paste(if(AuxTrait)"+",paste(paste("V_", included, sep = ""), collapse = "+")) else if(AuxTrait) "" else "1" ,sep="")
    lm1 <- lm(as.formula(paste(p1,p2,sep="")), data = df1)
    k1 <- length(lm1$coef)
    AICc1 <- AIC(lm1) + (2 * k1 * (k1 + 1)/(length(y) - k1 - 1))
    AICc2 <- rep(NA, ncol(object$u))
    for (i in candidates) {
      lm2 <- lm(as.formula(paste(p1,p2,"+ V_",i,sep = "")), data = df1)
      k2 <- length(lm2$coef)
      AICc2[i] <- AIC(lm2) + (2 * k2 * (k2 + 1)/(length(y) - k2 - 1))
    }
    if (min(AICc2, na.rm = TRUE) < AICc1) {
      included <- c(included,candidates[candidates==which.min(AICc2)])
      candidates <- candidates[candidates!=which.min(AICc2)]
    } else {
      lm1$AICc <- AICc1
      return(lm1)
    }
  }
}
#
lmforwardsequentialsidak <- function (y, x, object, alpha = 0.05) {
  AuxTrait <- !missing(x)## AuxTrait <- TRUE ## AuxTrait <- FALSE
  included <- numeric(0)
  candidates <- 1L:ncol(object$u)
  df1 <- if(AuxTrait) {x <- cbind(x) ; data.frame(x,object)} else object
  p1 <- if(AuxTrait) paste("y~",paste(colnames(x),collapse="+"),sep="") else "y~"
  while (TRUE) {
    p2 <- paste(if(length(included)) paste(if(AuxTrait)"+",paste(paste("V_", included, sep = ""), collapse = "+")) else if(AuxTrait) "" else "1" ,sep="")
    pval <- rep(NA, ncol(object$u))
    lm1 <- lm(as.formula(paste(p1,p2,sep="")), data = df1)
    for (i in candidates) {
      # i <- candidates[1L]
      lm2 <- lm(as.formula(paste(p1,p2,"+ V_",i,sep = "")), data = df1)
      aovcomp <- anova(lm1, lm2)
      pval[i] <- 1 - (1 - aovcomp[["Pr(>F)"]][2L])^(length(candidates) - length(included))
    }
    if (min(pval, na.rm = TRUE) < alpha) {
      included <- c(included, candidates[candidates == which.min(pval)])
      candidates <- candidates[candidates != which.min(pval)]
    } else {
      aovlm1 <- anova(lm1)
      lm1[["Familiwise"]] <- aovlm1[["Pr(>F)"]]
      idx <- match(colnames(object$u)[included],rownames(aovlm1))
      lm1[["Familiwise"]][idx] <- 1-(1-aovlm1[["Pr(>F)"]][idx])^(ncol(object$u):(ncol(object$u)-length(included)+1L))
      return(lm1)
    }
  }
}
#
