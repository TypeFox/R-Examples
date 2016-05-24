.reg1fitBasic <-
function(lm.out, dname="mydata", TotSS, digits.d=NULL, show.R=FALSE) {

  nm <- all.vars(lm.out$terms)  # names of vars in the model
  n.vars <- length(nm)
  n.pred <- n.vars - 1L
  n.obs <- nrow(lm.out$model)

  tx <- character(length = 0)

  # model fit
  sm <- summary(lm.out)

  if (is.null(options()$knitr.in.progress)) {
    tx[length(tx)+1] <- "Model Fit"
    tx[length(tx)+1] <- ""
  }

  se <- sm$sigma
  tx[length(tx)+1] <- paste("Standard deviation of residuals: ", .fmt(se,digits.d),
    "for", sm$df[2], "degrees of freedom")

  tcut <- -qt(0.025, df=sm$df[2])
  range <- 2*tcut*se
  #tx[length(tx)+1] <- "If normal, the approximate 95% range of residuals about each fitted"
  #tx[length(tx)+1] <- paste("  value is 2*t-cutoff*", .fmt(se,digits.d), 
    #", with a 95% interval t-cutoff of ", .fmt(tcut,3), sep="")
  #tx[length(tx)+1] <- paste("95% range of variation: ", .fmt(range,digits.d), sep="")

  # predicted residual sum of squares
  prs.terms <- residuals(lm.out)/(1 - lm.influence(lm.out)$hat)
  PRESS <- sum(prs.terms^2)
  RsqPRESS <- 1 - (PRESS / TotSS)

  if (n.pred > 0) {
    pvl <- 1-pf(sm$fstatistic[1],sm$fstatistic[2],sm$fstatistic[3])
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste("R-squared: ", .fmt(sm$r.squared,3), 
      "   Adjusted R-squared: ", .fmt(sm$adj.r.squared,3),
      "   PRESS R-squared: ", .fmt(RsqPRESS,3))
    tx[length(tx)+1] <- paste("\n", "Null hypothesis that all population slope coefficients are 0:\n", 
        "  F-statistic: ", .fmt(sm$fstatistic[1],3),
        "     df: ", sm$fstatistic[2], " and ", sm$fstatistic[3],
        "     p-value:", .fmt(pvl, 3, 7), sep="")
  }

  return(list(tx=tx, se=sm$sigma, range=range, Rsq=sm$r.squared,
    Rsqadj=sm$adj.r.squared, PRESS=PRESS, RsqPRESS=RsqPRESS))
 
}
