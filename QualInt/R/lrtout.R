lrtout <- function(betaout, nsbp, alpha) {
  
  over <- betaout[, 1] / betaout[, 2]
  oversq <- over ^ 2
  cc <- min(sum(oversq * (over > 0)), sum(oversq * (over < 0)))
  P_value <- alphaf(nsbp, cc)
  c_value <- cvaluef(nsbp, alpha)
  Powerv <- powerval(betaout, nsbp, test = "LRT", alpha)
  
  list(index = c_value, pvalue = P_value, power = Powerv)
  
}