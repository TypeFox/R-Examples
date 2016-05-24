lambda.reg <- function(object, columns){

  if (class(object) != "eiReg")
    stop("'object' must be output from 'ei.reg'")
  if (missing(columns) | length(columns) < 2)
    stop("'columns' requires at least two column names")

  coefs <- matrix(NA, nrow(object$coef), length(columns))
  rownames(coefs) <- rownames(object$coef)
  colnames(coefs) <- columns

  for(i in columns){
    coefs[,i] <- object$coef[,i]/apply(object$coef[,columns],1,sum)
  }

  se <- delta(object, columns)
  lambda.out <- list(call = match.call(), lambda = coefs, se = se)
  class(lambda.out) <- c("lambdaReg", "list")
  lambda.out
}
