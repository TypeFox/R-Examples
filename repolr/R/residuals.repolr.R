residuals.repolr <-
function(object, type = c("deviance", "pearson", "response"), ...){
  type <- match.arg(type)
  y <- object$y
  mu <- object$fitted.values
  res <- switch(type, response = y - mu, pearson = (y - mu)/sqrt(mu * (1 - mu)),
         deviance = ifelse(y == 0, - sqrt(2 * abs(log(1 - mu))), 
                    sqrt(2 * abs(log(mu)))))
  res
}
