fit.add <-
function(y, rates.frame, cov.frame, start)
{
  ## Modify covariate values
  rates.sum <- apply(rates.frame, 1, sum)
  cov.frame <- sweep(cov.frame, 1, rates.sum, "*")
  model.frame <- cbind(rates.frame, cov.frame)
  
  if (missing(start)) {
    glm.out <- fit.baseline(y, rates.frame)
    mu.inits <- predict(glm.out, type="response")
    glm.out <- glm(y~-1 + ., family=binomial(link=log),
                   data=model.frame, mustart=mu.inits, maxit=100)
  }
  else {
    glm.out <- glm(y~-1 + ., family=binomial(link=log),
                   data=model.frame, start=start, maxit=100)
  }
  return(list(rates=glm.out))
}
