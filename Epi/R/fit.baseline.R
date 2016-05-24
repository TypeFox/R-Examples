fit.baseline <-
function(y, rates.frame, start)
{
  if (missing(start)) {
    ## Get starting values from logistic regression model
    glm.out.init <- glm(y~., family=binomial, data=rates.frame)
    mu.init <- fitted(glm.out.init, type="response")
    glm(y~-1 + ., family=binomial(link=log), data=rates.frame, mustart=mu.init)
  }
  else {
    glm(y~-1 + ., family=binomial(link=log), data=rates.frame, start=start)
  }
}
