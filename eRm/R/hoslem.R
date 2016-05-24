hoslem <- function(object, groups.hl = 10, pi.hat)
{
# computes the Hosmer-Lemeshow test for objects of class "ppar"
# groups.hl ... number of groups for percentile splitting

  K <- dim(object$X)[2]
  N <- dim(object$X.ex)[1]
  
  #Pi <- pmat(object)                            #expected values
  if (length(object$pers.ex) > 0) {
    y <- as.vector(t(object$X[-object$pers.ex,]))   #observed values
  } else {
    y <- as.vector(t(object$X))
  }
  pi.hat <- as.vector(t(pi.hat))

  cutpoints <- quantile(pi.hat, probs = seq(0, 1, 1/groups.hl))                     #perzentiles
  groupvec <- cut(pi.hat, cutpoints, include.lowest = TRUE, labels = 1:groups.hl)   #recode ph.hat

  
  o.g <- tapply(y, groupvec, sum)               #number of 1-responses in group
  n.g <- table(groupvec)                        #number of responses in group
  pi.mean <- tapply(pi.hat, groupvec, mean)     #average response probabilites

  value <- sum((o.g - n.g*pi.mean)^2/(n.g *pi.mean*(1-pi.mean)))    #HM-test statistic
  df <- groups.hl - 2
  p.value <- 1 - pchisq(value, df)

  result <- list(value = value, df = df, p.value = p.value)
  result
}
