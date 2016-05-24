cor0.test <- function(x,y,rho0=0,alternative = c("two.sided", "less", "greater")) {
  n <- length(x)
  r <- cor(x,y)
  z <- log((1+r)/(1-r))/2
  Ez <- log((1+rho0)/(1-rho0))/2
  stat <- (z-Ez)/sqrt(1/(n-3))
  if (alternative[1] == "two.sided") pvalue <- 2*pnorm(abs(stat),lower.tail = FALSE)
  if (alternative[1] == "less") pvalue <- pnorm(stat,lower.tail = TRUE)
  if (alternative[1] == "greater") pvalue <- pnorm(stat,lower.tail = FALSE)
  
  return(list(statistic=stat,p.value=pvalue))
}

