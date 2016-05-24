twosample.cor.test <- function(x1,y1,x2,y2,alpha=0.05,alternative = c("two.sided", "less", "greater")) {
  n1 <- length(x1)
  n2 <- length(x2)
  r1 <- cor(x1,y1)
  r2 <- cor(x2,y2)
  z1 <- log((1+r1)/(1-r1))/2
  z2 <- log((1+r2)/(1-r2))/2
  stat <- (z1-z2)/sqrt(1/(n1-3)+1/(n2-3))
  if (alternative[1] == "two.sided") pvalue <- 2*pnorm(abs(stat),lower.tail = FALSE)
  if (alternative[1] == "less") pvalue <- pnorm(stat,lower.tail = TRUE)
  if (alternative[1] == "greater") pvalue <- pnorm(stat,lower.tail = FALSE)
  
  return(list(statistic=stat,p.value=pvalue))
}
