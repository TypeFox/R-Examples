ormidp.test <- function(a1, a0, b1, b0, or = 1){
  x <- matrix(c(a1,a0,b1,b0),2,2, byrow=TRUE)
  lteqtoa1 <- fisher.test(x,or=or,alternative="l")$p.val
  gteqtoa1 <- fisher.test(x,or=or,alternative="g")$p.val
  pval1 <- 0.5*(lteqtoa1-gteqtoa1+1)
  one.sided <- min(pval1, 1-pval1)
  two.sided <- 2*one.sided
  data.frame(one.sided=one.sided, two.sided=two.sided)
}
