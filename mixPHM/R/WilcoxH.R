# Computes the independence test following Wilcox (1997), Biometrical Journal
# It works also for non-normal distributions
# The null hypothesis is H0: R = I.

WilcoxH <- function(x, use = "pairwise.complete.obs")
{
  #x ... data set with dwell times

  p <- dim(x)[2]
  Rmat <- cor(x, use = use)

  #------- compute matrix with number of non-NA's ---------
  if (use == "pairwise.complete.obs") {
    n <- matrix(NA, p, p)                                 #n is matrix
    Rmat.ind <- cbind(rep(1:p, each = p), rep(1:p, p))
    n.i <- apply(Rmat.ind, 1, function(ind) {
                     pair <- x[,ind]
                     return(sum(!is.na(rowSums(pair))))
                   })
    n[Rmat.ind] <- n.i
  }
  if (use == "complete.obs") n <- sum(!is.na(rowSums(x))) #n integer
  if (use == "all.obs") n <- dim(x)[1]
 
  # Steiger-Hakstian test
  Tjk.full <- Rmat*sqrt((n-2)/(1-(Rmat^2)))
  Tjk <- Tjk.full[lower.tri(Tjk.full)]            #T_jk on p. 185
  SH <- sum(Tjk^2)                                #Steiger & Hakstian test
  df.SH <- length(Tjk)
  pval.SH <- 1 - pchisq(SH, df = df.SH)
  SHvec <- c(SH, df.SH, pval.SH)
  names(SHvec) <- c("SH","df","pval")

  # Wilcox Test
  eta <- (n-2)[lower.tri(n)]
  a <- eta - 0.5
  b <- 48*a^2
  cjk <- sqrt(a*log(1 + Tjk^2/eta))           
  term2 <- (cjk^3 + 3*cjk)/b
  term3 <- (4*cjk^7 + 33*cjk^5 + 240*cjk^3 + 855*cjk)/(10*b^2 + 8*b*cjk^4 + 1000*b)
  zjk <- cjk + term2 - term3
  WH <- sum(zjk^2)                               #Wilcox H
  df.WH <- df.SH
  pval.WH <- 1 - pchisq(WH, df = df.WH)
  WHvec <- c(WH, df.WH, pval.WH)
  names(WHvec) <- c("SH","df","pval")

  result <- list(Rmat = Rmat, SH.res = SHvec, WH.res = WHvec)
  class(result) <- "wilcoxh"
  result
}



