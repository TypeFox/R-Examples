Score.Test.Overall <-
function(Y, K1, K2, K3)
{
      n <- length(Y)
b <- mean(Y); sig2 <- var(Y)
      U <- t(Y-b)%*%(K1+K2+K3)%*%(Y-b)/(2*sig2) 
M <- (K1+K2+K3)
e <- TRACE(PROJECT(M))/2
c11 <- TT(PROJECT(K1), PROJECT(K1))
c12 <- TT(PROJECT(K1), PROJECT(K2))
c13 <- TT(PROJECT(K1), PROJECT(K3))
c22 <- TT(PROJECT(K2), PROJECT(K2))
c23 <- TT(PROJECT(K2), PROJECT(K3))
c33 <- TT(PROJECT(K3), PROJECT(K3))
COV <- matrix(c(c11, c12, c13, c12, c22, c23, c13, c23, c33), 3, 3)
Its <- c(TRACE(PROJECT(K1)), TRACE(PROJECT(K2)), TRACE(PROJECT(K3)))
correct_COV <- (COV-Its%*%t(Its)/(n-1))/2
Itt <- sum(correct_COV)
     k <- Itt/(2*e); v <- 2*e^2/Itt
      pvalue <- pchisq(U/k, df=v, lower.tail=FALSE)  
      object <- list(Score=U, p.value=pvalue, df=v, scale=k)
      class(object) <- "Score Test: tau1=tau2=tau3=0"
      return(object)
}
