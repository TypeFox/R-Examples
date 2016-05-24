chi.plot<-function(Y1, Y2, ...){
###################
if(length(Y1) != length(Y2)) stop("Y1 and Y2 must have equal length")
n <- length(Y1)
o1 <- outer(Y1, Y1, "-")
o2 <- outer(Y2, Y2, "-")
diag(o1) <- NA
diag(o2) <- NA
t1 <- o1 > 0
t2 <- o2 > 0
t3 <- t1 & t2

#####################
H <- apply(t3, 1, function(x){sum(x, na.rm = TRUE)})/(n - 1)
Z <- apply(t1, 1, function(x){sum(x, na.rm = TRUE)})/(n - 1)
G <- apply(t2, 1, function(x){sum(x, na.rm = TRUE)})/(n - 1)
S <- sign((Z - 0.5)*(G - 0.5))
chi <- (H - Z * G)/sqrt(Z * (1 - Z) * G * (1 - G))
lambda <- 4 * S * pmax((Z - 0.5)^2, (G - 0.5)^2)
thresh <- 4 * (1/(n - 1) - 0.5)^2

####################### 
plot(lambda[abs(lambda)<thresh], chi[abs(lambda)<thresh], ylim = c(-1,1),
xlab = expression(lambda), ylab = expression(chi), ...)
abline(h = 1.78/sqrt(n), lty = 2)
abline(h = -1.78/sqrt(n), lty = 2)
invisible(cbind(lambda = lambda, chi = chi))
}

