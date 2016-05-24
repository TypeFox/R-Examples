p2.est <-
function(x, y, n, j, k)
{
	xa <- x[1:j]
	ya <- y[1:j]
	jp1 <- j + 1
	xb <- x[jp1:k]
	xb2 <- x[jp1:k]^2
	yb <- y[jp1:k]
	kp1 <- k + 1
	xc <- x[kp1:n]
	yc <- y[kp1:n]
	g1 <- lm(ya ~ xa)
	g2 <- lm(yb ~ xb + xb2)
	g3 <- lm(yc ~ xc)
	beta<-c(g1$coef[1],g1$coef[2],g2$coef[1],g2$coef[2],g2$coef[3],g3$coef[1],g3$coef[2])
	s2 <- (sum((ya - g1$fit)^2) + sum((yc - g3$fit)^2))/(n - k + j)
	t2 <- sum((yb - g2$fit)^2)/(k - j)
	xj <- x[j]
	xk <- x[k]
	list(a0=beta[1],a1=beta[2],b0=beta[3],b1=beta[4],b2=beta[5],c0=beta[6],c1=beta[7], sigma2 = 
		s2, tau2 = t2, xj = xj, xk = xk)
}
