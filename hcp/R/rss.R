rss <- function(x, y, n, j0, k0)
{
	xa <- x[1:j0]
	ya <- y[1:j0]
	jp1 <- j0 + 1
	xb <- x[jp1:k0]
	yb <- y[jp1:k0]
	kp1 <- k0 + 1
	xc <- x[kp1:n]
	yc <- y[kp1:n]
	g1 <- lm(ya ~ xa)
	g2 <- lm(yb ~ xb)
	g3 <- lm(yc ~ xc)
	beta<-c(g1$coef[1],g1$coef[2],g2$coef[1],g2$coef[2],g3$coef[1],g3$coef[2])
	s1 <- sum((g1$res)^2)
	s2 <- sum((g2$res)^2)
	s3 <- sum((g3$res)^2)
	list(a0=beta[1],a1=beta[2],b0=beta[3],b1=beta[4],c0=beta[5],c1=beta[6],S1 = s1, S2 = s2, 
		S3 = s3)
}
