# Constrainted Method LQL Model
# Different variances for all three segments 


con2.search.D <- function(x, y, n, jlo, jhi, klo, khi,plot)
{
	fjk <- matrix(0, n, n)
	fxy <- matrix(0, (jhi - jlo + 1), (khi - klo + 1))
	
 ##Yulei's edit to avoid using for-loop
	jkgrid <- expand.grid(jlo:jhi, klo:khi)
	res <- data.frame(j = jkgrid[,1],
	                  k = jkgrid[,2],
	                  k.ll = apply(jkgrid, 1, con2.parmsFUN.D, x = x, 
	                               y = y, n = n))
	
	res.m <- matrix(res$k.ll, nrow = jhi-jlo+1, ncol = khi-klo+1)
	rownames(res.m) <- jlo:jhi
	colnames(res.m) <- klo:khi
  
	fxy <- res.m
 if (plot == "TRUE") {
   jx<-jlo:jhi
   ky<-klo:khi
   persp(jx, ky, fxy, xlab = "j", ylab = "k", zlab = "LL(x,y,j,k)")
   title("Log-likelihood Surface")
 }
	z <- findmax(fxy)
	jcrit <- z$imax + jlo - 1
	kcrit <- z$jmax + klo - 1
	list(jhat = jcrit, khat = kcrit, value = max(fxy))
}

con2.parmsFUN.D<-
  function(jk,x,y,n){
    j = jk[1]
    k = jk[2]
    a <- con2.parms.D(x,y,n,j,k,1,1,1)
nr <- nrow(a$theta)
est <- a$theta[nr,  ]
b <- con2.est.D(x[j], x[k], est)
s2 <- 1/b$eta1
t2 <- 1/b$eta2
u2 <- 1/b$eta3
return(p.ll.D(n, j, k, s2, t2, u2))
}

con2.parms.D <-
function(x,y,n,j0,k0,e10,e20,e30){

th <- matrix(0,100,8)

# Iteration 0

 th[1,1] <- e10
 th[1,2] <- e20
 th[1,3] <- e30
 bc <- beta2.calc.D(x,y,n,j0,k0,e10,e20,e30)
 th[1,4:8] <- bc$B

# Iterate to convergence (100 Iter max)

 for (iter in 2:100){
 m <- iter-1
 ec <- eta2.calc.D(x,y,n,j0,k0,th[m,4:8])
 th[iter,1] <- ec$eta1
 th[iter,2] <- ec$eta2
 th[iter,3] <- ec$eta3
 bc <- beta2.calc.D(x,y,n,j0,k0,ec$eta1,ec$eta2,ec$eta3)
 th[iter,4:8] <- bc$B
 theta <- th[1:iter,]
 #delta <- abs(th[iter,]-th[m,])
 delta <- abs(th[iter,]-th[m,])/th[m,]
 if( (delta[1]<.001) & (delta[2]<.001) & (delta[3]<.001)
   & (delta[4]<.001) & (delta[5]<.001) & (delta[6]<.001)
   & (delta[7]<.001) & (delta[8]<.001))
 break
 }
 list(theta=theta)
}

con2.est.D <-
function(xj, xk, est)
{
	eta1 <- est[1]
	eta2 <- est[2]
  eta3 <- est[3]
	a0 <- est[4]
	a1 <- est[5]
	b1 <- est[6]
	b2 <- est[7]
	c1 <- est[8]
	b0 <- a0 + (a1 - b1) * xj - b2 * xj^2
	c0 <- b0 + (b1 - c1) * xk + b2 * xk^2
	list(eta1 = eta1, eta2 = eta2, eta3 = eta3, a0 = a0, a1 = a1, b0 = b0, b1 = b1, b2 = b2, c0 = c0, c1
		 = c1)
}


con2.vals.D <-
function(x, y, n, j, k)
{
	a <- con.parms.D(x, y, n, j, k, 1, 1, 1)
	nr <- nrow(a$theta)
	est <- a$theta[nr,  ]
	b <- con.est.D(x[j], x[k], est)
	eta <- c(b$eta1, b$eta2, b$eta3)
  beta <- c(b$a0, b$a1, b$b0, b$b1, b$b2, b$c0, b$c1)
	tau <- c(x[j], x[k])
	list(eta = eta, beta = beta, tau = tau)
}


p.ll.D<-function(n, j, k, s2, t2, u2){
 q1 <- n * log(sqrt(2 * pi))
 q2 <- 0.5 * (j) * (1 + log(s2))
 q3 <- 0.5 * (k - j) * (1 + log(t2))
## Yulei's edit
 q4 <- 0.5 * (n - k) * (1 + log(u2))
 - (q1 + q2 + q3 + q4)
}

## Yulei's edit to avoid using for-loop
findmax <-function(a){
	maxa<-max(a)
  imax<- which(a==max(a),arr.ind=TRUE)[1]
  jmax<-which(a==max(a),arr.ind=TRUE)[2]
	list(imax = imax, jmax = jmax, value = maxa)
}

beta2.calc.D <-
function(x, y, n, j, k, e1, e2, e3)
{
	aa <- wmat2.D(x, y, n, j, k, e1, e2, e3)
	W <- aa$w
	bb <- rvec2.D(x, y, n, j, k, e1, e2, e3)
	R <- bb$r
	beta <- solve(W, R)
	list(B = beta)
}

eta2.calc.D <-
function(x, y, n, j, k, theta)
{
	jp1 <- j + 1
	kp1 <- k + 1
	a0 <- theta[1]
	a1 <- theta[2]
	b1 <- theta[3]
  b2 <- theta[4]
	c1 <- theta[5]
	b0 <- a0 + (a1 - b1) * x[j] - b2 * x[j]^2
	c0 <- b0 + (b1 - c1) * x[k] + b2 * x[k]^2
	rss1 <- sum((y[1:j] - a0 - a1 * x[1:j])^2)
	rss2 <- sum((y[jp1:k] - b0 - b1 * x[jp1:k] - b2 * x[jp1:k]^2)^2)
	rss3 <- sum((y[kp1:n] - c0 - c1 * x[kp1:n])^2)
	e1 <- j/rss1
	e2 <- (k - j)/rss2
  e3 <- (n - k)/rss3
	list(eta1 = e1, eta2 = e2, eta3 = e3)
}
wmat2.D<-
function(x, y, n, j, k, e1, e2, e3)
{
	W <- matrix(0, 5, 5)
	jp1 <- j + 1
	kp1 <- k + 1
	W[1, 1] <- e1 * j + e2 * (k - j) + e3 * (n - k)
	W[1, 2] <- e1 * sum(x[1:j]) + e3 * (n - k) * x[j] + e2 * (k - j) * x[j]
	W[1, 3] <- e3 * (n - k) * (x[k] - x[j]) + e2 * sum(x[jp1:k] - x[j])
  W[1, 4] <- e3 * (n - k) * (x[k]^2 - x[j]^2) + e2 * sum(x[jp1:k]^2 - x[j]^2)
	W[1, 5] <- e3 * sum(x[kp1:n] - x[k])
	
  W[2, 2] <- e1 * sum(x[1:j] * x[1:j]) + e3 * (n - k) * x[j] * x[j] + e2 * (k - j) *
 		x[j] * x[j]
	W[2, 3] <- e3 * (n - k) * x[j] * (x[k] - x[j]) + e2 * x[j] * sum(x[jp1:k] - x[
		j])
  W[2, 4] <- e3 * (n - k) * x[j] * (x[k]^2 - x[j]^2) + e2 * x[j] * sum(x[jp1:k]^2 - x[
		j]^2)
	W[2, 5] <- e3 * x[j] * sum(x[kp1:n] - x[k])
	
W[3, 3] <- e3 * (n - k) * (x[k] - x[j]) * (x[k] - x[j]) + e2 * sum((x[jp1:k] - 
		x[j]) * (x[jp1:k] - x[j]))
W[3, 4] <- e3 * (n - k) * (x[k] - x[j]) * (x[k]^2 - x[j]^2) + e2 * sum((x[jp1:k] - 
		x[j]) * (x[jp1:k]^2 - x[j]^2))
	W[3, 5] <- e3 * (x[k] - x[j]) * sum(x[kp1:n] - x[k])
       
W[4, 4] <- e3 * (n - k) * (x[k]^2 - x[j]^2) * (x[k]^2 - x[j]^2) + e2 * sum((x[jp1:k]^2 - 
		x[j]^2) * (x[jp1:k]^2 - x[j]^2))
        W[4, 5] <- e3 * (x[k]^2 - x[j]^2) * sum(x[kp1:n] - x[k])
	
W[5, 5] <- e3 * sum((x[kp1:n] - x[k]) * (x[kp1:n] - x[k]))
	
W[2, 1] <- W[1, 2]
	W[3, 1] <- W[1, 3]
	W[4, 1] <- W[1, 4]
        W[5, 1] <- W[1, 5]
	W[3, 2] <- W[2, 3]
	W[4, 2] <- W[2, 4]
        W[5, 2] <- W[2, 5]
	W[4, 3] <- W[3, 4]
        W[5, 3] <- W[3, 5]
        W[5, 4] <- W[4, 5]
	list(w = W)
}
rvec2.D<-
function(x, y, n, j, k, e1, e2, e3)
{
	R <- array(0, 5)
	jp1 <- j + 1
	kp1 <- k + 1
	y1j <- sum(y[1:j])
	yjk <- sum(y[jp1:k])
	ykn <- sum(y[kp1:n])
	xy1j <- sum(x[1:j] * y[1:j])
	xyjk <- sum(x[jp1:k] * y[jp1:k])
  x2yjk <- sum(x[jp1:k]^2 * y[jp1:k])
	xykn <- sum(x[kp1:n] * y[kp1:n])
	R[1] <- e1 * y1j + e3 * ykn + e2 * yjk
	R[2] <- e1 * xy1j + e3 * x[j] * ykn + e2 * x[j] * yjk
	R[3] <- e3 * (x[k] - x[j]) * ykn + e2 * (xyjk - x[j] * yjk)
  R[4] <- e3 * (x[k]^2 - x[j]^2) * ykn + e2 * (x2yjk - x[j]^2 * yjk) 
	R[5] <- e3 * (xykn - x[k] * ykn)
	list(r = R)
}

p2.est.D <- function(x, y, n, j, k)
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
  s2 <- sum((ya - g1$fit)^2)/j
  t2 <- sum((yb - g2$fit)^2)/(k - j)
  u2<-sum((yc - g3$fit)^2)/(n - k)
  xj <- x[j]
  xk <- x[k] 
  list(a0=beta[1],a1=beta[2],b0=beta[3],b1=beta[4],b2=beta[5],c0=beta[6],c1=beta[7], sigma2 = 
                s2, tau2 = t2, u2 = u2, xj = xj, xk = xk)
}


