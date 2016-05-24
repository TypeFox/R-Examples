# Constraint Method LLL Model
# Different variances for all three segments


con.search.D <- function(x, y, n, jlo, jhi, klo, khi,plot)
{
	fjk <- matrix(0, n, n)
	fxy <- matrix(0, (jhi - jlo + 1), (khi - klo + 1))	
## Yulei's edit to avoid using for-loop
	jkgrid <- expand.grid(jlo:jhi, klo:khi)
	res <- data.frame(j = jkgrid[,1],
	                  k = jkgrid[,2],
	                  k.ll = apply(jkgrid, 1, con.parmsFUN.D, x = x, 
	                               y = y, n = n))
	
	fxy <- matrix(res$k.ll, nrow = jhi-jlo+1, ncol = khi-klo+1)
	rownames(fxy) <- jlo:jhi
	colnames(fxy) <- klo:khi
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

con.parmsFUN.D <- function(jk, x, y, n){
  j = jk[1]
  k = jk[2]
  a <- con.parms.D(x,y,n,j,k,1,1,1)
  nr <- nrow(a$theta)
  est <- a$theta[nr,  ]
  b<-con.est.D(x[j],x[k],est)
  s2<-1/b$eta1
  t2<-1/b$eta2
  u2<-1/b$eta3
  return(p.ll.D(n, j, k, s2, t2, u2))
}

## Yulei's edit: e30, e3, and eta3 is the eta used for the third segment
con.parms.D <-
function(x,y,n,j0,k0,e10,e20,e30){

th <- matrix(0,100,7)

# Iteration 0

 th[1,1] <- e10
 th[1,2] <- e20
 th[1,3] <- e30
 bc <- beta.calc.D(x,y,n,j0,k0,e10,e20,e30)
 
th[1,4:7] <- bc$B

# Iterate to convergence (100 Iter max)

 for (iter in 2:100){
 m <- iter-1

 ec <- eta.calc.D(x,y,n,j0,k0,th[m,4:7])
 th[iter,1] <- ec$eta1
 th[iter,2] <- ec$eta2
 th[iter,3] <- ec$eta3

 bc <- beta.calc.D(x,y,n,j0,k0,ec$eta1,ec$eta2,ec$eta3)

 th[iter,4:7] <- bc$B
 theta <- th[1:iter,]
 #delta <- abs(th[iter,]-th[m,])
 delta <- abs(th[iter,]-th[m,])/th[m,]
 if( (delta[1]<.001) & (delta[2]<.001) & (delta[3]<.001)
   & (delta[4]<.001) & (delta[5]<.001) & (delta[6]<.001) &(delta[7]<.001))
 break
 }
list(theta=theta)
}

con.est.D <-
function(xj, xk, est)
{
	eta1 <- est[1]
	eta2 <- est[2]
        eta3 <- est[3]
	a0 <- est[4]
	a1 <- est[5]
	b1 <- est[6]
	c1 <- est[7]
	b0 <- a0 + (a1 - b1) * xj
	c0 <- b0 + (b1 - c1) * xk
	list(eta1 = eta1, eta2 = eta2,eta3 = eta3, a0 = a0, a1 = a1, b0 = b0, b1 = b1, c0 = c0, c1
		 = c1)
}

con.vals.D <-
function(x, y, n, j, k)
{
	a <- con.parms.D(x, y, n, j, k, 1, 1, 1)
	nr <- nrow(a$theta)
	est <- a$theta[nr,  ]
	b <- con.est.D(x[j], x[k], est)
## Yulei's edit
	eta <- c(b$eta1, b$eta2, b$eta3)
	beta <- c(b$a0, b$a1, b$b0, b$b1, b$c0, b$c1)
	tau <- c(x[j], x[k])
	list(eta = eta, beta = beta, tau = tau)
}

p.ll.D <-function(n, j, k, s2, t2, u2){
 q1 <- n * log(sqrt(2 * pi))
 q2 <- 0.5 * (j) * (1 + log(s2))
 q3 <- 0.5 * (k - j) * (1 + log(t2))
## Yulei's edit
q4<-0.5*(n-k)*(1+log(u2))
 - (q1 + q2 + q3 + q4)
}

##Yulei's edit to avoid using for-loop
findmax <-function(a){
	maxa<-max(a)
  imax<- which(a==max(a),arr.ind=TRUE)[1]
  jmax<-which(a==max(a),arr.ind=TRUE)[2]
	list(imax = imax, jmax = jmax, value = maxa)
}

##Yulei's edit
beta.calc.D <-
function(x, y, n, j, k, e1, e2,e3)
{
	aa <- wmat.D(x, y, n, j, k, e1, e2, e3)
	W <- aa$w
	bb <- rvec.D(x, y, n, j, k, e1, e2, e3)
	R <- bb$r
	beta <- solve(W, R)
	list(B = beta)
}

eta.calc.D <-
function(x, y, n, j, k, theta)
{
	jp1 <- j + 1
	kp1 <- k + 1
	a0 <- theta[1]
	a1 <- theta[2]
	b1 <- theta[3]
	c1 <- theta[4]
	b0 <- a0 + (a1 - b1) * x[j]
	c0 <- b0 + (b1 - c1) * x[k]
	rss1 <- sum((y[1:j] - a0 - a1 * x[1:j])^2)
	rss2 <- sum((y[jp1:k] - b0 - b1 * x[jp1:k])^2)
	rss3 <- sum((y[kp1:n] - c0 - c1 * x[kp1:n])^2)
	e1 <- j/rss1
	e2 <- (k - j)/rss2
## Yulei's edit
  e3 <- (n - k)/rss3
      list(eta1 = e1, eta2 = e2, eta3 = e3)
}

##Yulei's edit: e3 is eta used for the third segment
wmat.D <- 
function(x, y, n, j, k, e1, e2, e3)
{
	W <- matrix(0, 4, 4)
	jp1 <- j + 1
	kp1 <- k + 1
	W[1, 1] <- e1 * j + e2 * (k - j) + e3 * (n - k)
	W[1, 2] <- e1 * sum(x[1:j]) + e2 * (k - j) * x[j] + e3 * (n - k) * x[j]
	W[1, 3] <- e2 * sum(x[jp1:k] - x[j]) + e3 * (n - k) * (x[k] - x[j])
	W[1, 4] <- e3 * sum(x[kp1:n] - x[k])
	W[2, 2] <- e1 * sum(x[1:j] * x[1:j]) + e2 * (k - j) * x[j] * x[j]  + e3 * (n - k) * x[j] * x[j] 
	W[2, 3] <- e2 * x[j] * sum(x[jp1:k] - x[j]) + e3 * (n - k) * x[j] * (x[k] - x[j]) 
	W[2, 4] <- e3 * x[j] * sum(x[kp1:n] - x[k])
	W[3, 3] <- e2 * sum((x[jp1:k] - x[j]) * (x[jp1:k] - x[j])) + e3 * (n - k) * (x[k] - x[j]) * (x[k] - x[j])
	W[3, 4] <- e3 * (x[k] - x[j]) * sum(x[kp1:n] - x[k])
	W[4, 4] <- e3 * sum((x[kp1:n] - x[k]) * (x[kp1:n] - x[k]))
	W[2, 1] <- W[1, 2]
	W[3, 1] <- W[1, 3]
	W[4, 1] <- W[1, 4]
	W[3, 2] <- W[2, 3]
	W[4, 2] <- W[2, 4]
	W[4, 3] <- W[3, 4]
	list(w = W)
}
## Yulei's edit
rvec.D <-
function(x, y, n, j, k, e1, e2, e3)
{
	R <- array(0, 4)
	jp1 <- j + 1
	kp1 <- k + 1
	y1j <- sum(y[1:j])
	yjk <- sum(y[jp1:k])
	ykn <- sum(y[kp1:n])
	xy1j <- sum(x[1:j] * y[1:j])
	xyjk <- sum(x[jp1:k] * y[jp1:k])
	xykn <- sum(x[kp1:n] * y[kp1:n])
R[1] <- e1 * y1j + e3 * ykn + e2 * yjk
R[2] <- e1 * xy1j + e3 * x[j] * ykn + e2 * x[j] * yjk
R[3] <- e3 * (x[k] - x[j]) * ykn + e2 * (xyjk - x[j] * yjk)
R[4] <- e3 * (xykn - x[k] * ykn)
	list(r = R)
}