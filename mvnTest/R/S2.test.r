## function to calculate Y2, U2, S2 tests...

S2.test <- function(data,M=5,qqplot=FALSE){
## data - as n times p matrix or data frame, where p is dimension
## M - number of equiprobable intervals
## qqplot - idea taken from Korkmaz et al (2015)
   if (!is.data.frame(data) && !is.matrix(data)) stop('data supplied must be either of class \"data frame\" or \"matrix\"')
   if (dim(data)[2] < 2 || is.null(dim(data))) {stop("data dimesion has to be more than 1")}
   if (dim(data)[1] < 3) {stop("not enough data for assessing mvn")}
   data.name <- deparse(substitute(data))
   xp <- as.matrix(data)
   p <- dim(xp)[2]      
   n <- dim(xp)[1]
   b <- rep(0, M)  ## getting c_i - ending points of the intervals
   for (i in 1:(M-1))   ## other c[i]... 
        b[i] <- qchisq(i/M,df=p) 
   b[M] <- 1e20

  ## getting MLEs of initial data...
  s.mean <- colMeans(xp)
  s.cov <- (n-1)/n*cov(xp)
  ## Karhunen-Loeve transformation....
  eigen.s <- eigen(s.cov,symmetric=TRUE) # eigendecomposition
  y <- matrix(NA,n,p)
  for (i in 1:n)  y[i,] <- t(eigen.s$vectors)%*%xp[i,]
  ## getting MLEs of the transformed data...
 # s.mean <- colMeans(y)
  y.centered <- scale(y,scale=FALSE) # centering by subtracting the column means from the corresponding columns
  s.cov <- (n-1)/n*cov(y)
  s.cov.inv <- solve(s.cov) # inverse matrix of S (matrix of sample covariances)
  xSx <- rep(NA,n) # vector of (Xi-mu)'S^-1(Xi-mu)...
  for (i in 1:n) # xSx[i] <- t(y[i,]-s.mean)%*%s.cov.inv%*%(y[i,]-s.mean)
       xSx[i] <- t(y.centered[i,])%*%s.cov.inv%*%y.centered[i,]
  ## calculation of the cells frequencies...
  Nn <- rep(0,M)
  for (i in 1:n){
      j <-1
      done <- FALSE 
      while (j <= (M) & !done){
           if (xSx[i] > b[M-1]) {
                   Nn[M] <- Nn[M] + 1
                   done <- TRUE;
           }
           else if (xSx[i] > b[j]) {j <- j+1}
                else {
                       Nn[j] <- Nn[j] + 1
                       done <- TRUE;
                }
      }
  }
  Vn <- rep(NA,M)     # define M-vector of standardized cell frequencies 
  Vn <- (Nn-(n/M))/sqrt(n/M)
  m <- p+p*(p+1)/2    # dimension of the parameter vector theta
  d <- rep(NA,M)      # sub-elements of matrix B
  bp1 <- 1
  for (i in 0:(floor(p/2)-1))  bp1 <- bp1*(p-2*i)
  if (p%%2 != 0)    bp <- (2/pi)^0.5/bp1  # for p odd
   else        bp <- 1/bp1      # for p even
  d[1] <- (-b[1]^(p/2)*exp(-b[1]/2))*bp/2
  for (i in 2:M)  
    d[i] <- (b[i-1]^(p/2)*exp(-b[i-1]/2)- b[i]^(p/2)*exp(-b[i]/2))*bp/2  
  Vd <- (sum(Vn*d))^2 
  d2 <- sum(d^2)
  ChLeh <- sum(Vn^2)
  y2 <- ChLeh + (2*p*M*Vd)/(1-2*p*M*d2)
  u2 <- ChLeh -Vd/d2  # value of the DN statistics
  s2 <- y2 - u2
  if (qqplot) {   ## NOT checked yet....
            xp <- scale(xp, scale = FALSE)
            Sa <- cov(xp)
            D <- xp %*%solve(Sa)%*% t(xp)
            d <- diag(D)    
            r <- rank(d)  
            chi2q <- qchisq((r-0.5)/n,p)
            plot(d, chi2q, pch = 19, main = "Chi-Square Q-Q Plot", xlab = "Squared Mahalanobis Distance", ylab = "Chi-Square Quantile")
            abline(0, 1,lwd = 2, col = "black")
        }
  pv.y2 <- pchisq(y2, df=M-1, lower.tail=FALSE) ## p-value of Y2
  pv.s2 <- pchisq(s2, df=1, lower.tail=FALSE) ## p-value of S2
  pv.u2 <- pchisq(u2, df=M-2, lower.tail=FALSE) ## p-value of U2

  result <- new("S2",y2=y2, s2=s2, u2=u2, data.name=data.name,  p.value.y2=pv.y2, p.value.s2=pv.s2, p.value.u2=pv.u2)
  result
} ## S2.test


