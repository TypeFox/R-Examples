###################################################
### chunk number 1: 
###################################################
#line 6 "cap7.Rnw"
options(prompt="R> ")
options(width=80)


###################################################
### chunk number 2: grid eval=FALSE
###################################################
## #line 58 "cap7.Rnw"
## Smin <- 0
## Smax <- 60
## 
## T <- 10
## N <- 10
## M <- 6
## 
## Dt = T/N
## DS = (Smax-Smin)/M
## 
## t <- seq(0, T, by =Dt)
## S <- seq(Smin, Smax, by=DS)
## 
## plot(range(t),range(S), type="n",xlab="",ylab="",axes=FALSE)
## 
## axis(1, 0:10, c(0,"","","",expression(i*Delta*t), expression((i+1)*Delta*t),"","","","",expression(T)), padj=0.5)
## axis(2, S , c(expression(S[min]), "",expression((j-1)*Delta*S), expression(j*Delta*S), expression((j+1)*Delta*S),"",expression(S[max])))
## 
## abline(v = t, h = S, col = "darkgray", lty = "dotted")
## 
## for(i in t)
##  for(j in S)
##   points(i,j,pch=1,cex=.5)
## 
## text(t[5]-0.1,S[4], expression(C[i][j]),adj=1)


###################################################
### chunk number 3: gridPlot
###################################################
#line 88 "cap7.Rnw"
par(mar=c(3,3,1,1))
#line 58 "cap7.Rnw#from line#89#"
Smin <- 0
Smax <- 60

T <- 10
N <- 10
M <- 6

Dt = T/N
DS = (Smax-Smin)/M

t <- seq(0, T, by =Dt)
S <- seq(Smin, Smax, by=DS)

plot(range(t),range(S), type="n",xlab="",ylab="",axes=FALSE)

axis(1, 0:10, c(0,"","","",expression(i*Delta*t), expression((i+1)*Delta*t),"","","","",expression(T)), padj=0.5)
axis(2, S , c(expression(S[min]), "",expression((j-1)*Delta*S), expression(j*Delta*S), expression((j+1)*Delta*S),"",expression(S[max])))

abline(v = t, h = S, col = "darkgray", lty = "dotted")

for(i in t)
 for(j in S)
  points(i,j,pch=1,cex=.5)

text(t[5]-0.1,S[4], expression(C[i][j]),adj=1)
#line 90 "cap7.Rnw"


###################################################
### chunk number 4: grid2 eval=FALSE
###################################################
## #line 177 "cap7.Rnw"
## Smin <- 0
## Smax <- 60
## 
## T <- 10
## N <- 10
## M <- 6
## 
## Dt = T/N
## DS = (Smax-Smin)/M
## 
## t <- seq(0, T, by =Dt)
## S <- seq(Smin, Smax, by=DS)
## 
## plot(range(t),range(S), type="n",xlab="",ylab="",axes=FALSE)
## 
## axis(1, 0:10, c(0,"","","",expression(i*Delta*t), expression((i+1)*Delta*t),"","","","",expression(T)),padj=0.5)
## axis(2, S , c(expression(S[min]), "",expression((j-1)*Delta*S), expression(j*Delta*S), expression((j+1)*Delta*S),"",expression(S[max])))
## 
## abline(v = t, h = S, col = "darkgray", lty = "dotted")
## 
## for(i in t)
##  for(j in S)
##   points(i,j,pch=1,cex=.5)
## 
## 
## 
## text(t[5]-0.1,S[4], expression(C[i][j]),adj=1)
## 
## text(t[6]+0.1,S[5], expression(C[i+1][j+1]),adj=0)
## text(t[6]+0.1,S[4], expression(C[i+1][j]),adj=0)
## text(t[6]+0.1,S[3], expression(C[i+1][j-1]),adj=0)
## 
## 
## 
## arrows(t[6]-0.1,S[5]-0.5,t[5]+0.1,S[4]+0.5,length=0.05,lwd=2)
## arrows(t[6]-0.1,S[4],t[5]+0.1,S[4],length=0.05,lwd=2)
## arrows(t[6]-0.1,S[3]+0.5,t[5]+0.1,S[4]-0.5,length=0.05,lwd=2)
## 
## 
##  for(j in S[2:M])
##   points(t[N+1],j,pch=21,bg="red",cex=1.5)
## 
##  for(i in t)
##   points(i,S[1],pch=22,bg="green",cex=1.5)
## 
##  for(i in t)
##   points(i,S[M+1],pch=23,bg="blue",cex=1.5)


###################################################
### chunk number 5: gridPlot2
###################################################
#line 228 "cap7.Rnw"
par(mar=c(3,3,1,1))
#line 177 "cap7.Rnw#from line#229#"
Smin <- 0
Smax <- 60

T <- 10
N <- 10
M <- 6

Dt = T/N
DS = (Smax-Smin)/M

t <- seq(0, T, by =Dt)
S <- seq(Smin, Smax, by=DS)

plot(range(t),range(S), type="n",xlab="",ylab="",axes=FALSE)

axis(1, 0:10, c(0,"","","",expression(i*Delta*t), expression((i+1)*Delta*t),"","","","",expression(T)),padj=0.5)
axis(2, S , c(expression(S[min]), "",expression((j-1)*Delta*S), expression(j*Delta*S), expression((j+1)*Delta*S),"",expression(S[max])))

abline(v = t, h = S, col = "darkgray", lty = "dotted")

for(i in t)
 for(j in S)
  points(i,j,pch=1,cex=.5)



text(t[5]-0.1,S[4], expression(C[i][j]),adj=1)

text(t[6]+0.1,S[5], expression(C[i+1][j+1]),adj=0)
text(t[6]+0.1,S[4], expression(C[i+1][j]),adj=0)
text(t[6]+0.1,S[3], expression(C[i+1][j-1]),adj=0)



arrows(t[6]-0.1,S[5]-0.5,t[5]+0.1,S[4]+0.5,length=0.05,lwd=2)
arrows(t[6]-0.1,S[4],t[5]+0.1,S[4],length=0.05,lwd=2)
arrows(t[6]-0.1,S[3]+0.5,t[5]+0.1,S[4]-0.5,length=0.05,lwd=2)


 for(j in S[2:M])
  points(t[N+1],j,pch=21,bg="red",cex=1.5)

 for(i in t)
  points(i,S[1],pch=22,bg="green",cex=1.5)

 for(i in t)
  points(i,S[M+1],pch=23,bg="blue",cex=1.5)
#line 230 "cap7.Rnw"


###################################################
### chunk number 6: 
###################################################
#line 255 "cap7.Rnw"
AmericanPutExp <- function(Smin=0, Smax,  T=1, N=10, M=10, K, r=0.05, sigma=0.01){
 Dt = T/N 
 DS = (Smax-Smin)/M
 t <- seq(0, T, by =Dt) 
 S <- seq(Smin, Smax, by=DS)
 A <- function(j) (-0.5*r*j*Dt + 0.5*sigma^2*j^2*Dt)/(1+r*Dt) 
 B <- function(j) (1-sigma^2*j^2*Dt)/(1+r*Dt) 
 C <- function(j) (0.5*r*j*Dt + 0.5*sigma^2*j^2*Dt)/(1+r*Dt)
 P <- matrix(, M+1, N+1)
 colnames(P) <- round(t,2)
 rownames(P) <- round(rev(S),2)
 P[M+1, ] <- K   # C(,j=0) = K
 P[1,] <- 0   # C(,j=M) = 0
 P[,N+1] <- sapply(rev(S), function(x) max(K-x,0))
 optTime <- matrix(FALSE, M+1, N+1)
 optTime[M+1,] <- TRUE
 optTime[which(P[,N+1]>0),N+1] <- TRUE

 for(i in (N-1):0){
  for(j in 1:(M-1)){
   J <- M+1-j
   I <- i+1
   P[J,I] <- A(j)*P[J+1,I+1] + B(j)*P[J,I+1] + C(j)*P[J-1,I+1]
   if(P[J,I] < P[J,N+1])
    optTime[J,I] <- TRUE
  }
 }
 colnames(optTime) <- colnames(P)
 rownames(optTime) <- rownames(P)
 ans <- list(P=P, t=t, S=S, optTime=optTime,N=N,M=M)
 class(ans) <- "AmericanPut"
 return(invisible(ans))
}


###################################################
### chunk number 7: 
###################################################
#line 293 "cap7.Rnw"
plot.AmericanPut <- function( obj ){
 plot(range(obj$t),range(obj$S),type="n",axes=F,xlab="t", ylab="S")
 axis(1,obj$t,obj$t)
 axis(2,obj$S,obj$S)
 abline(v = obj$t, h = obj$S, col = "darkgray", lty = "dotted")
 for(i in 0:obj$N){
  for(j in 0:obj$M){
   J <- obj$M+1-j
   I <- i+1
   cl <- "grey"; 
   if(obj$optTime[J,I])
    cl <- "black"
   text(obj$t[i+1],obj$S[j+1], round(obj$P[J,I],2),cex=0.75, col=cl)
  }
 }
 DS <- mean(obj$S[1:2])
 y <- as.numeric(apply(obj$optTime,2, function(x) which(x)[1]))
 lines(obj$t, obj$S[obj$M+2-y]+DS, lty=2)
}


###################################################
### chunk number 8: 
###################################################
#line 315 "cap7.Rnw"
put <- AmericanPutExp(Smax = 60, sigma = 0.4, K = 30)
round(put$P,2)


###################################################
### chunk number 9: pput1 eval=FALSE
###################################################
## #line 320 "cap7.Rnw"
## plot(put)


###################################################
### chunk number 10: pput1Plot
###################################################
#line 324 "cap7.Rnw"
par(mar=c(3,3,1,1))
#line 320 "cap7.Rnw#from line#325#"
plot(put)
#line 326 "cap7.Rnw"


###################################################
### chunk number 11: 
###################################################
#line 330 "cap7.Rnw"
S0 <- 36
myval <- round(put$P[which(rownames(put$P)==S0),1],2)


###################################################
### chunk number 12: 
###################################################
#line 347 "cap7.Rnw"
put.bad <- AmericanPutExp(Smax = 60, sigma = 0.4, K = 30, M=15)
round(put.bad$P,2)


###################################################
### chunk number 13: 
###################################################
#line 420 "cap7.Rnw"
AmericanPutImp <- function( Smin=0, Smax,  T=1, N=10, M=10, K, r=0.05, sigma=0.01){
 Dt = T/N 
 DS = (Smax-Smin)/M
 t <- seq(0, T, by =Dt) 
 S <- seq(Smin, Smax, by=DS)

 A <- function(j) 0.5*r*j*Dt - 0.5*sigma^2*j^2*Dt 
 B <- function(j) 1+sigma^2*j^2*Dt+r*Dt
 C <- function(j) -0.5*r*j*Dt - 0.5*sigma^2*j^2*Dt

 a <- sapply(0:M, A)
 b <- sapply(0:M, B)
 c <- sapply(0:M, C)

 P <- matrix(, M+1, N+1)
 colnames(P) <- round(t,2)
 rownames(P) <- round(rev(S),2)
 
 P[M+1, ] <- K   # C(,j=0) = K
 P[1,] <- 0   # C(,j=M) = 0
 P[,N+1] <- sapply(rev(S), function(x) max(K-x,0))
 
 AA <- matrix(0, M-1, M-1)
 for(j in 1:(M-1)){
  if(j>1) AA[j,j-1] <- A(j)
  if(j<M) AA[j,j] <- B(j)
  if(j<M-1) AA[j,j+1] <- C(j)
 }

 optTime <- matrix(FALSE, M+1, N+1)
 for(i in (N-1):0){
  I <- i+1
  bb <- P[M:2,I+1]
  bb[1] <- bb[1]-A(1)*P[M+1-0,I+1]
  bb[M-1] <- bb[M-1]-C(M-1)*P[M+1-M,I+1] 
  P[M:2,I] <- solve(AA,bb)
  idx <- which(P[,I] < P[,N+1])
  P[idx,I] <- P[idx,N+1] 
  optTime[idx, I] <- TRUE
 }
 optTime[M+1,] <- TRUE 
 optTime[which(P[,N+1]>0),N+1] <- TRUE
 colnames(optTime) <- colnames(P)
 rownames(optTime) <- rownames(P)
 ans <- list(P=P, t=t, S=S, optTime=optTime,N=N,M=M)
 class(ans) <- "AmericanPut"
 return(invisible(ans))
}


###################################################
### chunk number 14: 
###################################################
#line 472 "cap7.Rnw"
put <- AmericanPutImp(Smax = 60, sigma = 0.4, K = 30)
round(put$P,2)


###################################################
### chunk number 15: pput2 eval=FALSE
###################################################
## #line 477 "cap7.Rnw"
## plot(put)


###################################################
### chunk number 16: pput2Plot
###################################################
#line 482 "cap7.Rnw"
par(mar=c(3,3,1,1))
#line 477 "cap7.Rnw#from line#483#"
plot(put)
#line 484 "cap7.Rnw"


###################################################
### chunk number 17: 
###################################################
#line 586 "cap7.Rnw"
require(fOptions)
T <- 1
sigma=0.4
r=0.05
S0 <- 36
K <- 30
BAWAmericanApproxOption("p", S=S0, X=K, Time=T, r=r, b=r, sigma=sigma)@price


###################################################
### chunk number 18: 
###################################################
#line 596 "cap7.Rnw"
put <- AmericanPutImp(0,100,T=T, K=K,r=r,sigma=sigma,M=100,N=100)
put$P["36",1]


###################################################
### chunk number 19: BG1
###################################################
#line 646 "cap7.Rnw"
par(mar=c(0,0,0,0))
plot(c(0,10), c(30,120), type = "n", xlab=" ", ylab="", axes=F)
  
 points(0,75)     
 points(5,100)
points(5,75)
points(5,50)

points(9,110)
points(9,100)
points(9,90)

points(9,85)
points(9,75)
points(9,65)

points(9,60)
points(9,50)
points(9,40)

segments(0,75,5,100)
segments(0,75,5,75)
segments(0,75,5,50)

segments(5,100,9,110)
segments(5,100,9,100)
segments(5,100,9,90)

segments(5,75,9,85)
segments(5,75,9,75)
segments(5,75,9,65)

segments(5,50,9,60)
segments(5,50,9,50)
segments(5,50,9,40)

text(0,79,expression(S[0]))     
 text(5,106,expression(S[1]^{1}))
text(5,81,expression(S[1]^{2}))
text(5,56,expression(S[1]^{3}))

text(9.5,110,expression(S[T]^{11}))
text(9.5,100,expression(S[T]^{12}))
text(9.5,90,expression(S[T]^{13}))

text(9.5,85,expression(S[T]^{21}))
text(9.5,75,expression(S[T]^{22}))
text(9.5,65,expression(S[T]^{23}))

text(9.5,60,expression(S[T]^{31}))
text(9.5,50,expression(S[T]^{32}))
text(9.5,40,expression(S[T]^{33}))


###################################################
### chunk number 20: 
###################################################
#line 715 "cap7.Rnw"
simTree <- function(b,d, S0, sigma, T, r){
 tot <- sum(b^(1:(d-1)))
 S <- numeric(tot+1) 
 S[1] <- S0
 dt <- T/d
 for(i in 0:(tot - b^(d-1))){
  for(j in 1:b){
   S[i*b+j +1] <- S[i+1]*exp((r - 0.5*sigma^2)*dt + sigma*sqrt(dt)*rnorm(1))
  }
 }
 S
}


###################################################
### chunk number 21: 
###################################################
#line 730 "cap7.Rnw"
upperBG <- function(S, b, d, f){
 tot <- sum(b^(1:(d-1)))
 start <- tot - b^(d-1) +1
 end <- tot +1
 P <- S
 P[start:end] <- f(S[start:end])
 tot1 <- sum(b^(1:(d-2)))
 for(i in tot1:0){
  m <- mean(P[i*b+1:b+1])
  v <- f(S[i+1])
  P[i+1] <- max(v,m)
 }
 P
}

lowerBG <- function(S, b, d, f){
 tot <- sum(b^(1:(d-1)))
 start <- tot - b^(d-1) +1
 end <- tot +1
 p <- S 
 p[start:end] <- f(S[start:end])
 tot1 <- sum(b^(1:(d-2)))

 m <- numeric(b)
 for(i in tot1:0){
  v <- f(S[i+1])
  for(j in 1:b){
   m[j] <- mean(p[i*b+(1:b)[-j]+1])
   m[j] <- ifelse( v>m[j], v, p[i*b+(1:b)[j]+1])
  }
  p[i+1] <- mean(m)
 }
 p
}


###################################################
### chunk number 22: 
###################################################
#line 767 "cap7.Rnw"
b <- 3
d <- 3
S0 <- 101
S <- c(101, 114, 50, 115, 74, 88, 102, 38, 47, 65, 88, 149, 116)
K <- 100
f <- function(x) sapply(x, function(x) max(x-K,0))


###################################################
### chunk number 23: 
###################################################
#line 777 "cap7.Rnw"
lowerBG(S, b,d,f)
upperBG(S, b,d,f)


###################################################
### chunk number 24: 
###################################################
#line 783 "cap7.Rnw"
set.seed(123)
b <- 3
d <- 3
K <- 100
f <- function(x) sapply(x, function(x) max(x-K,0))
T <- 1
r <- 0.05
sigma <- 0.4
S0 <- 101

S <- simTree(b,d, S0, sigma, T, r)
lowerBG(S, b,d,f)[1]
upperBG(S, b,d,f)[1]


###################################################
### chunk number 25: BG2
###################################################
#line 801 "cap7.Rnw"
par(mar=c(0,0,0,0))
plot(c(-1,10), c(30,120), type = "n", xlab=" ", ylab="", axes=F)
  
points(0,75)     
points(5,100)
points(5,75)
points(5,50)

points(9,110)
points(9,100)
points(9,90)

points(9,85)
points(9,75)
points(9,65)

points(9,60)
points(9,50)
points(9,40)

segments(0,75,5,100)
segments(0,75,5,75)
segments(0,75,5,50)

segments(5,100,9,110)
segments(5,100,9,100)
segments(5,100,9,90)

segments(5,75,9,85)
segments(5,75,9,75)
segments(5,75,9,65)

segments(5,50,9,60)
segments(5,50,9,50)
segments(5,50,9,40)

text(0,81,"101(1)[11.9]")     
text(5,106,"114(14)[14]")
text(5,81,"50(0)[0]")
text(5,56,"115(15)[21.7]")

text(9.6,110,"74(0)")
text(9.6,100,"88(0)")
text(9.6,90,"102(2)")

text(9.6,85,"38(0)")
text(9.6,75,"47(0)")
text(9.6,65,"65(0)")

text(9.6,60,"88(0)")
text(9.6,50,"149(49)")
text(9.6,40,"116(16)")

text(0,30, expression(t==0))
text(5,30, expression(t==t[1]))
text(9.5,30, expression(t==T))


###################################################
### chunk number 26: BG3
###################################################
#line 866 "cap7.Rnw"
par(mar=c(0,0,0,0))
plot(c(-1,10), c(30,120), type = "n", xlab=" ", ylab="", axes=F)
  
points(0,75)     
points(5,100)
points(5,75)
points(5,50)

points(9,110)
points(9,100)
points(9,90)

points(9,85)
points(9,75)
points(9,65)

points(9,60)
points(9,50)
points(9,40)

segments(0,75,5,100)
segments(0,75,5,75)
segments(0,75,5,50)

segments(5,100,9,110)
segments(5,100,9,100)
segments(5,100,9,90)

segments(5,75,9,85)
segments(5,75,9,75)
segments(5,75,9,65)

segments(5,50,9,60)
segments(5,50,9,50)
segments(5,50,9,40)

text(0,81,"101(1)[8.1]")     
text(5,106,"114(14)[14]")
text(5,81,"50(0)[0]")
text(5,56,"115(15)[10.3]")

text(9.6,110,"74(0)")
text(9.6,100,"88(0)")
text(9.6,90,"102(2)")

text(9.6,85,"38(0)")
text(9.6,75,"47(0)")
text(9.6,65,"65(0)")

text(9.6,60,"88(0)")
text(9.6,50,"149(49)")
text(9.6,40,"116(16)")

text(0,30, expression(t==0))
text(5,30, expression(t==t[1]))
text(9.5,30, expression(t==T))


###################################################
### chunk number 27: 
###################################################
#line 1079 "cap7.Rnw"
LSM <- function(n, d, S0, K, sigma, r, T){
 s0 <- S0/K
 dt <- T/d
 z <- rnorm(n)
 s.t <- s0*exp((r-1/2*sigma^2)*T+sigma*z*(T^0.5))
 s.t[(n+1):(2*n)] <- s0*exp((r-1/2*sigma^2)*T-sigma*z*(T^0.5))
 CC <- pmax(1-s.t, 0)
 payoffeu <- exp(-r*T)*(CC[1:n]+CC[(n+1):(2*n)])/2*K
 euprice <- mean(payoffeu)

 for(k in (d-1):1){
   z <- rnorm(n)
   mean <- (log(s0) + k*log(s.t[1:n]))/(k+1)
   vol <- (k*dt/(k+1))^0.5*z
   s.t.1 <- exp(mean+sigma*vol)
   mean <- (log(s0) + k*log( s.t[(n+1):(2*n)] )) / ( k + 1 )
   s.t.1[(n+1):(2*n)] <- exp(mean-sigma*vol)
   CE <- pmax(1-s.t.1,0)
   idx<-(1:(2*n))[CE>0]
   discountedCC<- CC[idx]*exp(-r*dt)
   basis1 <- exp(-s.t.1[idx]/2)
   basis2 <- basis1*(1-s.t.1[idx])
   basis3 <- basis1*(1-2*s.t.1[idx]+(s.t.1[idx]^2)/2)

   p <- lm(discountedCC ~ basis1+basis2+basis3)$coefficients
   estimatedCC <- p[1]+p[2]*basis1+p[3]*basis2+p[4]*basis3
   EF <- rep(0, 2*n)
   EF[idx] <- (CE[idx]>estimatedCC)
   CC <- (EF == 0)*CC*exp(-r*dt)+(EF == 1)*CE
   s.t <- s.t.1
  }

  payoff <- exp(-r*dt)*(CC[1:n]+CC[(n+1):(2*n)])/2
  usprice <- mean(payoff*K)
  error <- 1.96*sd(payoff*K)/sqrt(n)
  earlyex <- usprice-euprice
  data.frame(usprice, error, euprice)
}


###################################################
### chunk number 28: 
###################################################
#line 1120 "cap7.Rnw"
S0 <- 36
K <- 30
T <- 1
r <- 0.05
sigma <- 0.4
LSM(10000, 3, S0, K, sigma, r, T)
require(fOptions)
BSAmericanApproxOption("p", S0,K,T,r,r,sigma)@price
BAWAmericanApproxOption("p",S0, K, T,r,r, sigma)@price


