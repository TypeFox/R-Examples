mult.chart <-
function(type = c("chi","t2", "mewma", "mcusum", "mcusum2"),
                 x, Xmv, S, colm, alpha = 0.01, lambda = 0.1, k = 0.5,
 h = 5.5, phase = 1, method = "sw", ...){

type <- match.arg(type) 
 
# Notes: 
# 1. Multivariate Control Chart (mcc) performs the Chi-squared, T2, MEWMA and MCUSUM according
#### to Crosier(1988) and Pignatiello and Runger(1990)
#### Arguments used (optionally or not) by function: 
#### -chi-squared(x, Xmv, S, alpha, method)
#### -t2(x, Xmv, S, colm, alpha, phase, method)
#### -memwma(x, Xmv, S, lambda, method)
#### -mcusum(x, Xmv, S, k, h, method)
 
# 2. To perform the chart using a historical mean vector and covariance matrix
#### just enter this values such as Xmv, S. Moreover the number of samples (colm) is requiered
#### Then, the UCL used in the graph and in the decomposition is calculated 
#### based on the historical data.
#### 
     
 ###Variables
p <- ncol(x) # quality characteristics
m <- nrow(x) # number of samples or observations
if (class(x) == "matrix" || class(x) == "data.frame") (x <- array(data.matrix(x),c(m,p,1)))
n <- dim(x)[3] # observations or sample size 
 
if(!missing(Xmv))(phase <- 2)

x.jk <- matrix(0,m,p)
t2 <- matrix(0,m,1)

x.jk <- apply(x,1:2,mean)

if(missing(Xmv))(Xmv <- colMeans(x.jk)) 
if(missing(S))(S <- covariance(x,method = method))
if(missing(colm))(colm <- nrow(x))

if (type == "chi") { # Chi squared Control Chart
name <- paste("Chi-squared Control Chart")
for (ii in 1 : m){
t2[ii,1] <- n * t(x.jk[ii,] - Xmv) %*% solve(S) %*% (x.jk[ii,] - Xmv)
}
ucl <- qchisq(1 - alpha, p)
if(any(t2>ucl)){
cat("The following(s) point(s) fall outside the control limits" )
t3 <- which(t2 > ucl)
print( t3)
}


}
if (type == "t2") { # Hotelling Control Chart
name <- paste("Hotelling Control Chart")
for (ii in 1 : m){
t2[ii,1] <- n * t(x.jk[ii,] - Xmv) %*% solve(S) %*% (x.jk[ii,] - Xmv)
}

ifelse(n == 1, ifelse(phase == 1, 
 ucl <- ((colm - 1) ^ 2) / colm * qbeta(1 - alpha,p / 2,(((2 * (colm - 1) ^ 2) / (3 * colm - 4) - p - 1) / 2)),
 ucl<-((p * (colm + 1) * (colm - 1)) / ((colm ^ 2) - colm * p)) * qf(1 - alpha,p,colm - p)),
 ifelse(phase == 1, 
 ucl <- (p * (colm - 1) * (n - 1)) / (colm * n - colm - p + 1) * qf(1 - alpha,p,colm * n - colm - p + 1),
 ucl <- (p * (colm + 1) * (n - 1)) / (colm * n - colm - p + 1) * qf(1 - alpha,p,colm * n - colm - p + 1))
)
if(any(t2>ucl)){
cat("The following(s) point(s) fall outside of the control limits" )
t3 <- which(t2 > ucl)
print( t3)

#Decomposition of T2
for(ii in 1:length(t3)){
    
v=1
k=0 
      
for(i in 1:p) {
k <- k + factorial(p) / (factorial(i) * factorial(p - i)) # total number of permutations 
}
 
q <- matrix(0, k, p + 3)  # matrix of all the permutations

for(i in 1:p) { 
a <- t(combn(p, i))
   
 for(l in 1:nrow(a)) {
 for(j in 1:ncol(a)) { 
 q[v, j + 3] <- a[l, j]
} 
 v = v + 1
}
}
 
 for(i in 1:nrow(q)) { 
 b <- subset(q[i,4:ncol(q)], q[i,4:ncol(q)] > 0) 
 di<-length(b)
   
 if(length(b) > 1) {q[i,1] <- n * t(Xmv[b] - x.jk[t3[ii],][b]) %*% solve(S [b,b]) %*% (Xmv[b]-x.jk[t3[ii],][b])}  #this column is the t2
 else( q[i,1] <- n * (x.jk[t3[ii],][b] - Xmv[b]) ^ 2 / S [b,b])

 ifelse(n == 1, ifelse(phase == 1, 
 q[i,2] <- ((colm - 1) ^ 2) / colm * qbeta(1 - alpha,di/ 2,(((2 * (colm - 1) ^ 2) / (3 * colm - 4) - di- 1) / 2)),
 q[i,2]<-((di* (colm + 1) * (colm - 1)) / ((colm ^ 2) - colm * di)) * qf(1 - alpha,di,colm - di)),
 ifelse(phase == 1, 
 q[i,2] <- (di* (colm - 1) * (n - 1)) / (colm * n - colm - di+ 1) * qf(1 - alpha,di,colm * n - colm - di+ 1),
 q[i,2]<- (di* (colm + 1) * (n - 1)) / (colm * n - colm - di+ 1) * qf(1 - alpha,di,colm * n - colm - di+ 1))
)
  
 q[i,3] <- 1 - pf(q[i,1], di, colm - 1) #this column is of p-values
}

 colnames(q)<-c("t2 decomp","ucl","p-value",1:p)
 print(list("Decomposition of"=t3[ii]))
 print(round(q,4))
          
     }
 
 
        }


}

if (type == "mewma") { # MEWMA Control Chart
h4 <- matrix(c(8.6336,9.6476,10.083,10.3114,10.4405,10.5152,10.5581,10.5816,10.5932,10.814,
11.8961,12.3505,12.5845,12.7143,12.788,12.8297,12.8524,12.8635,12.7231,13.8641,14.3359,
14.576,14.7077,14.7818,14.8234,14.846,14.857,14.5363,15.7293,16.217,16.4629,16.5965,16.6711,
16.7127,16.7352,16.7463,16.2634,17.5038,18.0063,18.2578,18.3935,18.4687,18.5105,18.5331,
18.5442,17.9269,19.2113,19.7276,19.9845,20.1223,20.1982,20.2403,20.2631,20.2743,19.541,
20.8665,21.396,21.6581,21.798,21.8747,21.9171,21.9401,21.9515,21.1152,22.4796,23.0217,23.2887,
23.4307,23.5082,23.551,23.5742,23.5858,22.6565,24.0579,24.6119,24.8838,25.0278,25.1062,25.1493,
25.1728,25.1846),nrow = 9)

rownames(h4) <- c(seq(0.1,0.9, by  = 0.1))
colnames(h4) <- c(1:9)

z<-matrix(0, m, p)
m1 <- rownames(h4)
m2 <- colnames(h4)
l <- lambda * 10
ucl <- h4[m1[l], m2[p - 1]]

name <- paste("MEWMA Control Chart")

for (i in 1 : m){
 if(i==1){
 z[i,] <- lambda * (x.jk[i,] - Xmv)}
 else{ z[i,] <- lambda * (x.jk[i,] - Xmv) + (1 - lambda) * z[i - 1,]}
 weig <- S * (lambda * (1 - ((1 - lambda) ^ (2 * i))) / (2 - lambda))
 t2[i,1] <- t(z[i,]) %*% solve(weig) %*% z[i,]
}
}

if (type == "mcusum") { # MCUSUM Control Chart according to Crosier(1988)
name <- paste("MCUSUM Control Chart by Crosier (1988)")
ucl <- h
dif <- sweep(x.jk,2,Xmv)
s <- matrix(0,m,p)
ci <- matrix(0,m,1)

ci[1] <- sqrt(dif[1,] %*% solve((S / n)) %*% dif[1,])

if (ci[1] > k) { s[1,] <- (s[1,] + dif[1,]) * (1 - k / ci[1])}
else(s[1,] = matrix(0,ncol = p)) #compute s
  
for (i in 2:m){ 
ci[i,]=sqrt((s[i - 1,] + dif[i,]) %*% solve(S / n) %*% (s[i - 1,] + dif[i,])) 
if (ci[i] > k){ s[i,] = (s[i - 1,] + dif[i,]) * (1 - k / ci[i])} 
else {s[i,] = matrix(0,ncol = p)}
}

for (i in 1:m){
 t2[i]=sqrt(s[i,]%*%solve((S/n))%*%(s[i,]))
} 
}

if (type == "mcusum2") { # MCUSUM Control Chart according to Pignatiello and Runger(1990)
name <- paste("MCUSUM Control Chart by Pignatiello (1990)")
ucl <- h
dif <- sweep(x.jk,2,Xmv)

s <- matrix(0,m,p)
l <- matrix(0,m,1)


for (i in 1:m){
 if (i == 1){l[i,1] <- 1 }
 if (i > 1) {if (t2[i-1,1] > 0) {l[i,1] <- l[i-1,1] + 1} else{ l[i,1] <- 1}}
 if(i == ((i-l[i,1]+1))) {s[i,] <- dif[i,]} else{s[i,] <- colSums(dif[(i-l[i,1]+1):i,])}      
 t2[i,1] <- max(0, (t(s[i,]) %*% solve(S / n) %*% s[i,]) ^.5 - k * l[i,1])
}

}

t3 <- which(t2 > ucl)
par(mar = c(4,5,3,5))
plot(t2, ylim = c(0,1.1 * max(max(t2),ucl)), main = name,
xlab = "Sample", ylab = expression(T ^ 2), type = "o",las = 1) 

 
points(t3, t2[t3], col = 2) 
segments(0, ucl, m, ucl, col = 2)

mtext(paste(" UCL=", round(ucl, 2)),side = 4, at = ucl, las = 2)

outList = list (name,"ucl" = round(ucl,2), t2 = round(t2,2), Xmv = round(Xmv,2), covariance = signif(S,2))
 return(outList)

### Falta
#eliminar los ifelse de t2
#poner h4 como dataset
###

}
