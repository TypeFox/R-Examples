ddhaml <-
function(x, y, L, U) {
 N <- length(y)
 v <- c()
v1 <- c()
v2 <- c()
v3 <- c()
v11 <- c()
v12 <- c()
v13 <- c()
v22 <- c()
v23 <- c()
v33 <- c()
ddhamm <- mat.or.vec(3,3)
 v[1] <- x[3] / (1 - x[1] - x[2])
v1[1] <- x[3] / (1 - x[1] - x[2]) ^ 2
v2[1] <- x[3] / (1 - x[1] - x[2]) ^ 2
v3[1] <- 1 / (1 - x[1] - x[2])
u <- y[1] ^ 2
v11[1] <- -2 * x[3] / (1 - x[1] - x[2]) ^ 3
v12[1] <- -2 * x[3] / (1 - x[1] - x[2]) ^ 3
v13[1] <- 1 / (1 - x[1] - x[2]) ^ 2
v22[1] <- -2 * x[3] / (1 - x[1] - x[2]) ^ 3
v23[1] <- 1 / (1 - x[1] - x[2]) ^ 2
v33[1] <- 0
if (y[1] <= L){
temp1 <- -L*jNormPdf(L/sqrt(v[1]))/(jNormCdf(L/sqrt(v[1]))*sqrt(v[1])^3)
temp2 <- 1.5*L*jNormPdf(L/sqrt(v[1]))/(jNormCdf(L/sqrt(v[1]))*sqrt(v[1])^5)+
0-0.5*L^3/(v[1]^2)*jNormPdf(L/sqrt(v[1]))/(jNormCdf(L/sqrt(v[1]))*sqrt(v[1])^3)+
0-0.5*L^2/(sqrt(v[1])^3)*jNormPdf(L/sqrt(v[1]))^2/(jNormCdf(L/sqrt(v[1]))^2*sqrt(v[1])^3)
}
if (y[1] >= U){
temp1 <- U*jNormPdf(-U/sqrt(v[1]))/(jNormCdf(-U/sqrt(v[1]))*sqrt(v[1])^3)
temp2 <- -1.5*U*jNormPdf(-U/sqrt(v[1]))/(jNormCdf(-U/sqrt(v[1]))*sqrt(v[1])^5)+
0.5*U^3/(v[1]^2)*jNormPdf(-U/sqrt(v[1]))/(jNormCdf(-U/sqrt(v[1]))*sqrt(v[1])^3)+
0-0.5*U^2/(sqrt(v[1])^3)*jNormPdf(-U/sqrt(v[1]))^2/(jNormCdf(-U/sqrt(v[1]))^2*sqrt(v[1])^3)
}
if (y[1] > L & y[1] < U){
temp1 <- -1 / v[1] + u / v[1] ^ 2
temp2 <- 1 / v[1] ^ 2 - 2 * u / v[1] ^ 3
}
ddhamm[1, 1]<- temp1 * v11[1] + temp2 * v1[1] ^ 2
ddhamm[1, 2]<- temp1 * v12[1] + temp2 * v1[1] * v2[1]
ddhamm[1, 3]<- temp1 * v13[1] + temp2 * v1[1] * v3[1]
ddhamm[2, 2]<- temp1 * v22[1] + temp2 * v2[1] ^ 2
ddhamm[2, 1]<- ddhamm[1, 2]
ddhamm[2, 3]<- temp1 * v23[1] + temp2 * v2[1] * v3[1]
ddhamm[3, 3]<- temp1 * v33[1] + temp2 * v3[1] ^ 2
ddhamm[3, 1]<- ddhamm[1, 3]
ddhamm[3, 2]<- ddhamm[2, 3]
for (i in 2:N) {
  v[i] <- x[3] + x[1] * u + x[2] * v[i - 1]
 v1[i] <- u + x[2] * v1[i - 1]
 v2[i] <- x[2] * v2[i - 1] + v[i - 1]
 v3[i] <- 1 + x[2] * v3[i - 1]
v11[i] <- x[2] * v11[i - 1]
v12[i] <- x[2] * v12[i - 1] + v1[i - 1]
v13[i] <- 0
v22[i] <- x[2] * v22[i - 1] + 2 * v2[i - 1]
v23[i] <- x[2] * v23[i - 1] + v3[i - 1]
v33[i] <- x[2] * v33[i - 1]
u <- y[i] ^ 2
if (y[i] <= L){
temp1 <- -L*jNormPdf(L/sqrt(v[i]))/(jNormCdf(L/sqrt(v[i]))*sqrt(v[i])^3)
temp2 <- 1.5*L*jNormPdf(L/sqrt(v[i]))/(jNormCdf(L/sqrt(v[i]))*sqrt(v[i])^5)+
0-0.5*L^3/(v[i]^2)*jNormPdf(L/sqrt(v[i]))/(jNormCdf(L/sqrt(v[i]))*sqrt(v[i])^3)+
0-0.5*L^2/(sqrt(v[i])^3)*jNormPdf(L/sqrt(v[i]))^2/(jNormCdf(L/sqrt(v[i]))^2*sqrt(v[i])^3)
}
if (y[i] >= U){
temp1 <- U*jNormPdf(-U/sqrt(v[i]))/(jNormCdf(-U/sqrt(v[i]))*sqrt(v[i])^3)
temp2 <- -1.5*U*jNormPdf(-U/sqrt(v[i]))/(jNormCdf(-U/sqrt(v[i]))*sqrt(v[i])^5)+
0.5*U^3/(v[i]^2)*jNormPdf(-U/sqrt(v[i]))/(jNormCdf(-U/sqrt(v[i]))*sqrt(v[i])^3)+
0-0.5*U^2/(sqrt(v[i])^3)*jNormPdf(-U/sqrt(v[i]))^2/(jNormCdf(-U/sqrt(v[i]))^2*sqrt(v[i])^3)
}
if (y[i] > L & y[i] < U){
temp1 <- -1 / v[i] + u / v[i] ^ 2
temp2 <- 1 / v[i] ^ 2 - 2 * u / v[i] ^ 3
}
ddhamm[1, 1]<- ddhamm[1, 1]+ temp1 * v11[i] + temp2 * v1[i] ^ 2
ddhamm[1, 2]<- ddhamm[1, 2]+ temp1 * v12[i] + temp2 * v1[i] * v2[i]
ddhamm[1, 3]<- ddhamm[1, 3]+ temp1 * v13[i] + temp2 * v1[i] * v3[i]
ddhamm[2, 2]<- ddhamm[2, 2]+ temp1 * v22[i] + temp2 * v2[i] ^ 2
ddhamm[2, 1]<- ddhamm[1, 2]
ddhamm[2, 3]<- ddhamm[2, 3]+ temp1 * v23[i] + temp2 * v2[i] * v3[i]
ddhamm[3, 3]<- ddhamm[3, 3]+ temp1 * v33[i] + temp2 * v3[i] ^ 2
ddhamm[3, 1]<- ddhamm[1, 3]
ddhamm[3, 2]<- ddhamm[2, 3]
}
return(ddhamm)
}
