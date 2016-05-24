jNewRapGarchLim <-
function(x0, y, L, U){
b <- x0 #x0 <- c(0.1, 0.88, 0.000001) or c(0.05,0.85,mean(y^2)*0.1)
N <- length(x0)
threshold <- 10 ^ (-10)
maxloop <- 100
delta <- threshold + 1
numloop <- 1
while (abs(delta) > threshold & numloop < maxloop){
numloop <- numloop + 1
maxvalue <- haml(b, y, L, U)
Func <- dhaml(b, y, L, U)
Grad <- ddhaml(b, y, L, U)
Ginv <- solve(Grad)
Inov <- Func %*% Ginv
a <- b
stepp <- 1
if (Inov[1] > 0) {
stepp <- min(stepp, b[1] / Inov[1])
}
if (Inov[2] > 0) {
stepp <- min(stepp, b[2] / Inov[2])
}
if (Inov[3] > 0) {
stepp <- min(stepp, b[3] / Inov[3])
}
if (Inov[1] + Inov[2] < 0) {
stepp = min(stepp, (b[1] + b[2] - 0.99999) / (Inov[1] + Inov[2]))
}
Record <- 0.1
for (i in 1:10){
for (j in 1:N){
b[j] <- a[j] - Inov[j] *(i-1) * stepp / 10
}
tam <- haml(b, y, L, U)
if (tam > maxvalue) {
maxvalue <- tam
Record <- i
}
}
for (j in 1:N){
b[j] <- a[j] - Inov[j] * (Record - 0.1) * stepp / 10
}
delta <- abs(max(Inov)) * Record * stepp / 10
}
return(b)
}
