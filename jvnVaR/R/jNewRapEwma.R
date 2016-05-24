jNewRapEwma <-
function(x0, y){
b <- x0 #x0 <- 0.85
threshold <- 10 ^ (-10)
maxloop <- 100
delta <- threshold + 1
numloop <- 1
while (abs(delta) > threshold & numloop < maxloop){
numloop <- numloop + 1
maxvalue <- hame(b, y)
Func <- dhame(b, y)
Grad <- ddhame(b, y)
Ginv <- Grad^(-1)
Inov <- Func * Ginv
a <- b
stepp <- 1
if (Inov > 0) {# dk lambda >0
stepp <- min(stepp, b / Inov)
}
if (Inov < 0) {# dk lambda <1 
stepp = min(stepp, (b - 0.99999) / Inov)
}
Record <- 0.1
for (i in 1:10){
b <- a - Inov *(i-1) * stepp / 10
tam <- hame(b, y)
if (tam > maxvalue) {
maxvalue <- tam
Record <- i
}
}
b <- a - Inov * (Record - 0.1) * stepp / 10
delta <- abs(Inov) * Record * stepp / 10
}
return(b)
}
