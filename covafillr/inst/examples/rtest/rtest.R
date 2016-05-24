library(covafillr)

coord <- as.matrix(expand.grid(seq(-10,10,0.2),seq(-10,10,0.2)))
ftrue <- function(x)sum(x^3) + prod(x)
covObs <- apply(coord,1,function(x)ftrue(x) + rnorm(1,0,0.01))

## Test covafill
cf <- covafill(coord=coord,obs=covObs,h=c(1,1),p=2L)
cf2 <- covafill(coord=coord,obs=covObs)

val1 <- cf$predict(matrix(c(0,0),1,2))
xx <- matrix(runif(200,-10,10),100,2)
val2 <- cf$predict(xx)

## Test covatree

ct <- covatree(coord=coord,obs=covObs,h=c(1,1),p=2L, minLeft = 100)
val3 <- ct$predict(xx)



## Test with dim > 2
coord <- as.matrix(expand.grid(seq(-10,10,0.5),seq(-10,10,0.5), seq(-10,10,0.5)))
ftrue <- function(x)sum(x^3) + prod(x)
covObs <- apply(coord,1,function(x)ftrue(x) + rnorm(1,0,0.01))

cf <- covafill(coord=coord,obs=covObs,h=c(1,1),p=2L)
ct <- covatree(coord=coord,obs=covObs,h=c(1,1),p=2L, minLeft = 1000)

system.time(val1 <- cf$predict(matrix(c(0,0,0),1,3)))
system.time(val3 <- ct$predict(matrix(c(0,0,0),1,3)))

## These should give an error
ct <- covatree(coord=cbind(coord,1),obs=covObs,h=c(1,1,1),p=2L, minLeft = 1000)
ct <- covatree(coord=coord,obs=covObs,h=c(1,1,1),p=0L, minLeft = 1000)


## Cross validate - this takes a lot of time!


bw <- seq(0.5,5,0.5)
cv <- numeric(length(bw))
for(i in 1:length(bw)){
    cf$setBandwith(bw[i])
    res <- cf$residuals(0.1)
    cv[i] <- mean(res ^ 2)
}
    

plot(bw,cv)
