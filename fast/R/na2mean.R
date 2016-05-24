`na2mean` <-
function(x){
x1 <- c(x[-1], NA,NA)
x2 <- c(NA,x)
mean <- (x1 + x2) / 2
x[is.na(x)]<- mean[is.na(x)]
return(x)
}

