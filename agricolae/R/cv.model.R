`cv.model` <-
function(x) {
suma2 <- sum(x$residual^2)
gl <- x$df.residual
promedio<-mean(x$fitted.values)
return(sqrt(suma2/gl)*100/promedio)
}

