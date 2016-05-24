"vare36" <-
function(x)
 {
n <- x$n
rb <- rbar(x)
T <- sum(n,na.rm=TRUE)
k <- length (x$Rxy[!(is.na(x$Rxy))])
ve <- ((1-rb^2)^2*k)/T
return(ve)
}

