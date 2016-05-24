"SERHET" <-
function (x)
{
 N <- sum(x$n,na.rm=TRUE)
rb <- rbar(x)
k <- length (x$Rxy[!(is.na(x$Rxy))])
se <- sqrt((1-rb^2)^2/(N-k)+varRes(x)/k)
return(se)
}

