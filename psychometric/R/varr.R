"varr" <-
function(x)
 {
rxy <- x$Rxy
n <- x$n
rb <- rbar(x)
vr <- sum(n*(rxy-rb)^2,na.rm=TRUE)/sum(n,na.rm=TRUE)
return(vr)
}

