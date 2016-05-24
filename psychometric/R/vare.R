"vare" <-
function(x)
 {
n <- x$n
rb <- rbar(x)
ve <- sum(n*(1-rb^2)^2/(n-1),na.rm=TRUE)/sum(n,na.rm=TRUE)
return(ve)
}

