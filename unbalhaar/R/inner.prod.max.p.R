inner.prod.max.p <-
function(x, p = 0.8) {
ipi <- inner.prod.iter(x)
n <- length(ipi)
pip <- ipi[(1+floor((1-p)*n)):ceiling(p*n)]
return(med(which(abs(pip) == max(abs(pip))))+floor((1-p)*n))
}

