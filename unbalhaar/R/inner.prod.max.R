inner.prod.max <-
function(x) {
ipi <- inner.prod.iter(x)
return(med(which(abs(ipi) == max(abs(ipi)))))
}

