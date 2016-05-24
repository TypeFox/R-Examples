uh <-
function(x, criterion = inner.prod.max) {
n <- length(x)
sigma <- mad((x[2:n] - x[1:(n-1)])/sqrt(2))
x.buh <- best.unbal.haar(x, criterion)
x.buh.t <- hard.thresh(x.buh, sigma)
x.buh.t.r <- reconstr(x.buh.t)
return(x.buh.t.r)
}

