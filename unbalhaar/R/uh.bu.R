uh.bu <-
function(x, stretch = length(x)) {
sigma <- mad(diff(x)/sqrt(2))
x.d <- best.unbal.haar.bu(x, stretch)
x.d.t <- hard.thresh.bu(x.d, sigma)
x.d.t.r <- reconstr.bu(x.d.t)
return(x.d.t.r)
}

