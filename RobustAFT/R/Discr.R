"Discr" <-
function(yo,Po,cu){
# minimum ratio between Fn and Po (for argument > cu)
n  <- length(yo)
# if (cu > yo[n]) {cat("cu is too large ","\n"); return(list(j=NA))}
Fn <- (0:(n-1))/n; r  <- Fn/Po
j  <- (1:n)[yo > cu]
if (length(j)>0) {j <- j[1]; rj <- r[j:n]} else {j <- n+1; rj <- 1}
r.min <- min(rj,1)
if (r.min < 1) j.min <- (j:n)[rj==r.min] else j.min <- NULL
list(j=j,j.min=j.min,r.min=r.min)}

