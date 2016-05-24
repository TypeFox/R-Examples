"AveS2Pw" <-
function(v,args){
Beta <- args$Beta; n <- args$n
nr <- length(args$rr); nv <- length(v); A <- rep(0,nv)
for (i in 1:nv) {
rs   <- args$rr/v[i]
tmp  <- rs*exp(rs)*(1+rs)-rs
A[i] <- -(sum(tmp)/v[i])/n}; A}

