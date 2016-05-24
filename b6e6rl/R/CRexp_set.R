CRexp_set <- function(p,d){
y <- matrix(0,1,d+1)
y[1] <- 1
y[d] <- -d*p
y[d+1] <- d*p-1
r <- polyroot(rev(y))		
ri <- Im(r)
zrus <- zrus <- which(ri < -1e-7 | ri > 1e-7)
r <- r[-zrus]
r <- Re(r)
zrus <- which(r<0 | r>=1)
if(is.empty(zrus)){
r <- sort(r)
CR <- r[1]
}else{
r <- r[-zrus]
r <- sort(r)
CR <- r[1]
}
if(is.empty(r)){
CR <- runif(1)
}
return(CR)
}