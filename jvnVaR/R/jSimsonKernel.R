jSimsonKernel <-
function(n,v,v0,s,hh,k){
# compute integral by simson formula
h <- (v0-v)/(2*n)
g <- 0
for (i in 1:n){
   g <- g + h / 3 * (k(s,hh,v) + 4 * k(s,hh,v+h) + k(s,hh,v+h+h))
   v <- v + h + h
}
return(g)
}
