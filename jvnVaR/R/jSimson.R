jSimson <-
function(n,v,v0,k){# compute integral by simson formula
h <- (v0-v)/(2*n)
g <- 0
for (i in 1:n){
   g <- g + h / 3 * (k(v) + 4 * k(v+h) + k(v+h+h))
   v <- v + h + h
}
return(g)
}
