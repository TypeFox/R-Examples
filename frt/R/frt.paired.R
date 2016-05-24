frt.paired <-
function(x, y, alternative="two.sided"){
if (length(x) != length(y))
return("x and y should have the same length")
d <- y - x
n <- length(d)
alld <- NA
top <- 2^n - 1
for (i in 0:top){
v1 <- bin(i)
v2 <- c(v1,array(0,dim=(n-length(v1))))
alld[i+1] <- 0
for (j in 1:n){
alld[i+1] <- alld[i+1] + ifelse(v2[j]==1, d[j], -d[j])
}
alld[i+1] <- alld[i+1] / n
}
md <- mean(d)
pl <- length(alld[alld > md])
pg <- length(alld[alld < md])
if (substring(alternative,1,1) == "l") return(pl/2^n)
if (substring(alternative,1,1) == "g") return(pg / 2^n)
return(ifelse(pl < pg, 2*pl / 2^n, 2*pg / 2^n))
}

