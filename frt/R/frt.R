frt <-
function(x, y, alternative="two.sided"){
yobs <- mean(y) - mean(x)
y1 <- NA
allObs <- c(x,y)
lx <- length(x)
ly <- length(y)
k <- min(lx,ly)
for(j in 0:k){
rnd <- concat(comb(lx-j,j), comb(j,ly-j))
rnd <- matrix(rnd, ncol=lx+ly)
for (i in 1:dim(rnd)[1])
y1[length(y1)+1] <- mean(allObs[rnd[i,]==0]) - mean(allObs[rnd[i,]==1])
}
y1 <- y1[-1]
n <- factorial(lx + ly) / (factorial(lx)*factorial(ly))
pl <- length(y1[y1 > yobs])
pg <- length(y1[y1 < yobs])
if (substring(alternative,1,1) == "l") return(pl / n) 
if (substring(alternative,1,1) == "g") return(pg / n)
return(ifelse(pl < pg, 2*pl/n, 2*pg/n))
}

