maketall.fcn <-
function(ymtx)
{
ab <- prod(dim(ymtx))
y <- rep(NA,ab)
rows <- rep(NA,ab)
cols <- rep(NA,ab)
a <- nrow(ymtx)
b <- ncol(ymtx)
counter <- 0
for(i in 1:a)
for(j in 1:b)
{
counter <- counter+1
rows[counter] <- i
cols[counter] <- j
y[counter] <- ymtx[i,j]
}
list(y=y,rows=as.factor(rows),cols=as.factor(cols))
}
