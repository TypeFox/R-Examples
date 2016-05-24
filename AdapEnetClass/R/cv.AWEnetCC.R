cv.AWEnetCC <-
function(X, Y, delta, weight, kFold = 10, C, s, lambda2, AEnetCC=T)
{
n <- nrow(X) # number of samples
p <- ncol(X) # number of predictors
if(n != length(delta) || n != length(Y))
stop("dimensions of X, Y and censorship don't match!")
aft.cv.censorkfold <- function(n, delta, kfold, maxTry = 10)
{
for (i in 1:maxTry) # try at most maxTry times
{
cvFold <- split(sample(1:n), rep(1:kfold, length = n))
kfoldsuccess <- T
for (k in 1:kfold)
{
if (sum(delta[cvFold[[k]]] == 1) <= 1)
{
kfoldsuccess <- F
break
}
}
if (kfoldsuccess)
break
}
if(!kfoldsuccess)

stop(message="Too small number of samples in the given number of
folds! Try decreasing number of folds!")
return(cvFold)
}

cvFold <- aft.cv.censorkfold(n, delta, kFold)
Betas <- numeric(p)
cvscore <- numeric(n)
for (i in 1:kFold)
{
foldIndex <- cvFold[[i]]
if(AEnetCC)
beta<-AEnetCC.aft(as.matrix(X[-foldIndex, ]), Y[-foldIndex], delta[-foldIndex], weight, C, s, lambda2)
else 
beta<-WEnetCC.aft(as.matrix(X[-foldIndex, ]), Y[-foldIndex], delta[-foldIndex], weight, C, s, lambda2)

Betas <- cbind(Betas, beta)
cvscore[foldIndex] <- Y[foldIndex] - as.matrix(X[foldIndex, ]) %*% beta
}
kw <- aft.kmweight(Y, delta)$kmwts * delta
cvscore <- sum((cvscore^2) * kw) / 2
list(beta = Re(apply(Betas[,], 1, mean)), betavar =
Re(apply(Betas[,], 1, var)), cvscore = Re(cvscore) )
}
