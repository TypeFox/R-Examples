"empty.solution" <- function(pi, a, a0, lambda = lambda, x = x, y = y, K = K, Left = Left, Right = Right)
{
   index.plus <- which(y==1)
   index.minus <- which(y==-1)

   lplus <- intersect(Left, index.plus)
   lminus <- intersect(Left, index.minus)
   rplus <- intersect(Right, index.plus)
   rminus <- intersect(Right, index.minus)

   candidate <- union(rplus, lminus)
   nl <- length(Left)
   pi <- length(lplus)/nl

   a <- NULL
   a[lplus] <- 1-pi
   a[lminus] <- pi
   a[Right] <- 0
   Kscript <- outer(y,y) * K
   KscriptA <- Kscript %*% a

   can.a0 <- y[candidate] * (lambda - KscriptA[candidate])
   Elbow <- candidate[which.max(can.a0)]

   Left <- setdiff(Left, Elbow)
   Right <- setdiff(Right,Elbow)
   a0 <- max(can.a0)

obj = list(pi = pi, alpha = a, alpha0 = a0, Elbow = Elbow, Left = Left, Right = Right)    
obj
}
