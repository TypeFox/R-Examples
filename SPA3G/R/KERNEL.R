KERNEL <-
function(G, weight)
{
if (length(dim(G))==0)
{
size <- length(G)
k <- matrix(1, size, size)
     for (i in 1 : (size-1))
{ 
j <- seq(1, i, 1)
remain <- G[-j]
Ones <- matrix(1, length(remain), 1)
leading <- Ones*G[i]
D <- abs(remain-leading)
AM <- D
AM[AM==0] <- 4
AM[AM==2] <- 0
AM[AM==1] <- 2
AM[remain==1 & leading==1] <- 2
k[i, (i+1):size]  <- k[(i+1):size, i]<- AM*weight/sum(4*weight)
}
}
if (length(dim(G))>0)
{
      size <- nrow(G)
      k <- matrix(1, size, size)
      for (i in 1 : (size-1))
{ 
j <- seq(1, i, 1)
if (i<(size-1)) 
{remain=as.matrix(G[-j, ])}
if (i==(size-1)) 
{remain <- t(as.matrix(G[-j, ]))}
Ones <- matrix(1, nrow(remain), 1)
leading <- Ones%*%G[i, ]
D <- abs(remain-leading)
AM <- as.matrix(D)
AM[AM==0] <- 4
AM[AM==2] <- 0
AM[AM==1] <- 2
AM[remain==1 & leading==1] <- 2
k[i, (i+1):size]  <- k[(i+1):size, i]<- AM%*%weight/sum(4*weight)
}
}
return(k)
}
