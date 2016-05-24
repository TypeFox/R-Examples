classTAD <- function(xl, yl, xt, alpha=0)
{

yl<-as.matrix(yl)
ylf<-as.factor(yl)
classes<-levels(ylf)
no.classes <- length(classes)
xl<-as.matrix(xl);xt<-as.matrix(xt)
n <- c()
distances <- c(); pred<-c()

train.average <- list()
p <- list()
D <- list()
I <- list()
xi <- list()
for (i in 1:no.classes)
  {
     n[i] <- sum(ylf==classes[i])                # size of each learning group
     xi[[i]] <- xl[ylf==classes[i],]             # i-th group of the learning set
     train.average[[i]] <- apply(xi[[i]],2,mean) # Mean of the learning groups
     p[[i]] <- MBD(xi[[i]],plotting=FALSE)$MBD                  # MBD of elements in the i-th learning group
     D[[i]] <- sort(p[[i]])
     I[[i]] <- order(p[[i]])
  }  ### end FOR i

mmin <- min(n);                                  # Minimum size of the learning groups
mmin <- mmin-floor(alpha*mmin);                  # leave out a proportion of alfa*mmin points

trimmed.train <- list()
p<-list(); pm<-list(); 
distances <- c(); pred<-c()

for (i in 1:no.classes)
  {
   trimmed.train[[i]] <- xi[[i]][I[[i]][(n[i]-mmin+1):(n[i])],]        # Consider the mmin deepest points from training groups
   p[[i]]<- D[[i]] [(n[i]-mmin+1):(n[i])]
   pm[[i]] <- (p[[i]]/sum(p[[i]]))                                     #  Weights needed to compute TAD 
  }  ### end FOR i

for (j in 1:nrow(xt))
 {
   x <- xt[j,]
   for (k in 1:no.classes)
     {
       m <-nrow(trimmed.train[[k]])
       B <- (matrix(1,m,1)%*%x)-trimmed.train[[k]];
       dL <- apply(t(abs(B)),2,sum)
       D <- diag(as.vector(pm[[k]]))
       distances[k] <- sum(dL%*%D)

     }  ### end FOR k

   pred[j] <- as.numeric(classes[which.min(distances)])
 }  ### end FOR j

return(pred)
}
