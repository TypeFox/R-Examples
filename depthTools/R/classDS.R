classDS <- function(xl, yl, xt, alpha=0.2)
{
yl<-as.matrix(yl)
ylf<-as.factor(yl)
classes<-levels(ylf)
no.classes <- length(classes)
xl<-as.matrix(xl);xt<-as.matrix(xt)
n <- c(); m<-nrow(xt)
distances <- c(); pred<-c()
train.t.average <- list()

for (i in 1:no.classes)
   {
     n[i] <- sum(ylf==classes[i])      # size of each learning group
     train<-xl[ylf==classes[i],]
     train.t.average[[i]] <- tmean(train,alpha)$tm      # trimmed mean of each reference group
   }   ### end FOR i

for (j in 1:m)
  {
    for (k in 1:no.classes)
      {
        z <- xt[j,]-train.t.average[[k]]
        distances[k] <- t(z)%*%z
      }  ### end FOR k
    pred[j] <- as.numeric(classes[which.min(distances)])
  }   ### end FOR j

return(pred)
}


