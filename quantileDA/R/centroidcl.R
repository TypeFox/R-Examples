centroidcl <-
function(train,test,cl,cl.test=NULL)
{

## train is a matrix of observations of dimension n.train x p
## cl is a vector of length n with the classification labels (factor or numerical values are allowed)
## test is a matrix of dimension n.test x p
## cl.test (optional; if it is given the misclassification error in the test set is reported)

n.test<-dim(test)[1]
n.train<-dim(train)[1]
p<-dim(test)[2]
levels=unique(cl)
k<-length(levels) ## number of classes
dist.test<-matrix(0,n.test,k)
dist.train<-matrix(0,n.train,k)


for (j in 1:k) { if (p>1) medie<-colMeans(train[cl==levels[j],]) else medie<-mean(train[cl==levels[j],])
                dist.test[,j]<- rowSums((test-t(matrix(medie,p,n.test)))^2)
                dist.train[,j]<- rowSums((train-t(matrix(medie,p,n.train)))^2)
                }
cl.test2<-levels[apply(dist.test,1,which.min)]
cl.train2<-levels[apply(dist.train,1,which.min)]

centroid.train<-misc(cl.train2,cl)
if (!is.null(cl.test)) centroid.test<-misc(cl.test2,cl.test) else centroid.test<-NULL

out<-list(cl.train=cl.train2,cl.test=cl.test2,me.train=centroid.train,me.test=centroid.test)
return(out)
}
