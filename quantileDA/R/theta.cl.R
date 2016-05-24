theta.cl <-
function(train,test,cl,theta,cl.test=NULL)
{

## train is a matrix of observations of dimension n.train x p
## cl is a vector of length n.train with the classification labels
## test is a matrix of dimension n.test x p
## theta is a value for the quantiles; if theta<-0.5 the median classifiers can be obtained

n.train<-dim(train)[1]
n.test<-dim(test)[1]
p<-dim(train)[2]
levels=unique(cl)
k<-length(levels) ## number of classes
dist.train<-matrix(0,n.train,k)
dist.test<-matrix(0,n.test,k)




for (j in 1:k) {y<-train[cl==levels[j],]
                if (p>1) ym<-apply(y,2,quantile,theta) else ym<-quantile(y,theta)
                dist.test[,j]<-rowSums((theta+((1-2*theta)*(test<t(matrix(ym,p,n.test)))))*abs(test-t(matrix(ym,p,n.test))))
                dist.train[,j]<-rowSums((theta+((1-2*theta)*(train<t(matrix(ym,p,n.train)))))*abs(train-t(matrix(ym,p,n.train))))

}
cl.test2<-levels[apply(dist.test,1,which.min)]
cl.train2<-levels[apply(dist.train,1,which.min)]

quant.train<-misc(cl.train2,cl)
if (!is.null(cl.test)) quant.test<-misc(cl.test2,cl.test) else quant.test<-NULL

out<-list(cl.train=cl.train2,cl.test=cl.test2,me.train=quant.train,me.test=quant.test)
return(out)

}
