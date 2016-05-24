quantilecldiff <-
function(train,test,cl,theta=NULL,cl.test=NULL)
{

p<-dim(train)[2]
n.train<-dim(train)[1]
n.test<-dim(test)[1]
k<-max(as.double(as.factor(cl))) ## number of classes
theta.selected<-rep(0,p)

if (is.null(theta)) {theta<-seq(0,1,0.05)
                    theta<-theta[-length(theta)]
                    theta<-theta[-1]}

ymt<-matrix(0,p,k)

for (h in 1:p) {

train.rates<-rep(0,length(theta))

for (j in 1:length(theta)) {

train.h<-as.matrix(train[,h])

dist.train<-matrix(0,n.train,k)
for (i in 1:k) {y<-train.h[cl==i,]
                ym<-quantile(y,theta[j])
                dist.train[,i]<-rowSums((theta[j]+((1-2*theta[j])*(train.h<ym)))*abs(train.h-ym))
}

cl.train2<-apply(dist.train,1,which.min)
train.rates[j]<-misc(cl,cl.train2)       
}           
theta.choice<-theta[which.min(train.rates)]

theta.selected[h]<-theta.choice

ymt[h,]<- as.double(by(train[,h],cl,quantile,theta.selected[h]))

}

dist.train<-rep(0,k)
dist.test<-rep(0,k)
cl.test2<-rep(0,n.test)
cl.train2<-rep(0,n.train)


for (i in 1:n.test) { for (j in 1:k) dist.test[j]<-sum((theta.selected+((1-2*theta.selected)*(test[i,]<ymt[,j])))*abs(test[i,]-ymt[,j]))
                    cl.test2[i]<-which.min(dist.test)
                    }

for (i in 1:n.train) {for (j in 1:k)  dist.train[j]<-sum((theta.selected+((1-2*theta.selected)*(train[i,]<ymt[,j])))*abs(train[i,]-ymt[,j]))
                     cl.train2[i]<-which.min(dist.train)
                     }   


theta.choice<-mean(theta.selected)
if (!is.null(cl.test)) me.test<-misc(cl.test,cl.test2) else me.test<-NULL
me.train<-misc(cl,cl.train2)

out<-list(thetas=theta.selected,theta.choice=theta.choice,me.train=me.train,me.test=me.test,cl.train=cl.train2,cl.test=cl.test2)
return(out)
}
