quantileCV<-function(x,cl,nfold=min(table(cl)),folds=balanced.folds(cl,nfold),theta=NULL,seed=1,varying=FALSE)
{
	
	set.seed(seed)
	
	p<-ncol(x)
	
	c.median<-rep(0,nfold)
	c.centroid<-rep(0,nfold)
	c.quantile<-rep(0,nfold)
	c.quantile2<-rep(0,nfold) ##train
	
	if (is.null(theta)) {theta<-seq(0,1,0.02)
		theta<-theta[-length(theta)]
		theta<-theta[-1]}

	c.quantile.train<-matrix(0,length(theta),nfold)
	c.quantile.test<-matrix(0,length(theta),nfold)
	c.quantile.theta<-rep(0,nfold)
	
	for (h in 1:nfold) {
		
        print(paste(round(h/nfold*100,0),"%",sep=""))
        
        test<- x[folds[[h]],] 
        train<- x[-folds[[h]],]    
        cl.test<-cl[folds[[h]]]
        cl.train<-cl[-folds[[h]]]
        
# centroid classifier
        out.c<-centroidcl(train,test,cl.train,cl.test)
        c.centroid[h]<-out.c$me.test
		
		
# median classifiers
        out.m<-theta.cl(train,test,cl.train,0.5,cl.test)
        c.median[h]<-out.m$me.test
		
        
# quantile classifiers
        
        if (varying) out.q<-quantilecldiff(train,test,cl.train,theta,cl.test) else out.q<-quantilecl(train,test,cl.train,theta,cl.test)
        c.quantile[h]<-out.q$me.test
        c.quantile2[h]<-out.q$me.train
        if (!varying) c.quantile.test[,h]<-out.q$test.rates
        if (!varying) c.quantile.train[,h]<-out.q$train.rates
        c.quantile.theta[h]<-out.q$theta.choice
		
	}
	
	if (!varying) test.rates<-rowMeans(c.quantile.test) else test.rates<-NULL
	if (!varying) train.rates<-rowMeans(c.quantile.train) else train.rates<-NULL
	
	
	out<-list(folds=folds,test.rates=test.rates,train.rates=train.rates,thetas=theta,theta.choice=c.quantile.theta,me.test=c.quantile,me.train=c.quantile2,me.median=c.median,me.centroid=c.centroid)
#if (!varying) class(out)<-"quantileDA" else class(out)<-"quantiled"
	return(out)
	
}

