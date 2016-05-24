quantilecl<-function(train,test,cl,theta=NULL,cl.test=NULL,skew.correct="Galton")
{
  p<-ncol(train)
  if (p==1) skew.correct="none" 
  if (skew.correct=="Galton") {  
                                 n.train<-nrow(train)
                                 n.test<-nrow(test)
                                 levels<-unique(cl)
                                 k<-length(levels) ## number of classes
                                 
                                 skew<-matrix(0,k,p)
                                 for (i in 1:k) {skew[i,]<-apply(train[cl==levels[i],,drop=FALSE],2,galtonskew)
                                 }
                                 total.skew<-colSums(skew)
                                 total.skew<-ifelse(is.na(total.skew),1,total.skew)
                                 train<-train*t(matrix(sign(total.skew),p,n.train))
                                 test<-test*t(matrix(sign(total.skew),p,n.test))
  }
  
  if (skew.correct=="Kelley") { 
                                 n.train<-nrow(train)
                                 n.test<-nrow(test)
                                 levels<-unique(cl)
                                 k<-length(levels) ## number of classes
                                 
                                 skew<-matrix(0,k,p)
                                 for (i in 1:k) {skew[i,]<-apply(train[cl==levels[i],,drop=FALSE],2,kelleyskew)
                                 }
                                 total.skew<-colSums(skew)
                                 total.skew<-ifelse(is.na(total.skew),1,total.skew)                 
                                 train<-train*t(matrix(sign(total.skew),p,n.train))
                                 test<-test*t(matrix(sign(total.skew),p,n.test))
  }
  
  
  if (skew.correct=="Skewness") {  
                                   n.train<-nrow(train)
                                   n.test<-nrow(test)
                                   levels<-unique(cl)
                                   k<-length(levels) ## number of classes
                                   
                                   skew<-matrix(0,k,p)
                                   for (i in 1:k) {skew[i,]<-apply(train[cl==levels[i],,drop=FALSE],2,skewness)
                                                   skew[i,]<-sign(skew[i,])*(abs(skew[i,])^(1/3))
                                   }
                                   total.skew<-colSums(skew)
                                   total.skew<-ifelse(is.na(total.skew),1,total.skew)  
                                   train<-train*t(matrix(sign(total.skew),p,n.train))
                                   test<-test*t(matrix(sign(total.skew),p,n.test))
  }
  
  
	if (is.null(theta)) {theta<-seq(0,1,0.02)
		theta<-theta[-length(theta)]
		theta<-theta[-1]}
	
	test.rates<-rep(0,length(theta))
	train.rates<-rep(0,length(theta))
	
	
	for (j in 1:length(theta)) {
		
		out<-theta.cl(train,test,cl,theta[j],cl.test=cl.test)
		if (!is.null(cl.test)) test.rates[j]<-misc(cl.test,out$cl.test)    
		train.rates[j]<-misc(cl,out$cl.train)       
	}     
  
  minimi<-(train.rates==min(train.rates))
  count<-sum(minimi)
  if (count>1) index.theta.choice<-as.double(names(which.min((predict(lm(train.rates~theta+I(theta^2))))[minimi]))) else index.theta.choice<-which.min(train.rates)
  
  
	theta.choice<-theta[which.min(train.rates)]
	me.test<-test.rates[which.min(train.rates)]
	me.train<-train.rates[which.min(train.rates)]
	out.q<-theta.cl(train,test,cl,theta.choice,cl.test)	
	out<-list(train.rates=train.rates,test.rates=test.rates,thetas=theta,
            theta.choice=theta.choice,me.train=me.train,me.test=me.test,
            train=train,test=test,cl.train=out.q$cl.train,cl.test=out.q$cl.test,
            cl.train.0=cl,cl.test.0=cl.test)
	class(out)<-"quantileDA"
	return(out)
}