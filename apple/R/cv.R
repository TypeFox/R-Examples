cv.apple<-function(X,y, family='binomial', penalty='LASSO', gamma, K=10, alpha=0, seed=1, cha.poi=1, eps=1e-15, lambda.min.ratio, max.iter=100,  max.num, n.lambda=100){
	n<-dim(X)[1]
	p<-dim(X)[2]
	
	if(penalty=='MCP'){
		g=apple(X=X,y=y, gamma=gamma, lambda.min.ratio=lambda.min.ratio, n.lambda=n.lambda, family=family, max.num=max.num, max.iter=max.iter, penalty=penalty, eps=eps, cha.poi=cha.poi)
		}
		
	if(penalty=='LASSO'){
		g=apple(X=X, y=y, lambda.min.ratio=lambda.min.ratio, n.lambda=n.lambda, family=family, max.num=max.num, max.iter=max.iter, penalty=penalty, eps=eps, cha.poi=cha.poi)
		}
	
	ori.lam<-as.vector(g$lambda)
	#print(length(ori.lam))
	step<-length(ori.lam)
	est.beta<-rbind(g$beta,g$a0)
	
	set.seed(seed)
	di<-sample((1:n),n,replace=F)
	n.fold <- floor(n/K)
	ind<-NULL
	for(j in 1:K){
		ind[[j]]<-(di[((j-1)*n.fold+1):(j*n.fold)])
		}
	ind[[K]] <- di[((K-1)*n.fold+1):n]
	loss1<-matrix(0,nrow=K,ncol=step)
	loss2<-matrix(0,nrow=K,ncol=step)
	loss.list<-matrix(0,nrow=K,ncol=step)
	cv.aaa<-NULL
	
    	
    	for(j in 1:K){
    		cat('fold ',j,' of ',K,'\n')
    		subX<-cbind(X[ind[[j]],],rep(1,length(ind[[j]])))
    		deleteX<-X[-ind[[j]],]
		    suby<-y[ind[[j]]]
		    deletey<-y[-ind[[j]]]
		    if(penalty=='MCP'){
		    	inter.list=apple(X=deleteX, y=deletey,gamma=gamma, lambda.min.ratio=lambda.min.ratio, n.lambda=n.lambda,  family=family,max.num=max.num, max.iter=max.iter, lam.list=ori.lam, penalty=penalty,cha.poi=cha.poi, eps=eps)
		    	}
		    if(penalty=='LASSO'){
		    	inter.list=apple(X=deleteX, y=deletey, lambda.min.ratio=lambda.min.ratio, n.lambda=n.lambda, family=family, max.num=max.num, max.iter=max.iter, lam.list=ori.lam, penalty=penalty,cha.poi=cha.poi, eps=eps)
		    	}
		    
		    inter.list$est.beta=rbind(as.matrix(inter.list$beta),inter.list$a0)
		    cv.bbb=subX%*%(inter.list$est.beta)
		    cv.bb=(exp(-cv.bbb))
		    
		    if(family=='poisson'){
	    		yhat=1/cv.bb
	    		}
	    	if(family=='binomial'){
	    		yhat=1/(1+cv.bb)
	    		}
	    	cv.aaa[j]<-dim(yhat)[2]
	    	loss1[j,(step-cv.aaa[j]+1):step]=loss(y.hat=yhat, y.sub=suby,family=family)
	    	for(ll in (step-cv.aaa[j]+1):step){
  	            loss2[j,ll]<-sum(inter.list$est.beta[-(p+1),ll-(step-cv.aaa[j])]!=0)*log(n.fold)/(n.fold)	  	        }
  	        loss.list[j,(step-cv.aaa[j]+1):step]<-(1-alpha)*loss1[j,(step-cv.aaa[j]+1):step]+alpha*loss2[j,(step-cv.aaa[j]+1):step]
		    
	    #end iteration of computing and analysing the sub samples
	    }
	    
	    cv.aa<-max(cv.aaa)
	    cv<-apply(loss.list,2,function(x){mean(x[x!=0])})
	    cv.loc=which.min(cv[(step-cv.aa+1):step])+step-cv.aa 
	    cv.est.beta=est.beta[,cv.loc]
	    
	    ebic.loc=g$ebic.loc
	    ebic.est.beta=est.beta[,ebic.loc]
	    
	    family=family
	    
	    val=list(cv=cv, cv.loc=cv.loc,  lambda=ori.lam, beta=est.beta[-(p+1),], a0=est.beta[(p+1),], ebic.loc=ebic.loc, cv.beta=cv.est.beta[-(p+1)], ebic.beta=ebic.est.beta[-(p+1)],
cv.a0=cv.est.beta[p+1], ebic.a0=ebic.est.beta[p+1], family=family)

        class(val)='apple'
        
        return(val) 	    	
 	    	
	}


	
	

	