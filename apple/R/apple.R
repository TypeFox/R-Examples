apple<-function(X, y, family='binomial', penalty='LASSO', gamma,  cha.poi=1, eps=1e-15, lam.list, lambda.min.ratio, max.iter=100, max.num, n.lambda=100){
	
	if(penalty=='MCP' & missing(gamma)){
		gamma=1.3
	}
	if(penalty=='LASSO' & !missing(gamma)){
		stop('gamma is not needed')
	}
	
	n=dim(X)[1]
	p=dim(X)[2]
	
	if(missing(max.num)) max.num=p+1
	else max.num=min(max.num, p+1)
	
	if(missing(lambda.min.ratio)) lambda.min.ratio=ifelse(n<p, 0.01, 0.0001)
	
	################prepare variables
	lambda=NULL
	est.beta=matrix(0,nrow=(p+1),ncol=1)
	index=NULL
	y.mean=NULL
	sparse.check=NULL
	delta=NULL
	lik.list=NULL
	ebic.list=NULL
	adjust.step=100
	
	################standardize design matrix
	meanx=apply(X,2,mean)
	normx=sqrt(apply((t(X)-meanx)^2, 1, sum)/n)
	X=scale(X, meanx, normx)
	
	###########add 'one' column to the design matrix
	X=cbind(X,rep(1,n))
	
	#############lambda list
	if(family=='binomial') max.lam=max(abs((1/n)*(t(X[,-(p+1)])%*%(y-1/2))))
	if(family=='poisson') max.lam=max(abs((1/n)*(t(X[,-(p+1)])%*%(y-1))))
	if(missing(lam.list)){
		min.lam=max.lam*lambda.min.ratio
		lam.list=exp(seq(from=log(max.lam), to=log(min.lam), length.out=n.lambda))
	}
	else{
		lam.list=lam.list[lam.list<=max.lam+.Machine$double.eps]
		n.lambda=min(n.lambda, length(lam.list))
	}
	
	
	##############main part
	for(k in 1:(n.lambda-1)){
		#print(k)
		
		#############used to update and correct
		new.beta=est.beta[,k]
		
		#######if the est. model is sparse enough, use apple
		#######otherwise, use CD
		index=(1:p)[new.beta!=0]
		if(sqrt(n)>cha.poi*(length(index)+1)) sparse.check=TRUE
		else sparse.check=FALSE
		
		#######decrease lambda
		delta=lam.list[k+1]-lam.list[k]
				
		#######if it's sparse enough, use our method
		
			
		    break.ind=FALSE
		    
		    ########detect active set
		    eta=exp(-X%*%new.beta)
		    if(family=='poisson') y.mean=1/eta
		    if(family=='binomial') y.mean=1/(1+eta)
		    corr=(1/n)*t(X[,-(p+1)])%*%(y-y.mean)
		    new.index=(1:p)[abs(corr)>=lam.list[k]-1e-14]
		    index=unique(sort(c(index, new.index)))
		    if(sum(index!=0)>max.num) break
		    
		    #########update by linear approx.
		    lik.cov=lik.cov(X=X, y.mean=y.mean, family=family, index=index, n=n, p=p)
            if(penalty=='LASSO'){
		    	sgn=sgn(X=X,y=y,y.mean=y.mean,n=n,p=p)
			    deri=deri(cov=lik.cov, sgn=sgn, index=index)
			    new.beta[index]=new.beta[index]+delta*deri
			}
		    if(penalty=='MCP'){
		    	least.eigen=least.eigen(cov=lik.cov)
		    	sgn=sgn(X=X, y=y, y.mean=y.mean, n=n, p=p, beta=new.beta, gamma=gamma, lambda=lam.list[k])
		    	deri=deri(cov=lik.cov, sgn=sgn, index=index, beta=new.beta, gamma=gamma, lambda=lam.list[k], least.eigen=least.eigen)
		    	new.beta[index]=new.beta[index]+delta*deri
		    	
		    	lik.cov=lik.cov(X=X, y.mean=y.mean, family=family, index=index, n=n, p=p)
		    	least.eigen=least.eigen(cov=lik.cov)
		    	W=W.matrix(X=X, beta=new.beta, family=family)
		    	y.app=y.app(index=index, p=p, beta=new.beta, X=X, y=y, family=family)
		    	
		    }
		    
		    
		    ##########correct in each step
		    if(sparse.check==TRUE){
		    for(j in 1:max.iter){
		    
		        
		        eta=as.vector(exp(-X%*%new.beta))
		        if(sum(eta)==0 | sum(is.na(eta))>0){
		        	break.ind=TRUE
		         	break
		        }
		        if(family=='poisson') y.mean=1/eta
		        if(family=='binomial') y.mean=1/(1+eta)
		        
		        
		        #####newton-raphson
		        if(penalty=='LASSO'){
		            NR.deri=NR.deri(X=X, y=y, lambda=lam.list[k+1], n=n, index=index, y.mean=y.mean)
		            if(max(abs(NR.deri))<eps) break
		            NR.secderi=NR.secderi(X=X, y.mean=y.mean, n=n, p=p, index=index, family=family)	
		            new.beta[c(index,p+1)]=new.beta[c(index,p+1)]-ginv(NR.secderi)%*%NR.deri
		             
		        }
		        if(penalty=='MCP'){
		            sgn=sgn(X=X, y=y, y.mean=y.mean, n=n, p=p, beta=new.beta, gamma=gamma, lambda=lam.list[k+1])
		            mcp.NR.deri=mcp.NR.deri(X=X, y=y, n=n, index=index, beta=new.beta, sgn=sgn, y.app=y.app,lambda=lam.list[k+1],gamma=gamma, W=W, p=p)
		            if(max(abs(mcp.NR.deri))<eps){
		            	break
		            } 
		            mcp.divide=ifelse(j>=adjust.step,runif(1,0.5,1),1) 
		            mcp.NR.secderi=mcp.NR.secderi(X=X, index=index, beta=new.beta, lambda=lam.list[k], gamma=gamma, n=n, W=W,least.eigen = least.eigen, p=p)	
		            new.beta[c(index,p+1)]=new.beta[c(index,p+1)]-ginv(mcp.NR.secderi)%*%mcp.NR.deri*mcp.divide
		        	 
		        }
		        
		    #end correction	
		    }
		    }
		    
		    if(sparse.check==FALSE){
                  lam = lam.list[k+1]
                  break.ind = FALSE
                        
                  for(ite in 1:max.iter){
                      
                      if(penalty=='LASSO'){
			    	      if(ite==1){
			                  cd.temp1=cd_lasso1(X=X, y=y, beta=new.beta, epsilon=eps, max.iter=1, lambda=lam, family=family, bInd=break.ind)		
			    	          new.beta=cd.temp1$beta
			    	          break.ind=cd.temp1$bInd	
		        	          index=(1:(p+1))[new.beta!=0]
		                  }
		                  if(sum(index!=0)>max.num){
		                      break.ind=TRUE
		                      break
		                  } 	
		            
			              #CD
			              if(ite > 1){
			              cd.temp2=cd_lasso1(X=X, y=y, beta=new.beta, epsilon=eps, max.iter=max.iter, lambda=lam, family=family, bInd=break.ind)
			              new.beta=cd.temp2$beta
			              break.ind=cd.temp2$bInd
			        
				          ########redetect active set
			              check.beta=new.beta
			              cd.temp3=cd_lasso1(X=X, y=y, beta=check.beta, epsilon=eps, max.iter=1, lambda=lam, family=family, bInd=break.ind)
			              check.beta=cd.temp3$beta
				          break.ind=cd.temp3$bInd
				     
				          check.index=(1:(p+1))[check.beta!=0]
				          temp.index=unique(sort(c(index, check.index)))
				          if(length(temp.index)==length(index)) break
				          else{
		                      index=temp.index	
		                  } 
			              }
			          #END LASSO 
			          }
			         
			          if(penalty=='MCP'){
				      
				          if(ite==1){
		        	
		        	          ########detect active set
		                      mcp.temp1=cd_mcp1(X=X, y=y, beta=new.beta, epsilon=eps, max.iter=1, lambda=lam, gamma=gamma, family=family, bInd=break.ind)
		                      new.beta=mcp.temp1$beta
		                      break.ind=mcp.temp1$bInd
		                      index=(1:(p+1))[new.beta!=0]
		                  }
		         	
		                  if(sum(index!=0)>max.num){
		                      break.ind=TRUE
		                      break
		                  }
		              
		                  #CD
		                  if(ite > 1){
			              mcp.temp2=cd_mcp1(X=X, y=y, beta=new.beta, epsilon=eps, max.iter=max.iter, lambda=lam, gamma=gamma, family=family, bInd=break.ind)
		                  new.beta=mcp.temp2$beta
		                  break.ind=mcp.temp2$bInd
		            
			              ########redetect active set
			              check.beta=new.beta
			              mcp.temp3=cd_mcp1(X=X, y=y, beta=check.beta, epsilon=eps, max.iter=1, lambda=lam, gamma=gamma, family=family, bInd=break.ind)
		                  check.beta=mcp.temp3$beta
		                  break.ind=mcp.temp3$bInd
		            
		                  check.index=(1:(p+1))[check.beta!=0]
		                  temp.index=unique(sort(c(index, check.index)))
		                  
		                  if(length(temp.index)==length(index)) break
		                  else{
		                      index=temp.index
		                  }
		          	      }
			          #END MCP
			          }
			
			      #end check-cd-recheck
			      }
                  
		    #end case if sparse enough
		    }
		

		
		lik.list[k+1]=loglik(X=X, y=y, beta=new.beta, family=family)
		ebic.list[k+1]=ebic(loglik=lik.list[k+1], beta=new.beta[-(p+1)], n=n, p=p)
		if(sum(new.beta[-(p+1)]!=0)==0){
			new.beta[(p+1)] = f.intercept(y, family)
		}
		
		est.beta=cbind(est.beta, new.beta)
		
		#break criteria
		if(break.ind==TRUE){
			print('1Warning: Algorithm failed to converge for all values of lambda')
			break
		}
		
		if(family=='binomial'){
			prob=1/(1+exp(-X%*%new.beta))
			if(max(prob)>1-1e-10 | min(prob)<1e-10){
				print('2Warning: Algorithm failed to converge for all values of lambda')
				break
			}
		}
		
		if(family=='poisson'){
			current.lambda=exp(X%*%new.beta)
			if(max(current.lambda)>1/eps | min(current.lambda)<eps){
				print('2Warning: Algorithm failed to converge for all values of lambda')
				break
			}
		}
	
	#end the outer loop: decrease in lambda	
	}
	
	ebic.loc=which.min(ebic.list[-1])
	
	est.beta[-(p+1),]=est.beta[-(p+1),]/normx
	est.beta[(p+1),]=est.beta[(p+1),]-crossprod(meanx,est.beta[-(p+1),])
	colnames=c()
	for(i in 1:dim(est.beta)[2]){
		colnames[i] = paste('s',i,sep='')
	}
	colnames(est.beta) = colnames
	
	lambda=lam.list[1:(dim(est.beta)[2]-1)]
	val=list(a0=est.beta[(p+1),-1], beta=est.beta[-(p+1),-1], lambda=lambda, ebic=ebic.list[-1], ebic.loc=ebic.loc, family=family)
	class(val)='apple'
	return(val)
}
