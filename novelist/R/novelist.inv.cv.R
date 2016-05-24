######################################################################################
################## NOVELIST estimator for precision matrix   ########################
################## Based on cross validation delta and lambda ########################


novelist.inv.cv<-function(x,data=TRUE,is.cov=TRUE,lambda=seq(0,1,by=0.05),delta=NULL,CV=TRUE,CV.inv.cov=TRUE,Th.method=softt,rep.cv=50)  
  
{

	if (data==FALSE)

    	{
	
		if (CV==TRUE)

		{

			print('Error:Cross validation can only be conducted when data is given.')

			break
		}


		if (is.cov==TRUE)

		{
				
			cov.inv.novel.cv<-novelist.assign.inv(x,lambda,delta,Th.method)$inv.cov.novel
	
			cor.inv.novel.cv<-novelist.assign.inv(x,lambda,delta,Th.method)$inv.cor.novel

					
		}

		if (is.cov==FALSE)

		{

			cov.inv.novel.cv<-NA

			cor.inv.novel.cv<-novelist.assign.inv(x,lambda,delta,Th.method)$inv.cor.novel
					
		}

	}

	if (data==TRUE)

	{
    	 
     		n=dim(x)[1]
  		
		p=dim(x)[2]

     		cov.hat<-cov(x)
  
		cor.hat<-cor(x)

    		if (CV==FALSE)

		{
			
			cov.inv.novel.cv <-novelist.assign.inv(cov.hat,lambda,delta,Th.method)$inv.cov.novel

			cor.inv.novel.cv <-novelist.assign.inv(cov.hat,lambda,delta,Th.method)$inv.cor.novel
			
		}

		if (CV==TRUE)

		{
			
			if (is.null(delta)==FALSE)

			{

				print('Error:no value can be assigned for delta if the method is cross validation.')

				break
			}
	
 			l<-length(lambda) 
  				
			spec.inv<-rep(0,l)
  
  			position.inv<-matrix(0,1,1)
    
   			cov.soft<-Th.method(cor.hat,lambda)
    
  			for (repcv in 1:rep.cv)
    
  			{
    
  				print(repcv)
				
				r.cv<-rep(0,l)

    				sample.cv<-sample(c(1:n),n/2)
    
    				x.a<-x[sample.cv,]
    
    				x.b<-x[-sample.cv,]
    
    				cov.hat.a<-cov(x.a)
    
    				novelist.a<-delta.star(x.a, lambda)
    
    				cov.novelist.a<-novelist.a$covariance.novelist.candidates
    
				cor.novelist.a<-novelist.a$correlation.novelist.candidates

    				cov.hat.b<-cov(x.b)

    				cor.hat.b<-cor(x.b)

    
    				for (i in 1:l)
      
    				{
      
					e<-eigen(cov.novelist.a[,,i])$values
      
      					if(e[p]>0&(e[1]/e[p])<10^6)       
        
      					{
        
        					if (2*p<n)             # inverse for n>2p

        					{
        	
	   						if (CV.inv.cov==TRUE)

	   						{

	      							spec.inv[i]<-norm(solve(cov.novelist.a[,,i])-solve(cov.hat.b),"2") 

	   						}

	   						if (CV.inv.cov==FALSE)

	  						{
           
	      							spec.inv[i]<-norm(solve(cor.novelist.a[,,i])-solve(cor.hat.b),"2") 

	   						}
        
        					}
        
       						if (2*p>=n)   # precision for n<2p
       	
       						{
       
	   						if (CV.inv.cov==TRUE)

	   						{

	      							spec.inv[i]<-norm(solve(cov.novelist.a[,,i])%*%cov.hat.b-diag(p),"2") 

	   						}

	   						if (CV.inv.cov==FALSE)

	   						{
           
	      							spec.inv[i]<-norm(solve(cor.novelist.a[,,i])%*%cor.hat.b-diag(p),"2") 

	   						}

     						}
       
        					r.cv[i]<-r.cv[i]+1
        
      					}
      
				}
    
				position.inv<-rbind(position.inv,as.matrix(which((spec.inv/r.cv)==min((spec.inv/r.cv)[!is.na(spec.inv/r.cv)]))))    
  			
			}
  
  			position.inv<-position.inv[-1]
  
  			m.po<-round(mean(position.inv))
    
  			delta.cv<-delta.star(x,lambda)$delta.star[,2][m.po]  
  
  			cor.inv.novel.cv<-solve((1-delta.cv)*cor.hat+delta.cv*(cov.soft[,,m.po]))
  
 			cov.inv.novel.cv<-solve(sqrt(diag(diag(cov.hat)))%*%cor.inv.novel.cv%*%sqrt(diag(diag(cov.hat))))
  	
		}

	}

	if(CV==TRUE)

	{

		return(list('inv.cov.novel'=cov.inv.novel.cv,'inv.cor.novel'= cor.inv.novel.cv, 'lambda.star'=lambda[m.po],'delta.star'=delta.cv)) 

	}

	if (CV==FALSE)

	{

		return(list('inv.cov.novel'=cov.inv.novel.cv,'inv.cor.novel'= cor.inv.novel.cv)) 


	}

}
 