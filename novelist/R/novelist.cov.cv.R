######################################################################################
################## Novelist estimator for covariance matrix   ########################
################## Based on cross validation delta and lambda ########################

#' @export novelist.cov.cv

novelist.cov.cv<-function(x,data=TRUE,is.cov=TRUE,lambda=seq(0,1,by=0.05),delta=NULL,CV=TRUE,CV.cov=TRUE,Th.method=softt,rep.cv=50)  
  
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

			cov.novel.cv<-novelist.assign(x,lambda,delta,Th.method)$cov.novel
	
			cor.novel.cv<-novelist.assign(x,lambda,delta,Th.method)$cor.novel

					
		}

		if (is.cov==FALSE)

		{

			cov.novel.cv<-NA

			cor.novel.cv<-novelist.assign(x,lambda,delta,Th.method)$cor.novel
					
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
			
			cov.novel.cv<-novelist.assign(cov.hat,lambda,delta,Th.method)$cov.novel

			cor.novel.cv<-novelist.assign(cov.hat,lambda,delta,Th.method)$cor.novel
			
		}

		if (CV==TRUE)

		{
			
			if (is.null(delta)==FALSE)

			{

				print('Error:no value can be assigned for delta if the method is cross validation.')

				break
			}
	
 			l<-length(lambda) 
  				
			spec<-rep(0,l)
  
  			position<-matrix(0,1,1)
    
   			cov.soft<-Th.method(cor.hat,lambda)
    
  			for (repcv in 1:rep.cv)
    
  			{
    
  				print(repcv)
    
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
      
					if (CV.cov==TRUE)

					{

						spec[i]<-norm(cov.novelist.a[,,i]-cov.hat.b,"2") 

					}

					if (CV.cov==FALSE)

					{

						spec[i]<-norm(cor.novelist.a[,,i]-cor.hat.b,"2") 

					}
  
      
    				}
    
    					position<-rbind(position,as.matrix(which(spec==min(spec))))
    
  			}
  
  			position<-position[-1]
  
  			m.po<-round(mean(position))
    
  			delta.cv<-delta.star(x,lambda)$delta.star[,2][m.po]  
  
  			cor.novel.cv<-(1-delta.cv)*cor.hat+delta.cv*(cov.soft[,,m.po])
  
 			cov.novel.cv=sqrt(diag(diag(cov.hat)))%*%cor.novel.cv%*%sqrt(diag(diag(cov.hat)))
  	
		}

	}

	if(CV==TRUE)

	{

		return(list('cov.novel'=cov.novel.cv,'cor.novel'=cor.novel.cv, 'lambda.star'=lambda[m.po],'delta.star'=delta.cv)) 

	}

	if (CV==FALSE)

	{

 		return(list('cov.novel'=cov.novel.cv,'cor.novel'=cor.novel.cv)) 

	}

}
 