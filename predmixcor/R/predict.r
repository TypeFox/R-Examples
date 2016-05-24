Predict.mix <- function(test,spl,start)
{
	iter <- length(spl$r )
	phitimesr <- spl$phi[start:iter,]
	
	phitimesr[,1] <- phitimesr[,1] * spl$r1[start:iter]
	phitimesr[,-1] <-  phitimesr[,-1] * spl$r[start:iter]
        rplusN2 <- rplusN1 <- matrix(0,iter-start+1,ncol(spl$phi))
	rplusN1[,1] <- spl$r1[start:iter] + spl$N1[start:iter]
	rplusN1[,-1] <- spl$r[start:iter] + spl$N1[start:iter]
	rplusN2[,1] <- spl$r1[start:iter] + spl$N2[start:iter]
	rplusN2[,-1] <- spl$r[start:iter] + spl$N2[start:iter]
	theta1 <- (phitimesr + spl$I1[start:iter,]) / rplusN1
	theta2 <- (phitimesr + spl$I2[start:iter,]) / rplusN2
	
	PC = ( spl$N1[start:iter] + 1 ) / 
	     ( spl$N1[start:iter] + spl$N2[start:iter] + 2 )	
	pred = rep( 0,nrow(test) )
	
	for( i in 1:nrow(test) )
	{
		den <- num <- rep( 0, iter - start + 1 )
		
		for( t in 1:length( num ) )
		{
		   num[t] = LOG.add.exp(	   
	             	sum( dbinom( c(1,test[i,-1]),1,theta1[t,],log=TRUE) ) +
		     	log( PC[t] ) ,
                     	sum( dbinom( c(1,test[i,-1]),1,theta2[t,],log=TRUE) ) +
                     	log( 1-PC[t] )
			)
		   den[t] = LOG.add.exp(	   
	             	sum( dbinom( c(0,test[i,-1]),1,theta1[t,],log=TRUE) ) +
		     	log( PC[t] ) ,
                     	sum( dbinom( c(0,test[i,-1]),1,theta2[t,],log=TRUE) ) +
                     	log( 1-PC[t] )
			)			
		}
		
		pred[i] = 1/(1 + exp( LOG.sum.exp(den) - LOG.sum.exp(num)))
		
	}
		
	pred	
}
