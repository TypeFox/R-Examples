tft <-
function(x,G,pig,dhfgs78,p,mug,sigmainv,n,sigma,univar,delta,gauss){
	log_num <- matrix(0,n,G)
	if(gauss){
		if(univar){
			for(g in 1:G){
				log_num[,g] <- log(pig[g])-(p/2)*log(2*pi)-(1/2)*log(sigma[,,g])-
				(1/2)*delta[,g] 
			}
		}
		else{
			for(g in 1:G){
				log_num[,g] <- log(pig[g])-(p/2)*log(2*pi)-(1/2)*log(det(sigma[,,g]))-
				(1/2)*delta[,g] 
			}
		}
	}
	else{
		if(univar){
			for(g in 1:G){
				log_num[,g]<-log(pig[g])+lgamma((dhfgs78[g]+1)/2)-(1/2)*log(sigma[,,g])-
				((p/2)*(log(pi)+log(dhfgs78[g]))+lgamma(dhfgs78[g]/2)+((dhfgs78[g]+p)/2)*
				 (log(1+ delta[,g]/dhfgs78[g])))
			}
		}
		else{
			for(g in 1:G){
				log_num[,g]<-log(pig[g])+lgamma((dhfgs78[g]+p)/2)-(1/2)*log(det(sigma[,,g]))-
				((p/2)*(log(pi)+log(dhfgs78[g]))+lgamma(dhfgs78[g]/2)+((dhfgs78[g]+p)/2)*(log(1+
				delta[,g]/dhfgs78[g])))
			}
		}
	}
	log_num
}

