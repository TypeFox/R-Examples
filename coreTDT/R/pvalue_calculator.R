##comput the LRT for all x1,x2 combination

loglr_comp = function(x1,x2,n1,n2,theta1 = 0.25,theta2 = 0.5){
	
	lr1 = x1*log(theta1) +(n1-x1)*log(1-theta1) + x2*log(theta2) +(n2-x2)*log(1-theta2)
	if(x1 >= (theta1*n1) & x2>=(theta2*n2)){
		
		theta1_mle = ifelse(n1==0,0,x1/n1)
		theta2_mle = ifelse(n2==0,0,x2/n2)
	}else if(x1 >= (theta1*n1) & x2 < (theta2*n2)){
		theta1_mle = ifelse(n1==0,0,x1/n1)
		theta2_mle = theta2
	}else if(x1 < (theta1*n1) & x2 >= (theta2*n2)){
		theta1_mle = theta1
		theta2_mle = ifelse(n2==0,0,x2/n2)
					
	}else{		
		theta1_mle = theta1
		theta2_mle = theta2
	}	 
	
	if(theta1_mle==0 | theta1_mle==1 | theta2_mle ==0| theta2_mle==1){
		lr2 = log(theta1_mle^x1*(1-theta1_mle)^(n1-x1)*theta2_mle^x2*(1-theta2_mle)^(n2-x2) )
		lr = lr1 - lr2
	}else{
		lr2 = x1*log(theta1_mle) +(n1-x1)*log(1-theta1_mle) + x2*log(theta2_mle) +(n2-x2)*log(1-theta2_mle)	
		lr = lr1 - lr2
	}
	
	
	
	return(lr)	
}


loglr_comp_2side = function(x1,x2,n1,n2,theta1 = 0.25,theta2 = 0.5){
	
	lr1 = x1*log(theta1) +(n1-x1)*log(1-theta1) + x2*log(theta2) +(n2-x2)*log(1-theta2)
			
	theta1_mle = ifelse(n1==0,0,x1/n1)
	theta2_mle = ifelse(n2==0,0,x2/n2)
	
	if(theta1_mle==0 | theta1_mle==1 | theta2_mle ==0| theta2_mle==1){
		lr2 = log(theta1_mle^x1*(1-theta1_mle)^(n1-x1)*theta2_mle^x2*(1-theta2_mle)^(n2-x2) )
		lr = lr1 - lr2
		
	}else{
		lr2 = x1*log(theta1_mle) +(n1-x1)*log(1-theta1_mle) + x2*log(theta2_mle) +(n2-x2)*log(1-theta2_mle)	
		lr = lr1 - lr2
	}
	
	
	
	return(lr)	
}



pvalue_calculator = function(y1,y2,n1,n2,theta1 = 0.25,theta2 = 0.5){
	prob_null = array(0,dim=c((n1+1),(n2+1)))
	lrt_all = array(0,dim=c((n1+1),(n2+1)))
	lrt2_all = array(0,dim=c((n1+1),(n2+1)))
	for(i in c(0:n1)){
		for(j in c(0:n2)){
			prob_null[(i+1),(j+1)] = dbinom(i,n1,theta1) *dbinom(j,n2,theta2) 	
			lrt_all[(i+1),(j+1)] = loglr_comp(i,j,n1=n1,n2=n2,theta1=theta1,theta2=theta2)
			lrt2_all[(i+1),(j+1)] = loglr_comp_2side(i,j,n1=n1,n2=n2,theta1=theta1,theta2=theta2)
		
		}
	}
	
	pr = prob_null[(y1+1),(y2+1)]
	lrt = lrt_all[(y1+1),(y2+1)]
	lrt2 = lrt2_all[(y1+1),(y2+1)]
	pvalue_pr = sum(c(prob_null)[which(c(prob_null) <= pr)])
		 	
	pvalue_lr = sum(c(prob_null)[which(c(lrt_all) <= lrt)])
	pvalue_lr2 = sum(c(prob_null)[which(c(lrt2_all) <=lrt2)])

	return(list("pvalue_pr" = pvalue_pr,"pvalue_lr" = pvalue_lr,"pvalue_lr2"= pvalue_lr2))

}
