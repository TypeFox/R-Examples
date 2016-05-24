envelopes <-
function(iterations=100,countdata,dimention){



	## Randobly arrange the counts;



	GSpermprev<-array();

	for (i in 1:iterations){
	
		perm<-matrix(sample(as.vector(countdata),size=2^(dimention*2),replace=FALSE),nrow=2^dimention,byrow=TRUE);
		GSperm<-iterate(perm,dimention);

		
		GSpermprev<-cbind(GSpermprev,GSperm[,3])
	


	}


	ret.val<-cbind(apply(GSpermprev[,-1],1,quantile,probs=c(.05)),apply(GSpermprev[,-1],1,quantile,probs=c(.95)))


}

