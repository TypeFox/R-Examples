entropy.joint <-function(pmXY,q,iicc){

HXY<-sapply(pmXY,function(x){
		sapply(x,function(y){
			   y<-matrix(data=y,nrow=16 ,ncol=1)
			   h<-switch(iicc$classentropy, "Shannon"=.C("entropyShannon",prob=as.double(y),lengthprob=as.double(length(y)),H=double(1)),"Renyi"=.C("entropyRenyi",prob=as.double(y),lengthprob=as.double(length(y)),q=as.double(q),H=double(1)))
               h$H
			   })
	})
	return(HXY)
}

