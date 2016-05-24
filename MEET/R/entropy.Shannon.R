entropy.Shannon <-function(wind){
	
	H<-apply(wind,2,function(y){
        n<-length(y)
        out<-.C("entropyShannon",prob=as.double(y),lengthprob=as.double(n), H=double(1))
        out$H
             })
		return(H)			 
	}
