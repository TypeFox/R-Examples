tab.summary <-
function(res,daic=2,show.rate=FALSE){
	Ks<-unique(res[,1])
	
		temps<-res[which(res[,1]==Ks[1]),]
		tempminaic<-min(temps[,length(temps[1,])-1])
		targetrow<-(temps[,length(temps[1,])-1])<tempminaic+daic
		tobind<-temps[targetrow,]
		out<-tobind
	
	
	if (length(Ks)>1){

	for (i in 2:length(Ks)){
		temps<-res[which(res[,1]==Ks[i]),]
		tempminaic<-min(temps[,length(temps[1,])-1])
		targetrow<-(temps[,length(temps[1,])-1])<tempminaic+daic
		tobind<-temps[targetrow,]
		out<-rbind(out,tobind)
		}
		}
		

				if(show.rate!=FALSE){
					keep<- (dim(out)[2]-4-length(Ks))  :dim(out)[2]
					out<-out[,c(1,2,keep)]}
				if(show.rate==FALSE){
					keep<-(dim(out)[2]-4):dim(out)[2]
					out<-out[,c(1,2,keep)]}
				
		return(out)
	}

