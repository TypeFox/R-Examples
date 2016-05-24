print.minDNF<-function(x,which=0,...){
	if(is.character(x)){
		cat("Minimum DNF",if(which>0) paste("",which),":\n",sep="")
		cat("   ",x[1],"\n")
		for(i in 2:length(x))
			cat(" | ",x[i],"\n")
	}
	else{
		cat("Because of cyclic covering, there are",length(x),"solutions:\n")
		for(i in 1:length(x)){
			cat("\n")
			class(x[[i]])<-"minDNF"
			print(x[[i]],which=i)
		}
	}
}

