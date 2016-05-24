
##################################################################################################################
#####################################      Control functions     #################################################
##################################################################################################################

hbstmControl <- function(time=TRUE,timerem=FALSE,seed=-1){
    list(
		time = as.logical(time),
		timerem = as.logical(timerem),
        seed = as.numeric(seed)
	)
}


controlspatlags=function(lengthSpatPar,spatlags){
	text1="The length of the vector "
	text2=" exceeds the maximum possible length. The function assign the length of the vector to the maximum length "
	ifelse(test=lengthSpatPar[1]<spatlags[1],yes=print(paste(text1,"alpha",text2,lengthSpatPar[1],sep="")),no=lengthSpatPar[1]<-spatlags[1])
	ifelse(test=lengthSpatPar[2]<spatlags[2],yes=print(paste(text1,"beta",text2,lengthSpatPar[2],sep="")),no=lengthSpatPar[2]<-spatlags[2])
	ifelse(test=lengthSpatPar[3]<spatlags[3],yes=print(paste(text1,"phi",text2,lengthSpatPar[3],sep="")),no=lengthSpatPar[3]<-spatlags[3])
	ifelse(test=lengthSpatPar[4]<spatlags[4],yes=print(paste(text1,"theta",text2,lengthSpatPar[4],sep="")),no=lengthSpatPar[4]<-spatlags[4])
	return(lengthSpatPar)
}


controlfun=function(Zt,K,newGrid,reglag,spatlags,posterior){
	if(dim(Zt)[1]!=dim(K)[1])	stop(" The number of spatial points in Zt differs from K \n")
	if(dim(K)[2]!=dim(newGrid)[1])	stop(" The number of spatial points in newGrid differs from K \n")
	if(dim(newGrid)[2]!=2)	stop(" newGrid must have two columns, first one with the Longitude values and second column with the Latitude \n")
	if(max(reglag)>=dim(Zt)[2]) stop(" The maximum temporal lag cannot be bigger than the maximum time \n")
	if(length(spatlags)!=4) stop("The length of the vector spatlags must be 4, containing the spatial lags of the four main directions \n")
	if((posterior!="mean")&(posterior!="median")) stop("The attribute 'posterior' is not 'mean' or 'median' \n")
}


estYt.check=function(S,T,nIter,nBurn,val){
	if(as.logical(Sys.info()["sysname"] == "Windows")){
		memory2=object.size(matrix(val,nrow=S,ncol=T))/1e6
		if(memory2*(nIter-nBurn)>0.9*memory.limit()){
			stop(paste("The memory required to save the fitted parameters ",nIter-nBurn," iterations requires more than the 90% of the maximum RAM memory.Please, try with less iterations.\n For your case, we suggest to save less than ",round((0.9*memory.limit())/memory2,0)," iterations\n",sep=""))
		}
	}
}

save.check=function(save,parameters,nIter,nBurn,S,T,val){
	
	if(length(save)>1) stop("The length of the argument 'save' has to be 1\n")
	
	opt=c("Parameters","Mu","Mt","Xt")
	if(sum(which(save==opt))==0){
		stop("The argument 'save' has to be between the options 'Parameters', 'Mu', 'Mt' or 'Xt'\n")
	}
	
	if(as.logical(Sys.info()["sysname"] == "Windows")){
		if(save=="Parameters"){
			memory=object.size(parameters)/1e6
		}else{
			if(save=="Mu"){
				memory=object.size(parameters@Mu)/1e6
			}else{
				if(save=="Mt"){
					memory=object.size(parameters@Mt)/1e6
				}else{
					if(save=="Xt"){
						memory=object.size(parameters@Xt)/1e6
					}
				}
			}
		}	
		
		memory2=object.size(matrix(val,nrow=S,ncol=T))/1e6
		if((memory+memory2)*(nIter-nBurn)>0.9*memory.limit()){
			stop(paste("The memory required to save the MCMC sample parameters ",nIter-nBurn," iterations requires more than the 90% of the maximum RAM memory.Please, try with less iterations.\n For your case, we suggest to save less than ",round((0.9*memory.limit())/(memory+memory2),0)," iterations\n",sep=""))
		}
	}
	if(nBurn>=nIter) stop("nBurn must be less than nIter.\n")
}


component.check=function(component,clas){
	if(sum(which(component==c("Parameters","Mu","Mt","Xt")))==0){
		stop("The argument 'component' has to be between the options 'Parameters', 'Mu', 'Mt' or 'Xt'\n")
	}else{
		if(clas!="Parameters"){
			if(component!=clas) stop("The argument 'component' is different from the class of the saved values of the argument 'object'\n")
		}
	}
}