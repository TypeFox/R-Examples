dim.thetaMSAR <-
function(x) {
	thetadim=list()
    if(!is.thetaMSAR(x)){thetadim=NULL}
    else{thetadim[[1]]=attributes(x)$n_par
    	thetadim[[2]]=attributes(x)$order
    	thetadim[[3]]=attributes(x)$NbRegimes
    	thetadim[[4]]=attributes(x)$NbComp
    	names(thetadim)=c("Parameters","Order","NbRegimes","Dimension")
    	for(i in 1:4){
    		names(thetadim[[i]])=""}
    	}
   return(thetadim)
    	
}
