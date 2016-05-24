long2matrices <- function(id,time=NULL,X=NULL,Y){

# preliminaries
	idu = unique(id)
	n = length(idu)
	if(is.null(time)){
		Tv = table(id)
		TT = max(Tv)
	}else{
		TT = max(time)	
	}
	if(!is.null(X)){
		X = as.matrix(X)
		nx = ncol(X)
	} 
	Y = as.matrix(Y)
	ny = ncol(Y)
# create matrices
	if(!is.null(X)) XX = array(NA,c(n,TT,nx)) else XX = NULL
	YY = array(NA,c(n,TT,ny))
	for(i in 1:n){
		ind = which(id==idu[i])
		if(is.null(time)){
			for(t in 1:TT){
				if(!is.null(X)) XX[i,t,] = X[ind[t],] 
				YY[i,t,] = Y[ind[t],]
			}
		}else{
			tmp = 0
			for(t in time[ind]){
				tmp=tmp+1
				if(!is.null(X)) XX[i,t,] = X[ind[tmp],] 
				YY[i,t,] = Y[ind[tmp],]
			}	
		}		
	}
# output
	out = list(XX=XX,YY=YY)
	out
  
}