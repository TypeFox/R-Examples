`knncatimpute` <-
function(x,dist=NULL,nn=3,weights=TRUE){
	n.row<-nrow(x)
	if(nn<1)
		stop("nn must be at least 1.")
	if(nn>=n.row)
		stop("nn must be smaller than the number of rows of x.")
	check4Monomorphism(x)
	if(is.null(dist))
		dist<-smc(x,dist=TRUE)
	else
		dist<-checkDist(x,dist)
	tmp<-colSums(is.na(x))
	na.cols<-which(tmp>0)
	if(length(na.cols)==0)
		stop("There are no missing values in x.")
	rownames(dist)<-colnames(dist)<-1:nrow(x)
	for(i in na.cols){
		ids<-which(is.na(x[,i]))
		for(j in ids){
			neighbors<-sort(dist[j,-ids])[1:nn]
			ids.nn<-as.numeric(names(neighbors))
			if(round(neighbors[1],8)==0)
				x[j,i]<-x[ids.nn[1],i]
			else{
				if(weights)
					x[j,i]<-weightMode(neighbors^-1,x[ids.nn,i]) 
				else 
					x[j,i]<-modeDist(x[ids.nn,i])
			}
		}
	} 
	x
}

