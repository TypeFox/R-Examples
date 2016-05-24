`gknn` <-
function(data,cl,newdata,nn=5,distance=NULL,use.weights=FALSE,...){
	if(!is.matrix(data) | !is.matrix(newdata))
		stop("Both data and newdata must be matrix objects.")
	if(length(cl)!=nrow(data))
		stop("Length of cl must be equal to the number of rows of data.")
	r.cl<-range(cl)
	le.cl<-length(unique(cl))
	if(le.cl>10)
		stop("cl contains more than 10 different values.")
	if(r.cl[1]!=1 | r.cl[2]!=le.cl)
		stop("cl must consist of integers between 1 and ",le.cl,".")
	if(any(!cl%in%1:le.cl))
		stop("cl must consist of integers between 1 and ",le.cl,".") 
	if(nn<1)
		stop("nn must be at least 1.")
	if(nn>nrow(data))
		stop("nn must be smaller than or equal to the number of rows of data.")
	if(is.null(distance))
		distance<-getDistance(data)
	distance<-match.arg(distance,c("smc","cohen","pcc","euclidean","maximum","manhattan",
		"canberra","minkowski"))
	if(distance%in%c("smc","cohen","pcc"))
		checkX1X2(data,newdata,impute=FALSE)
	distance<-paste(distance,"2Mats",sep="")
	FUN<-match.fun(distance)
	mat.dist<-FUN(data,newdata,...)
	rownames(mat.dist)<-1:nrow(mat.dist)
	preds<-numeric(nrow(newdata))
	colS<-colSums(mat.dist<0.000001)
	ids1<-which(colS>0)
	ids2<-which(colS==0)
	for(i in ids1){
		tmp.ids<-which(mat.dist[,i]<0.000001)
		preds[i]<-modeDist(cl[tmp.ids])
	}
	for(i in ids2){
		tmp.dist<-sort(mat.dist[,i])[1:nn]
		tmp.ids<-as.numeric(names(tmp.dist))
		preds[i]<-if(use.weights) weightMode(tmp.dist^-1,cl[tmp.ids]) else modeDist(cl[tmp.ids])
	}
	preds
}

