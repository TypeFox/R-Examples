dist.GDM<-function(x,method="GDM1",weightsType="equal",weights=NULL)
{
	if(is.data.frame(x)) x<-as.matrix(x)
	if (method !="GDM1" && method !="GDM2" && method !="GDM3")
	{
		print("Bad GDM method")
		print("Choose 'GDM1' for metric data or 'GDM2' for ordinal data")
		stop()
	}
  	if(weightsType=="equal" || is.null(weights)){
 		 weights<-array(1,c(ncol(x)))
	}
  	if(weightsType=="different2"){
    if(sum(weights)!=ncol(x) || sum(weights<0)!=0) {
    stop("for wegithsType='different2' weights should satisfy: each weight takes value from [0; m] and sum of weights ' m (m-num. of variables)")
	    }
	}
	if(weightsType=="different1"){
    if(sum(weights)!=1 || sum(weights<0)!=0){
     stop("for wegithsType='different1' weights should satisfy: each weight takes value from [0; 1] and sum of weights eguals one")
     	}
	}
	method_int<- switch(method,
			"GDM1"=1,
			"GDM2"=2,
			"GDM3"=3)

	nr=nrow(x)
	t<-.C("fngdm",as.double(x),as.integer(nrow(x)),as.integer(ncol(x)),
	as.integer(method_int),as.double(weights),wynik=double(nrow(x)*nrow(x)),PACKAGE="clusterSim")$wynik
	wynik<-matrix(nrow=nr,ncol=nr,dimnames=names(x))
	for (i in 1:nr)
	for (j in 1:nr)
	{
		wynik[i,j]=t[(i-1)*nr+j]
		wynik[j,i]=t[(j-1)*nr+i]
	}
	as.dist(wynik)
}
GDM<-function(x,method="GDM1",weightsType="equal",weights=NULL){dist.GDM(x,method=method,weightsType,weights)}
GDM1<-function(x,weightsType="equal",weights=NULL){dist.GDM(x,"GDM1",weightsType,weights)}
GDM2<-function(x,weightsType="equal",weights=NULL){dist.GDM(x,"GDM2",weightsType,weights)}
