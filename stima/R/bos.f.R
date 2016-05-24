bos.f <-
function(x,nv,col,mb,crit="f2")
{
	## x=data,nv=nodevec,mb=minbucket
	##anova(lm(data),lm(newdata))$F Value[2]
	##(y2rsq-y1rsq)*(n1-n2)/(1-y2rsq)
	critn<--1
	if(crit=="F-value") {critn<-2}
	if(crit=="R2change") {critn<-1}
	if(crit=="f2") {critn<-0}
	m<-dim(x)[1]
	n<-dim(x)[2]
	res<-numeric(2)
  	res<-.Fortran("rs_bos",res=as.double(res),x=as.matrix(x),
	     as.integer(nv),as.double(col),as.integer(m),as.integer(n),
	     as.integer(mb),as.integer(critn),package="stima")$res
	return(res)
}
