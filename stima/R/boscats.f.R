boscats.f <-
function(data,nv,col,mb,crit="f2")
{
	## nv=nodevec,mb=minbucket
	##anova(lm(data),lm(newdata))$F Value[2]
	##(y2rsq-y1rsq)*(n1-n2)/(1-y2rsq)
	critn<--1
	if(crit=="F-value") {critn<-2}
	if(crit=="R2change") {critn<-1}
	if(crit=="f2") {critn<-0}
	m<-dim(data)[1]
	n<-dim(data)[2]
	res<-numeric(2)
  	res<-.Fortran("rs_boscats",res=as.double(res),data=as.matrix(data),
	     as.integer(nv),as.double(col),as.integer(m),as.integer(n),
	     as.integer(mb),as.integer(critn),package="stima")$res
	object<-res
}
