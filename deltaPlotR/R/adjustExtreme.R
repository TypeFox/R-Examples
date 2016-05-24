adjustExtreme<-function(data=NULL,group=NULL,focal.name=NULL,prop,method="constraint",const.range=c(0.001,0.999),nrAdd=1){
res<-prop
if (!is.null(data) & !is.null(group) & !is.null(focal.name) & method=="add"){
meth<-"add"
for (i in 1:nrow(res)){
if (res[i,1]==0 | res[i,1]==1){
x<-data[,i][!is.na(data[,i]) & group!=focal.name]
res[i,1]<-(sum(x)+nrAdd)/(length(x)+nrAdd*2)
}
if (res[i,2]==0 | res[i,2]==1){
x<-data[,i][!is.na(data[,i]) & group==focal.name]
res[i,2]<-(sum(x)+nrAdd)/(length(x)+nrAdd*2)
}
}
}
else{
meth<-"constraint"
res[res<const.range[1]]<-const.range[1]
res[res>const.range[2]]<-const.range[2]
}
return(list(adj.prop=res,method=meth,range=const.range,nrAdd=nrAdd))
}
