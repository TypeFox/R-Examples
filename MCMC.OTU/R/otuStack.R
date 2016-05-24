otuStack <-
function(data,count.columns,condition.columns){
	summ=apply(data[,count.columns],1,sum)
	data$summ=as.integer(summ)
	row.names(data)=NULL
	dds=stack(data[,c(count.columns,length(data[1,]))])
	dds=cbind(dds,data[,condition.columns])
	names(dds)[1:2]=c("count","otu")
	for(co in 3:length(dds[1,])) { dds[,co]=as.factor(dds[,co]) }
	dds$otu=relevel(dds$otu,ref="summ")
	return(dds)
}
