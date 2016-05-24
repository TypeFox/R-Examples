extractFeatures.T <-
function(output.ratio,feat.names){
load(paste(output.ratio$dataname,".Rdata"))
if(output.ratio$pvalue==FALSE){
data = 1 - data
}

dim1=dim(data)[1]
lists = dim(data)[2]


#Decision rules:
max.T = max(output.ratio$ratios)
threshold.max = output.ratio$h[output.ratio$ratios==max.T]

#function Table
table=function(threshold){
temp=matrix(0,dim1,lists)
for(i in 1:dim1){
for(j in 1:lists){
if(data[i,j]<= threshold){temp[i,j]<-1}
}
}
return(temp)
}

#Table
temp<-table(threshold.max)

name="Names"
for(i in 1:ncol(data)){
name=c(name,paste("List",as.character(i)))
}

table.max <- data[apply(temp,1,sum)==lists,]
names.max <- feat.names[apply(temp,1,sum)==lists]

if(output.ratio$pvalue==FALSE){
table.max <- 1-table.max
}

if(is.vector(table.max)==TRUE){
table.max=as.matrix(table.max)
table.max=t(table.max)
}

table.max <- data.frame(Names=names.max,RankingStat = table.max)
colnames(table.max)<-name

write.csv(table.max,"featuresT.csv")

return(list(max = table.max))
}

