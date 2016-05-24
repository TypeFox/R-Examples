mergeNodes <-
function(dis,na.rm.row.col=FALSE,save.distance=FALSE, save.distance.name="Merged_Distance.txt"){
if(length(which(dis==0))==nrow(dis)) stop("Matrix with no off-diagonal zeros")
DIS<-dis

if(length(which(is.na(dis)))!=0 & na.rm.row.col==FALSE) stop("NA values found")
{
if(length(which(is.na(dis)))!=0 & na.rm.row.col==TRUE)
{
dis<-as.matrix(dis)
repeat{
conNA<-c()
for (i in 1:nrow(dis))
conNA<-c(conNA,length(which(is.na(dis[i,]))))
Out<-sort(which(conNA==sort(conNA,decreasing=T)[1]),decreasing=T)[1]
dis<-dis[-Out,-Out]
if(nrow(dis)==0) stop ("The algorithm could not find a matrix without NA values")
if(length(which(is.na(dis)))==0) break
}
}

repes<-list()
length(repes)<-nrow(dis)
names(repes)<-colnames(dis)
for (zero1 in 1:nrow(dis))
{
kk<-c()
	for (zero2 in 1:nrow(dis))
	{
	if(dis[zero1,zero2]==0)
	kk<-c(kk,zero1,zero2)
	}
	repes[[zero1]]<-c(colnames(DIS)[kk])
	
}

for (i in 1:length(repes))
{
repes[[i]]<-sort(unique(repes[[i]]))
if(length(repes[[i]])==1)
repes[[i]]<-NA
}

groups<-unique(repes)
groups<-groups[which(is.na(groups)==FALSE)]

ToRemove<-c()
dis3<-dis
for (i in 1:length(groups))
{
colus<-match(groups[[i]],colnames(dis))
newC<-as.matrix(rowMeans(dis3[,colus]))
colnames(newC)<-paste(groups[[i]],collapse="-")
NEW<-cbind(dis3,newC)
NEW<-rbind(NEW,c(newC,0))
dis3<-NEW
row.names(NEW)<-colnames(NEW)
ToRemove<-c(ToRemove,colus)
}
if(save.distance==TRUE) write.table(file=save.distance.name,dis)
if(length(which(NEW==0))==nrow(NEW)) print("Matrix free of off-diagonal zeros not found")
NEW[-ToRemove,-ToRemove]
}
}
