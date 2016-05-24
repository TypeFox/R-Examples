"GRmarking"=function(res,lev,col.neg="lightblue",col.pos="pink",colour="black"){
#res:objet de type mark_chi
#lev:choix du niveau à représenter
for (m in 1:length(res)){
if(length(res[[m]][[lev]])==1)
  print(paste("Level",lev,": nothing to plot for the category",names(res)[m]))
if(length(res[[m]][[lev]])>1){
col<-rep(0)
dev.new()
res2=res[[m]][[lev]][[2]][,5]
for (i in 1:length(res2)){
		if (res2[i]<0)
			col<-c(col,col.neg)
		else col<-c(col,col.pos)
}
col<-col[-1]
par(cex=0.7)
par(yaxt="n")
coord=barplot(res2,beside=TRUE,horiz=TRUE,las=2,col=col,main=paste(c("Category ",names(res)[m]),collapse=": "),cex.lab=5,legend.text=FALSE)
if (min(res2)>0 | max(res2)<0)
text(x=(res2/2),y=coord,labels=rownames(res[[m]][[lev]][[2]]),adj=0.5,col=colour)
else
text(x=0,y=coord,labels=rownames(res[[m]][[lev]][[2]]),adj=0.5,col=colour)
mtext(paste("Level : ",lev),side=3)}
}}
