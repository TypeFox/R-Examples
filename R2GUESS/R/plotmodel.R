plotmodel <-
function(x,threshold.model=20){
if(is.null(x$label.Y)) Pheno <- paste("Y",1:x$q,sep="",collapse="_") else Pheno <- paste("Y: ",paste(x$label.Y,collapse="/"),sep="") 

if(threshold.model<1) {
ListSelect <- which(x$BestModels$postProb>threshold.model) }else{
ListSelect <-1:threshold.model
}

n.model <- length(ListSelect)
list.best.model <- strsplit(x$BestModels$modelName[ListSelect]," ")
distance <- matrix(NA,ncol=n.model,nrow=n.model)
for( i in 1:(n.model-1)){
for (j in (i+1):n.model){
distance[i,j] <- length(intersect(list.best.model[[i]],list.best.model[[j]]))/min(length(list.best.model[[i]]),length(list.best.model[[j]]))
}}
breaks=seq(0,1,by=0.1)
par(mar=c(5,5,5,7))
toto <-substitute(expression(paste(" Proximity measure:   ",frac("#"*(M[i]*intersect()*M[j]),"MIN"(M[i],M[j])))),list(Pheno=Pheno))
toto <-substitute(paste(" Proximity measure:   ",frac("#"*(M[i]*intersect()*M[j]),"MIN"(M[i],M[j]))),list(Pheno=Pheno))

break.seq <- seq(0,1,by=0.1)
xlab <- paste("Model for ",Pheno,sep="")
image(x=1:length(ListSelect),y=1:length(ListSelect),t(distance),zlim=c(0,1),axes=FALSE,xlab=xlab,ylab="Model",main=toto,col=rev(heat.colors(10)),breaks=break.seq)
par(las=2)
par(cex.axis=0.5)
axis(2, at = 1:length(ListSelect),paste("M",1:length(ListSelect),sep=""))
par(cex.axis=0.5)
axis(1, at = 1:length(ListSelect),paste("M",1:length(ListSelect),sep=""))
box()
image.plot(x=1:length(ListSelect),y=1:length(ListSelect),t(distance),nlevel=10,zlim=c(0,1),xlab=xlab,ylab="Model",main=toto,col=rev(heat.colors(10)),horizontal = FALSE,legend.only=TRUE,lab.breaks =break.seq,breaks=break.seq)


}
