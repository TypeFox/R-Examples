hca.plot <-
function(g,o,method="correlation",link="ward",colored=palette(),border=NA,code=colnames(o),cex.code=1,breaks=round(nrow(oreihe)/4),cutcolors=colorpanel(breaks,low="green",mid="black",high="red")){
require(amap)
require(gplots)

if(is.data.frame(o)){
      if(!identical(rownames(o),colnames(g))){
        stop("Colnames of g are not the same as rownames of o")}
          }
if(!is.data.frame(o)){
  o<-data.frame(o,row.names=colnames(g))} ## vector in data.frame

   classes<-unlist(lapply(unclass(o),class))
   if(all(classes%in%c("factor","numeric","integer"))==F){stop("o can only contain factors and numeric")}
      
fit<-hcluster(t(g),method=method,link=link)
dend1<-as.dendrogram(fit)

Reihenfolge<-colnames(g)[fit$order]
oreihe<-o[Reihenfolge,]
if(!is.data.frame(oreihe)){
  oreihe<-data.frame(oreihe)}  ## vector in data.frame

    def.par <- par(no.readonly = TRUE)
layout(c(1,2),heights=c(1,dim(oreihe)[2]/10))
par(mar=c(0,6,0,1))
plot(dend1,xlab="",ylab="",cex=0.5,main="",leaflab="none")
par(mar=c(0,6+(15/nrow(oreihe)),0,1+(15/nrow(oreihe))))
plot(1,type="n",axes=F,ylab="",xlab="",xlim=c(1,nrow(oreihe)),ylim=c(0,dim(oreihe)[2]))
for(j in 1:dim(oreihe)[2]){
for (i in 1:nrow(oreihe)){
  if (classes[j]=="factor"){
  rect(i-0.5,j-1,i+0.5,j-0.1,col=colored[as.numeric(oreihe[,j])[i]],border=border)
  }
  if (classes[j]%in%c("numeric","integer")){
  cutted<-cut(oreihe[,j],breaks)
  cutcolors<-cutcolors
  rect(i-0.5,j-1,i+0.5,j-0.1,col=cutcolors[as.numeric(cutted)[i]],border=border)
  }
}
}
axis(side=2,at=(1:ncol(oreihe))-0.5,labels=code,tick=F,pos=0,las=2,cex.axis=cex.code)
    par(def.par)
}

