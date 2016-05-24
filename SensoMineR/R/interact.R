"interact" <- function(donnee,col.p=1,col.j,firstvar,lastvar=ncol(donnee)){


############################################################################
plotinteract<-function(tab,cex=1.1,xlegend=ncol(tab)-5,ylegend=max(tab),xlab=NULL,ylab=NULL,main=NULL){
  x <- as.factor(1:ncol(tab))
  miny <- min(tab,na.rm=TRUE)
  maxy <- max(tab,na.rm=TRUE)
  plot(as.integer(x),apply(tab,2,mean),type="n",main=main,xlab=xlab,ylab=ylab,ylim=c(miny,maxy),xlim=c(1,ncol(tab)),cex=0.8)
  abline(v = x, lty = "dotted")
  abline(h = 0)
  for (i in 1:nrow(tab)) points(x,tab[i,],col=i,cex=cex,pch=20)
  legend("topright",legend=rownames(tab),col=1:nrow(tab),pch=rep(20,nrow(tab)),cex=0.8,bg="white")
}
############################################################################
old.contr = options()$contrasts
options(contrasts=c("contr.sum", "contr.sum"))
for (j in 1 :(firstvar-1))  donnee[,j] <- as.factor(donnee[,j])
nbprod <- length(levels(donnee[,col.p]))
nbjuge <- length(levels(donnee[,col.j]))
tab<-array(0,dim=c(nbprod,nbjuge,lastvar-firstvar+1))

for (varendo in firstvar:lastvar) {
  aux <- summary.lm(aov( donnee[,varendo]~donnee[,col.p]+donnee[,col.j]+donnee[,col.p]:donnee[,col.j], data = donnee, na.action =na.exclude))$coef
  for (k in 1:(nbjuge-1)) tab[1:(nbprod-1),k,varendo-firstvar+1] <- aux[((nbprod+nbjuge-1)+(k-1)*(nbprod-1)+1):((nbprod+nbjuge-1)+k*(nbprod-1)),1]
  tab[,nbjuge,varendo-firstvar+1] <- - apply(tab[,,varendo-firstvar+1],1,sum)
  tab[nbprod,,varendo-firstvar+1] <- - apply(tab[,,varendo-firstvar+1],2,sum)
}

dimnames(tab) = list(levels(donnee[,col.p]),levels(donnee[,col.j]),labels(donnee)[[2]][firstvar:lastvar])

for (k in 1:dim(tab)[[3]]){
  plotinteract(tab[,,k],main=colnames(donnee)[firstvar+k-1],xlab=colnames(donnee)[col.j],ylab=paste(colnames(donnee)[col.p],"-",colnames(donnee)[col.j],"interaction coefficients"))
  if (k != dim(tab)[[3]]) dev.new()
}
barrow(t(apply(tab^2,c(1,3),sum)/matrix(rep(apply(tab^2,3,sum),nrow(tab)),byrow=TRUE,nrow=nrow(tab))))

## Make a graph to visualize the panelist which contribute the product-panelist interaction for each descriptor
barrow(t(  apply(tab^2,c(2,3),sum) /matrix(rep(apply(tab^2,3,sum),ncol(tab)),byrow=TRUE,nrow=ncol(tab))))
return(tab)
options(contrasts=old.contr)
}
