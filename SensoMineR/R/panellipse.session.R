"panellipse.session" <- function(donnee,col.p,col.j,col.s,firstvar,lastvar=ncol(donnee),alpha=0.05,coord=c(1,2),scale.unit=TRUE,nbsimul=500,nbchoix=NULL,level.search.desc=0.2,centerbypanelist=TRUE,scalebypanelist=FALSE,name.panelist=FALSE,variability.variable=FALSE,cex=1,color=NULL){

hotelling <- function(d1,d2,n1=nrow(d1),n2=nrow(d2)){
    k <- ncol(d1)
#    n1 <- nrow(d1)
#    n2 <- nrow(d2)
    xbar1 <- apply(d1,2,mean)
    xbar2 <- apply(d2,2,mean)
    dbar <- xbar2-xbar1
    v <- ((n1-1)*var(d1)+(n2-1)*var(d2))/(n1+n2-2)
    t2 <- n1*n2*dbar%*%solve(v)%*%dbar/(n1+n2)
    f <- (n1+n2-k-1)*t2/((n1+n2-2)*k)
    return(pf(f,k,n1+n2-k-1,lower.tail=FALSE))
}

if (length(color)==0) color = c("black","red","green3","blue",
  "cyan","magenta","darkgray","darkgoldenrod","darkgreen","violet","turquoise","orange","lightpink","lavender","yellow","lightgreen","lightgrey",
  "lightblue","darkkhaki", "darkmagenta","darkolivegreen","lightcyan", "darkorange",
  "darkorchid","darkred","darksalmon","darkseagreen","darkslateblue","darkslategray","darkslategrey",
  "darkturquoise","darkviolet", "lightgray","lightsalmon","lightyellow", "maroon")
for (j in 1:(firstvar-1)) donnee[,j]=as.factor(donnee[,j])
labseance=levels(as.factor(donnee[,col.s]))
nbseance <- length(labseance)
labprod = levels(as.factor(donnee[,col.p]))
nbprod <- length(levels(donnee[,col.p]))
nbjuge <- length(levels(donnee[,col.j]))

  donnee <- search.desc(donnee,col.j=col.j,col.p=col.p,firstvar=firstvar,lastvar=lastvar,level=level.search.desc)
  if (nbseance <2) print("This procedure is not adapted, there is only one session")
  oo=order(donnee[,col.j])
  donnee<- donnee[oo,]
  oo=order(donnee[,col.p])
  donnee<- donnee[oo,]
  oo=order(donnee[,col.s])
  donnee<- donnee[oo,]

  don <- cbind.data.frame(donnee[donnee[,col.s]==labseance[1],col.j],donnee[donnee[,col.s]==labseance[1],col.p])

  for (seance in 1:nbseance)  don <- cbind.data.frame(don,data.frame(donnee[donnee[,col.s]==labseance[seance],firstvar:ncol(donnee)],row.names=paste(donnee[donnee[,col.s]==labseance[seance],col.p],donnee[donnee[,col.s]==labseance[seance],col.j],sep=".")))
  colnames(don) <- colnames(donnee)[c(col.j,col.p,rep(firstvar:ncol(donnee),nbseance))]
  colnames(don) <- paste(colnames(don),c("","",rep(paste(".S",1:nbseance,sep=""),rep(ncol(donnee)-firstvar+1,nbseance))),sep="")
  bb=panellipse(don,group=c(rep(ncol(donnee)-firstvar+1,nbseance)),name.group=c(paste("S",1:nbseance,sep="")),col.j=1,col.p=2,firstvar=3,alpha=alpha,coord=coord,scale.unit=scale.unit,nbsimul=nbsimul,nbchoix=nbchoix,level.search.desc=1,centerbypanelist=centerbypanelist,scalebypanelist=scalebypanelist,name.panelist=name.panelist,variability.variable=variability.variable,cex=cex,color=color)
  legend("bottomleft",legend=paste(colnames(donnee)[col.s],1:nbseance,sep=" "),lty=1:nbseance,cex=0.8,bg="white")

  mat = list(bysession=bb$hotelling$bygroup,byproduct=bb$hotelling$byproduct)
  dimnames(mat$bysession)=list(labprod,labprod,c(paste(colnames(donnee)[col.s],1:nbseance,sep=" "),"global"))
  dimnames(mat$byproduct)=list(paste(colnames(donnee)[col.s],1:nbseance,sep=" "),paste(colnames(donnee)[col.s],1:nbseance,sep=" "),labprod)

  aa <- matrix(0,ncol(donnee)-firstvar+1,2)
  res.average=averagetable(don,formul=as.formula(paste("~",colnames(don)[2])),firstvar=3)
  for (j in 1:(ncol(donnee)-firstvar+1)) aa[j,] <- as.matrix(PCA(res.average[,(ncol(donnee)-firstvar+1)*(0:(nbseance-1))+j],graph=FALSE)$eig[1,1:2])
  rownames(aa) <- colnames(donnee[,firstvar:ncol(donnee)])
  colnames(aa) <- c("eig1","Reproductibility")
  res <- list()
  res$bysession =don
  res$eig =bb$eig
  res$coordinates =bb$coordinates
  if (variability.variable) res$correl =bb$correl
  res$hotelling =mat
  res$variability=aa
  return(res)
}
