"panellipse" <- function(donnee,col.p,col.j,firstvar,lastvar=ncol(donnee),alpha=0.05,coord=c(1,2),scale.unit=TRUE,nbsimul=500,nbchoix=NULL,group=NULL,name.group=NULL,level.search.desc=0.2,centerbypanelist=TRUE,scalebypanelist=FALSE,name.panelist=FALSE,variability.variable=TRUE,cex=1,color=NULL){

hotelling <- function(d1,d2,n1=nrow(d1),n2=nrow(d2)){
    k <- ncol(d1)
    xbar1 <- apply(d1,2,mean)
    xbar2 <- apply(d2,2,mean)
    dbar <- xbar2-xbar1
    if (n1+n2<3) return(NA)
    v <- ((n1-1)*var(d1)+(n2-1)*var(d2))/(n1+n2-2)
    if (sum(v^2) < 1/10^10) return (NA)
    else t2 <- n1*n2*dbar%*%solve(v)%*%dbar/(n1+n2)
    f <- (n1+n2-k-1)*t2/((n1+n2-2)*k)
    return(pf(f,k,n1+n2-k-1,lower.tail=FALSE))
}

variab.variable <- function(donnee,echantillon,mfa=FALSE,coord=c(1,2),scale.unit=TRUE,centerbypanelist=TRUE,scalebypanelist=TRUE,color=color){

  for (j in 1 :2)  donnee[,j] <- as.factor(donnee[,j])
  nbjuge <- length(levels(donnee[,1]))
  nbprod <- length(levels(donnee[,2]))
  nbdesc=ncol(donnee)-2
  nbcoord=max(coord)

  oo <- order(donnee[,2])
  donnee <- donnee[oo,]
  oo <- order(donnee[,1])
  donnee <- donnee[oo,]
  tab=scalebypanelist(donnee,col.j=1,col.p=2,firstvar=3,center=centerbypanelist,scale=scalebypanelist)

  tab.moy <- as.matrix(tab[1:nbprod,-(1:2)])
  tab.byjudge <- array(0,dim=c(nbprod,nbdesc,nbjuge))
  for (j in 1:nbjuge) tab.byjudge[,,j] <- as.matrix(tab[(j*nbprod+1):((j+1)*nbprod),-(1:2)])
  correl = array(NA,dim=c(nbdesc,nbdesc,nbsimul))
  res = array(NA,dim=c(nbsimul,ncol(tab.moy),nbcoord))
  for (k in 1:nbsimul){
     Xb = apply(tab.byjudge[,,echantillon[k,]],c(1,2),mean)
     correl[,,k] = cor(Xb)
     if (!mfa){
       resAF <-PCA(cbind(tab.moy,Xb),quanti.sup=(ncol(tab.moy)+1):(2*ncol(tab.moy)),graph=FALSE,ncp=nbcoord,scale.unit=scale.unit)
       res[k,,] = as.matrix(resAF$quanti.sup$coord)
     }
     if (mfa){
       if (scale.unit) resAF <- MFA(cbind(tab.moy,Xb),group=c(group,group),type=rep("s",2*length(group)),num.group.sup=(length(group)+1):(2*length(group)),graph=FALSE,ncp=nbcoord)
       if (!scale.unit) resAF <- MFA(cbind(tab.moy,Xb),group=c(group,group),type=rep("c",2*length(group)),num.group.sup=(length(group)+1):(2*length(group)),graph=FALSE,ncp=nbcoord)
       res[k,,] = as.matrix(resAF$quanti.var.sup$cor)
     }
     if (k==1){
        dev.new(width=12,height=8)
		plot(resAF,choix="var", invisible=c("quanti.sup","sup"),new.plot = FALSE)
        legend("topleft",legend=colnames(tab.moy),fill=color[1:ncol(tab.moy)],cex=0.7)
     }
     points(res[k,,coord[1]],res[k,,coord[2]],col=color[1:ncol(tab.moy)],pch=15,cex=0.3)
  }

  mini=maxi=matrix(0,ncol(tab.moy),ncol(tab.moy))
  for (i in 1:ncol(tab.moy)){
    for (j in 1:ncol(tab.moy)){
      mini[i,j] = correl[i,j,order(correl[i,j,])[round(nbsimul*0.025,0)]]
      maxi[i,j] = correl[i,j,order(correl[i,j,])[round(nbsimul*0.975,0)]]
  }}
  colnames(mini)=rownames(mini)=colnames(maxi)=rownames(maxi)=colnames(tab.moy)
  res <- list()
  res$moy <- cor(tab.moy)
  res$mini <- mini
  res$maxi <- maxi
  return(res)
}

  if (length(color)==0) color = c("black","red","green3","blue",
    "cyan","magenta","darkgray","darkgoldenrod","darkgreen","violet","turquoise","orange","lightpink","lavender","yellow","lightgreen","lightgrey",
    "lightblue","darkkhaki", "darkmagenta","darkolivegreen","lightcyan", "darkorange",
    "darkorchid","darkred","darksalmon","darkseagreen","darkslateblue","darkslategray","darkslategrey",
    "darkturquoise","darkviolet", "lightgray","lightsalmon","lightyellow", "maroon")

  don.interesting <- search.desc(donnee,col.j=col.j,col.p=col.p,firstvar=firstvar,lastvar=lastvar,level=level.search.desc)
  don.interesting <- don.interesting[,c(col.j,col.p,firstvar:ncol(don.interesting))]
if (!is.null(group)){
  group.aux = sum((colnames(donnee[,firstvar:lastvar])%in%colnames(don.interesting[,3:ncol(don.interesting)]))[1:group[1]])
  for (j in 2:length(group)) group.aux = c(group.aux,sum((colnames(donnee[,firstvar:lastvar])%in%colnames(don.interesting[,3:ncol(don.interesting)]))[(sum(group[1:(j-1)])+1):sum(group[1:j])]))
  group.old = group
  group = group.aux
}
  axe <- construct.axes(don.interesting,group=group,name.group=name.group,coord=coord,scale.unit=scale.unit,centerbypanelist=centerbypanelist,scalebypanelist=scalebypanelist)
  labprod = axe$moyen[axe$moyen[,ncol(axe$moyen)]==0,ncol(axe$moyen)-1]
  nbprod = length(labprod)
  nbjuge = nlevels(as.factor(don.interesting[,1]))
  if (is.null(nbchoix)) nbchoix = nbjuge
if (length(group)<2) {
  mat = matrix(NA,nbprod,nbprod)
    aa = axe$moyen[-(1:nbprod),]
    for (i in 1:(nbprod-1)){
      for (j in i:nbprod) mat[i,j] = mat[j,i] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa[aa[,ncol(aa)-1]==labprod[j],coord],nbchoix,nbchoix)
    }
    diag(mat)=1
    colnames(mat)=rownames(mat)=labprod
}

if (length(group)>1) {
  mat = array(NA,dim=c(nbprod,nbprod,length(group)+1))
  for (k in 1:length(group)){
    aa = cbind.data.frame(axe$partiel[-(1:length(labprod)),max(coord)*(k-1)+(1:max(coord))],axe$partiel[-(1:length(labprod)),(ncol(axe$partiel)-1):ncol(axe$partiel)])
    for (i in 1:(nbprod-1)){
      for (j in i:nbprod) mat[i,j,k] = mat[j,i,k] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa[aa[,ncol(aa)-1]==labprod[j],coord],nbchoix,nbchoix)
    }
    for (i in 1:nbprod) mat[i,i,k]=1
  }
  aa=axe$moyen[-(1:length(labprod)),]
  for (i in 1:(length(labprod)-1)){
    for (j in (i+1):length(labprod)){
      if (length(nbchoix)==0) mat[i,j,length(group)+1] = mat[j,i,length(group)+1] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa[aa[,ncol(aa)-1]==labprod[j],coord])
      if (length(nbchoix)!=0) mat[i,j,length(group)+1] = mat[j,i,length(group)+1] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa[aa[,ncol(aa)-1]==labprod[j],coord],nbchoix,nbchoix)
  }}
  for (i in 1:nbprod) mat[i,i,length(group)+1]=1
  dimnames(mat)=list(labprod,labprod,c(paste("Group",1:length(group),sep=" "),"global"))
  mat2 = array(NA,dim=c(length(group),length(group),nbprod))
  for (i in 1:nbprod){
    for (k in 1:(length(group)-1)){
      aa = cbind.data.frame(axe$partiel[-(1:length(labprod)),max(coord)*(k-1)+(1:max(coord))],axe$partiel[-(1:length(labprod)),(ncol(axe$partiel)-1):ncol(axe$partiel)])
      for (kk in (k+1):length(group)){
        aa2 = cbind.data.frame(axe$partiel[-(1:length(labprod)),max(coord)*(kk-1)+(1:max(coord))],axe$partiel[-(1:length(labprod)),(ncol(axe$partiel)-1):ncol(axe$partiel)])
        mat2[k,kk,i] = mat2[kk,k,i] = hotelling(aa[aa[,ncol(aa)-1]==labprod[i],coord],aa2[aa2[,ncol(aa2)-1]==labprod[i],coord],nbchoix,nbchoix)
      }}
    for (k in 1:length(group)) mat2[k,k,i]=1
  }
  dimnames(mat2)=list(c(paste("Group",1:length(group),sep=" ")),c(paste("Group",1:length(group),sep=" ")),labprod)
}

  if (variability.variable==FALSE) dev.new()
	plotpanelist(axe$moyen,coord=coord,eig=signif(axe$eig,4),color=color,name=name.panelist,cex=cex)
  long.group <- length(group)
  if (long.group==0) long.group <- 1
  simul <- simulation(axe,nbgroup=long.group,nbchoix=nbchoix,nbsimul=nbsimul)
  if (variability.variable) auxil <- variab.variable(don.interesting,simul$sample,mfa=(!is.null(group)),coord=coord,scale.unit=scale.unit,centerbypanelist=centerbypanelist,scalebypanelist=scalebypanelist,color=color)
dev.new()
  plotellipse(simul,alpha=alpha,coord=coord,eig=signif(axe$eig,4),color=color,cex=cex)  
  res <- list()
  res$eig= axe[[length(names(axe))]]
  res$coordinates= axe[-length(names(axe))]
  if (nbchoix!=1){
    if (length(group)<2) res$hotelling=mat
    if (length(group)>1) res$hotelling=list(bygroup=mat,byproduct=mat2)
  }
  if (variability.variable) res$correl <- auxil
  return(res)
}
