"ardi" <- function(donnee,col.p,col.j,firstvar,lastvar=ncol(donnee),nbval=10,center=TRUE,scale=FALSE){

########################################################################
searchvalue <- function(mat,nbval){
  aa <- rev(order(mat))[1:nbval]
  bb <- (aa[1]-0.1)%/%nrow(mat)
  res <- data.frame(colnames(mat)[bb+1],rownames(mat)[aa[1]-nrow(mat)*bb])
  for (i in 1:2) res[,i]=as.character(res[,i])
  for (i in 2:nbval){
    bb <- (aa[i]-0.1)%/%nrow(mat)
    res <- rbind.data.frame(res,c(colnames(mat)[bb+1],rownames(mat)[aa[i]-nrow(mat)*bb]))
  }
  res <- cbind.data.frame(res,mat[aa])
}
########################################################################

for (j in 1:(firstvar-1))  donnee[,j] <- as.factor(donnee[,j])
donnee <- donnee[,1:lastvar]
nbprod <- length(levels(donnee[,col.p]))
nbpanelist <- length(levels(donnee[,col.j]))
labprod <- levels(donnee[,col.p])
labpanelist <- levels(donnee[,col.j])
nbdesc <- lastvar-firstvar+1

  moy.aux <- scalebypanelist(donnee,col.j=col.j,col.p=col.p,firstvar=firstvar,center=center,scale=scale)
###########
  moy <- matrix(apply(cbind(moy.aux[1:nbprod,-(1:2)]),2,mean),nrow=nrow(cbind(moy.aux[-(1:nbprod),-(1:2)])),ncol=ncol(cbind(moy.aux[-(1:nbprod),-(1:2)])),byrow=TRUE)
  ecart <- matrix(apply(cbind(moy.aux[1:nbprod,-(1:2)]),2,var),nrow=nrow(cbind(moy.aux[-(1:nbprod),-(1:2)])),ncol=ncol(cbind(moy.aux[-(1:nbprod),-(1:2)])),byrow=TRUE)*(nbprod-1)/nbprod
  ecart <- sqrt(ecart)
  tab <- moy.aux[-(1:nbprod),-(1:2),]-moy / sqrt(ecart)
  tab <- cbind.data.frame(moy.aux[-(1:nbprod),1:2],abs(tab))  
  colnames(tab)=colnames(moy.aux)
  auxpanelist <- matrix(0,nbpanelist,nbdesc)
  auxprod <- matrix(0,nbprod,nbdesc)
  for (j in 1:nbpanelist) auxpanelist[j,] <- apply(cbind(tab[tab[,1]==labpanelist[j],-(1:2)]),2,mean)
  for (j in 1:nbprod) auxprod[j,] <- apply(cbind(tab[tab[,2]==labprod[j],-(1:2)]),2,mean)
  dimnames(auxpanelist) <- list(labpanelist,colnames(tab)[-(1:2)])
  dimnames(auxprod) <- list(labprod,colnames(tab)[-(1:2)])
  respanelist <- searchvalue(auxpanelist,nbval=min(nbval,nbpanelist*nbdesc))
  colnames(respanelist) <- c("Descriptor","Panelist","Distance to the mean")
  resproduct <- searchvalue(auxprod,nbval=min(nbval,nbprod*nbdesc))
  colnames(resproduct) <- c("Descriptor","Product","Distance to the mean")
  tabinteract <- as.matrix(tab[,-(1:2)])
  dimnames(tabinteract) <- list(paste(tab[,2],"-",tab[,1]),colnames(tab)[-(1:2)])
  resinteract <- searchvalue(tabinteract,nbval=min(nbval,nbprod*nbdesc*nbpanelist))
  colnames(resinteract) <- c("Descriptor","Product-Panelist","Distance to the mean")
  res <- list()
  res$tab <- tab
  res$panelist <- respanelist
  res$product <- resproduct
  res$combination <- resinteract
  return(res)
}
