"construct.axes" <- function(matrice,coord=c(1,2),scale.unit=TRUE,group=NULL,name.group=NULL,centerbypanelist=FALSE,scalebypanelist=FALSE,method="coeff"){

  nbcoord=max(coord)
  oo <- order(matrice[,2])
  matrice <- matrice[oo,]
  oo <- order(matrice[,1])
  matrice <- matrice[oo,]

  nbjuge <- nlevels(matrice[,1])
  if (0%in% summary(matrice[,1])) nbjuge <- nbjuge-1
  nbprod <- length(levels(matrice[,2]))
  nbdesc <- ncol(matrice)-2

  moy.aux=scalebypanelist(matrice,col.j=1,col.p=2,firstvar=3,center=centerbypanelist,scale=scalebypanelist,method=method)
  rownames(moy.aux) <- paste("i",1:nrow(moy.aux),sep="")
  rownames(moy.aux)[1:nbprod] <- as.character(moy.aux[1:nbprod,2])
  ###AF with active data the averages for all the panelist 
  axe <- list()
  if (is.null(group)){
    res.af <- PCA(moy.aux[,-c(1,2)],ind.sup = (nbprod+1):nrow(moy.aux), scale.unit = scale.unit, ncp = nbcoord,graph=FALSE)
    axe$moyen <- data.frame(rbind(res.af$ind$coord,res.af$ind.sup$coord),as.factor(moy.aux[,2]),as.factor(moy.aux[,1]))
    plot(res.af,choix="var",axes=coord)
  }
  else {
    if (scale.unit) res.af <- MFA(moy.aux[,-c(1,2)],ind.sup = (nbprod+1):nrow(moy.aux), group = group, name.group = name.group, type = rep("s",length(group)), ncp = nbcoord,graph=FALSE)
    else res.af <- MFA(moy.aux[,-c(1,2)],ind.sup = (nbprod+1):nrow(moy.aux), group = group, name.group = name.group, type = rep("c",length(group)), ncp = nbcoord,graph=FALSE)
    axe$moyen <- data.frame(rbind(res.af$ind$coord,res.af$ind.sup$coord),as.factor(moy.aux[,2]),as.factor(moy.aux[,1]))
    axe$partiel <- data.frame(rbind(t(matrix(t(as.matrix(res.af$ind$coord.partiel)),nrow=nbcoord*length(group),byrow=FALSE)),t(matrix(t(as.matrix(res.af$ind.sup$coord.partiel)),nrow=nbcoord*length(group),byrow=FALSE))),as.factor(moy.aux[,2]),as.factor(moy.aux[,1]))
    dimnames(axe$partiel)[2][[1]][(dim(axe$partiel)[2]-1):dim(axe$partiel)[2]] <- c("Product","Panelist")
    for (i in 1:length(group)) dimnames(axe$partiel)[2][[1]][((i-1)*nbcoord+1):(i*nbcoord)]<-paste("Dim", 1:nbcoord, sep = "",name.group[i])
    plot(res.af,choix="var",habillage="group",axes=coord)
  }
  plot(res.af,choix="ind",invisible="ind.sup",axes=coord)
  dev.new()
	plot(res.af,choix="var",axes=coord)
  dimnames(axe$moyen)[2][[1]]<-c (paste("Dim", 1:nbcoord, sep = ""),"Product","Panelist")
  axe$eig = res.af$eig
  return(axe)
}
