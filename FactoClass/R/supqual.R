#-------------------------------------------------
# Funcion para calcular coordenadas y ayudas para
# variables suplementarias en ACP
# Campo Elias Pardo Julio 31 de 2015
#------------------------------------------------
supqual<-function(du,qual){
  # pendiente control de entrada
  if (class(du)[1]=="acm") method="mca"
  if (class(du)[1]=="pca") method="pca"
  neje<-paste("Axis",1:du$nf,sep="")
if (method=="pca")  {
  acp<-du  # cambio de nombre 
  Ysupcat<-data.frame(qual) # cambio de nombre
  redo.dudi(acp,acp$rank)->reacp # para calcular distancias
  sup<-NULL
  nq<-ncol(Ysupcat) # numero de variables categoricas
  nrow(Ysupcat)->n  # numero de filas (individuos)
  njs<-NULL
  # ciclo de calculo para cada variable
    for (i in 1:nq){
      nl <- length(attributes(Ysupcat[,i])$levels)  # numero de categorias
      wcat<-centroids(acp$li,Ysupcat[,i])$weights
      sup$wcat<-c(sup$wcat,wcat)
      d2 <- rowSums(centroids(reacp$li,Ysupcat[,i])$centroids^2)
      sup$dis2<-c(sup$dis2,d2) # distancias al cuadrado
      centroids(acp$li,Ysupcat[,i])$centroids->Fsupcat #coordenadas
      sup$coor<-rbind(sup$coor,Fsupcat)
      # numero de individuos en las categorias
      nj<-centroids(acp$li,Ysupcat[,i])$weights*n
      njs<-c(njs,nj)
      # raices de valores propios
      sqrtValP<-matrix(1,nl,1)%*%sqrt(acp$eig[1:acp$nf])
      # valores test
      tv <- sqrt(nj*(n-1)/(n-nj))*Fsupcat/sqrtValP
      sup$tv <-rbind(sup$tv,tv)
      cos2<-Fsupcat^2/d2
      sup$cos2<-rbind(sup$cos2,cos2)
    }  
  # relaciones de correlacion
  ss<-ncol(qual)
  scr<-NULL
  bin<-1; bfi<-0
  encat<-njs/n * sup$coor * sup$coor
  for (q in 1:ss) {
    bfi<-bfi+ length(attributes(qual[,q])$levels)
    scr<-rbind(scr,colSums(encat[bin:bfi,]))
    bin <- bfi+1
  }
  rownames(scr)<-colnames(qual)
  colnames(scr)<-neje
  sup$scr<-scr
  return(sup)
} # if method == "pca"
  #======================= mca
  if (method == "mca"){
    acm<-du
    n<-nrow(acm$tab); s<-ncol(acm$tab)
    qual<-data.frame(qual)
    Zs<-acm.disjonctif(qual)
    njs<-colSums(Zs)
    # calculo de coordenadas con fomula de transicion
    coor<-diag(1/njs)%*%t(Zs)%*%as.matrix(acm$li)%*%diag(1/sqrt(acm$eig[1:acm$nf]))
    # valores test  
    vt<-sqrt(njs*(n-1)/(n-njs))*coor
    # distancias al origen
    dis2<-n/njs-1
    names(dis2)<-colnames(Zs)
    # cosenos cuadrados
    cos2<-diag(1/dis2)%*%(coor*coor)
    # relaciones de correlacion
    ss<-ncol(qual)
    scr<-NULL
    bin<-1; bfi<-0
    encat<-njs/n * coor * coor
    for (q in 1:ss) {
      bfi<-bfi+ length(attributes(qual[,q])$levels)
      scr<-rbind(scr,colSums(encat[bin:bfi,]))
      bin <- bfi+1
    }
    #salida
    colnames(coor)<-colnames(vt)<-colnames(cos2)<-colnames(scr)<-neje
    rownames(coor)<-rownames(vt)<-rownames(cos2)<-colnames(Zs)
    rownames(scr)<-colnames(qual)
    sup<-list(ncat=njs,dis2=dis2,coor=coor,tv=vt,cos2=cos2,scr=scr)
    return(sup)
  }
}
