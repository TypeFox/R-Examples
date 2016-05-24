"simulation" <- function(axeAFM,nbchoix=NULL,nbgroup=1,nbsimul=500){

##################################################################################
"simulate.judgement" <- function(axeACP,echantillon){
  axeACP<-as.data.frame(axeACP)
  nbcoord<-ncol(axeACP)-2
  axeACP[,(ncol(axeACP)-1)]<-as.factor(axeACP[,(ncol(axeACP)-1)])
  axeACP[,ncol(axeACP)]<-as.factor(axeACP[,ncol(axeACP)])
  tab.a.diviser <- calculate.average(axeACP)
  tab <- tab.a.diviser$tab
  moy <- tab.a.diviser$moy
  
  nbjuge<-nrow(tab)
  nbprod<-ncol(tab)
  nbsimul <- nrow(echantillon)
  nbchoix <- ncol(echantillon)
  res<-array(0,dim=c(nbsimul,nbprod,nbcoord))
  if (nbchoix!=1) for (j in 1:nbsimul) res[j,,] <- apply(tab[echantillon[j,],,],c(2,3),mean)
  else for (j in 1:nbsimul) res[j,,] <- tab[echantillon[j,1],,]

  aux1 <- aux2 <- NULL
  for (k in 1:nbcoord){
    aux1 <- cbind(aux1,c(tab[,,k]))
    aux2 <- cbind(aux2,c(res[,,k]))
  }
  aux1 <- cbind.data.frame(aux1,rep(axeACP[axeACP[,ncol(axeACP)]==0,ncol(axeACP)-1],rep(nbjuge,nbprod)))
  aux2 <- cbind.data.frame(aux2,rep(axeACP[axeACP[,ncol(axeACP)]==0,ncol(axeACP)-1],rep(nbsimul,nbprod)))
  colnames(aux1)=colnames(aux2)=colnames(cbind.data.frame(t(moy),axeACP[axeACP[,ncol(axeACP)]==0,ncol(axeACP)-1]))
  donnee<-as.data.frame(rbind(cbind.data.frame(t(moy),axeACP[axeACP[,ncol(axeACP)]==0,ncol(axeACP)-1]),aux1,aux2))
  return(donnee)
}
##################################################################################
  nbcoord= ncol(axeAFM$moyen)-2
  nbjuge <- length(levels(axeAFM$moyen[,ncol(axeAFM$moyen)]))-1
  nbprod <- length(levels(axeAFM$moyen[,ncol(axeAFM$moyen)-1]))
  if (length(nbchoix)==0) nbchoix <- nbjuge
#  print(paste("Number of panelists in the virtual panel: ",nbchoix,"  (in the real panel: ",nbjuge,")."))
  echantillon <- matrix(0,nbsimul,nbchoix)
  for (j in 1:nbsimul) echantillon[j,] <- sample(nbjuge,nbchoix,replace=TRUE)

  if (nbgroup==1) simul.moy<-simulate.judgement(axeAFM$moyen,echantillon)
  simulres <- list()
  if (nbgroup>1){
    simul.analyse.partiel<-simulate.judgement(axeAFM$partiel,echantillon)
    nom.des.prod<-simul.analyse.partiel[,ncol(simul.analyse.partiel)]
    simul.analyse.partiel[,ncol(simul.analyse.partiel)] <- as.integer(simul.analyse.partiel[,ncol(simul.analyse.partiel)])
    simul.analyse.partiel <- as.matrix(simul.analyse.partiel)
    simul.moy<-matrix(0,nrow(simul.analyse.partiel),nbcoord) 
    for (i in 1:nbcoord) simul.moy[,i]<-apply(simul.analyse.partiel[,i+nbcoord*(0:(nbgroup-1))],1,mean)
    simul.moy <- cbind.data.frame(simul.moy,as.factor(nom.des.prod))
    ### pour les coord des prod
    simul.partiel <- matrix(0,0,nbcoord)
    for (i in 1:nbgroup)  simul.partiel <- rbind(simul.partiel,simul.analyse.partiel[1:nbprod,((i-1)*nbcoord+1):(i*nbcoord)])
    for (i in 1:nbgroup)  simul.partiel <- rbind(simul.partiel,simul.analyse.partiel[(nbprod+1):(nbprod*(nbchoix+1)),((i-1)*nbcoord+1):(i*nbcoord)])
    for (i in 1:nbgroup)  simul.partiel <- rbind(simul.partiel,simul.analyse.partiel[(nbprod*(nbchoix+1)+1):nrow(simul.analyse.partiel),((i-1)*nbcoord+1):(i*nbcoord)])
    ### pour les noms des prod
    simul.nom.partiel <- NULL
    for (i in 1:nbgroup)  simul.nom.partiel <- rbind(simul.nom.partiel,t(t(as.vector(paste(nom.des.prod[1:nbprod],gsub("Comp1","",colnames(axeAFM$partiel)[((i-1)*nbcoord+1)]),sep="")))))
    for (i in 1:nbgroup)  simul.nom.partiel <- rbind(simul.nom.partiel,t(t(as.vector(paste(nom.des.prod[(nbprod+1):(nbprod*(nbchoix+1))],gsub("Comp1","",colnames(axeAFM$partiel)[((i-1)*nbcoord+1)]),sep="")))))
    for (i in 1:nbgroup)  simul.nom.partiel <- rbind(simul.nom.partiel,t(t(as.vector(paste(nom.des.prod[(nbprod*(nbchoix+1)+1):nrow(simul.analyse.partiel)],gsub("Comp1","",colnames(axeAFM$partiel)[((i-1)*nbcoord+1)]),sep="")))))
    rownames(simul.partiel)<-NULL
    simul.partiel <- cbind.data.frame(simul.partiel,simul.nom.partiel)
    simulres$partiel$P <- simul.partiel[1:(nbprod*nbgroup),]
    simulres$partiel$PJ <- simul.partiel[(nbprod*nbgroup+1):(nbgroup*nbprod*(1+nbjuge)),]
    simulres$partiel$simul <- simul.partiel[(nbgroup*nbprod*(1+nbjuge)+1):nrow(simul.partiel),]
  }
  simulres$moy$P <- simul.moy[1:nbprod,]
  simulres$moy$PJ <- simul.moy[(nbprod+1):(nbprod*(1+nbjuge)),]
  simulres$moy$simul <- simul.moy[(1+nbprod*(nbjuge+1)):nrow(simul.moy),]
  simulres$sample <- echantillon
  return(simulres)
}
