#Etapa 1

DefCap <- function(cau="character")
{
  d<-length(cau)
  cat("~~~ Lectura de ",d, " c\u00F3digos de defunci\u00F3n      ~~~ \n")
  Capitulo<-rep(0,d)
  for (i in 1:d)
  {
    ca<-substr(cau[i],1,3)
    if ("A00"<=ca & ca<="B99") {Capitulo[i]<-1}
    else if ("C00"<=ca & ca<="D48") {Capitulo[i]<-2}
    else if ("D50"<=ca & ca<="D89") {Capitulo[i]<-3}
    else if ("E00"<=ca & ca<="E90") {Capitulo[i]<-4}
    else if ("F00"<=ca & ca<="F99") {Capitulo[i]<-5}
    else if ("G00"<=ca & ca<="G99") {Capitulo[i]<-6}
    else if ("H00"<=ca & ca<="H59") {Capitulo[i]<-7}
    else if ("H60"<=ca & ca<="H95") {Capitulo[i]<-8}
    else if ("I00"<=ca & ca<="I99") {Capitulo[i]<-9}
    else if ("J00"<=ca & ca<="J99") {Capitulo[i]<-10}
    else if ("K00"<=ca & ca<="K93") {Capitulo[i]<-11}
    else if ("L00"<=ca & ca<="L99") {Capitulo[i]<-12}
    else if ("M00"<=ca & ca<="M99") {Capitulo[i]<-13}
    else if ("N00"<=ca & ca<="N99") {Capitulo[i]<-14}
    else if ("O00"<=ca & ca<="O99") {Capitulo[i]<-15}
    else if ("P00"<=ca & ca<="P96") {Capitulo[i]<-16}
    else if ("Q00"<=ca & ca<="Q99") {Capitulo[i]<-17}
    else if ("R00"<=ca & ca<="R99") {Capitulo[i]<-18}
    else if ("S00"<=ca & ca<="T98") {Capitulo[i]<-19}
    else if ("V01"<=ca & ca<="Y98") {Capitulo[i]<-20}
    else if ("Z00"<=ca & ca<="Z99") {Capitulo[i]<-21}
    else if ("U00"<=ca & ca<="U99") {Capitulo[i]<-22}
    else Capitulo[i]<-23
  }
  Error<-sum(Capitulo==23)
  cat("~~~ Finalizo lectura de los c\u00F3digos con",Error, "Errores ~~~ \n")
  list(Capitulo,Error)
}

OptFact<-function(V1,Ns,n)
  {
  cat("~~~ Inicia etapa de muestra                       ~~~ \n")
  t<-0;pr<-0
  while (pr==0){t<-t+1
                if (Ns[t]>0){
                  NT<-c(1:Ns[t])
                  Capit<-as.integer(rep(V1[t],Ns[t]))
                  FactorExp<-rep(Ns[t]/n[t],Ns[t])
                  EnMuestra<-as.integer(srswor(n[t],Ns[t]))
                  MuestraGr<-data.frame(NT,Capit,FactorExp,EnMuestra)
                  cat(t," ")
                  Ac<-Ns[t]
                  pr<-1}}
  for (i in (t+1):20) {Ac[i]<-Ac[i-1]+Ns[i]}
  for (i in (t+1):20)
    {
      if (Ns[i]>0){
        NT<-c((Ac[i-1]+1):Ac[i])
        N<-Ns[i];m<-n[i]
        Capit<-as.integer(rep(V1[i],N))
        FactorExp<-rep(N/m,N)
        EnMuestra<-as.integer(srswor(m,N))
        Mu<-data.frame(NT,Capit,FactorExp,EnMuestra)
        MuestraGr<-rbind(MuestraGr,Mu)
        cat(i," ")}
    }
  cat("~~ Fin ~~~ \n")
  return(MuestraGr)
  }

Etapa01<-function(IDm,CausaD,Cap,Es,ps,n,Grande){
  cat("~~~ Probando dimenciones                          ~~~ \n")
  if(length(IDm)!=length(CausaD)){
    stop ("No se puede seguir, numero de ID y Casos de defunci\u00F3n no coinciden")}
  IDm <- as.integer(IDm)
  CausaD <- as.character(CausaD)
  CapA<-DefCap(cau=CausaD)
  CapAut<-as.integer(CapA[[1]])
  codcap<-as.data.frame(cbind(IDm,CausaD,CapAut))
  cat("~~~ Ordenando Datos                               ~~~ \n")
  codcapor<-codcap[order(CausaD),]
  Id<-1:dim(codcapor)[1]
  codord<-cbind(codcapor,Id)
  Tem<-codord[codord$CapAut==1,]
  PN<-c(1,dim(Tem)[1])
  for (i in 2:20){
    Tem<-codord[codord$CapAut==i,]
    PN<-rbind(PN,c(i,dim(Tem)[1]))}
  if(missing(Cap)){Cap <- 1:20}
  else{Cap<-as.integer(Cap)}
  if(missing(Es)){Es <- rep(.03,20)}
  else{Es <- Es}
  if(missing(ps)){ps <- rep(.5,20)}
  else{ps <- ps}
  if(missing(n))
  {
    Ns<-PN[,2]
    z<-rep(1.96, 20)
    S<-z^2*ps*(1-ps)*Ns
    I<-Es^2*(Ns-1)+z^2*ps*(1-ps)
    n<-as.integer(S/I)
  }
  else{n <- as.integer(n)}
  MuestraGr<-OptFact(Cap,PN[,2],n)
  Muestra<-merge(codord, MuestraGr, by.x="Id", by.y="NT")
  if(missing(Grande)){
  Mu<-Muestra
  Muestra<-Mu[Mu$EnMuestra==1,]}
  return(Muestra)
}