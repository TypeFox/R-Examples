#______________________________________________________________
#Etapa 1
#_____________________________________________________________
#####
Ordenar<-function(IDm,CausaD)
{
  IDm <- as.integer(IDm)
  CausaD <- as.character(CausaD)
  CapA<-DefCap(cau=CausaD)
  CapAut<-as.integer(CapA[[1]])
  codcap<-as.data.frame(cbind(IDm,CausaD,CapAut))
  cat("~~~ Ordenando Datos                               ~~~ \n")
  codcapor<-codcap[order(codcap$CausaD),]
  Id<-1:dim(codcapor)[1]
  IDo <- as.integer(as.character(codcapor[,1]))
  codord<-cbind(codcapor,Id,IDo)
  Tem<-codord[codord$CapAut==1,]
  PN<-c(1,dim(Tem)[1])
  for (i in 2:20){
    Tem<-codord[codord$CapAut==i,]
    PN<-rbind(PN,c(i,dim(Tem)[1]))}
  list(codord,PN)
}


OptFact<-function(V1,Ns,n)
  {
  cat("~~~ Inicia etapa de muestra                       ~~~ \n")
  #set.seed(2013)
  set.seed(as.integer(format(Sys.time(), "%Y")))
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

optn<-function(N,p,E)
    {
      z<-1.96
      S<-z^2*p*(1-p)*N
      I<-E^2*(N-1)+z^2*p*(1-p)
      n<-as.integer(S/I)
      return(n)
    }

#####
#_____________________________________________________________
#Etapa 2
#_______________________________________________________________
#####
Revic<-function(CAUSADEF,RECODCBD,RECODCBD2){
dimn<-length(CAUSADEF)
#dimn
Valor3<-rep(0,dimn)
Valor4<-rep(0,dimn)
Bien3<-rep(0,dimn)
Bien<-rep(0,dimn)
Rev<-rep(0,dimn)
Caus<-as.character(CAUSADEF)
Rec1<-as.character(RECODCBD)
Rec2<-as.character(RECODCBD2)

for (i in 1:dimn)
{
  Ca<-substr(Caus[i],1,3)
  R1<-substr(Rec1[i],1,3)
  R2<-substr(Rec2[i],1,3)
  if (Ca==R1 & R1==R2) {Valor3[i]<-1;Bien3[i]<-1}
  else if (Ca==R1 & R1!=R2) {Valor3[i]<-2}
  else if (Ca==R2 & R1!=R2) {Valor3[i]<-3}
  else if (Ca!=R1 & R1==R2) {Valor3[i]<-4}
  else if ((Ca!=R1) & (R1!=R2 & Ca!=R2)) {Valor3[i]<-5}
  else Valor[i]<-6
}

for (i in 1:dimn)
{
  Ca<-substr(Caus[i],1,4)
  if (substr(Caus[i],4,4)=="X") {Ca<-substr(Caus[i],1,3)}
  R1<-substr(Rec1[i],1,4)
  if (substr(Rec1[i],4,4)=="X") {R1<-substr(Rec1[i],1,3)}
  R2<-substr(Rec2[i],1,4)
  if (substr(Rec2[i],4,4)=="X") {R2<-substr(Rec2[i],1,3)}
  if (Ca==R1 & R1==R2) {Valor4[i]<-1;Bien[i]<-1}
  else if (Ca==R1 & R1!=R2) {Valor4[i]<-2;Rev[i]<-1}
  else if (Ca==R2 & R1!=R2) {Valor4[i]<-3;Rev[i]<-1}
  else if (Ca!=R1 & R1==R2) {Valor4[i]<-4}
  else if ((Ca!=R1) & (R1!=R2 & Ca!=R2)) {Valor4[i]<-5;Rev[i]<-1}
  else Valor[i]<-6
}

Valor3<-as.integer(Valor3)
Valor4<-as.integer(Valor4)
Bien3<-as.integer(Bien3)
Bien<-as.integer(Bien)
Rev<-as.integer(Rev)
Cap1<-DefCap(CAUSADEF)
Cap<-DefCap(RECODCBD)
CapAut<-as.integer(Cap1[[1]])
ManualD<-as.integer(Cap[[1]])
Etapa4<-cbind(Caus,Rec1,Rec2,CapAut,ManualD,Valor3,Valor4,Bien3,Bien,Rev)
#Etapa4<-cbind(CapAut,ManualD,Valor3,Valor4,Bien3,Bien,Rev)
return(Etapa4)
}


Frecu<-function(TabFrec){
  d<-dim(TabFrec)
  #cat("TabFrec",d)
  PosDif<-matrix(NA, d[1], 50)
  for (i in 1:d[1]){
    PosDif[i,1]<-i
    PosDif[i,2]<-sum(TabFrec[i,])
    p<-3
    for (j in 1:d[2]) {
      if (TabFrec[i,j]!=0) {p<-p+1
                            PosDif[i,p]<-j
                            #cat("Valor de j ",j,"\n")
      }
    }
    PosDif[i,3]<-p-3
  }
  ColFin<-max(PosDif[,3])+2
  TabFrec2<-matrix(NA, d[1]*2, ColFin)
  rown<-rownames(TabFrec)
  coln<-colnames(TabFrec)
  TabFrec2[1,1]<-rown[1]
  for (j in 1:PosDif[1,3]){
  TabFrec2[1,j+1]<-coln[PosDif[1,j+3]]
  TabFrec2[2,j+1]<-TabFrec[1,PosDif[1,j+3]]
  TabFrec2[1,ColFin]<-PosDif[1,2]
  TabFrec2[2,ColFin]<-PosDif[1,2]
  }
  for (i in 2:d[1]){
  TabFrec2[i+(i-1),1]<-rown[i]
  for (j in 1:PosDif[i,3]){
    TabFrec2[i+(i-1),j+1]<-coln[PosDif[i,j+3]]
    TabFrec2[i+(i-1)+1,j+1]<-TabFrec[i,PosDif[i,j+3]]
    }
  TabFrec2[i+(i-1),ColFin]<-PosDif[i,2]
  TabFrec2[i+(i-1)+1,ColFin]<-PosDif[i,2]
  }
return(TabFrec2)
}

InterVal<-function(CapAutBien,Pob,Error){
  #INTERVALOS DE CONFIANZA
  Muestra<-rep(0,21)
  Bien<-Muestra;BienT<-Muestra;Poblacion<-Muestra;P<-Muestra;Pn<-Muestra;FactorExp<-Muestra
  
  for (i in 1:20)
  {TemCap<-CapAutBien[CapAutBien[,1]==i,,drop = FALSE]
   #TemCap<-as.matrix(TemCap)
   Dim_Cap<-dim(TemCap)
   Muestra[i]<-Dim_Cap[1]
   #cat("Muestra ",i," es ",Muestra[i]," dim ",Dim_Cap[2], " \n")
   if (Dim_Cap[2]==1){
     Bien[i]<-sum(TemCap[2,1])
     Muestra[i]<-1
     #cat("Fue 1")
   } else {
     Bien[i]<-sum(TemCap[,2])
     #cat("No es 1 ")
   }
   FactorExp[i]<-Pob[i]/Muestra[i]
   Pn[i]<-(Bien[i]/Muestra[i])#*100
   if (Dim_Cap[2]==1){
     BienT[i]<-sum(TemCap[2,1])*FactorExp[i]
   } else {
     BienT[i]<-sum(TemCap[,2])*FactorExp[i]
   }
   
   Poblacion[i]<-Pob[i]
   P[i]<-(BienT[i]/Pob[i])#*100
  }

  i<-21
  Muestra[i]<-nrow(CapAutBien)
  Bien[i]<-sum(CapAutBien[,2])
  Pn[i]<-(Bien[i]/Muestra[i])#*100
  Poblacion[i]<-sum(Pob)
  BienT[i]<-sum(BienT,na.rm = T)
  P[i]<-(BienT[i]/Poblacion[i])#*100
  Cap<-c(1:21)
  
  psup<-(qbeta(1-Error/2,Bien+.5,Muestra-Bien+.5))#*100
  pinf<-(qbeta(Error/2,Bien+.5,Muestra-Bien+.5))#*100
  psup2<-P+(psup-Pn)
  pinf2<-P-(Pn-pinf)
  pinf3<-P-(1-(Muestra/Poblacion))^(1/2)*(P-pinf2)
  psup3<-P+(1-(Muestra/Poblacion))^(1/2)*(psup2-P)
  Int.Conf<-cbind(Cap,Poblacion,BienT,Muestra,Bien,Pn,P,pinf3,psup3)
  Int.Conf<-Int.Conf[!is.na(Int.Conf[,3]),]
  return(Int.Conf)
}

#___________________________________________________________________
#                  Etapa 4 y 5
#__________________________________________________________________

RevicE4<-function(CAUSADEF,RECODCBD,RECODCBD2,COD_SEL){
dimn<-length(CAUSADEF)
#dimn
Fin3<-rep(0,dimn)
Fin4<-rep(0,dimn)
CAUSADEF<-as.character(CAUSADEF)
CauFin<-as.character(CAUSADEF)
COD_SEL<-as.character(COD_SEL)

for (i in 1:dimn)
{
  Ca<-substr(CAUSADEF[i],1,3)
  CF<-substr(COD_SEL[i],1,3)
  if (is.na(CF)==FALSE) {if (Ca==CF){Fin3[i]<-1}}
}

for (i in 1:dimn)
{
  Ca<-substr(CAUSADEF[i],1,4)
  if (substr(CAUSADEF[i],4,4)=="X") {Ca<-substr(CAUSADEF[i],1,3)}
  CF<-substr(COD_SEL[i],1,4)
  if (is.na(CF)==FALSE){CauFin[i]<-CF}
  if (is.na(CF)==FALSE){if (substr(COD_SEL[i],4,4)=="X"){CF<-substr(COD_SEL[i],1,3)}
                        else if (Ca==CF){Fin4[i]<-1}}
}

Fin3<-as.integer(Fin3)
Fin4<-as.integer(Fin4)

Valor3<-rep(0,dimn)
Valor4<-rep(0,dimn)
Bien3<-rep(0,dimn)
Bien4<-rep(0,dimn)
Rev<-rep(0,dimn)
#CAUSADEF<-as.character(CAUSADEF)
RECODCBD<-as.character(RECODCBD)
RECODCBD2<-as.character(RECODCBD2)

for (i in 1:dimn)
{
  Ca<-substr(CAUSADEF[i],1,3)
  R1<-substr(RECODCBD[i],1,3)
  R2<-substr(RECODCBD2[i],1,3)
  if (Ca==R1 & R1==R2) {Valor3[i]<-1;Bien3[i]<-1}
  else if (Ca==R1 & R1!=R2) {Valor3[i]<-2}
  else if (Ca==R2 & R1!=R2) {Valor3[i]<-3}
  else if (Ca!=R1 & R1==R2) {Valor3[i]<-4}
  else if (Ca!=R1 & R1!=R2) {Valor3[i]<-5}
  else Valor[i]<-6
}

for (i in 1:dimn)
{
  Ca<-substr(CAUSADEF[i],1,4)
  if (substr(CAUSADEF[i],4,4)=="X") {Ca<-substr(CAUSADEF[i],1,3)}
  R1<-substr(RECODCBD[i],1,4)
  if (substr(RECODCBD[i],4,4)=="X") {R1<-substr(RECODCBD[i],1,3)}
  R2<-substr(RECODCBD2[i],1,4)
  if (substr(RECODCBD2[i],4,4)=="X") {R2<-substr(RECODCBD2[i],1,3)}
  if (Ca==R1 & R1==R2) {Valor4[i]<-1;Bien4[i]<-1}
  else if (Ca==R1 & R1!=R2) {Valor4[i]<-2;Rev[i]<-1}
  else if (Ca==R2 & R1!=R2) {Valor4[i]<-3;Rev[i]<-1}
  else if (Ca!=R1 & R1==R2) {Valor4[i]<-4;CauFin[i]<-R1}
  else if (Ca!=R1 & R1!=R2) {Valor4[i]<-5;Rev[i]<-1}
  else Valor[i]<-6
}

Valor3<-as.integer(Valor3)
Valor4<-as.integer(Valor4)
Bien3<-as.integer(Bien3)
Bien4<-as.integer(Bien4)
Rev2<-as.integer(Rev)
Cap1<-DefCap(CAUSADEF)
Cap<-DefCap(RECODCBD)
CapAut<-as.integer(Cap1[[1]])
ManualD<-as.integer(Cap[[1]])

Cap<-DefCap(CauFin)
CapFin<-as.integer(Cap[[1]])

BienFin4<-rep(0,dimn)
BienFin4[CauFin==CAUSADEF]<-1
BienFin4<-as.integer(BienFin4)
BienFin3<-rep(0,dimn)
BienFin3[substr(CauFin,1,3)==substr(CAUSADEF,1,3)]<-1
BienFin3<-as.integer(BienFin3)

#Etapa5rev<-cbind(CAUSADEF,RECODCBD,RECODCBD2,COD_SEL,CauFin,Valor3,Valor4,BienFin3,BienFin4,Rev2,Bien3,Bien4,Fin3,Fin4,CapAut,ManualD,CapFin)

Fin3c1<-rep(0,dimn)
Fin4c1<-rep(0,dimn)
RECODCBD<-as.character(RECODCBD)
CodFin<-CauFin

for (i in 1:dimn)
{
  Ca<-substr(RECODCBD[i],1,3)
  CF<-substr(CodFin[i],1,3)
  if (Ca==CF) {Fin3c1[i]<-1} 
  else {Fin3c1[i]<-0}
}

for (i in 1:dimn)
{
  Ca<-substr(RECODCBD[i],1,4)
  if (substr(RECODCBD[i],4,4)=="X") {Ca<-substr(RECODCBD[i],1,3)}
  CF<-substr(CodFin[i],1,4)
  if (substr(CodFin[i],4,4)=="X") {CF<-substr(CodFin[i],1,3)}
  if (Ca==CF) {Fin4c1[i]<-1}
  else {Fin4c1[i]<-0}
}

Fin3c1<-as.integer(Fin3c1)
Fin4c1<-as.integer(Fin4c1)

Fin3c2<-rep(0,dimn)
Fin4c2<-rep(0,dimn)
RECODCBD2<-as.character(RECODCBD2)

for (i in 1:dimn)
{
  Ca<-substr(RECODCBD2[i],1,3)
  CF<-substr(CodFin[i],1,3)
  if (Ca==CF) {Fin3c2[i]<-1}
  else {Fin3c2[i]<-0}
}

for (i in 1:dimn)
{
  Ca<-substr(RECODCBD2[i],1,4)
  if (substr(RECODCBD2[i],4,4)=="X") {Ca<-substr(RECODCBD2[i],1,3)}
  CF<-substr(CodFin[i],1,4)
  if (substr(CodFin[i],4,4)=="X") {CF<-substr(CodFin[i],1,3)}
  if (Ca==CF) {Fin4c2[i]<-1} 
  else {Fin4c2[i]<-0}
}

Fin3c2<-as.integer(Fin3c2)
Fin4c2<-as.integer(Fin4c2)

#####################################
##Para compararar los codificadores manuales C1 y C2 con el 
##codificador automatico

Aut3c1<-rep(0,dimn)
Aut4c1<-rep(0,dimn)
RECODCBD<-as.character(RECODCBD)
CAUSADEF<-as.character(CAUSADEF)

for (i in 1:dimn)
{
  Ca<-substr(RECODCBD[i],1,3)
  CF<-substr(CAUSADEF[i],1,3)
  if (Ca==CF) {Aut3c1[i]<-1} 
  else {Aut3c1[i]<-0}
}

for (i in 1:dimn)
{
  Ca<-substr(RECODCBD[i],1,4)
  if (substr(RECODCBD[i],4,4)=="X") {Ca<-substr(RECODCBD[i],1,3)}
  CF<-substr(CAUSADEF[i],1,4)
  if (substr(CAUSADEF[i],4,4)=="X") {CF<-substr(CAUSADEF[i],1,3)}
  if (Ca==CF) {Aut4c1[i]<-1} 
  else {Aut4c1[i]<-0}
}

Aut3c1<-as.integer(Aut3c1)
Aut4c1<-as.integer(Aut4c1)

Aut3c2<-rep(0,dimn)
Aut4c2<-rep(0,dimn)
RECODCBD2<-as.character(RECODCBD2)

for (i in 1:dimn)
{
  Ca<-substr(RECODCBD2[i],1,3)
  CF<-substr(CAUSADEF[i],1,3)
  if (Ca==CF) {Aut3c2[i]<-1}
  else {Aut3c2[i]<-0}
}

for (i in 1:dimn)
{
  Ca<-substr(RECODCBD2[i],1,4)
  if (substr(RECODCBD2[i],4,4)=="X") {Ca<-substr(RECODCBD2[i],1,3)}
  CF<-substr(CAUSADEF[i],1,4)
  if (substr(CAUSADEF[i],4,4)=="X") {CF<-substr(CAUSADEF[i],1,3)}
  if (Ca==CF) {Aut4c2[i]<-1}
  else {Aut4c2[i]<-0}
}

Aut3c2<-as.integer(Aut3c2)
Aut4c2<-as.integer(Aut4c2)
Etapa5rev<-cbind(CAUSADEF,RECODCBD,RECODCBD2,COD_SEL,CauFin,Valor3,Valor4,BienFin3,BienFin4,Rev2,Bien3,Bien4,Fin3,Fin4,CapAut,ManualD,CapFin,Fin3c1,Fin4c1,Fin3c2,Fin4c2,Aut3c1,Aut4c1,Aut3c2,Aut4c2)
return(Etapa5rev)
}

PonerFact<-function(Base,Ns){
  #Base_ord<-data.frame(Base[order(Base[,1]),])
  Base_ord<-Base[order(Base[,1]),]
  cat("Dim orde",dim(Base_ord),"\n")
  
  #i<-1 #Esto es para el primer capitulo
  #n<-0
  n<-rep(0,20)
  Tempo<-subset(Base_ord,as.integer(Base_ord[,15])==1)#nrow(Base_ord[Base_ord[,1]==i,])
  n[1]<-nrow(Tempo)
  cat("Dim n 1 ",dim(n),"Valor Ns ",Ns[1],"Valor n ",n[1], "\n")
  Capit<-as.integer(rep(1,n[1]))
  FactorExp<-rep(Ns[1]/n[1],n[1])
  #cat("Capit ",Capit," Fsa ",FactorExp, "\n")
  CapFact<-data.frame(Capit,FactorExp)
  #cat("Dim CapFact 1",dim(CapFact),"\n")
  
  Ac<-Ns[1]
  for (i in 2:20) {Ac[i]<-Ac[i-1]+Ns[i]}
  
  for (i in 2:20){
    n[i]<-nrow(Base_ord[as.integer(Base_ord[,15])==i,])
    if (n[i]!=0){
      N<-Ns[i];m<-n[i]
      Capit<-as.integer(rep(i,m))
      FactorExp<-rep(N/m,m)
      Mu<-data.frame(Capit,FactorExp)
      cat("lugar",i, "dim Mu ",dim(Mu),"\n")
      CapFact<-rbind(CapFact,Mu)}}
  
  BaseFacExp<-data.frame(Base_ord,CapFact)
  cat("dim ",dim(BaseFacExp))
  return(BaseFacExp)
}

#Crear tablas de ponderados
TabPon<-function(CapAut_CapFin){
  #DimT<-dim(Tabla)
  #Muestra<-rep(0,DimT[1])
  #FactorExp<-rep(0,DimT[1])
  #FactorExp<-rep(0,20)
  #Tabla<-as.matrix(Tabla)
  #MatPon<-matrix(DimT[1], DimT[2])
  #Pob<-Pob[Pob!=0]
  #cat("Tabla dim",DimT,"Pob",Pob)
  #for (i in 1:DimT[1]){
  #  Muestra[i]<-length(PaMue[PaMue==i])
  #  FactorExp[i]<-Pob[i]/Muestra[i]
  #  cat("Muesttra",Muestra,"Factor",FactorExp)
  #  for (j in 1:DimT[2]){MatPon[i,j]<-as.integer(Tabla[i,j])*FactorExp[i]}
  #}
  #CapAut_CapFin<-as.data.frame(CapAut_CapFin)
  #MatPon<-matrix(nrow = 20, ncol = 20)
  #for (i in 1:20)
  #{
  #  Tempo<-subset(CapAut_CapFin,CapAut_CapFin[,1]==i)
  #  #Tempo<-subset(CapAut_CapFin,CapAut==i)
  #  FactorExp<-as.integer(Pob[i]/nrow(Tempo))
  #  cat("Tempo 1",Tempo[,1],"Tempo 2",Tempo[,2], "factor",FactorExp," i",i,"\n")
  #  for (j in 1:20)
  #  {line<-subset(Tempo,Tempo[,2]==j)
  #    cont<-nrow(line)
  #   cat("line",line,"  Cont",cont," j",j,"\n")
  #    MatPon[i,j]<-cont*FactorExp}
  #}
  
  
  MatPon<-matrix(nrow = 20, ncol = 20)
  for (i in 1:20)
  {
    Tempo<-CapAut_CapFin[CapAut_CapFin[,1]==i,]
    #Tempo<-subset(CapAut_CapFin,CapAut_CapFin[,1]==i)
    cat("  i",i, "dim Tempo",dim(Tempo),"head",head(Tempo), "\n")
    for (j in 1:20)
    {line<-subset(Tempo,Tempo[,2]==j)
     #cat("j",j, "line",line[,3],"\n")
      MatPon[i,j]<-sum(line[,3],na.rm = T)}
  }
  cat("MatPon",MatPon)
  return(MatPon)
}


SalvarDatos <- function(Dat_Def,Cod_fin,Cod_Tex) {#
  Dat_Def[,"COD_SEL"]<-Cod_fin
  Dat_Def[,"Cap_Tex"]<-Cod_Tex
  saveRDS(Dat_Def, "Dat_Def.rds")
  }

CargarDatos <- function() {
  Dat<-readRDS("Dat_Def.rds")
  Dat
}

CapDef <- function(cau="character")
{
    ca<-substr(cau,1,3)
    if ("A00"<=ca & ca<="B99") {Capitulo<-1
    Cau_Def_Tex<-"Cap\u00edtulo I Enfermedades infecciosas y parasitarias"}
    else if ("C00"<=ca & ca<="D48") {Capitulo<-2
    Cau_Def_Tex<-"Cap\u00edtulo II Neoplasias"}
    else if ("D50"<=ca & ca<="D89") {Capitulo<-3
    Cau_Def_Tex<-"Cap\u00edtulo III Enfermedades de la sangre y de los \u00f3rganos hematopoy\u00e9ticos"}
    else if ("E00"<=ca & ca<="E90") {Capitulo<-4
    Cau_Def_Tex<-"Cap\u00edtulo IV Enfermedades endocrinas, nutricionales y metab\u00f3licas."}
    else if ("F00"<=ca & ca<="F99") {Capitulo<-5
    Cau_Def_Tex<-"Cap\u00edtulo V Transtornos mentales y del comportamiento."}
    else if ("G00"<=ca & ca<="G99") {Capitulo<-6
    Cau_Def_Tex<-"Cap\u00edtulo VI Enfermedades del sistema nervioso."}
    else if ("H00"<=ca & ca<="H59") {Capitulo<-7
    Cau_Def_Tex<-"Cap\u00edtulo VII Enfermedades del ojo y sus anexos."}
    else if ("H60"<=ca & ca<="H95") {Capitulo<-8
    Cau_Def_Tex<-"Cap\u00edtulo VIII Enfermedades del o\u00eddo y de la ap\u00f3fisis mastoides"}
    else if ("I00"<=ca & ca<="I99") {Capitulo<-9
    Cau_Def_Tex<-"Cap\u00edtulo IX Enfermedades del sistema cardiocirculatorio"}
    else if ("J00"<=ca & ca<="J99") {Capitulo<-10
    Cau_Def_Tex<-"Cap\u00edtulo X Enfermedades del sistema respiratorio"
    }
    else if ("K00"<=ca & ca<="K93") {Capitulo<-11
    Cau_Def_Tex<-"Cap\u00edtulo XI Enfermedades del sistema digestivo"
    }
    else if ("L00"<=ca & ca<="L99") {Capitulo<-12
    Cau_Def_Tex<-"Cap\u00edtulo XII Enfermedades de la piel y tejido subcut\u00e1neo"
    }
    else if ("M00"<=ca & ca<="M99") {Capitulo<-13
    Cau_Def_Tex<-"Cap\u00edtulo XIII Enfermedades del sistema osteomuscular y del tejido conjuntivo."
    }
    else if ("N00"<=ca & ca<="N99") {Capitulo<-14
    Cau_Def_Tex<-"Cap\u00edtulo XIV Enfermedades del sistema genitourinario"
    }
    else if ("O00"<=ca & ca<="O99") {Capitulo<-15
    Cau_Def_Tex<-"Cap\u00edtulo XV Enfermedades del embarazo, parto y puerperio"
    }
    else if ("P00"<=ca & ca<="P96") {Capitulo<-16
    Cau_Def_Tex<-"Cap\u00edtulo XVI Enfermedades del feto y reci\u00e9n nacido"
    }
    else if ("Q00"<=ca & ca<="Q99") {Capitulo<-17
    Cau_Def_Tex<-"Cap\u00edtulo XVII Enfermedades cong\u00e9nitas, malformaciones y alteraciones cromos\u00f3micas"
    }
    else if ("R00"<=ca & ca<="R99") {Capitulo<-18
    Cau_Def_Tex<-"Cap\u00edtulo XVIII S\u00edntomas y observaciones cl\u00ednicas o de laboratorio anormales no clasificados en otras parte"
    }
    else if ("S00"<=ca & ca<="T98") {Capitulo<-19
    Cau_Def_Tex<-"Cap\u00edtulo XIX Lesiones, heridas, intoxicaciones y otros factores externos"
    }
    else if ("V01"<=ca & ca<="Y98") {Capitulo<-20
    Cau_Def_Tex<-"Cap\u00edtulo XX Causas externas de mortalidad y morbilidad"
    }
    else if ("Z00"<=ca & ca<="Z99") {Capitulo<-21
    Cau_Def_Tex<-"Cap\u00edtulo XXI Factores que influyen en el estado de salud y contacto con los servicios de salud."
    }
    else if ("U00"<=ca & ca<="U99") {Capitulo<-22
    Cau_Def_Tex<-"Cap\u00edtulo XXII C\u00f3digos para prop\u00f3sitos especiales"
    }
    else Capitulo<-23
    return(list(Capitulo,Cau_Def_Tex))
#return(Capitulo,Cau_Def_Tex)
}

