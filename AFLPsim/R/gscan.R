gscan <-
function(mat, type=c("F1","BxA","BxB"), method=c("bal&gar-ca","gagnaire")){
if (method=="bal&gar-ca" && type=="F1"){
pa<-apply(mat$PA,2,mean)
pb<-apply(mat$PB,2,mean)
pf1<-apply(mat$F1,2,mean)
Na<-dim(mat$PA)[1]
Nb<-dim(mat$PB)[1]
Nf1<-dim(mat$F1)[1]
Nmarker<-dim(mat$PA)[2]

#Calcula los intervalos de confianza de la binomial
p.L <- function(x, n, alpha) {  
  if (x == 0)   0
  else qbeta(alpha, x, n - x + 1)
  }
p.U <- function(x,n, alpha) {  
  if (x == n)    1
  else qbeta(1 - alpha, x + 1, n - x)
  }

#Calcula las coordenadas x e y de las esquinas de la sabana

xU<-sapply(pa*Na,p.U,n=Na, alpha=0.02532)
xL<-sapply(pa*Na,p.L,n=Na, alpha=0.02532)
yL<-sapply(pb*Nb,p.L,n=Nb, alpha=0.02532)
yU<-sapply(pb*Nb,p.U,n=Nb, alpha=0.02532)

#Modelo de frecuencias esperadas en F1
model<-function(x,y){
  a<-1-sqrt(1-x)
  b<-1-sqrt(1-y)
  a+b-a*b}

#Calcula el z para cada esquina
zUU<-model(xU, yU)
zUL<-model(xU, yL)
zLL<-model(xL, yL)
zLU<-model(xL, yU)

# Calcula la diferencia entre el valor de z estimado y el teorico 

DUU<-pf1-zUU
DLL<-pf1-zLL
DLU<-pf1-zLU
DUL<-pf1-zUL

#Crea una dataframe con todos los valores anteriores para su uso en la funcion

datos<-cbind(pf1*Nf1, pf1,zLL,zUU,zLU,zUL,DLL,DUU,DLU,DUL )
datos<-as.data.frame(datos)

colnames(datos)<- c("x","z","zLL","zUU","zLU","zUL","dLL","dUU","dLU","dUL")

#Funcion  que calcula la p de outlier en un fragmento basada en un test binomial eligiendo la esquina de la sabana mas cercana al valor estimado
zbinomial<-function(x,z,zLL,zUU,zLU,zUL,dLL,dUU,dLU,dUL) {
  if (sign(dLL)== sign(dUU) && sign(dUU) == sign(dLU) && sign(dLU) == sign(dUL) && z > zLL)
  {maxz<-max(zLL,zUU,zLU,zUL); bi<-binom.test(x,n=Nf1, p=maxz);print(bi$p.value)}
  else if (sign(dLL)== sign(dUU) && sign(dUU) == sign(dLU) && sign(dLU) == sign(dUL) && z < zLL)
  {minz<-min(zLL,zUU,zLU,zUL); bi<-binom.test(x,n=Nf1, p=minz);print(bi$p.value)}
  else {1}}


#Loop para correr la funcion zbinomial para cada uno de los fragmentos
marker <-Nmarker #numero de markers#
pvaluebinomf1<-numeric(length(marker))
for(i in 1:marker){
  pvaluebinomf1[i]<-zbinomial(datos$x[i],datos$z[i],datos$zLL[i],datos$zUU[i],datos$zLU[i],datos$zUL[i],datos$dLL[i], datos$dUU[i],datos$dLU[i], datos$dUL[i])}

#Ajuste de los p valores con el False Discovery Rate
fdrf1<-p.adjust(pvaluebinomf1,method="fdr")
fdrf1<-as.data.frame(fdrf1)
rownames(fdrf1)<- 1:marker
#Crea un vector con los loci outiler (alpha=0.05) tras el ajuste
sigf1<-as.numeric(rownames(fdrf1)[fdrf1<0.05])
result<-list(Pvalues=fdrf1,Outliers=sigf1)
return(result)
}

if (method=="gagnaire" && type=="F1"){
pa<-apply(mat$PA,2,mean)
pb<-apply(mat$PB,2,mean)
pf1<-apply(mat$F1,2,mean)
Na<-dim(mat$PA)[1]
Nb<-dim(mat$PB)[1]
Nf1<-dim(mat$F1)[1]
Nmarker<-dim(mat$PA)[2]

#Modelo de frecuencias esperadas en F1
model<-function(x,y){
  a<-1-sqrt(1-x)
  b<-1-sqrt(1-y)
  a+b-a*b}

#Frecuencias esperadas
z<-model(pa,pb)

marker <-Nmarker #numero de markers#
pvaluebinomf1<-numeric(length(marker))
for(i in 1:marker){
  pvaluebinomf1[i]<-(binom.test(pf1[i]*Nf1,n=Nf1, p=z[i]))$p.value
  }

fdrf1<-p.adjust(pvaluebinomf1,method="fdr")
fdrf1<-as.data.frame(fdrf1)
rownames(fdrf1)<- 1:marker

#Crea un vector con los loci outiler (alpha=0.05) tras el ajuste
sigf1<-as.numeric(rownames(fdrf1)[fdrf1<0.05])
result<-list(Pvalues=fdrf1,Outliers=sigf1)
return(result)
}


if (method=="bal&gar-ca" && type =="BxA"){
pa<-apply(mat$PA,2,mean)
pb<-apply(mat$PB,2,mean)
pf1<-apply(mat$BxA,2,mean)
Na<-dim(mat$PA)[1]
Nb<-dim(mat$PB)[1]
Nf1<-dim(mat$BxA)[1]
Nmarker<-dim(mat$PA)[2]


#Calcula los intervalos de confianza de la binomial
p.L <- function(x, n, alpha) {  
  if (x == 0)   0
  else qbeta(alpha, x, n - x + 1)
  }
p.U <- function(x,n, alpha) {  
  if (x == n)    1
  else qbeta(1 - alpha, x + 1, n - x)
  }

#Calcula las coordenadas x e y de las esquinas de la sabana

xU<-sapply(pa*Na,p.U,n=Na, alpha=0.02532)
xL<-sapply(pa*Na,p.L,n=Na, alpha=0.02532)
yL<-sapply(pb*Nb,p.L,n=Nb, alpha=0.02532)
yU<-sapply(pb*Nb,p.U,n=Nb, alpha=0.02532)

#Modelo de frecuencias esperadas en Bxa
model<-function(x,y){
  a<-1-sqrt(1-x)
  b<-1-sqrt(1-y)
  ((3*a)+b-(a^2)-(a*b))/2}

#Calcula el z para cada esquina
zUU<-model(xU, yU)
zUL<-model(xU, yL)
zLL<-model(xL, yL)
zLU<-model(xL, yU)

# Calcula la diferencia entre el valor de z estimado y el teorico 

DUU<-pf1-zUU
DLL<-pf1-zLL
DLU<-pf1-zLU
DUL<-pf1-zUL

#Crea una dataframe con todos los valores anteriores para su uso en la funcion

datos<-cbind(pf1*Nf1, pf1,zLL,zUU,zLU,zUL,DLL,DUU,DLU,DUL )
datos<-as.data.frame(datos)

colnames(datos)<- c("x","z","zLL","zUU","zLU","zUL","dLL","dUU","dLU","dUL")

#Funcion  que calcula la p de outlier en un fragmento basada en un test binomial eligiendo la esquina de la sabana mas cercana al valor estimado
zbinomial<-function(x,z,zLL,zUU,zLU,zUL,dLL,dUU,dLU,dUL) {
  if (sign(dLL)== sign(dUU) && sign(dUU) == sign(dLU) && sign(dLU) == sign(dUL) && z > zLL)
  {maxz<-max(zLL,zUU,zLU,zUL); bi<-binom.test(x,n=Nf1, p=maxz);print(bi$p.value)}
  else if (sign(dLL)== sign(dUU) && sign(dUU) == sign(dLU) && sign(dLU) == sign(dUL) && z < zLL)
  {minz<-min(zLL,zUU,zLU,zUL); bi<-binom.test(x,n=Nf1, p=minz);print(bi$p.value)}
  else {1}}


#Loop para correr la funcion zbinomial para cada uno de los fragmentos
marker <-Nmarker #numero de markers#
pvaluebinomf1<-numeric(length(marker))
for(i in 1:marker){
  pvaluebinomf1[i]<-zbinomial(datos$x[i],datos$z[i],datos$zLL[i],datos$zUU[i],datos$zLU[i],datos$zUL[i],datos$dLL[i], datos$dUU[i],datos$dLU[i], datos$dUL[i])}

#Ajuste de los p valores con el False Discovery Rate
fdrf1<-p.adjust(pvaluebinomf1,method="fdr")
fdrf1<-as.data.frame(fdrf1)
rownames(fdrf1)<- 1:marker

#Crea un vector con los loci outiler (alpha=0.05) tras el ajuste
sigbx<-as.numeric(rownames(fdrf1)[fdrf1<0.05])
result<-list(Pvalues=fdrf1,Outliers=sigbx)
return(result)
}
if (method=="gagnaire" && type=="BxA"){
pa<-apply(mat$PA,2,mean)
pb<-apply(mat$PB,2,mean)
pf1<-apply(mat$BxA,2,mean)
Na<-dim(mat$PA)[1]
Nb<-dim(mat$PB)[1]
Nf1<-dim(mat$BxA)[1]
Nmarker<-dim(mat$PA)[2]

#Modelo de frecuencias esperadas en BxA
model<-function(x,y){
  a<-1-sqrt(1-x)
  b<-1-sqrt(1-y)
  ((3*a)+b-(a^2)-(a*b))/2}

#Frecuencias esperadas
z<-model(pa,pb)

marker <-Nmarker #numero de markers#
pvaluebinomf1<-numeric(length(marker))
for(i in 1:marker){
  pvaluebinomf1[i]<-(binom.test(pf1[i]*Nf1,n=Nf1, p=z[i]))$p.value
  }

fdrf1<-p.adjust(pvaluebinomf1,method="fdr")
fdrf1<-as.data.frame(fdrf1)
rownames(fdrf1)<- 1:marker
#Crea un vector con los loci outiler (alpha=0.05) tras el ajuste
sigbx<-as.numeric(rownames(fdrf1)[fdrf1<0.05])
result<-list(Pvalues=fdrf1,Outliers=sigbx)
return(result)
}

if (method=="bal&gar-ca" && type =="BxB"){
pa<-apply(mat$PA,2,mean)
pb<-apply(mat$PB,2,mean)
pf1<-apply(mat$BxB,2,mean)
Na<-dim(mat$PA)[1]
Nb<-dim(mat$PB)[1]
Nf1<-dim(mat$BxB)[1]
Nmarker<-dim(mat$PA)[2]


#Calcula los intervalos de confianza de la binomial
p.L <- function(x, n, alpha) {  
  if (x == 0)   0
  else qbeta(alpha, x, n - x + 1)
  }
p.U <- function(x,n, alpha) {  
  if (x == n)    1
  else qbeta(1 - alpha, x + 1, n - x)
  }

#Calcula las coordenadas x e y de las esquinas de la sabana

xU<-sapply(pa*Na,p.U,n=Na, alpha=0.02532)
xL<-sapply(pa*Na,p.L,n=Na, alpha=0.02532)
yL<-sapply(pb*Nb,p.L,n=Nb, alpha=0.02532)
yU<-sapply(pb*Nb,p.U,n=Nb, alpha=0.02532)

#Modelo de frecuencias esperadas en Bxa
model<-function(x,y){
  a<-1-sqrt(1-x)
  b<-1-sqrt(1-y)
  ((3*b)+a-(b^2)-(a*b))/2}

#Calcula el z para cada esquina
zUU<-model(xU, yU)
zUL<-model(xU, yL)
zLL<-model(xL, yL)
zLU<-model(xL, yU)

# Calcula la diferencia entre el valor de z estimado y el teorico 

DUU<-pf1-zUU
DLL<-pf1-zLL
DLU<-pf1-zLU
DUL<-pf1-zUL

#Crea una dataframe con todos los valores anteriores para su uso en la funcion

datos<-cbind(pf1*Nf1, pf1,zLL,zUU,zLU,zUL,DLL,DUU,DLU,DUL )
datos<-as.data.frame(datos)

colnames(datos)<- c("x","z","zLL","zUU","zLU","zUL","dLL","dUU","dLU","dUL")

#Funcion  que calcula la p de outlier en un fragmento basada en un test binomial eligiendo la esquina de la sabana mas cercana al valor estimado
zbinomial<-function(x,z,zLL,zUU,zLU,zUL,dLL,dUU,dLU,dUL) {
  if (sign(dLL)== sign(dUU) && sign(dUU) == sign(dLU) && sign(dLU) == sign(dUL) && z > zLL)
  {maxz<-max(zLL,zUU,zLU,zUL); bi<-binom.test(x,n=Nf1, p=maxz);print(bi$p.value)}
  else if (sign(dLL)== sign(dUU) && sign(dUU) == sign(dLU) && sign(dLU) == sign(dUL) && z < zLL)
  {minz<-min(zLL,zUU,zLU,zUL); bi<-binom.test(x,n=Nf1, p=minz);print(bi$p.value)}
  else {1}}


#Loop para correr la funcion zbinomial para cada uno de los fragmentos
marker <-Nmarker #numero de markers#
pvaluebinomf1<-numeric(length(marker))
for(i in 1:marker){
  pvaluebinomf1[i]<-zbinomial(datos$x[i],datos$z[i],datos$zLL[i],datos$zUU[i],datos$zLU[i],datos$zUL[i],datos$dLL[i], datos$dUU[i],datos$dLU[i], datos$dUL[i])}

#Ajuste de los p valores con el False Discovery Rate
fdrf1<-p.adjust(pvaluebinomf1,method="fdr")
fdrf1<-as.data.frame(fdrf1)
rownames(fdrf1)<- 1:marker
#Crea un vector con los loci outiler (alpha=0.05) tras el ajuste
sigbx<-as.numeric(rownames(fdrf1)[fdrf1<0.05])
result<-list(Pvalues=fdrf1,Outliers=sigbx)
return(result)
}
if (method=="gagnaire" && type=="BxB"){
pa<-apply(mat$PA,2,mean)
pb<-apply(mat$PB,2,mean)
pf1<-apply(mat$BxA,2,mean)
Na<-dim(mat$PA)[1]
Nb<-dim(mat$PB)[1]
Nf1<-dim(mat$BxA)[1]
Nmarker<-dim(mat$PA)[2]

#Modelo de frecuencias esperadas en BxA
model<-function(x,y){
  a<-1-sqrt(1-x)
  b<-1-sqrt(1-y)
  ((3*b)+a-(b^2)-(a*b))/2}

#Frecuencias esperadas
z<-model(pa,pb)

marker <-Nmarker #numero de markers#
pvaluebinomf1<-numeric(length(marker))
for(i in 1:marker){
  pvaluebinomf1[i]<-(binom.test(pf1[i]*Nf1,n=Nf1, p=z[i]))$p.value
  }

fdrf1<-p.adjust(pvaluebinomf1,method="fdr")
fdrf1<-as.data.frame(fdrf1)
rownames(fdrf1)<- 1:marker

#Crea un vector con los loci outiler (alpha=0.05) tras el ajuste
sigf1<-as.numeric(rownames(fdrf1)[fdrf1<0.05])
result<-list(Pvalues=fdrf1,Outliers=sigf1)
return(result)}
}
