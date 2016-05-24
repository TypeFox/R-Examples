OverlapTaxa<-function(data, Level, digits =2,
file1 = "Overlap among taxa.csv", file2 = "Mean overlap among taxa.csv",
na = "NA", dec = ",",row.names = FALSE,fileEncoding = ""){

datos<-na.exclude(data)


if (missing(digits)) digits<-2 else digits<-digits
if (missing(file1)) file1= "Overlap among taxa.csv" else file1 = file1
if (missing(file2)) file2= "Mean overlap among taxa.csv" else file2 = file2
if (missing(na)) na="NA" else na=na
if (missing(dec)) dec="," else dec=dec
if (missing(row.names)) row.names=FALSE else row.names=row.names
if (missing(fileEncoding)) fileEncoding="" else fileEncoding=fileEncoding


datosD<-subset(datos[,1:5], !duplicated(datos$Species))
datosD<-datosD[order(datosD[,1],datosD[,2],datosD[,3],datosD[,4],datosD[,5]) , ]
a<-dim(datosD)

datosJ<-subset(datos[,1:5], !duplicated(datos[, Level]))
datosJ<-datosJ[order(datosJ[,1],datosJ[,2],datosJ[,3],datosJ[,4],datosJ[,5]) , ]
hj<-dim(datosJ)

datosF<-datos[1,1:7]
colnames(datosF)<-c(colnames(datos[1:4]),"Species1","Species2","Overlap")
gg<-1
ww<-0
ZZ<-matrix(c("gg","gg",1,2,"","","",""), nrow=4)
ZZ[1,1]<-Level

com2<-0
runmax<-0
run1<-1

for (u in 1:hj[1]){
pp<-1

begin.time<-Sys.time() 
begin.times <- format(begin.time, "%b %d, %Y at %X")

Fa<-subset(datos,datos[,Level] %in% datosJ[u,Level])
Fa<-subset(Fa[,1:5], !duplicated(Fa$Species))

jj<-0

a<-dim(Fa)

if(u>=2){
ZZ[1,2]<-as.character(datosJ[u,Level])
ZZ[2,1]<-u
ZZ[2,2]<- paste("of",hj[1])


if(a[1]<160){
if(a[1]>=2) com<-factorial(a[1])/(2*factorial(a[1]-2)) else com<-1
run1<-run1*com
}
else{
com<-0
run1<-0
}


run1[is.na(run1)]<-0
run1[is.null(run1)]<-0
com[is.na(com)]<-0
com[is.null(com)]<-0

if(run1>=3600){
ZZ[4,2]<-"remaining hours in this taxon...."
}
else{
if(run1<=60) ZZ[4,2]<-"remaining seconds in this taxon...." else ZZ[4,2]<-"remaining minutes in this taxon...."
}

if(com==0){
ZZ[4,2]<-"It is not possible to estimate remaining time...." 
minutes<-""
}
else{
if(run1>=3600){
minutes<-run1/3600
}
else{
if(run1<=60) minutes<-run1 else minutes<-run1/60 
}
ZZ[4,2]<-ZZ[4,2]
}

ZZ[3,1]<-end.times

if(com==0) minutes<-"" else minutes<-round(minutes, digits=1)

ZZ[4,1]<-minutes

print(paste(u, " of ", hj[1]))
print(datosJ[u,Level])
print(end.times)
print(paste(minutes,ZZ[4,2]))

write.table(ZZ,"Inf.txt", row.names=FALSE,col.names=FALSE)

}
else{
print(paste(u, " of ", hj[1]))
print(datosJ[u,Level])
print(begin.times)
ZZ[1,2]<-as.character(datosJ[u,Level])
ZZ[2,1]<-u
ZZ[2,2]<- paste("of",hj[1])
ZZ[3,1]<- begin.times
ZZ[3,2]<- ""
ZZ[4,1]<- ""
ZZ[4,2]<- ""
write.table(ZZ,"Inf.txt", row.names=FALSE,col.names=FALSE)
}

for (z in 1:a[1]){
jj<-jj+1
if(jj<=a[1]){
for (h in 1:a[1]){
A<-subset(datos,datos[,5] %in% Fa[z,5])
dimA<-dim(A)

sp<-match(Fa[z,"Species"], Fa[h,"Species"], nomatch=0)

if(sp==0){

if(Fa[z,Level]==Fa[h,Level]){



B<-subset(datos,datos[,5] %in% Fa[h,5])
dimB<-dim(B)


###################
#Overlap

AA<-cbind(A$Pixel.X,A$Pixel.Y)
AA<-unique(AA)
DA<-dim(AA)

BB<-cbind(B$Pixel.X,B$Pixel.Y)
BB<-unique(BB)
DB<-dim(BB)

CC<-rbind(AA,BB)
DUP<-duplicated(CC)
DC<-length(DUP[DUP==TRUE])

OA<-DC*100/DA[1]

OA<-round(OA,digits=digits)

datosH<-cbind(Fa[z,1:5],Fa[h,5],OA)
if(pp==1){
datosK<-cbind(Fa[z,1:5],Fa[h,5],OA)
pp<-2
ww<-1
}
else{
datosK<-rbind(datosK,datosH)
}

}
else{
tt<-1
}

}
else{
tt<-1
}


}

}
else{
tt<-1
}
}

if(ww==1){
if(gg==1){
datosF<-datosK
gg<-2
ww<-1
}
else{
datosF<-rbind(datosF,datosK)
}
}
else{
ww<-0
}

end.time<-Sys.time() 
end.times <- format(end.time, "%b %d, %Y at %X")
run.time<-difftime(end.time,begin.time,units="secs")
run<-as.numeric(run.time)

if(a[1]<160){
if(a[1]>=2) com<-factorial(a[1])/(2*factorial(a[1]-2)) else com<-1
run1<-run/com
}
else{
com<-0
run1<-0
}

if(com>com2) runmax<-run1 else runmax<-runmax
if(com>com2) run1<-run1 else run1<-runmax
if(com>com2) com2<-com else com2<-com2
}



colnames(datosF)<-c(colnames(datos[1:4]),"Species1","Species2","Overlap")
datosF<-datosF[order(datosF[,1],datosF[,2],datosF[,3],datosF[,4],datosF[,5],datosF[,6]) , ]
datosF<-unique(datosF)

datosM<-aggregate(datosF[,7],by=list(datosF[,Level]),mean)
datosM<-datosM[order(datosM[,1]) , ]

datosSD<-aggregate(datosF[,7],by=list(datosF[, Level]),sd)
datosSD<-datosSD[order(datosSD[,1]) , ]

datosL<-aggregate(datosD[,5],by=list(datosD[,Level]),length)
colnames(datosL)<-c(Level, "Number")
datosL<-datosL[which(datosL$Number>1),]
datosL<-datosL[order(datosL[,1]) , ]


datosFM<-cbind(datosM,datosSD[,2], datosL[,2])

colnames(datosFM)<-c(Level, "Mean", "SD", "Species")


#Output file

if(dec=="."){
write.csv(x=datosF,file = file1, fileEncoding = fileEncoding,row.names=row.names,na=na)
write.csv(x = datosFM, file = file2, fileEncoding = fileEncoding, row.names=row.names,na=na)
}
else{
write.csv2(x = datosF,file = file1, fileEncoding = fileEncoding, row.names=row.names,na=na)
write.csv2(x = datosFM, file = file2, fileEncoding = fileEncoding, row.names=row.names,na=na)
}
ZZ<-matrix(c("END",""), nrow=1)
write.table(ZZ,"Inf.txt", row.names=FALSE,col.names=FALSE)
rm(datos)
rm(datosD)
rm(datosF)
rm(datosFM)
rm(datosH)
rm(datosJ)
rm(datosK)
rm(datosL)
rm(datosM)
rm(datosSD)
}

