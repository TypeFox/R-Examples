Niche<-function(data, variables, Level=NULL, Taxon=NULL, cor=TRUE,
 d.main=0.5, xlab="Polar coordinate X in pixels",ylab="Polar coordinate Y in pixels",
cex.labS=1.5, font.lab=1, main="", colramp = IDPcolorRamp, cex.main = 2, font.main=2, nlab.xaxis = 5,
nlab.yaxis = 5, minL.axis = 3, las = 1, border = FALSE, tcl = -0.3, boxplot=TRUE, outline=FALSE,
color=NULL, range = 1.5, width = NULL, varwidth = FALSE, plot = TRUE, pars = list(boxwex = 0.8,
staplewex = 0.5, outwex = 0.5), cex.boxplot=1.5, cex.labB=1.5, namesB, family="serif", line=1,
file1 = "List of species.csv", file2 = "Environmental variables.csv", file3 = "Polar coordinates.csv",
na = "NA", dec = ",", row.names = FALSE, fileEncoding = ""){

if (requireNamespace("IDPmisc", quietly = TRUE)) {
IDPmisc::iplot
IDPmisc::IDPcolorRamp
}
else{
## do something else not involving IDPmisc.
}



if (missing(cor)) cor=TRUE else cor=cor


datos<-na.exclude(data)

selection<-subset(datos, select=variables)

dimdt<-dim(selection)

if(cor==TRUE){
if(dimdt[2]>3){
datosT2<-selection
datosT3<-selection[,1]

names<-names(selection)[1]

for(zz in 1:(dimdt[2]-2)){
cor<-cor(datosT2)
corT1<-cor[,-1]
datosT2<-datosT2[,-1]

df<-which.max(cor[!cor[,1]==1,1])

if(zz==(dimdt[2]-2)){
dimdim<-dim(datosT3)
cor1<-cor(datosT3[,dimdim[2]], datosT2[,1])
cor2<-cor(datosT3[,dimdim[2]], datosT2[,2])
if(cor1>cor2){
datosT3<-data.frame(datosT3, datosT2[,1],datosT2[,2])
names(datosT3)<-c(names,names(datosT2))
}
else{
datosT3<-data.frame(datosT3, datosT2[,2],datosT2[,1])
kk<-datosT2
kk<-data.frame(datosT2[,2],datosT2[,1])
names(kk)<-c(names(datosT2)[2],names(datosT2)[1])
names(datosT3)<-c(names,names(kk))
} 

}
else{
datosT3<-data.frame(datosT3, datosT2[,df])
names(datosT3)<-c(names,names(datosT2)[df])
}

names<-names(datosT3)


datosP<-datosT2[,-df]


datosZZZ<-datosT2
datosT2<-data.frame(datosT2[,names(datosT2)[df]], datosP)
names(datosT2)<-c(names(datosZZZ)[df], colnames(datosP))
}

selection<-datosT3
}
}#end cor

datos<-data.frame(datos[,1:7],selection)


if (missing(d.main)) d.main=0.5 else d.main=d.main
if (missing(xlab)) xlab="Polar coordinate X in pixels" else xlab=xlab
if (missing(ylab)) ylab="Polar coordinate Y in pixels" else ylab=ylab
if (missing(cex.labS)) cex.labS=1.5 else cex.labS=cex.labS
if (missing(cex.labB)) cex.labB=1.5 else cex.labB=cex.labB
if (missing(font.lab)) font.lab=1 else font.lab=font.lab
if (missing(main)) main="" else main=main
if (missing(cex.main)) cex.main=2 else cex.main=cex.main
if (missing(font.main)) font.main=2 else font.main=font.main
if (missing(nlab.xaxis)) nlab.xaxis=5 else nlab.xaxis=nlab.xaxis
if (missing(nlab.yaxis)) nlab.yaxis=5 else nlab.yaxis=nlab.yaxis
if (missing(minL.axis)) minL.axis=3 else minL.axis=minL.axis
if (missing(las)) las=1 else las=las
if (missing(border)) border=FALSE else border=border
if (missing(tcl)) tcl=-0.3 else tcl=tcl
if (missing(boxplot)) boxplot=TRUE else boxplot=boxplot
if (missing(outline)) outline=FALSE else outline=outline
if (missing(range)) range=1.5 else range=range
if (missing(width)) width=NULL else width=width
if (missing(varwidth)) varwidth=FALSE else varwidth=varwidth
if (missing(plot)) plot=TRUE else plot=plot
if (missing(pars)) pars=list(boxwex = 0.8, staplewex = 0.5, outwex = 0.5) else pars=pars
if (missing(cex.boxplot)) cex.boxplot=1.5 else cex.boxplot=cex.boxplot
if (missing(family)) family="serif" else family=family
if (missing(line)) line=1 else line=line
if (missing(file1)) file1= "List of species.csv" else file1 = file1
if (missing(file2)) file2= "Environmental variables.csv" else file2 = file2
if (missing(file3)) file3= "Polar coordinates.csv" else file3 = file3
if (missing(na)) na="NA" else na=na
if (missing(dec)) dec="," else dec=dec
if (missing(row.names)) row.names=FALSE else row.names=row.names
if (missing(fileEncoding)) fileEncoding="" else fileEncoding=fileEncoding
if (missing(Level)) Level=NULL else Level=Level
if (missing(Taxon)) Taxon=NULL else Taxon=Taxon

lentih<-length(variables)
if(lentih>12){
boxplot=FALSE
}


ZZ<-matrix(c("","","","","","","",""), nrow=4)

begin.time<-Sys.time() 
begin.times <- format(begin.time, "%b %d, %Y at %X") 
ZZ[2,1]<-"Estimating the niche...."
ZZ[3,1]<-begin.times
write.table(ZZ,"Inf.txt", row.names=FALSE,col.names=FALSE)


#Standardization 1 to -1

a<-dim(datos)

datosE<-datos[,1:7]

for (z in 8:a[2]){
matrixE<-matrix(c(-1, 1, min(datos[,z],na.rm=TRUE),max(datos[,z],na.rm=TRUE)), nrow = 2 , ncol = 2)
reg<-lm(matrixE[,1]~matrixE[,2])
datosC<-reg$coefficients[1]+datos[,z]*reg$coefficients[2]
datosE<-cbind(datosE,datosC)
}

colnames(datosE)<-colnames(datos)



#Estimation of polar coordinates

angle<-360/((a[2]-7)*2)

datosX<-datosE[,1:7]
h<-0
for (z in 8:a[2]){
h<-h+1
datosC<- ifelse(datosE[,z] <=0, abs(datosE[,z])*cos(angle*h+180), abs(datosE[,z])*cos(angle*h)) 
datosX<-cbind(datosX,datosC)
}
par(mfrow = c(1,1))
pixels<-30
pxx<-dev.size("px")
pixs<-(-1.596756639+0.007875013*pxx[1])



XX<-apply(datosX[,8:a[2]],1,sum)
RX<-(max(XX)-min(XX))
XXX<-round(((XX-min(XX))*pixels/RX),digits=0)


datosY<-datosE[,1:7]
h<-0
for (z in 8:a[2]){
h<-h+1
datosC<- ifelse(datosE[,z] <=0, abs(datosE[,z])*sin(angle*h+180), abs(datosE[,z])*sin(angle*h)) 
datosY<-cbind(datosY,datosC)
}


YY<-apply(datosY[,8:a[2]],1,sum)
RY<-(max(YY)-min(YY))
YYY<-round(((YY-min(YY))*pixels/RY),digits=0)


datosF<-cbind(datosE[,1:7],XXX,YYY,XX,YY)
colnames(datosF)<-c(colnames(datosE[1:7]),"Pixel.X","Pixel.Y","X","Y")

datosF<-merge(datosF,datos)
b<-dim(datosF)

datosFF<-aggregate(datosF[,10:b[2]],by=list(datosF[,1],datosF[,2],
datosF[,3],datosF[,4],datosF[,5],datosF$"Pixel.X",datosF$"Pixel.Y"),mean)

colnames(datosFF)<-c(colnames(datosE[1:5]),"Pixel.X","Pixel.Y",colnames(datosF[,10:b[2]]))

rm(datos)
rm(datosY)
rm(datosC)
rm(datosE)
rm(datosF)

#Plot of polar coordinates

if(!is.null(Level)){
datosFF<-subset(datosFF,datosFF[,Level] %in% Taxon)
}
else{
datosFF<-datosFF
}


datosFFF<-rbind(datosFF[,6:7],c(0,0),c(pixels,pixels))

dimFF<-dim(datosFF)

for(kkk in 1:dimFF[2]){
if(colnames(datosFF)[kkk]=="Latitude.1"){
colnames(datosFF)[kkk]<-"Latitude"
}
if(colnames(datosFF)[kkk]=="Longitude.1"){
colnames(datosFF)[kkk]<-"Longitude"
}
}

write.table(15,"Pointsize.dat", row.names=FALSE,col.names=FALSE)



devact<-dev.cur()

if(devact==3){
dev.off()
}

iplot(x=datosFFF$Pixel.X,y=datosFFF$Pixel.Y, pixs=pixs, xlab=xlab, 
ylab=ylab, cex.lab=cex.labS, font.lab=font.lab,colramp = colramp, cex = 1,
legend = TRUE, d.legend = 1, nlab.xaxis = nlab.xaxis, nlab.yaxis = nlab.yaxis,
minL.axis = minL.axis, las = las, border = border,oma = c(5,4,1,0)+0.1, tcl= tcl, family=family)

polygon(x = c(-1,-1,4,4), y = c(-1,4,4,-1), col = "white", border = NA)
polygon(x = c(29,29,31,31), y = c(29,31,31,29), col = "white", border = NA)

if(main==""){
pp<-1
}
else{
mtext(main, 3, line=d.main, cex=cex.main, font=font.main)
}

maxx<-max(datosFF$Pixel.X, na.rm=TRUE)
minx<-min(datosFF$Pixel.X, na.rm=TRUE)
maxy<-max(datosFF$Pixel.Y, na.rm=TRUE)
miny<-min(datosFF$Pixel.Y, na.rm=TRUE)

matriz<-datosFF

datos3<-matriz[(matriz$Pixel.X>minx)&(matriz$Pixel.X<maxx),]
datos3<-datos3[(datos3$Pixel.Y>miny)&(datos3$Pixel.Y<maxy),]
datos3<-na.exclude(datos3)
names<-colnames(datos3)


b<-dim(datosFF)
m<-b[2]-9


if(boxplot==TRUE){



#Identify the coordinates
point<-""

begin.time<-Sys.time() 
begin.times <- format(begin.time, "%b %d, %Y at %X") 
ZZ[2,1]<-"Select the pixels just by clicking four times on the graph when a cross appears...."
ZZ[3,1]<-begin.times
write.table(ZZ,"Inf.txt", row.names=FALSE,col.names=FALSE)

point<-as.data.frame(locator(4))
hh<-length(point)
if (hh==0){
maxx<-max(datosFF$Pixel.X, na.rm=TRUE)
minx<-min(datosFF$Pixel.X, na.rm=TRUE)
maxy<-max(datosFF$Pixel.Y, na.rm=TRUE)
miny<-min(datosFF$Pixel.Y, na.rm=TRUE)
}
else{
maxx<-ceiling(max(point$x, na.rm=TRUE))
minx<-floor(min(point$x, na.rm=TRUE))
maxy<-ceiling(max(point$y, na.rm=TRUE))
miny<-floor(min(point$y, na.rm=TRUE))
}


par(mfrow = c(1,1))

remove(datosFFF)

matriz<-datosFF

datos3<-matriz[(matriz$Pixel.X>minx)&(matriz$Pixel.X<maxx),]
datos3<-datos3[(datos3$Pixel.Y>miny)&(datos3$Pixel.Y<maxy),]
datos3<-na.exclude(datos3)
names<-colnames(datos3)


b<-dim(datosFF)
m<-b[2]-9

#Boxplot with the range of the environmental variables of the polar coordinates selected


BB<-0

par(mfcol=c(1,m),oma=c(3,5,1,1))

for (z in 10:b[2]){
BB<-BB+1
h<-z/b[2]

if(!is.null(color)) col=color[z-7] else col=hsv(h = h, s = 1, v = 1, alpha=1)


boxplot(x=datos3[,z],outline= outline, col=col,
range = range, width = width, varwidth = varwidth, plot = plot,
pars = pars, horizontal = FALSE, cex.axis=cex.boxplot)
if (missing(namesB)) textB<-names[z] else textB<-namesB[BB]
mtext(text=textB, side=1, line=line,outer = FALSE,at = NA,adj = NA, padj = NA, cex = cex.labB, col = NA, font = font.lab, family=family, las=las)
}
}
else{
par(mfrow = c(1,1))
}


datos4<-subset(datos3[,1:5], !duplicated(datos3$Species))
y<-dim(datos3)
Env<-summary(datos3[,c(10:y[2])])



#Output files

if(dec=="."){
write.csv(x=datos4,file = file1, fileEncoding = fileEncoding,
row.names=row.names,na=na)
write.csv(x=Env,file = file2, fileEncoding = fileEncoding,
row.names=row.names,na=na)
write.csv(x=datosFF,file = file3, fileEncoding = fileEncoding,
row.names=row.names,na=na)
}
else{
write.csv2(x=datos4,file = file1, fileEncoding = fileEncoding,
row.names=row.names,na=na)
write.csv2(x=Env,file = file2, fileEncoding = fileEncoding,
row.names=row.names,na=na)
write.csv2(x=datosFF,file = file3, fileEncoding = fileEncoding,
row.names=row.names,na=na)
}
rm(datosFF)
rm(Env)
rm(datos3)
rm(datos4)
rm(datosX)
rm(matriz)
rm(XX)
rm(XXX)
rm(YY)
rm(YYY)
}
