`stat.freq` <-
function(histogram) {
xx <- histogram$mids
yy <- histogram$counts
media <- sum(yy * xx)/sum(yy)
variancia <- sum(yy * (xx - media)^2)/(sum(yy) - 1)
zz <- histogram$breaks
z <- length(zz)
#Localiza los puntos de un vector
names(yy)<-1:length(yy)
id<-as.numeric(names(yy[max(yy)==yy]))
x1 <- xx[1] - (zz[2] - zz[1])
x2 <- xx[z - 1] + (zz[z] - zz[z - 1])
zz<-c(x1,zz,x2)
yy<-c(0,yy,0)
z<- z+2
names(yy)<-1:length(yy)
# Calculo de la mediana
total<-sum(yy)
suma<-0
i<-0
while (suma<=total/2) {
i<-i+1
suma<-suma+yy[i]
}
tic   <- zz[i+1]-zz[i]
mediana <- zz[i] +(total/2 - suma+yy[i])*tic/yy[i]
mediana<-as.numeric(mediana)
# Calculo de las modas
id<-as.numeric(names(yy[max(yy)==yy]))
modas<-length(id)
clases <- rep(0,2*modas)

for (i in 1:modas) {
j<-id[i]
k<-2*i-1
clases[k] <- zz[j]
clases[k+1]<-zz[j+1]
}
dim(clases)<-c(2,modas)
clases<-t(clases)
colnames(clases)<-c("[-","-]")
mode<-rep(0,modas)
for (i in 1:modas) {
j<-id[i]
delta1<- yy[j]-yy[j-1]
delta2<- yy[j]-yy[j+1]
tic   <- zz[j+1]-zz[j]
mode[i]<-zz[j]+delta1*tic/(delta1+delta2)
}
Moda<-cbind(clases,mode)
lista<-list(variance=variancia,mean=media,median=mediana,mode=Moda)
return(lista)
}

