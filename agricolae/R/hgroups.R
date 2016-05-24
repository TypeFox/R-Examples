`hgroups` <-
function(hmerge) {
cc<-hmerge
n1<-nrow(cc)
clases= paste( -cc[1,1], -cc[1,2],sep="-")
clases <-as.matrix(clases)
for(i in 2:n1) {
if ((cc[i,1] < 0) & (cc[i,2] < 0) ) {
clases <-rbind(clases,paste( -cc[i,1], -cc[i,2],sep="-") )
}
if ((cc[i,1]) < 0 & (cc[i,2] > 0)) {
clave<- cc[i,2]
parte<-clases[clave,]
clases <-rbind(clases,paste( -cc[i,1], parte,sep="-") )
}
if ((cc[i,1]) > 0 & (cc[i,2] < 0)) {
clave<- cc[i,1]
parte<-clases[clave,]
clases <-rbind(clases,paste( parte,-cc[i,2],sep="-") )
}
if ((cc[i,1]) > 0 & (cc[i,2] > 0)) {
clave1<- cc[i,1]
clave2<- cc[i,2]
parte1<-clases[clave1,]
parte2<-clases[clave2,]
clases <-rbind(clases,paste( parte1,parte2,sep="-") )
}
}

for(i in 1:n1) {
x<-strsplit(clases[i,1],split="-")[[1]]
x
unique(x)
y<-as.numeric(unique(x))
y<-sort(y)
y<-as.character(y)
ny<- length(y)
s<-y[1]
for(j in 2:ny) s<-paste(s,y[j],sep="-")
clases[i,1]<- s
}
return(clases)
}

