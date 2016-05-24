plot.flu <-
function(x,...){
  opar<-par(mfrow=c(1,2))
  # Graph 1
  semanas<-dim(x$i.data)[1]
  anios<-dim(x$i.data)[2]
  datos.graf<-x$moving.epidemics
  colnames(datos.graf)<-names(x$i.data)
  lab.graf<-(1:semanas)+x$mean.start[2]-x$mean.start[1]
  lab.graf[lab.graf>52]<-(lab.graf-52)[lab.graf>52]
  lab.graf[lab.graf<1]<-(lab.graf+52)[lab.graf<1]
  tipos<-rep(1,anios)
  anchos<-rep(2,anios)
  colores<-c(rgb(runif(anios-1),runif(anios-1),runif(anios-1)),"#FF0000")
  limite.superior<-c(1.05*max.fix.na(datos.graf))
  matplot(1:semanas,datos.graf,type="l",sub=paste("Weekly incidence rates of the seasons matching their relative position in the model."),lty=tipos,lwd=anchos,col=colores,xlim=c(1,semanas),xlab="Week",ylab="Rate",font.axis=1,font.lab=1,font.main=2,font.sub=1, ylim=c(0,limite.superior),xaxt="n")
  axis(1,at=1:semanas,labels=as.character(lab.graf))
  i.temporada<-x$mean.start[1]
  f.temporada<-x$mean.start[1]+x$mean.length-1
  abline(v=c(i.temporada-0.5,f.temporada+0.5),col=c("#00C000","#FFB401"),lty=c(2,2))
  ya<-limite.superior*0.975
  if ((i.temporada-1)<=(semanas-f.temporada)){
    xa<-semanas*0.99
   	legend(x="topright",xjust=1,legend=names(x$i.data),lty=tipos,lwd=anchos,col=colores,cex=0.75)
  }else{
    xa<-semanas*0.01
    legend(x="topleft",legend=names(x$i.data),lty=tipos,lwd=anchos,col=colores,cex=0.75)
  }
  # Graph 2
  lineas.basicas<-rbind(matrix(rep(x$pre.post.intervals[1,1:3],i.temporada-1),ncol=3,byrow=T),
  matrix(rep(NA,3*(f.temporada-i.temporada+1)),ncol=3),
  matrix(rep(x$pre.post.intervals[2,1:3],semanas-f.temporada),ncol=3,byrow=T))
  n.niveles<-dim(x$epi.intervals)[1]
  limites.niveles<-as.vector(x$epi.intervals[1:n.niveles,4])
  nombres.niveles<-as.character(x$epi.intervals[1:n.niveles,1])
  limites.niveles[limites.niveles<0]<-0
  limites.niveles.mas<-array(dim=c(semanas,n.niveles))
  for (i in i.temporada:f.temporada) limites.niveles.mas[i,]<-limites.niveles
  datos.graf<-as.data.frame(cbind(lineas.basicas[,3],x$typ.curve[,2],limites.niveles.mas))
  etiquetas<-paste(round(limites.niveles,2)," (",as.numeric(nombres.niveles)*100,"%) ",sep="")
  names(datos.graf)<-c("Threshold",etiquetas)
  tipos<-c(1,2,rep(2,times=n.niveles))
  anchos<-c(3,2,rep(2,times=n.niveles))
  colores<-c("#FF0000","#800080","#C0C0FF","#8080FF","#4040FF","#0000C0","#000080")
  limite.superior<-c(0,1.1*max.fix.na(datos.graf))
  matplot(1:semanas,datos.graf,type="l",sub=paste("MEM Threshold"),lty=tipos,lwd=anchos,
    col=colores,xlim=c(1,semanas),xlab="Week",ylab="Rate",font.axis=1,font.lab=1,font.main=2,
    font.sub=1, ylim=limite.superior,xaxt="n")
  axis(1,at=1:semanas,labels=as.character(lab.graf))
  xa<-c(i.temporada/2,1+rep(f.temporada,n.niveles))
  ya<-c(lineas.basicas[1,3],limites.niveles)
  texto<-c(paste("Threshold: ",round(lineas.basicas[1,3],2),sep=""),etiquetas)
  posiciones<-c(3,rep(4,n.niveles))
  tamanios<-c(0.75,rep(0.75,n.niveles))
  colores<-c("#FF0000",colores[3:(3+n.niveles-1)])
  text(xa,ya,texto,pos=posiciones,col=colores,cex=tamanios)
  par(opar)
}
