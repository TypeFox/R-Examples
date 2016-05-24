#---------------------------------------------------------------------------------------------
# Agrega las clases a un plano factorial
# Campo Elías Pardo, Pedro del Campo
# Octubre 11-06, modificada octubre 30/08
# requiere ade4, utiliza s.class de ade4
#PARAMETROS 
#  Los mismos parametros de planfac pero entra objeto FactoClass y parámetros nclus para nombre de
# las clases y cex.clus para tamaño de clases, col.row es un vector para color de las clases
# cstar es de s.class 0 para no poner radios a los individuos de una clase
#  infaxes: lugar para imprimir información de ejes: "out","in","no" ("out")
#---------------------------------------------------------------------------------------------
plotFactoClass <- function(FC,x=1,y=2,xlim=NULL,ylim=NULL,rotx=FALSE,roty=FALSE,roweti=row.names(dudi$li),
                        coleti=row.names(dudi$co),titre=NULL,axislabel=TRUE,
                        col.row=1:FC$k,col.col="blue",cex=0.8,cex.row=0.8,cex.col=0.8,
                        all.point=TRUE,Trow=TRUE,Tcol=TRUE,cframe=1.2,ucal=0,
                        cex.global=1,infaxes="out",
                        nclus=paste("cl", 1:FC$k, sep=""),cex.clu=cex.row,cstar=1)
{

dudi <- FC$dudi

 col.r <- numeric(nrow(dudi$li))
 for(i in 1:FC$k){ col.r[ as.numeric(FC$cluster) == i ]  <-  col.row[i]   }
 
 names(col.r) <-  names(FC$cluster)
colpunto <- "black"

plot.dudi(dudi,ex=x,ey=y,xlim=xlim,ylim=ylim,main=titre,rotx=rotx,roty=roty,roweti=roweti,
     coleti=coleti,axislabel=axislabel,                        	col.row=colpunto,col.col=col.col,cex=cex,cex.row=cex.row,cex.col=cex.col,
                        all.point=all.point,Trow=Trow,Tcol=Tcol,cframe=cframe,ucal=ucal,
                        cex.global=cex.global,infaxes=infaxes)
# grafica de las clases
     if (rotx) rotx<- -1 else rotx<-1
     if (roty) roty<- -1 else roty<-1
     corli <- cbind(rotx*dudi$li[,x],roty*dudi$li[,y])
     
     s.class(corli,fac=FC$cluste, wt=dudi$lw,add.plot = TRUE,
              cstar = cstar, cellipse = 0,
              label=nclus,clabel=cex.clu,col=col.row)
}
#---------------------------------------------------------------------------------------------


