#' biplot for the dataset
#' 
#' Loading plot of the variables from a Principal Components Analysis.
#' scores of the observations are surimposed
#'
#' @param X the data matrix
#' @param sX TRUE/FALSE : standardization or not of the columns X (TRUE by default)
#' @param axeh component number for the horizontal axis 
#' @param axev component number for the vertical axis 
#' @param cex.lab : magnification to be used for labels (1 by default)
#' 
#' @export
#' 

data_biplot <-
function(X,sX=TRUE,axeh=1,axev=2,cex.lab=1) {
 

  
  X<- scale(X, center=T, scale=sX)
  p <- dim(X)[2] 
  n <- dim(X)[1]
   
   
  # PCA of X
  A<-max(axeh,axev)
  if (n<p) {
     reseig<-eigen(X%*%t(X)/n)
     valp<-100*reseig$values[1:A]/sum(reseig$values)
     coordvar<-t(X)%*%reseig$vectors[,1:A]/sqrt(n)
     coordind<-reseig$vectors[,1:A]%*%diag(sqrt(reseig$values[1:A]))
  } else {
     reseig<-eigen(t(X)%*%X/(n))
     valp<-100*reseig$values[1:A]/sum(reseig$values)
     coordvar<-reseig$vectors[,1:A]%*%diag(sqrt(reseig$values[1:A]))
     coordind<-X%*%reseig$vectors[,1:A]/sqrt(n)
  }
 #re-orientation of the PC so that the maximal nb of var have positive coordinate along this axe
  for (a in 1:A) {
      if (sign(mean(coordvar[,a]))==(-1)) {
         coordvar[,a]=coordvar[,a]*(-1)
         coordind[,a]=coordind[,a]*(-1)
      }
  }
    
par(pty="s")

### PCA biplot
# scaling factor
vp<-(coordind[,axeh]^2+coordind[,axev]^2)
lp=max(vp)
vv<-(coordvar[,axeh]^2+coordvar[,axev]^2)
lv=max(vv)
f=sqrt(lp/lv)
#plot
plot(c(coordvar[,axeh]*f,coordind[,axeh]),c(coordvar[,axev]*f,coordind[,axev]),type="n",
     xlab=paste("Dim ",axeh,"(",round(valp[axeh],2),"%)"),
     ylab=paste("Dim ",axev,"(",round(valp[axev],2),"%)"),
     main="PCA biplot")
arrows(0,0,coordvar[,axeh]*f,coordvar[,axev]*f,length=0.1,angle=10,lwd=0.5,col="gray")
posi=rep(1,n)
posi[which(coordind[,axeh]>max(c(coordvar[,axeh]*f,coordind[,axeh]))*0.8)]=2
posi[which(coordind[,axeh]<min(c(coordvar[,axeh]*f,coordind[,axeh]))*0.8)]=4
posi[which(coordind[,axev]>max(c(coordvar[,axev]*f,coordind[,axev]))*0.8)]=1
posi[which(coordind[,axev]<min(c(coordvar[,axev]*f,coordind[,axev]))*0.8)]=3
text(coordind[,axeh],coordind[,axev],labels=rownames(X),pos=posi,cex=cex.lab) 
#autoLab(coordind[,axeh],coordind[,axev],labels=rownames(X),cex=1.5) 
abline(h=0,v=0)

par(pty="m")
}
