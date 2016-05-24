#'@title  Representation of the variables and their group membership
#' 
#' @description
#' Loading plot of the variables from a Principal Components Analysis. The group membership of the variables is superimposed.
#'
#' @param resclv results of CLV(), CLV_kmeans() or LCLV()
#' @param K the number of groups in the partition (already defined if CLV_kmeans is used)
#' @param axeh component number for the horizontal axis 
#' @param axev component number for the vertical axis 
#' @param label = TRUE :the column names in X are used as labels / = FALSE: no labels (by default)
#' @param cex.lab : magnification to be used for labels (1 by default)
#' @param v_colors default NULL. If missing colors are given, by default
#' @param v_symbol =TRUE : symbols are given isntead of colors for the identification of the groups/ =FALSE: no symbol (by default). 
#' @param beside =TRUE : a plot per cluster of variables, side-by-side/ =FALSE :an unique plot with all the variables with the identification of their group membership (by default). 
#' 
#' @examples data(apples_sh)
#' resclvX <- CLV(X = apples_sh$senso, method = 1, sX = TRUE)
#' plot_var(resclvX, K = 4, axeh = 1, axev = 2)
#' 
#' @export
#' 
plot_var <-
function(resclv,K=NULL,axeh=1,axev=2,label=FALSE,cex.lab=1,v_colors=NULL,v_symbol=FALSE,beside=FALSE) {
  if (!inherits(resclv, c("clv","lclv"))) 
    stop("non convenient objects")
  
  X<-resclv$param$X
  
  if(is.null(resclv$param$K)) { 
    if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
    clusters<-resclv[[K]]$clusters[2,]
} else {
    clusters<-resclv$clusters[2,] 
    K<-resclv$param$K
 }
  
  X<- scale(X, center=T, scale=resclv$param$sX)
  p <- dim(X)[2] 
  n <- dim(X)[1]
   
  if (is.null(v_colors)) {v_colors <- c("blue","red","green","black",
                                        "purple","orange","yellow","tomato","pink",
                                        "gold","cyan","turquoise","violet","green",
                                        "khaki","maroon","chocolate","burlywood")}
  if(v_symbol) {v_colors<-rep("black",20)}
  
  
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

# ### PCA biplot
# # scaling factor
# vp<-(coordind[,axeh]^2+coordind[,axev]^2)
# lp=max(vp)
# vv<-(coordvar[,axeh]^2+coordvar[,axev]^2)
# lv=max(vv)
# f=sqrt(lp/lv)
# #plot
# plot(c(coordvar[,axeh]*f,coordind[,axeh]),c(coordvar[,axev]*f,coordind[,axev]),type="n",
#      xlab=paste("Dim ",axeh,"(",round(valp[axeh],2),"%)"),
#      ylab=paste("Dim ",axev,"(",round(valp[axev],2),"%)"),
#      main="PCA biplot")
# arrows(0,0,coordvar[,axeh]*f,coordvar[,axev]*f,length=0.1,angle=10,lwd=0.5,col="gray")
# posi=rep(1,n)
# posi[which(coordind[,axeh]>max(c(coordvar[,axeh]*f,coordind[,axeh]))*0.8)]=2
# posi[which(coordind[,axeh]<min(c(coordvar[,axeh]*f,coordind[,axeh]))*0.8)]=4
# posi[which(coordind[,axev]>max(c(coordvar[,axev]*f,coordind[,axev]))*0.8)]=1
# posi[which(coordind[,axev]<min(c(coordvar[,axev]*f,coordind[,axev]))*0.8)]=3
# text(coordind[,axeh],coordind[,axev],labels=rownames(X),pos=posi) 
# #autoLab(coordind[,axeh],coordind[,axev],labels=rownames(X),cex=1.5) 
# abline(h=0,v=0)
# #####



if (resclv$param$sX) {
  xmin=-1
  xmax=1
  ymin=-1
  ymax=1
} else {
  xmin=min(min(coordvar[,axeh]),-max(coordvar[,axeh]))
  xmax=max(max(coordvar[,axeh]),-min(coordvar[,axeh]))
  ymin=min(min(coordvar[,axev]),-max(coordvar[,axev]))
  ymax=max(max(coordvar[,axev]),-min(coordvar[,axev]))
}

  
  clean_var<-NULL
  if(resclv$param$strategy=="kplusone")  clean_var<-which(clusters==0)
  if(resclv$param$strategy=="sparselv"){
    names0=c()
    sloading = resclv$sloading
    for(k in 1:K) names0 = c(names0,dimnames(sloading[[k]])[[1]])
    names1 = names0[which(unlist(sloading)==0)]
    names2= dimnames(resclv$clusters)[[2]]
    clean_var = sort(match(names1,names2))
    names(clean_var) = colnames(X)[clean_var]
  }
  
  

 
  colpart<-NULL
  symbpart<-NULL
  symbpart<-20
  for (j in 1:p) {
    if (j %in% clean_var) {
        colpart[j]<-"gray"
        if(v_symbol) symbpart[j]="."
    }else{
        colpart[j]<-v_colors[clusters[j]]
        if(v_symbol) symbpart[j]=clusters[j]
    }
  }
  
if (beside==FALSE) {
  plot(coordvar[,c(axeh,axev)],col=colpart,pch=symbpart,cex.axis=0.7,cex.lab=0.7,
       xlab=paste("Dim ",axeh," (",round(valp[axeh],2),"%)"), 
       ylab=paste("Dim ",axev," (",round(valp[axev],2),"%)"), xlim=c(xmin,xmax), ylim=c(ymin,ymax))
  arrows(rep(0,p), rep(0,p), coordvar[,axeh], coordvar[,axev], col= colpart, angle=0)
  if(label) {
    for (j in 1:p) {
      if(coordvar[j,axev]>0) text(coordvar[j,axeh], coordvar[j,axev],colnames(X)[j],pos=3,cex=cex.lab,col=colpart[j])
      if(coordvar[j,axev]<0) text(coordvar[j,axeh], coordvar[j,axev],colnames(X)[j],pos=1,cex=cex.lab,col=colpart[j])
    }
  }
  abline(h=0,v=0)
  if (resclv$param$sX) symbols(0, 0, circles = 1, inches = FALSE, add = TRUE) 
  ncol=ceiling(K/4)
  if(v_symbol) {
    legend("topleft",paste("G",1:K,sep=""),col=v_colors[1:K],title="Groups", lty="solid", pch=1:K, seg.len=0.6,cex=0.7,ncol=ncol)     
  } else {
    legend("topleft",paste("G",1:K,sep=""),col=v_colors[1:K],title="Groups", lty="solid",  seg.len=0.6,cex=0.7,ncol=ncol) 
  }
}


if (beside==TRUE) {
  if (length(clean_var)==0) KK=K
  if (length(clean_var)>0) KK=K+1
  NL=ceiling(KK/3)
  par(mfrow=c(NL,3))
  case=levels(as.factor(colpart))
  for(k in 1:K) {
    who=which(clusters==k)
    plot(coordvar[who,c(axeh,axev)],col=colpart[who],pch=symbpart[who],cex.axis=0.7,cex.lab=0.7,
       xlab=paste("Dim ",axeh," (",round(valp[axeh],2),"%)"), 
       ylab=paste("Dim ",axev," (",round(valp[axev],2),"%)"), xlim=c(xmin,xmax), ylim=c(ymin,ymax))
   arrows(rep(0,p), rep(0,p), coordvar[who,axeh], coordvar[who,axev], col= colpart[who], angle=0)
   if(label) {
    for (j in who) {
      if(coordvar[j,axev]>0) text(coordvar[j,axeh], coordvar[j,axev],colnames(X)[j],pos=3,cex=cex.lab,col=colpart[j])
      if(coordvar[j,axev]<0) text(coordvar[j,axeh], coordvar[j,axev],colnames(X)[j],pos=1,cex=cex.lab,col=colpart[j])
    }
   }
  abline(h=0,v=0)
  if (resclv$param$sX) symbols(0, 0, circles = 1, inches = FALSE, add = TRUE) 
  if(v_symbol) {
    legend("topleft",paste("G",k,sep=""),col=v_colors[k], lty="solid", pch=k, seg.len=0.6,cex=1)     
  } else {
    legend("topleft",paste("G",k,sep=""),col=v_colors[k], lty="solid",  seg.len=0.6,cex=1) 
  }
  }
  # for the "others"
  if (length(clean_var)>0) {
    who=clean_var
    plot(coordvar[who,c(axeh,axev)],col=colpart[who],pch=symbpart[who],cex.axis=0.7,cex.lab=0.7,
       xlab=paste("Dim ",axeh," (",round(valp[axeh],2),"%)"), 
       ylab=paste("Dim ",axev," (",round(valp[axev],2),"%)"), xlim=c(xmin,xmax), ylim=c(ymin,ymax))
    arrows(rep(0,p), rep(0,p), coordvar[who,axeh], coordvar[who,axev], col= colpart[who], angle=0)
    if(label) {
     for (j in who) {
       if(coordvar[j,axev]>0) text(coordvar[j,axeh], coordvar[j,axev],colnames(X)[j],pos=3,cex=cex.lab,col=colpart[j])
       if(coordvar[j,axev]<0) text(coordvar[j,axeh], coordvar[j,axev],colnames(X)[j],pos=1,cex=cex.lab,col=colpart[j])
     }
    }
    abline(h=0,v=0)
    if (resclv$param$sX) symbols(0, 0, circles = 1, inches = FALSE, add = TRUE) 
    if(v_symbol) {
      legend("topleft","others",col="gray", lty="solid", pch=k, seg.len=0.6,cex=1)     
    } else {
      legend("topleft","others",col="gray", lty="solid",  seg.len=0.6,cex=1) 
    }
  }  
}

par(mfrow=c(1,1))
par(pty="m")
}
