#16:

#1:
faces<-function(xy=rbind(1:3,5:3,3:5,5:7),which.row,fill=FALSE,nrow,ncol,
                scale=TRUE,byrow=FALSE,main,labels){
  
#21:
spline<-function(a,y,m=200,plot=FALSE){
    n<-length(a)
  h<-diff(a)
  dy<-diff(y)
  sigma<-dy/h
  lambda<-h[-1]/(hh<-h[-1]+h[-length(h)])
  mu<-1-lambda
  d<-6*diff(sigma)/hh
  tri.mat<-2*diag(n-2)
  tri.mat[2+  (0:(n-4))*(n-1)] <-mu[-1]
  tri.mat[    (1:(n-3))*(n-1)] <-lambda[-(n-2)]
  M<-c(0,solve(tri.mat)%*%d,0)
  x<-seq(from=a[1],to=a[n],length=m)
  anz.kl <- hist(x,breaks=a,plot=FALSE)$counts
  adj<-function(i) i-1
  i<-rep(1:(n-1),anz.kl)+1
  S.x<-  M[i-1]*(a[i]-x          )^3 / (6*h[adj(i)])  +
         M[i]  *(x        -a[i-1])^3 / (6*h[adj(i)])  +
         (y[i-1] - M[i-1]*h[adj(i)]^2 /6) * (a[i]-x)/ h[adj(i)] +
        (y[i]   - M[i]  *h[adj(i)]^2 /6) * (x-a[i-1]) / h[adj(i)]
  if(plot){ plot(x,S.x,type="l"); points(a,y)    }
  return(cbind(x,S.x))
}

#:21

  
#4:
n.char<-15
xy<-rbind(xy)
if(byrow) xy<-t(xy)
if(!missing(which.row)&& all(  !is.na(match(which.row,1:dim(xy)[2]))  ))
       xy<-xy[,which.row,drop=FALSE]
mm<-dim(xy)[2];  n<-dim(xy)[1]
xnames<-dimnames(xy)[[1]]
if(is.null(xnames)) xnames<-as.character(1:n)
if(!missing(labels)) xnames<-labels
if(scale){
   xy<-apply(xy,2,function(x){
           x<-x-min(x); x<-if(max(x)>0) 2*x/max(x)-1 else x })
} else xy[]<-pmin(pmax(-1,xy),1)
xy<-rbind(xy);n.c<-dim(xy)[2]
xy<-xy[,(h<-rep(1:mm,ceiling(n.char/mm))),drop=FALSE]
if(fill) xy[,-(1:n.c)]<-0

#:4

  
#5:
face.orig<-list(
      eye  =rbind(c(12,0),c(19,8),c(30,8),c(37,0),c(30,-8),c(19,-8),c(12,0))
     ,iris =rbind(c(20,0),c(24,4),c(29,0),c(24,-5),c(20,0))
     ,lipso=rbind(c(0,-47),c( 7,-49),lipsiend=c( 16,-53),c( 7,-60),c(0,-62))
     ,lipsi=rbind(c(7,-54),c(0,-54))                  # add lipsiend
     ,nose =rbind(c(0,-6),c(3,-16),c(6,-30),c(0,-31))
     ,shape =rbind(c(0,44),c(29,40),c(51,22),hairend=c(54,11),earsta=c(52,-4),
                  earend=c(46,-36),c(38,-61),c(25,-83),c(0,-89))
     ,ear  =rbind(c(60,-11),c(57,-30))                # add earsta,earend
     ,hair =rbind(hair1=c(72,12),hair2=c(64,50),c(36,74),c(0,79)) # add hairend
)
lipso.refl.ind<-4:1
lipsi.refl.ind<-1
nose.refl.ind<-3:1
hair.refl.ind<-3:1
shape.refl.ind<-8:1
shape.xnotnull<-2:8
nose.xnotnull<-2:3

#:5

  
#2:
nr<-n^0.5; nc<-n^0.5
if(!missing(nrow)) nr<-nrow
if(!missing(ncol)) nc<-ncol
opar<-par(mfrow=c(ceiling(c(nr,nc))),oma=rep(6,4), mar=rep(.7,4))
on.exit(par(opar))

#:2

  
#6:
for(ind in 1:n){
  
#7:
factors<-xy[ind,]
face <- face.orig

#:7

  
#9:
m<-mean(face$lipso[,2])
face$lipso[,2]<-m+(face$lipso[,2]-m)*(1+0.7*factors[4])
face$lipsi[,2]<-m+(face$lipsi[,2]-m)*(1+0.7*factors[4])
face$lipso[,1]<-face$lipso[,1]*(1+0.7*factors[5])
face$lipsi[,1]<-face$lipsi[,1]*(1+0.7*factors[5])
face$lipso["lipsiend",2]<-face$lipso["lipsiend",2]+20*factors[6]

#:9

  
#10:
m<-mean(face$eye[,2])
face$eye[,2] <-m+(face$eye[,2] -m)*(1+0.7*factors[7])
face$iris[,2]<-m+(face$iris[,2]-m)*(1+0.7*factors[7])
m<-mean(face$eye[,1])
face$eye[,1] <-m+(face$eye[,1] -m)*(1+0.7*factors[8])
face$iris[,1]<-m+(face$iris[,1]-m)*(1+0.7*factors[8])


#:10

  
#11:
m<-min(face$hair[,2])
face$hair[,2]<-m+(face$hair[,2]-m)*(1+0.2*factors[9])
m<-0
face$hair[,1]<-m+(face$hair[,1]-m)*(1+0.2*factors[10])
m<-0
face$hair[c("hair1","hair2"),2]<-face$hair[c("hair1","hair2"),2]+50*factors[11]

#:11

  
#12:
m<-mean(face$nose[,2])
face$nose[,2]<-m+(face$nose[,2]-m)*(1+0.7*factors[12])
face$nose[nose.xnotnull,1]<-face$nose[nose.xnotnull,1]*(1+factors[13])

#:12

  
#13:
m<-mean(face$shape[c("earsta","earend"),1])
face$ear[,1]<-m+(face$ear[,1]-m)* (1+0.7*factors[14])
m<-min(face$ear[,2])
face$ear[,2]<-m+(face$ear[,2]-m)* (1+0.7*factors[15])

#:13

  
#8:
face<-lapply(face,function(x){ x[,2]<-x[,2]*(1+0.2*factors[1]);x})
face<-lapply(face,function(x){ x[,1]<-x[,1]*(1+0.2*factors[2]);x})
face<-lapply(face,function(x){ x[,1]<-ifelse(x[,1]>0,
                                        ifelse(x[,2] > -30, x[,1],
                  pmax(0,x[,1]+(x[,2]+50)*0.2*sin(1.5*(-factors[3])))),0);x})
#face$shape[,2]<-face$shape[,2]*(1+0.2*factors[1])
#face$shape[,1]<-face$shape[,1]*(1+0.2*factors[2])
#face$shape[,1]<-face$shape[,1]<-ifelse(face$shape[,1]>0,
#   ifelse(face$shape[,2] > -30, face$shape[,1],
#      pmax(0,face$shape[,1]+(face$shape[,2]+50)*0.2*sin(1.5*(-factors[3])))),0)


#:8

  
#14:
invert<-function(x) cbind(-x[,1],x[,2])
face.obj<-list(
     eyer=face$eye
    ,eyel=invert(face$eye)
    ,irisr=face$iris
    ,irisl=invert(face$iris)
    ,lipso=rbind(face$lipso,invert(face$lipso[lipso.refl.ind,]))
    ,lipsi=rbind(face$lipso["lipsiend",],face$lipsi,
                 invert(face$lipsi[lipsi.refl.ind,,drop=FALSE]),
                 invert(face$lipso["lipsiend",,drop=FALSE]))
    ,earr=rbind(face$shape["earsta",],face$ear,face$shape["earend",])
    ,earl=invert(rbind(face$shape["earsta",],face$ear,face$shape["earend",]))
    ,nose=rbind(face$nose,invert(face$nose[nose.refl.ind,]))
    ,hair=rbind(face$shape["hairend",],face$hair,invert(face$hair[hair.refl.ind,]),
                invert(face$shape["hairend",,drop=FALSE]))
    ,shape=rbind(face$shape,invert(face$shape[shape.refl.ind,]))
)

#:14

  
#15:
plot(1,type="n",xlim=c(-105,105)*1.1, axes=FALSE,
     ylab="",ylim=c(-105,105)*1.3)
title(xnames[ind])
for(ind in seq(face.obj)) {
       x <-face.obj[[ind]][,1]; y<-face.obj[[ind]][,2]
       xx<-spline(1:length(x),x,40,FALSE)[,2]
       yy<-spline(1:length(y),y,40,FALSE)[,2]
       lines(xx,yy)
}

#:15

}

#:6

  
#3:
if(!missing(main)){
  par(opar);par(mfrow=c(1,1))
  mtext(main, 3, 3, TRUE, 0.5)
  title(main)
}

#:3

}

#:1


#:16

