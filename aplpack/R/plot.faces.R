plot.faces<-function(x,x.pos,y.pos,face.type = 1, 
                     width=1,height=1,labels,
                     ncolors=20,
                     col.nose=rainbow(ncolors),                   # nose
                     col.eyes=rainbow(ncolors,start=0.6,end=0.85),# eyes
                     col.hair=terrain.colors(ncolors),            # hair
                     col.face=heat.colors(ncolors),               # face
                     col.lips=rainbow(ncolors,start=0.0,end=0.2), # lips
                     col.ears=rainbow(ncolors,start=0.0,end=0.2), # ears

                     ...){
  if(missing(x)) return("no face.list object in call")
  face.list<-x$faces; face.data<-x$xy
  if(class(face.list)!="faces") {
      if(!is.list(face.list) || !any(names(face.list[[1]])=="lipso") )
        return("input not of class faces")
  }
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

  n<-length(face.list)
  if(missing(x.pos)){
     co<-ro<-ceiling(n^0.5)
     plot((0:ro)+.5,(0:co)+.5,type="n",xlab="",ylab="",axes=FALSE)
     m<-matrix(1,ro,co); x.pos<-col(m); y.pos<-(1+ro)-row(m)
  } 
  if(!missing(labels)) names(face.list)<-labels 
  fac.x<-width/1.1/210; fac.y<-height/1.3/210
  xtrans<-function(x){x.pos[j]+fac.x*x};  ytrans<-function(y){y.pos[j]+fac.y*y}
  for(j in seq(face.list)){
    face.obj<-face.list[[j]]; factors<-face.data[,j]
    f<-1+(ncolors-1)*(factors+1)/2 # translate factors into color numbers
    for(obj.ind in seq(face.obj)[c(10:11,1:9)]) {
       x <-face.obj[[obj.ind]][,1]; y<-face.obj[[obj.ind]][,2]
       xx<-spline(1:length(x),x,40,FALSE)[,2]
       yy<-spline(1:length(y),y,40,FALSE)[,2]
       lines(xtrans(xx),ytrans(yy),...)
       if(face.type>0){ if(obj.ind==10) 
                          polygon(xtrans(xx),ytrans(yy),col=col.hair[ceiling(mean(f[9:11]))],xpd=NA) # hair
                        if(obj.ind==11){ 
                          polygon(xtrans(xx),ytrans(yy),col=col.face[ceiling(mean(f[1:2 ]))],xpd=NA) # face
                          if(face.type==2){
                              # beard
                              for(zzz in seq(hhh<-max(face.obj[[8]][,1]),-hhh,length=30)){
                                hrx<-rnorm(8,zzz,2); hry<-0:7*-3*rnorm(1,3)+abs(hrx)^2/150
                                hry<-min(face.obj[[9]][,2])+hry
                                lines(xtrans(hrx),ytrans(hry),lwd=5,col="#eeeeee",xpd=NA)
                              }
                              ind<-which.max(xx); wx<-xx[ind]; ind<-which.max(yy); wy<-yy[ind]
                              # edge of hat
                              wxh<-wx<-seq(-wx,wx,length=20); wyh<-wy<-wy-(wx-mean(wx))^2/250+runif(20)*3
                              lines(xtrans(wxh),ytrans(wyh)); wx<-c(wx,rev(wx)); wy<-c(wy-10,rev(wy)+20)
                              wmxy1<-wmxy0<-c(min(wx),min(wy)+20)
                              wmxy2<-wmxy3<-c(runif(1,wmxy0[1],-wmxy0[1]), wy[1]+100)
                              wmxy1[2]<-0.5*(wmxy0[2]+wmxy3[2]); wmxy2[1]<-0.5*(wmxy2[1]+wmxy0[1])
                              npxy<-20; pxy<-seq(0,1,length=npxy)
                              gew<-outer(pxy,0:3,"^")*outer(1-pxy,3:0,"^")*
                                   matrix(c(1,3,3,1),npxy,4,byrow=TRUE)
                              wxl<-wmxy0[1]*gew[,1]+wmxy1[1]*gew[,2]+wmxy2[1]*gew[,3]+wmxy3[1]*gew[,4]
                              wyl<-wmxy0[2]*gew[,1]+wmxy1[2]*gew[,2]+wmxy2[2]*gew[,3]+wmxy3[2]*gew[,4]
                              lines(xtrans(wxl),ytrans(wyl),col="green")
                              wmxy1[1]<- wmxy0[1]<- -wmxy0[1]
                              wmxy1[2]<-0.5*(wmxy0[2]+wmxy3[2]); wmxy2[1]<-0.5*(wmxy2[1]+wmxy0[1])
                              wxr<-wmxy0[1]*gew[,1]+wmxy1[1]*gew[,2]+wmxy2[1]*gew[,3]+wmxy3[1]*gew[,4]
                              wyr<-wmxy0[2]*gew[,1]+wmxy1[2]*gew[,2]+wmxy2[2]*gew[,3]+wmxy3[2]*gew[,4]
                              points(xtrans(wmxy3[1]),ytrans(wmxy3[2]),pch=19,cex=2,col="#ffffff",xpd=NA)
                              points(xtrans(wmxy3[1]),ytrans(wmxy3[2]),pch=11,cex=2.53,col="red",xpd=NA)
                              polygon(xtrans(c(wxl,rev(wxr))),ytrans(c(wyl,rev(wyr))),col="red",xpd=NA) # hat
                              polygon(xtrans(wx),ytrans(wy),col="#ffffff",xpd=NA) # edge of hat
                          }

                        }
                        xx<-xtrans(xx); yy<-ytrans(yy)
                        if(obj.ind %in% 1:2) polygon(xx,yy,col="#eeeeee") # eyes without iris
                        if(obj.ind %in% 3:4) polygon(xx,yy,col=col.eyes[ceiling(mean(f[7:8 ]))],xpd=NA) # eyes:iris
                        if(obj.ind %in% 9)   polygon(xx,yy,col=col.nose[ceiling(mean(f[12:13]))],xpd=NA)# nose
                        if(obj.ind %in% 5:6) polygon(xx,yy,col=col.lips[ceiling(mean(f[1:3]))],xpd=NA)  # lips
                        if(obj.ind %in% 7:8) polygon(xx,yy,col=col.ears[ceiling(mean(f[14:15]))],xpd=NA)# ears
 }
    }
    lab<-names(face.list)[j]
    text(x.pos[j],y.pos[j]-0.5*height,lab,xpd=NA)
  }
}

