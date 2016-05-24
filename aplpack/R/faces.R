faces<-function(xy,which.row,fill=FALSE,face.type=1,
                nrow.plot,ncol.plot,scale=TRUE,byrow=FALSE,main,
                labels,print.info = TRUE,na.rm = FALSE,
                ncolors=20,
                col.nose=rainbow(ncolors),                   # nose
                col.eyes=rainbow(ncolors,start=0.6,end=0.85),# eyes
                col.hair=terrain.colors(ncolors),            # hair
                col.face=heat.colors(ncolors),               # face
                col.lips=rainbow(ncolors,start=0.0,end=0.2), # lips
                col.ears=rainbow(ncolors,start=0.0,end=0.2), # ears

                plot.faces=TRUE){  # 070831 pwolf
  if((demo<-missing(xy))){
    xy<-rbind(
              c(1,3,5),c(3,5,7),
              c(1,5,3),c(3,7,5),
              c(3,1,5),c(5,3,7),
              c(3,5,1),c(5,7,3),
              c(5,1,3),c(7,3,5),
              c(5,3,1),c(7,5,3),
              c(1,1,1),c(4,4,4),c(7,7,7)
    )
    labels<-apply(xy,1,function(x) paste(x,collapse="-"))
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

  n.char<-15
  xy<-rbind(xy)
  if(byrow) xy<-t(xy)
  if(any(is.na(xy))){
    if(na.rm){ 
      xy<-xy[!apply(is.na(xy),1,any),,drop=FALSE]
      if(nrow(xy)<3) {print("not enough data points"); return()}
      print("Warning: NA elements have been removed!!")
     }else{
      xy.means<-colMeans(xy,na.rm=TRUE)
      for(j in 1:length(xy[1,])) xy[is.na(xy[,j]),j]<-xy.means[j]
      print("Warning: NA elements have been exchanged by mean values!!")
    }  
  }
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
  # expand input matrix xy by replication of cols
  xy<-xy[,(rows.orig<-h<-rep(1:mm,ceiling(n.char/mm))),drop=FALSE]
  if(fill) xy[,-(1:n.c)]<-0

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

  nr<-n^0.5; nc<-n^0.5
  if(!missing(nrow.plot)) nr<-nrow.plot
  if(!missing(ncol.plot)) nc<-ncol.plot
  if(plot.faces){
    opar<-par(mfrow=c(ceiling(c(nr,nc))),oma=rep(6,4), mar=rep(.7,4))
    on.exit(par(opar))
  }

  face.list<-list()
  for(ind in 1:n){
    factors<-xy[ind,]
    face<-face.orig

    m<-mean(face$lipso[,2])
    face$lipso[,2]<-m+(face$lipso[,2]-m)*(1+0.7*factors[4])
    face$lipsi[,2]<-m+(face$lipsi[,2]-m)*(1+0.7*factors[4])
    face$lipso[,1]<-face$lipso[,1]*(1+0.7*factors[5])
    face$lipsi[,1]<-face$lipsi[,1]*(1+0.7*factors[5])
    face$lipso["lipsiend",2]<-face$lipso["lipsiend",2]+20*factors[6]

    m<-mean(face$eye[,2])
    face$eye[,2] <-m+(face$eye[,2] -m)*(1+0.7*factors[7])
    face$iris[,2]<-m+(face$iris[,2]-m)*(1+0.7*factors[7])
    m<-mean(face$eye[,1])
    face$eye[,1] <-m+(face$eye[,1] -m)*(1+0.7*factors[8])
    face$iris[,1]<-m+(face$iris[,1]-m)*(1+0.7*factors[8])


    m<-min(face$hair[,2])
    face$hair[,2]<-m+(face$hair[,2]-m)*(1+0.2*factors[9])
    m<-0
    face$hair[,1]<-m+(face$hair[,1]-m)*(1+0.2*factors[10])
    m<-0
    face$hair[c("hair1","hair2"),2]<-face$hair[c("hair1","hair2"),2]+50*factors[11]

    m<-mean(face$nose[,2])
    face$nose[,2]<-m+(face$nose[,2]-m)*(1+0.7*factors[12])
    face$nose[nose.xnotnull,1]<-face$nose[nose.xnotnull,1]*(1+factors[13])

    m<-mean(face$shape[c("earsta","earend"),1])
    face$ear[,1]<-m+(face$ear[,1]-m)* (1+0.7*factors[14])
    m<-min(face$ear[,2])
    face$ear[,2]<-m+(face$ear[,2]-m)* (1+0.7*factors[15])

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
    face.list<-c(face.list,list(face.obj))

    if(plot.faces){
      plot(1,type="n",xlim=c(-105,105)*1.1, axes=FALSE,
           ylab="",ylim=c(-105,105)*1.3)
      title(xnames[ind])
      f<-1+(ncolors-1)*(factors+1)/2 # translate factors into color numbers
      xtrans<-function(x){x};  ytrans<-function(y){y}
      for(obj.ind in seq(face.obj)[c(10:11,1:9)]) {
        x <-face.obj[[obj.ind]][,1]; y<-face.obj[[obj.ind]][,2]
        xx<-spline(1:length(x),x,40,FALSE)[,2]
        yy<-spline(1:length(y),y,40,FALSE)[,2]
        if(plot.faces){ 
          lines(xx,yy)
          if(face.type>0){
            if(obj.ind==10) 
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
      }
    }

  }

  if(plot.faces&&!missing(main)){
    par(opar);par(mfrow=c(1,1))
    mtext(main, 3, 3, TRUE, 0.5)
    title(main)
  }

  info<-c(
    "var1"="height of face   ",
    "var2"="width of face    ",
    "var3"="structure of face",
    "var4"="height of mouth  ",
    "var5"="width of mouth   ",
    "var6"="smiling          ",
    "var7"="height of eyes   ",
    "var8"="width of eyes    ",
    "var9"="height of hair   ",
    "var10"="width of hair   ",
    "var11"="style of hair   ",
    "var12"="height of nose  ",
    "var13"="width of nose   ",
    "var14"="width of ear    ",
    "var15"="height of ear   ")
  var.names<-dimnames(xy)[[2]]
  if(0==length(var.names)) var.names<-paste("Var",rows.orig,sep="")
  info<-cbind("modified item"=info,"Var"=var.names[1:length(info)])
  rownames(info)<-rep("",15)
  if(print.info){
    cat("effect of variables:\n")
    print(info)
  }
  if(demo&&plot.faces) {
    plot(1:15,1:15,type="n",axes=FALSE,bty="n")
    text(rep(1,15),15:1,adj=0,apply(info,1,function(x) 
        paste(x,collapse="   -   ")),cex=0.7)
  }

  names(face.list)<-xnames
  out<-list(faces=face.list,info=info,xy=t(xy))
  class(out)<-"faces"
  invisible(out)

}

