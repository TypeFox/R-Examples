shape2d<-function(shtseq, gnum=500, type="radius", type2="slice", 
gnum2=1000, ngrid=30, norma=FALSE, xmax=10, modelim=2, exmalim=NULL,
maxnum=NULL)
{
# type "proba"    type2 "boundary"
lkm<-length(shtseq$level)
d<-length(shtseq$pcf$N)

if (type2=="slice"){

  if (type=="radius") x<-shtseq$level else x<-matrix(0,lkm,1)

  td<-shtseq$shtseq[[1]]
  if (type=="proba") td$volume<-td$proba
  td<-treedisc(td,shtseq$pcf,ngrid=ngrid)
  xy<-lst2xy(td,gnum=gnum)
  ylen<-length(xy$x)
  ystep<-1/(ylen-1)
  y<-seq(0,1,ystep)   #matrix(0,xlen,1)
  z<-matrix(0,length(x),length(y))
  delineator<-matrix(0,10*length(x),d)
  delinrun<-1
  delineatorlevel<-matrix(0,10*length(x),1)
  
  delineator.redu<-matrix(0,4*length(x),d)
  dr.redu<-1
  delineatorlevel.redu<-matrix(0,4*length(x),1)

  for (i in 1:lkm){
     td<-shtseq$shtseq[[i]]

     if (type=="proba"){ 
         tdvolume<-td$volume
         td$volume<-td$proba
         indi<-lkm-i+1
         voluu<-max(tdvolume)  #[1]  #root=1
         if (norma) x[indi]<-(voluu/volball(1,d))^(1/d)
         else x[indi]<-voluu 
     }
     else indi<-i
     td<-treedisc(td,shtseq$pcf,ngrid=ngrid)
     if (length(td$parent)==1) ynew<-0
     else{
        xy<-lst2xy(td,gnum=gnum)   #ma<-matchxy(xy$x,xy$y,y)

        ## normalize
        volu<-xy$x[length(xy$x)]-xy$x[1]
        int<-0
        step<-xy$x[2]-xy$x[1]
        for (j in 1:length(xy$x)){
            int<-int+step*xy$y[j]
        }
        if (norma){
            normavolu<-(volu/volball(1,d))^(1/d)
            b<-volu*normavolu/int
        }
        else b<-volu^2/int
        ynew<-b*xy$y
        ## end normalize

        # location
        ml<-moodilkm(td$parent)
        mc<-t(td$center[,ml$modloc])  #modecent(td)
        modenum<-dim(mc)[1]
        delineator[delinrun:(delinrun+modenum-1),]<-mc     
        delineatorlevel[delinrun:(delinrun+modenum-1)]<-x[indi] 
        delinrun<-delinrun+modenum

        if (modenum>modelim){
            prunum<-modenum-modelim
            pru<-prunemodes(td,prunum,exmalim,num=maxnum)
        }
        else pru<-td 
        ml<-moodilkm(pru$parent)
        mc<-t(pru$center[,ml$modloc])  #modecent(td)
        modenum<-dim(mc)[1]
        delineator.redu[dr.redu:(dr.redu+modenum-1),]<-mc     
        delineatorlevel.redu[dr.redu:(dr.redu+modenum-1)]<-x[indi] 
        dr.redu<-dr.redu+modenum
     }
     z[indi,]<-ynew   
  }
  delineator<-delineator[1:(delinrun-1),]
  delineatorlevel<-delineatorlevel[1:(delinrun-1)]

  delineator.redu<-delineator.redu[1:(dr.redu-1),]
  delineatorlevel.redu<-delineatorlevel.redu[1:(dr.redu-1)]
}

else{ #type2=="boundary"

if (is.null(xmax)){
    td<-shtseq$shtseq[[1]]
    if (type=="proba") td$volume<-td$proba
    xmax<-max(td$volume)
}

ymax<-xmax
step<-2*xmax/(gnum-1)
x<-seq(-xmax,xmax,step)
y<-x
z<-matrix(0,length(x),length(y))

for (i in 1:lkm){
  td<-shtseq$shtseq[[i]]
  if (type=="proba") td$volume<-td$proba
  xy<-lst2xy(td,gnum=gnum2,type=type)  

  ## normalize
  volu<-xy$x[length(xy$x)]-xy$x[1]
  int<-0
  step<-xy$x[2]-xy$x[1]
  for (j in 1:length(xy$x)){
       int<-int+step*xy$y[j]
  }
  b<-volu^2/int
  ynew<-b*xy$y
  ## end normalize

  for (j in 1:length(x)){
      for (k in 1:length(y)){
          len<-sqrt(x[j]^2+y[k]^2)
          xn<-x[j]/len
          yn<-y[k]/len
          th2<-atan(xn/yn)
          if (yn<0) th2<-atan(xn/yn)+pi else if (xn<0) th2<-atan(xn/yn)+2*pi
          propo<-th2/(2*pi) 
          dirind<-max(1,round( propo*length(xy$x) ))
          rho<-ynew[dirind]
          if (len<=rho) z[j,k]<-shtseq$level[i]
      }
  }
}

}

return(list(x=x,y=y,z=z,type=type,type2=type2,norma=norma,
            delineator=delineator,delineatorlevel=delineatorlevel,
            delineator.redu=delineator.redu,
            delineatorlevel.redu=delineatorlevel.redu))
}


