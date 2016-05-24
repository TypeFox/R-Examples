boxplot2D<-function( xy, add.to.plot=TRUE, box.size=10, box.shift=0, angle=0,
                angle.type="0", tukey.style=TRUE, coef.out=1.5, coef.h.out=3, design="sl",
                na.rm=FALSE, ... ){
##########################################################################
# boxplot for scatterplots, pw 03/05                                     #
# xy:          2-col matrix                                              #
# add.to.plot: if TRUE => plot(xy,...)                                   #
# box.size:    height of box in mm                                       #
# box.shift:   shift of boxplot in mm                                    #
# angle:       direction of projection in units defined by angle.type    #
# angle.type:  "0"=     angle    , # angle in (0,2*pi)                   #
#              "1"=2*pi*angle/12,  # clock-like                          #
#              "2"=2*pi*angle/360, # angle in (0,360)                    #
#              "3"=arctan(angle)   # by fraction: delta.y/delta.x        #
# tukey.style: if TRUE outliers are defined as described in Tukey (1977): EDA
# coef.out=1.5  # outliers are defined outside coef.out*boxwidth
# coef.h.out=3   # heavy outliers are defined outside coef.h.out*boxwidth
# design:      if "sl" then parallelogram else box                       #
##########################################################################
if(any(is.na(xy))){
  if(na.rm){ xy<-xy[!apply(is.na(xy),1,any),,drop=FALSE]
    print("Warning: NAs elements have been removed!!")
  }else{
    xy.means<-colMeans(xy,na.rm=TRUE)
    for(j in 1:ncol(xy)) xy[is.na(xy[,j]),j]<-xy.means[j]
    print("Warning: NAs elements have been exchanged by mean values!!")
  }  
}
if(!add.to.plot) plot(xy,...)

if(is.numeric(angle.type)){
  angle.type<-as.character(angle.type)
  if(all(angle.type!=c("0","1","2","3"))) angle.type<-"0"
}
w <- switch(angle.type, "0"=     angle, "1"=2*pi*(3-angle)/12,
                        "2"=2*pi*angle/360, "3"=atan(angle) )
TM <- matrix(c(cos(w),sin(w),-sin(w),cos(w)),2,2)

xyt<- xy %*% TM
##ermittle empirsche Quantile##=
z  <- xyt[,1]
if(tukey.style){
  z.stats1 <- boxplot.stats(z,coef=coef.out)
  z.stats2 <- boxplot.stats(z,coef=coef.h.out)
  z<-c(min(z),z.stats1$stats,max(z))
  names(z) <- c("min","fence","hinge","median","hinge","fence","max")
  outlier.heavy <- z.stats2$out
  outlier <- z.stats1$out
  outlier <- outlier[! outlier %in% outlier.heavy ]
}else{
  z  <-     c(min(z),quantile(z,c(0.10, 0.25, 0.5, 0.75, 0.90)),max(z))
  names(z)<-c("min",              ".1", ".25",".5",".75",".9", "max")
}
xy.q <- cbind(z,median(xyt[,2])) %*% t(TM)
if(tukey.style){
  xy.out<-if(0<length(outlier)) cbind(outlier,median(xyt[,2])) %*% t(TM) else numeric(0)
  xy.out.heavy<-if(0<length(outlier.heavy))
                                     cbind(outlier.heavy,median(xyt[,2])) %*% t(TM) else numeric(0)
}

# d. steht fuer delta.:
d.xy <-apply(xy.q[c(1,7),],           2,diff)
d.usr<-apply(matrix(par()$usr,2,2),2,diff)
d.pmm<-par()$pin*25.4
if(design=="sl"){
    # d.xy.s.mm: bzgl. Welt-coor senkrechten Vektor in mm-coor
    d.xy.s.mm <- (c(-1,1)*rev(d.xy))/d.usr*d.pmm
    # d.xy.0.s: bzgl. mm-coor senkrechten Einheits-Vektor
    d.xy.mm.s <- c(-1,1)*rev( d.xy/d.usr*d.pmm )
    d.xy.0.s  <- d.xy.mm.s/(d.xy.mm.s%*%d.xy.mm.s)^0.5
    # Streckung von d.xy.s.mm so, dass Projektion auf d.xy.0.s von Laenge 1
    d.xy.s.mm.proj1 <- d.xy.s.mm / (d.xy.s.mm %*% d.xy.0.s)
    # Boxsize: in mm
    real.step.mm <- box.size * d.xy.s.mm.proj1
    # Darstellung in Welt-coor
    real.step <- real.step.mm/d.pmm*d.usr

}else{
    # Einheitsvektor zu d.xy in mm-coor
    d.xy.mm   <- d.xy/d.usr*d.pmm
    d.xy.mm.ev<- d.xy.mm/(d.xy.mm%*%d.xy.mm)^0.5
    # Konstruktion eines Boxhoehenvektors der Laenge box.size in mm
    mm.step   <- box.size*c(-1,1)*rev(d.xy.mm.ev)
    # Darstellung in Welt-coor
    real.step <- mm.step/d.pmm*d.usr # real.step fuer Boxdickenvektor

}

# Boxshift: in mm
real.shift<-box.shift*real.step/box.size
xy.q.orig <- xy.q
xy.q      <-cbind(xy.q[,1]+real.shift[1], xy.q[,2]+real.shift[2])
if(tukey.style){
  if(0<length(xy.out)) xy.out<-cbind(xy.out[,1]+real.shift[1], xy.out[,2]+real.shift[2])
  if(0<length(xy.out.heavy))
            xy.out.heavy<-cbind(xy.out.heavy[,1]+real.shift[1], xy.out.heavy[,2]+real.shift[2])
}


# Boxplotkonstruktion
hs<-real.step/2
xyxy <-c(
          #min:xy.q[1,1]-hs[1],xy.q[1,2]-hs[2],xy.q[1,1]+hs[1],xy.q[1,2]+hs[2],
          xy.q[2,1]      ,xy.q[2,2]      ,xy.q[3,1]      ,xy.q[3,2]      ,
          xy.q[3,1]-hs[1],xy.q[3,2]-hs[2],xy.q[5,1]-hs[1],xy.q[5,2]-hs[2],
          xy.q[3,1]-hs[1],xy.q[3,2]-hs[2],xy.q[3,1]+hs[1],xy.q[3,2]+hs[2],
          xy.q[4,1]-hs[1],xy.q[4,2]-hs[2],xy.q[4,1]+hs[1],xy.q[4,2]+hs[2],
          xy.q[5,1]-hs[1],xy.q[5,2]-hs[2],xy.q[5,1]+hs[1],xy.q[5,2]+hs[2],
          xy.q[3,1]+hs[1],xy.q[3,2]+hs[2],xy.q[5,1]+hs[1],xy.q[5,2]+hs[2],
          xy.q[6,1]      ,xy.q[6,2]      ,xy.q[5,1]      ,xy.q[5,2]
          #max ,xy.q[7,1]-hs[1],xy.q[7,2]-hs[2],xy.q[7,1]+hs[1],xy.q[7,2]+hs[2]
        )
xyxy<-matrix(xyxy,length(xyxy)/4,4,TRUE)
segments(xyxy[,1],xyxy[,2],xyxy[,3],xyxy[,4])
if(tukey.style){
  if(0<length(xy.out)) points(xy.out[,1],xy.out[,2],pch=1,cex=1.5)
  if(0<length(xy.out)) points(xy.out[,1],xy.out[,2],pch=18)
  if(0<length(xy.out.heavy)) points(xy.out.heavy[,1],xy.out.heavy[,2],pch=19,cex=1.5)
}else{
  points(xy.q[c(1,7),1],xy.q[c(1,7),2],pch=19)
}
} # end of boxplot2D

