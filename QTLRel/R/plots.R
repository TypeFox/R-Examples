
plot.scanOne<- function(x,...){
   xTmp<- list(...)
   if(is.null(xTmp$cex)) xTmp$cex<- 0.5
   if(is.null(xTmp$main)) xTmp$main<- ""
   if(is.null(xTmp$xlab)) xlab<- "Chromosome"
   cv<- xTmp$cv

   lrt<- data.frame(y=x$p)
   if(is.null(x$chr)){
      if(is.null(xTmp$gmap)) stop("need: gmap...")
      gmap<- xTmp$gmap
      idx<- match(names(x$p),gmap$snp)
      lrt$chr<- gmap$chr[idx]
      lrt$dist<- gmap$dist[idx]
   }else{
      lrt$chr=x$chr
      lrt$dist=x$dist
   }
   if(is.element("None",class(x))){
      lrt$y<- x$p/(2*log(10))
      plot.lrt(lrt,cv,cex=xTmp$cex,main=xTmp$main,ylim=xTmp$ylim,xlab=xlab,ylab="LOD")
   }else{
      lrt$y<- -log10(x$p)
      plot.lrt(lrt,cv,cex=xTmp$cex,main=xTmp$main,ylim=xTmp$ylim,
         xlab=xlab,ylab=expression(paste(-log[10],"(p-value)")))
   }
}

plot.lrt<- function(lrt,cv,...){
# lrt: data.frame(y,chr,dist,...)
   lrt$chr<- reorder(factor(lrt$chr))
   lrt<- lrt[order(lrt$chr,lrt$dist),]

   chr<- unique(lrt$chr)
      nchr<- length(chr)
   for(i in 1:nchr){
      idx<- lrt$chr==chr[i]
      lrt$dist[idx]<- diff(c(0,lrt$dist[idx]))
   }
   lrt$dist<- cumsum(lrt$dist)

   plot(range(lrt$dist),range(lrt$y),type="n",xaxt="n",...)

   col<- NULL
   for(i in 1:nchr){
      idx<- lrt$chr==chr[i]
      col<- c(col,rep(i%%2 + 3,sum(idx)))
   }
   lines(lrt$dist,lrt$y,type="p",col=col,...)
   idx<- c(TRUE,lrt$chr[-length(lrt$chr)]!=lrt$chr[-1])
   min.p<- lrt$dist[idx]
   min.p<- c(min.p,max(lrt$dist))

   if(!missing(cv)) abline(h=cv,col=2,lty=4,...)
   abline(v=min.p,lty=3,lwd=0.1)
   pos<- (min.p[-1]+min.p[-c(nchr+1)])/2
   axis(1,at=pos,labels=chr,tick=T,las=2)
}


plotit<- function(lrt,cv,bychr=FALSE,chr.labels=TRUE,
   type="p",lty=NULL,col=NULL,pch=NULL,cex=NULL,...){
# lrt: data.frame(y,chr,dist,group,...)
   hh<- NULL
   if(!missing(cv) && !is.null(cv)) hh<- cv

   if(bychr){
      groups<- NULL
      if(!is.null(lrt$group)){
         groups=lrt$group
      }

      if(length(unique(lrt$chr))>1){
         xyplot(y~dist|chr,data=lrt,
            groups=groups,
            panel=function(x,y,...){
               panel.xyplot(x,y,...)
               if(!is.null(hh)) panel.abline(h=hh,...)
               if(!is.null(lrt$group)) panel.superpose(x,y,...)
            },
            type=type,
            lty=lty,
            col=col,
            pch=if(!is.null(pch)) pch else 1,
            cex=cex,
            ...
         )
      }else{
         xyplot(y~dist,data=lrt,
            groups=groups,
            panel=function(x,y,...){
               panel.xyplot(x,y,...)
               if(!is.null(hh)) panel.abline(h=hh,...)
               if(!is.null(lrt$group)) panel.superpose(x,y,...)
            },
            type=type,
            lty=lty,
            col=col,
            pch=if(!is.null(pch)) pch else 1,
            cex=cex,
            ...
         )
      }
   }else{
      if(is.null(lrt$group)) lrt$group<- 1
      lrt<- as.data.frame(lrt)
      lrt$chr<- reorder(factor(lrt$chr))
      lrt$group<- reorder(factor(lrt$group))
      lrt<- lrt[order(lrt$group,lrt$chr,lrt$dist),]

      groups<- lrt$group
         groups<- sort(unique(groups))
      ngr<- length(groups)
      if(ngr>1) cat("  Groups:",as.character(groups),"\n")
      chr<- unique(lrt$chr)
         nchr<- length(chr)
      for(i in 1:ngr){
         idx<- lrt$group==groups[i]
         lrt0<- lrt[idx,]
         chr0<- unique(lrt0$chr)
            nchr0<- length(chr0)
         for(i in 1:nchr0){
            idx0<- lrt0$chr==chr0[i]
            lrt0$dist[idx0]<- diff(c(0,lrt0$dist[idx0]))
         }
         lrt0$dist<- cumsum(lrt0$dist)
         lrt$dist[idx]<- lrt0$dist
      }
      if(chr.labels){
         plot(range(lrt$dist),range(lrt$y),type="n",xaxt="n",...)
      }else{
         plot(range(lrt$dist),range(lrt$y),type="n",...)
      }

      if(ngr>1){
         if(is.null(lty)){
            lty<- 1:5
            lty<- rep(lty,ngr)[1:ngr]
         }
         if(is.null(col)){
            col<- 1:ngr
         }
         pch<- rep(pch,ngr); pch<- pch[1:ngr]
         cex<- rep(cex,ngr); cex<- cex[1:ngr]
      }else{
         if(is.null(col)) col<- 1
      }
      if(length(col) < length(lrt$dist))
         colTmp<- col[match(lrt$group, groups)]
      else colTmp<- col

      min.p<- matrix(NA,nrow=ngr,ncol=nchr)
      for(g in 1:ngr){
         idx0<- lrt$group==groups[g]
         lrt0<- lrt[idx0,]
         col0<- colTmp[idx0]
         for(i in 1:nchr){
            idx<- lrt0$chr==chr[i]
            lines(lrt0$dist[idx],lrt0$y[idx],type=type,lty=lty[g],col=col0[idx],pch=pch[g],cex=cex[g],...)
         }
         idx<- c(TRUE,lrt0$chr[-length(lrt0$chr)]!=lrt0$chr[-1])
         min.p[g,]<- lrt0$dist[idx]
      }
      min.p<- apply(min.p,2,min)
      min.p<- c(min.p,max(lrt$dist))

      if(!missing(cv)) abline(h=cv,lty=lty,col=col)
      abline(v=min.p,lty=3,lwd=0.1)
      if(chr.labels){
         pos<- (min.p[-1]+min.p[-c(nchr+1)])/2
         axis(1,at=pos,labels=chr,tick=F)
      }
   }
}


plot.scanTwo<- function(x,...){
# x: object of scanTwo
# a genetic map 'gmap' is needed
   lst<- list(...)
   if(is.null(lst$gmap))
      stop("need a genetic map 'gmap'.")
   qqint<- x
   gmap<- lst$gmap
   v<- as.matrix(qqint)
   rv<- range(v,na.rm=TRUE)
      rv<- seq(floor(rv[1]),ceiling(rv[2]),by=0.5)
   dst<- gmap$dist
      dst<- diff(dst)
      dst[dst<0]<- 0
      dst[dst==0]<- 1e-5
      dst<- c(0,dst)
      dst<- cumsum(dst)
   chrs<- unique(gmap$chr)
   xat<- NULL
   for(chr in chrs){
      xat<- c(xat,mean(range(dst[gmap$chr==chr])))
   }

   scale<- max(dst)/(max(rv)-min(rv))*0.75
   image(x=dst,y=dst,z=v,axes=F,xlab=lst$xlab,xlim=c(-2,max(dst)),
      ylim=c(0,max(dst)+2),ylab=lst$ylab,
      main=lst$main,col=terrain.colors(12))
   image(x=max(dst)*c(19/20,1),y=(rv-min(rv))*scale,z=matrix(rv,nrow=1),
      col=terrain.colors(12),add=TRUE)
   axis(4,at=(rv-min(rv))*scale,labels=rv,pos=max(dst)*19.5/20,tick=FALSE)
   axis(2,at=xat,labels=chrs,tick=F,pos=max(dst)*0.025)
   axis(3,at=xat,labels=chrs,tick=F,pos=max(dst)*0.975)
   mtext("Chromosomes",2,line=1.5)
   mtext("Chromosomes",3,line=1.25,at=0)
   col<- 0
   for(chr in chrs){
      col<- col+1
      lines(c(-3,-3),range(dst[gmap$chr==chr]),col=col,lwd=5)
      lines(range(dst[gmap$chr==chr]),max(dst)-c(-3,-3),
      col=col,lwd=5)
   }
}

