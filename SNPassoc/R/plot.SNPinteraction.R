`plot.SNPinteraction` <-
function(x, main.tit, ...)
{
 
 control<-apply(x,1,function(x) sum(is.na(x))==length(x))
 x.OK<-x[!control,!control]
 
 if (!is.null(attr(x,"gen.info"))){
        genInfo<-attr(x,"gen.info")
        o <- order(genInfo[, 2], genInfo[, 3])
        label.SNPs <- as.character(genInfo[o, 1])
        label.SNPs <- label.SNPs[label.SNPs%in%dimnames(x.OK)[[1]]]
        orderSNPs.ok<-match(label.SNPs, dimnames(x.OK)[[1]])
        x.OK <- x.OK[orderSNPs.ok,orderSNPs.ok ]
        genInfo <- genInfo[genInfo[,1]%in%label.SNPs,]
    }
    else {
        label.SNPs <- dimnames(x.OK)[[1]]
    }
 

 old.xpd <- par("xpd")
 old.las <- par("las")
 old.mfrow <- par("mfrow")
 
 par(xpd=NA)
 
 m <- matrix(1:2, 1, 2)
 layout(m, widths=c(4.5, 1))

 on.exit(par(xpd = old.xpd, mfrow = old.mfrow, las = old.las))

# Other palettes:
# mypaletteOld<-brewer.pal(9,"Greens")
# mypaletteOld<-c("#F7FCF5", "#E5F5E0", "#C7E9C0","#A1D99B", "#74C476", "#41AB5D", "#238B45", "#006D2C", "#00441B")
# mypaletteOld<-brewer.pal(9,"Reds")
#  mypaletteOld<- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D")
# This is used:  mypaletteOld<-brewer.pal(9,"YlGn")

 mypaletteOld <- c("#FFFFE5", "#F7FCB9", "#D9F0A3", "#ADDD8E", "#78C679", "#41AB5D", "#238443", "#006837", "#004529")

 mypalette<-mypaletteOld[c(9,6,4,3,3,2,2,1,1)]

 pvalCut<-c(0,0.001,0.01,0.05,0.1,0.2,0.3,0.5,0.7,1)

 image(1:nrow(x.OK),1:ncol(x.OK),x.OK,col=mypalette,breaks=pvalCut,
     axes=FALSE,xlab="",ylab="")
 
 axis(1,at=c(1:nrow(x.OK)),labels=label.SNPs,las=3,cex.axis=0.7,col="darkgreen")
 axis(2,at=c(1:nrow(x.OK)),labels=label.SNPs,las=1,cex.axis=0.7,col="darkgreen")
 

 if (missing(main.tit))
  main.tit<-paste("SNPs interactions --",attr(x,"model"),"model")

 title(main.tit,line=3)

 if (!is.null(attr(x,"gen.info"))) 
        n.snps <- table(genInfo[, 2])
 else n.snps <- nrow(x.OK)

 
 a <- c(0.5, cumsum(n.snps) + 0.5)

 b <- par("usr")
 segments(a, b[3], a, b[4] + diff(b[3:4]) * 0.02, col="darkblue",lwd=2)
 segments(b[3], a, b[4]+diff(b[3:4]) * 0.02, a, col="darkblue",lwd=2)

 abline(coef=c(0,1),xpd=FALSE,col="yellow")

 if(!is.null(attr(x,"gen.info")))
  {
   a <- par("usr")
   wh <- cumsum(c(0.5, n.snps))
   names.geno<-unique(genInfo[,2])
   n.gen<-length(names.geno)

   for (i in 1:n.gen)
    { 
      text(mean(wh[i + c(0, 1)]), a[4] + (a[4] - a[3]) * 0.025, names.geno[i],srt=45,cex=0.8,adj=0.2)
      text(a[4] + (a[4] - a[3]) * 0.025, mean(wh[i + c(0, 1)]), names.geno[i],srt=45,cex=0.8,adj=0.5)
    }
  }  

 
 image(0.5,1:10,matrix(pvalCut,nrow=1,ncol=10),col=rev(mypalette),breaks=pvalCut,axes=FALSE,
          xlab="",ylab="")
 marcas<-c(0.5,3.5,4.5,5.5,7.5,8.5,9.5,10.5)
 axis(2,marcas,rev(c(0,0.001,0.01,0.05,0.1,0.2,0.3,1)),pos=0.5)
 text(30,5.5,"pvalues",srt=90)

}

