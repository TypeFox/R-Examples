#24-12-2013 MRC-Epid JHZ

h2 <- read.table("nshd.dat",header=FALSE, sep="\t",col.names=c("y","id","h2","se"),as.is=TRUE)
h2 <- h2[with(h2,order(y,id)),]
h2 <- within(h2, {
      label <- c(rep("bmi",12),rep("ht",12),rep("wt",12))
      v <- se^2
      w <- 1/v
      age0 <- c(2,4,6,7,11,15,20,26,36,43,53,63)
      age <- c(2,4.25,6.04,7.03,10.90,14.69,20,26.22,36.29,43.48,53.46,63.33)
      age2 <- age^2
      age3 <- age^3
})
library(metafor)
ma <-rma(h2,v,mods=~age+age2+age3,data=subset(h2,y=="bmi"))
ma
coef(ma)
diag(vcov(ma))
library(mgcv)
tiff("nshd.tiff",height=11.69,width=8.26,units="in",res=1200,compress="lzw")
l <- -7.75
with(subset(h2,y=="bmi"),{
   fit <<- gam(h2~s(age,bs="cr"),weights=w)
   plot(fit,se=TRUE,rug=FALSE,ylab=expression(h^2), xlab="", axes=FALSE, 
        main="Heritability estimates of BMI by age",shift=coef(fit)[1])
   axis(1,at=(0:7)*10,line=l)
   rug(age,line=l)
   mtext("age",1,line=l+3)
   axis(2,at=(0:5)*0.05)
   print(summary(fit))
   p <- predict(fit,se.fit=TRUE)
   print(p)
})
dev.off()
