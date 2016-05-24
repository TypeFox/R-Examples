NicheOverlap<-function(data, Level1, Taxon1, Level2=Level1, Taxon2,
colA=hsv(h=0,s=1,v=1, alpha=0.4),colB=hsv(h=0.7,s=1,v=1, alpha=0.4),
xlab="Polar coordinate X in pixel", ylab="Polar coordinate Y in pixels",
cex=1.57, cex.lab=1.5,font.lab=1, main="",cex.main = 2, font.main=2,
family="serif", digits =2, xlegend="topleft",ylegend=NULL, pch=15,
bty="n",text.font=3, cex.legend=1.2, ncol=1, x.intersp = 1, y.intersp = 1,
legend=TRUE){

datos<-na.exclude(data)



if (missing(Level2)) Level2<-Level1 else Level2<-Level2
if (missing(colA)) colA<-hsv(h=0,s=1,v=1, alpha=0.4) else colA<-colA
if (missing(colB)) colB<-hsv(h=0.7,s=1,v=1, alpha=0.4) else colB<-colB
if (missing(xlab)) xlab<-"Polar coordinate X in pixels" else xlab<-xlab
if (missing(ylab)) ylab<-"Polar coordinate Y in pixels" else ylab<-ylab
if (missing(cex)) cex<-1.57 else cex<-cex
if (missing(cex.lab)) cex.lab<-1.5 else cex.lab<-cex.lab
if (missing(font.lab)) font.lab<-1 else font.lab<-font.lab
if (missing(main)) main<-"" else main<-main
if (missing(cex.main)) cex.main<-2 else cex.main<-cex.main
if (missing(font.main)) font.main<-2 else font.main<-font.main
if (missing(family)) family<-"serif" else family<-family
if (missing(digits)) digits<-2 else digits<-digits
if (missing(xlegend)) xlegend<-"topleft" else xlegend<-xlegend
if (missing(ylegend)) ylegend<-NULL else ylegend<-ylegend
if (missing(pch)) pch<-15 else pch<-pch
if (missing(bty)) bty<-"n" else bty<-bty
if (missing(text.font)) text.font<-3 else text.font<-text.font
if (missing(cex.legend)) cex.legend<-1.2 else cex.legend<-cex.legend
if (missing(ncol)) ncol<-1 else ncol<-ncol
if (missing(x.intersp)) x.intersp<-1 else x.intersp<-x.intersp
if (missing(y.intersp)) y.intersp<-1 else y.intersp<-y.intersp
if (missing(legend)) legend<-TRUE else legend<-legend

A<-subset(datos,datos[,Level1] %in% Taxon1)
B<-subset(datos,datos[,Level2] %in% Taxon2)


###################
#Overlap

AA<-cbind(A$Pixel.X,A$Pixel.Y)
AA<-unique(AA)
DA<-dim(AA)

BB<-cbind(B$Pixel.X,B$Pixel.Y)
BB<-unique(BB)
DB<-dim(BB)

CC<-rbind(AA,BB)
DUP<-duplicated(CC)
DC<-length(DUP[DUP==TRUE])


OA<-DC*100/DA[1]
OB<-DC*100/DB[1]

OA<-round(OA,digits=digits)
OB<-round(OB,digits=digits)

####################

par(family=family)

#Plot
par(mfrow = c(1,1))
plot(A$Pixel.X,A$Pixel.Y, pch=pch, cex=cex, col=colA, xlim=c(0,30), ylim=c(0,30),
xlab=xlab, ylab=ylab, cex.lab=cex.lab, font.lab=font.lab, main=main, cex.main=cex.main,
font.main=font.main)
points(B$Pixel.X,B$Pixel.Y, pch=pch, cex=cex, col=colB)


if(legend==TRUE){
legend(x=xlegend,y=ylegend, legend=c(paste(Taxon1,OA,"%"),paste(Taxon2,OB,"%")), col=c(colA,colB), pch=pch, bty=bty,
text.font=text.font, cex=cex.legend, x.intersp=x.intersp, y.intersp=y.intersp, ncol=ncol)
}
else{
www<-1
}

}
