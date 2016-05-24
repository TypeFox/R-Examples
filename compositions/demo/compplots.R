if(FALSE) {
if(dev.cur() <= 1) get(getOption("device"))()

opar <- par(ask = interactive() &&
            (.Device %in% c("X11", "GTK", "gnome", "windows","quartz")))


cdata    <- rDirichlet.rcomp(12,c(Cu=3,Co=2,Na=6,S=0.1))
XXlanduse  <- factor(c("Field","Field","Field","Field","Field",
                     "Wood","Wood","Wood","Wood","Wood","House","House"))
                     
barplot(cdata)

barplot(mean.acomp(cdata))

pie(mean.acomp(cdata))

plot.acomp(cdata[,1:3])

plot.acomp(cdata,margin="margin",pca=T)

plot.acomp(cdata,margin="cmst",pca=T)

boxplot.acomp(cdata,XXlanduse,notch=T)

boxplot.acomp(cdata,log=F)

qqnorm.acomp(cdata,alpha=0.05)

biplot(princomp(cdata))

par(opar)
}

