##
## Computing and ploting a bpca object with 'scatterplot3d' package - 3d
##

oask <- devAskNewPage(dev.interactive(orNone=TRUE))

bp <- bpca(gabriel1971,
           d=1:3)

plot(bp,
     var.factor=3)

# Exploring the object 'bp' created by the function 'bpca'
class(bp)
names(bp)
str(bp)

summary(bp)
bp$call
bp$eigenval
bp$eigenvec
bp$numb
bp$import
bp$coord
bp$coord$obj
bp$coord$var
bp$var.rb
bp$var.rd

# Additional graphical parameters (nonsense)
plot(bpca(gabriel1971,
          d=1:3,
          meth='jk'),
     main='gabriel1971 - jk',
     sub='The graphical parameters are working fine!',
     var.factor=6,
     var.pch='+',
     var.cex=.6,
     var.col='green4',
     obj.pch='*',
     obj.cex=.8,
     obj.col=1:8,
     ref.lty='solid',
     ref.col='red',
     angle=70)

devAskNewPage(oask)

