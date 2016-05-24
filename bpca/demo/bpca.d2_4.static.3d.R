##
## Computing and ploting a bpca object with arbitrary choice of the first eigenvalue - 3d
##

oask <- devAskNewPage(dev.interactive(orNone=TRUE))

bp <- bpca(gabriel1971,
           d=2:4)

plot(bp,
     var.factor=2,
     xlim=c(-2,2),
     ylim=c(-2,2),
     zlim=c(-2,2))

# Exploring the object 'bp' created by the function 'bpca'
class(bp)
names(bp)
str(bp)

summary(bp)
bp$call
bp$eigenval
bp$eigenvec
bp$number
bp$import
bp$coord
bp$var.rb
bp$var.rd

# Changing the angle between x (PC2) and y (PC3) axis
plot(bp,
     var.factor=2,
     angle=65,
     xlim=c(-2,2),
     ylim=c(-2,2),
     zlim=c(-2,2))

devAskNewPage(oask)

