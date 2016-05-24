##
## Diagnostic of iris representation with 'var.rd' parameter - 2d and 3d
##

oask <- devAskNewPage(dev.interactive(orNone=TRUE))

bp1 <- bpca(iris[-5],
            var.rb=TRUE,
            var.rd=TRUE,
            limit=3)

plot(bp1,
     var.factor=.3,
     obj.names=FALSE,
     obj.pch=c('+', '-', '*')[unclass(iris$Species)],
     obj.col=c('red', 'green3', 'blue')[unclass(iris$Species)],
     obj.cex=1)

bp1$var.rd
bp1$eigenvec

# Graphical diagnostic
plot(bpca(iris[-5],
          d=3:4),
     var.factor=.6,
     obj.names=FALSE,
     obj.pch=c('+', '-', '*')[unclass(iris$Species)],
     obj.col=c('red', 'green3', 'blue')[unclass(iris$Species)],
     obj.cex=1)

# Interpretation:
# Sepal.length followed by Petal.Width contains information in dimensions
# (PC3) - the PC3 is, essentially, a contrast among both) that wasn't fully
# contemplated by the biplot reduction (PC1 and PC2) .
# Therefore, between all variables, they have a "poor" representation by a 2d
# biplot.

bp2 <- bpca(iris[-5],
            d=1:3,
            var.rb=TRUE,
            var.rd=TRUE,
            limit=2)

plot(bp2,
     obj.names=FALSE, 
     obj.pch=c('+', '-', '*')[unclass(iris$Species)],
     obj.col=c('red', 'green3', 'blue')[unclass(iris$Species)],
     obj.cex=1)

bp2$var.rd

bp2$eigenvec

round(bp2$var.rb,
      2)

round(cor(iris[-5]),
      2)

# Good representation of all variables with a 3d biplot!

devAskNewPage(oask)

