##
## Comparative example from YAN, W & KANG, M.S. GGE biplot analysis:
## a graphical tool for breeders, geneticists, and agronomists.
##

oask <- devAskNewPage(dev.interactive(orNone=TRUE))

bp <- bpca(t(gge2003),
           var.rb=TRUE)

as.dist(bp$var.rb)

op <- par(no.readonly=TRUE)

par(mfrow=c(1,2))

plot(bpca(gge2003,
          var.pos=2),
     main='Columns as variables \n (var.pos=2)',
     var.col=1,
     obj.col=c(2:4, 2),
     obj.cex=.8)

plot(bpca(gge2003,
          var.pos=1),
     main='Rows as variables \n (var.pos=1)',
     var.col=1,
     obj.col=2:4,
     obj.cex=.8)

par(op)

devAskNewPage(oask)

