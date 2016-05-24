##
## Example of 'var.rd=TRUE' parameter as a measure of the quality of the biplot - 2d
## Mainly recommended to large datasets.
##

oask <- devAskNewPage(dev.interactive(orNone=TRUE))

bp <- bpca(gabriel1971,
           meth='hj',
           var.rb=TRUE,
           var.rd=TRUE,
           limit=3)

bp$var.rd

# RUR followed by CRISTIAN contains information in dimensions that
# weren't contemplated by the biplot reduction (PC3).
# Between all, RUR followed by CRISTIAN, variables are bad represented by a 2d
# biplot.

plot(qbpca(gabriel1971,
           bp))

# Graphical visualization of the importance of the variables not contemplated
# in the prior (2d) reduction
plot(bpca(gabriel1971,
          meth='hj',
          d=3:4),
     main='hj',
     xlim=c(-1,1),
     ylim=c(-1,1))

devAskNewPage(oask)

