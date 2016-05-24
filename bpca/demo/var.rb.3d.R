##
## Example of 'var.rb=TRUE' parameter as a measure of the quality of the biplot - 3d
##

oask <- devAskNewPage(dev.interactive(orNone=TRUE))

## Differences between methods of factorization
# SQRT
bp1 <- bpca(gabriel1971,
            meth='sqrt',
            d=1:3,
            var.rb=TRUE)

qbp1 <- qbpca(gabriel1971,
              bp1)

plot(qbp1,
     main='sqrt - 3d \n (poor)')

# JK
bp2 <- bpca(gabriel1971,
            meth='jk',
            d=1:3,
            var.rb=TRUE)

qbp2 <- qbpca(gabriel1971,
              bp2)

plot(qbp2,
     main='jk - 3d \n (very poor)')

# GH
bp3 <- bpca(gabriel1971,
            meth='gh',
            d=1:3,
            var.rb=TRUE)

qbp3 <- qbpca(gabriel1971,
              bp3)

plot(qbp3,
     main='gh - 3d \n (whow!)')

# HJ
bp4 <- bpca(gabriel1971,
            meth='hj',
            d=1:3,
            var.rb=TRUE)

qbp4 <- qbpca(gabriel1971,
              bp4)

plot(qbp4,
     main='hj - 3d \n (whow!)')

devAskNewPage(oask)

