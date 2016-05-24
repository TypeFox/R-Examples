p <- splom(~iris[1:4],data=iris, groups=Species,
    pch=c(21,23,25), 
    cex=0.8,
    alpha=0.7,
    col='gray50',
    fill=trellis.par.get('superpose.symbol')$col
    )
print(p)
