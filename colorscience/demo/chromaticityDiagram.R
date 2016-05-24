par(mfrow=c(1,3),mar=c(1,1,1,1),oma=c(1, 0, 0,0))
xl<-yl<-c(0,.8)

# CIE 1931 xy Chromaticity Diagram
chromaticity.diagram.color.fill(xlim=xl,ylim=yl)
par(new=TRUE)
Maxwell.triangle(xlim=xl,ylim=yl)
par(new=TRUE)
chromaticity.diagram (xlim=xl,ylim=yl,main='CIE 1931 xy')

# CIE 1976 u'v' Chromaticity Diagram
chromaticity.diagram.color.fill(conversionFunction=CIE1931XYZ2CIE1976uv, xlim=xl,ylim=yl,xlab="u'",ylab="v'")
par(new=TRUE)
Maxwell.triangle(conversionFunction=CIE1931XYZ2CIE1976uv, xlim=xl,ylim=yl,xlab="u'",ylab="v'")
par(new=TRUE)
chromaticity.diagram (conversionFunction=CIE1931XYZ2CIE1976uv, xlim=xl,ylim=yl,xlab="u'",ylab="v'",main="CIE 1976 u'v'")

# CIE 1960 uv Chromaticity Diagram
chromaticity.diagram.color.fill(conversionFunction=CIE1931XYZ2CIE1960uv, xlim=xl,ylim=yl,xlab="u",ylab="v")
par(new=TRUE)
Maxwell.triangle(conversionFunction=CIE1931XYZ2CIE1960uv, xlim=xl,ylim=yl,xlab="u",ylab="v")
par(new=TRUE)
chromaticity.diagram (conversionFunction=CIE1931XYZ2CIE1960uv, xlim=xl,ylim=yl,xlab="u",ylab="v",main='CIE 1960 uv')
