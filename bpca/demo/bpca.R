##
## Computing and ploting a bpca object with 'graphics' package - 2d
##

oask <- devAskNewPage(dev.interactive(orNone=TRUE))

bp <- bpca(gabriel1971)

plot(bp, 
     var.factor=2)

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
          meth='sqrt'),
     main='gabriel1971 - sqrt',
     sub='The graphical parameters are working fine!',
     var.factor=2,
     var.cex=.6,
     var.col=rainbow(9),
     var.pch='v',
     obj.pch='o',
     obj.cex=.5,
     obj.col=rainbow(8),
     obj.pos=1,
     obj.offset=.5)

##
## Computing and ploting a bpca object with 'scatterplot3d' package - 3d
##

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

##
## Computing and ploting a bpca object with 'obj.identify=TRUE' parameter - 2d
##

bp <- bpca(gabriel1971)

# Normal labels
if(dev.interactive()) {
  plot(bp,
       obj.names=FALSE,
       obj.identify=TRUE)
}  

# Alternative labels
if(dev.interactive()) {
  plot(bp,
       obj.names=FALSE,
       obj.labels=c('toi', 'kit', 'bat', 'ele', 'wat', 'rad', 'tv', 'ref'),
       obj.identify=TRUE)
}       

##
## Computing and ploting a bpca object with 'obj.identify=TRUE' parameter - 3d
##

bp <- bpca(gabriel1971,
           d=1:3)

# Normal labels
if(dev.interactive()) {
  plot(bp,
       obj.names=FALSE,
       obj.identify=TRUE)
}  

# Alternative labels
if(dev.interactive()) {
  plot(bp,
       obj.names=FALSE,
       obj.labels=c('toi', 'kit', 'bat', 'ele', 'wat', 'rad', 'tv', 'ref'),
       obj.identify=T)
}

##
## Computes: vector variable lengths, angles between vector variables and
## variable correlations from dataframe or matrix objects (n x p)
## n = rows (objects)
## p = columns (variables)
##

dt <- dt.tools(iris,
               var.pos=2) # No numeric columns are removed in 'dt.tools'

# Exploring the object 'bp' created by the function 'var.tools'
class(dt)
names(dt)
str(dt)

dt$length
dt$angle
dt$r
dt

# Checking the determinations
(iris.tools <- round(dt.tools(iris[-5],
                              center=2)$r,
                     5))

(iris.obsv  <- round(cor(iris[-5]),
                     5))

all(iris.tools == iris.obsv)

##
## Grouping objects with different symbols and colors - 2d and 3d
##

# 2d
plot(bpca(iris[-5]),
     var.factor=.3,
     var.cex=.7,
     obj.names=FALSE,
     obj.cex=1.5,
     obj.col=c('red', 'green3', 'blue')[unclass(iris$Species)],
     obj.pch=c('+', '*', '-')[unclass(iris$Species)])

# 3d static
plot(bpca(iris[-5],
          d=1:3),
     var.factor=.2,
     var.color=c('blue', 'red'),
     var.cex=1,
     obj.names=FALSE,
     obj.cex=1,
     obj.col=c('red', 'green3', 'blue')[unclass(iris$Species)],
     obj.pch=c('+', '*', '-')[unclass(iris$Species)])

##
## Example of 'var.rb=TRUE' parameter as a measure of the quality of the biplot - 2d
##

## Differences between methods of factorization
# SQRT
bp1 <- bpca(gabriel1971,
            meth='sqrt',
            var.rb=TRUE)

qbp1 <- qbpca(gabriel1971,
              bp1)

plot(qbp1, main='sqrt - 2d \n (poor)')

# JK
bp2 <- bpca(gabriel1971,
            meth='jk',
            var.rb=TRUE)

qbp2 <- qbpca(gabriel1971, bp2)

plot(qbp2,
     main='jk - 2d \n (very poor)')

# GH
bp3 <- bpca(gabriel1971,
            meth='gh',
            var.rb=TRUE)

qbp3 <- qbpca(gabriel1971,
              bp3)

plot(qbp3,
     main='gh - 2d \n (good)')

# HJ
bp4 <- bpca(gabriel1971,
            meth='hj',
            var.rb=TRUE)

qbp4 <- qbpca(gabriel1971,
              bp4)

plot(qbp4,
     main='hj - 2d \n (good)')

##
## Example of 'var.rb=TRUE' parameter as a measure of the quality of the biplot - 3d
##

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

##
## Example of 'var.rd=TRUE' parameter as a measure of the quality of the biplot - 2d
## Mainly recommended for large datasets.
##

bp <- bpca(gabriel1971,
           meth='hj',
           var.rb=TRUE, 
           var.rd=TRUE, 
           limit=3)

bp$var.rd

# RUR followed by CRISTIAN contains information in dimensions that
# wasn't contemplated by the biplot reduction (PC3).
# Between all, RUR followed by CRISTIAN, variables are bad represented by a 2d
# biplot.

# Graphical visualization of the importance of the variables not contemplated
# in the reduction
plot(bpca(gabriel1971,
          meth='hj',
          d=3:4),
     main='hj',
     xlim=c(-1,1),
     ylim=c(-1,1),
     zlim=c(-1,1))

##
## New options ploting
##
data(ontario)

plot(bpca(ontario))

## Labels for all objects
(obj.lab <- paste('g',
                  1:18,
                  sep=''))

# Giving obj.labels
plot(bpca(ontario),
    obj.labels=obj.lab) 

# Evaluate an object (1 is the default)
plot(bpca(ontario),
     type='eo',
     obj.cex=1)

plot(bpca(ontario),
     type='eo',
     obj.id=7,
     obj.cex=1)

# Giving obj.labels
plot(bpca(ontario),
     type='eo',
     obj.labels=obj.lab,
     obj.id=7,
     obj.cex=1)

# The same as above
plot(bpca(ontario),
     type='eo',
     obj.labels=obj.lab,
     obj.id='g7',
     obj.cex=1)

# Evaluate a variable (1 is the default)
plot(bpca(ontario),
     type='ev',
     var.pos=2,
     var.cex=1)

plot(bpca(ontario),
     type='ev',
     var.id='E7',
     obj.labels=obj.lab,
     var.pos=1,
     var.cex=1)

# A complete plot
cl <- 1:3

plot(bpca(iris[-5]),
     type='ev',
     var.id=1,
     var.fac=.3,
     obj.names=FALSE,
     obj.col=cl[unclass(iris$Species)])

legend('topleft',
       legend=levels(iris$Species),
       text.col=cl,
       pch=19,
       col=cl,
       cex=.9,
       box.lty=0)   

# Compare two objects (1 and 2 are the default)
plot(bpca(ontario),
     type='co')

plot(bpca(ontario),
     type='co',
     obj.labels=obj.lab)

plot(bpca(ontario),
     type='co',
     obj.labels=obj.lab,
     obj.id=13:14)

plot(bpca(ontario),
     type='co',
     obj.labels=obj.lab,
     obj.id=c('g7', 'g13'))

# Compare two variables
plot(bpca(ontario),
     type='cv')

# Which won where/what
plot(bpca(ontario),
     type='ww')

# Discrimitiveness vs. representativeness
plot(bpca(ontario),
     type='dv')

# Means vs. stability
plot(bpca(ontario),
     type='ms')

# Rank objects with ref. to the ideal variable 
plot(bpca(ontario),
     type='ro')

# Rank variables with ref. to the ideal object
plot(bpca(ontario),
     type='rv')

plot(bpca(iris[-5]),
     type='eo',
     obj.id=42,
     obj.cex=1)

plot(bpca(iris[-5]),
     type='ev',
     var.id='Sepal.Width')

plot(bpca(iris[-5]),
     type='ev',
     var.id='Sepal.Width',
     var.factor=.3)

devAskNewPage(oask)

