## ----loading, results='hide', echo=FALSE, warning=FALSE, message=FALSE----
library(SOMbrero)

## ----dataGeneration------------------------------------------------------
set.seed(4031719)
the.data <- data.frame("x1"=runif(500), "x2"=runif(500))
plot(the.data, pch=19)

## ----dataTrain-----------------------------------------------------------
set.seed(593)
# run the SOM algorithm with 10 intermediate backups and 2000 iterations
my.som <- trainSOM(x.data=the.data, dimension=c(5,5), nb.save=10, maxit=2000, 
                   scaling="none", radius.type="letremy")

## ----energy--------------------------------------------------------------
plot(my.som, what="energy")

## ----hitmapObs-----------------------------------------------------------
plot(my.som, what="obs", type="hitmap")

## ----clusteredData, echo=FALSE, cache=TRUE-------------------------------
# prepare a vector of colors
my.colors <- rainbow(prod(my.som$parameters$the.grid$dim))[my.som$clustering]

# points depicted with the same color are in the same final cluster
plot(my.som$data[,1], my.som$data[,2], col=my.colors, pch=19, xlab="x1", 
     ylab="x2", main="Data according to final clustering")

## ----colorProto, fig.width=5, fig.height=2.5-----------------------------
par(mfrow=c(1,2))
plot(my.som, what="prototypes", type="color", var=1, main="prototypes - x1")
plot(my.som, what="prototypes", type="color", var=2, main="prototypes - x2")

## ----colorObs, fig.width=5, fig.height=2.5-------------------------------
par(mfrow=c(1,2))
plot(my.som, what="obs", type="color", var=1, main="obs mean values - x1")
plot(my.som, what="obs", type="color", var=2, main="obs mean values - x2")

## ----protoEvoluation, fig.width=15, fig.height=6, echo=FALSE-------------
# find out the prototypes to be linked
tra <- NULL
for (i in c(1,6,11,16,21)){
  tra <- c(tra, i, i+1, i+1, i+2, i+2, i+3, i+3, i+4)
}
for (i in c(1:5)){
  tra <- c(tra, i, i+5, i+5, i+10, i+10, i+15, i+15, i+20)
}
tmp <- matrix(tra, ncol=2, byrow=TRUE)
# plot the prototypes
par(mfrow=c(2, 5),mar=c(3,2,2,1))
invisible(sapply(1:my.som$parameters$nb.save, function(ind){
  plot(my.som$backup$prototypes[[ind]][,1], my.som$backup$prototypes[[ind]][,2],
       xlab="", ylab="", main=c("iteration ", my.som$backup$steps[ind]))
  for (i in 1:nrow(tmp)){
    segments(x0=my.som$backup$prototypes[[ind]][tmp[i,1],1], 
             y0=my.som$backup$prototypes[[ind]][tmp[i,1],2],
             x1=my.som$backup$prototypes[[ind]][tmp[i,2],1], 
             y1=my.som$backup$prototypes[[ind]][tmp[i,2],2], 
             col="red", pch=19)
  }
}))

## ----irisTrain, cache=TRUE-----------------------------------------------
set.seed(255)
# run the SOM algorithm with verbose set to TRUE
iris.som <- trainSOM(x.data=iris[,1:4], verbose=TRUE, nb.save=5)
iris.som

## ----energyIris----------------------------------------------------------
plot(iris.som, what="energy")

## ----irisClusters--------------------------------------------------------
iris.som$clustering
table(iris.som$clustering)

## ----irisHitmap----------------------------------------------------------
plot(iris.som, what="obs", type="hitmap")

## ----irisSummary---------------------------------------------------------
summary(iris.som)

## ----irisPred1-----------------------------------------------------------
# call predict.somRes
predict(iris.som, iris[1,1:4])
# check the result of the final clustering with the SOM algorithm
iris.som$clustering[1]

## ----irisGraphOP---------------------------------------------------------
par(mfrow=c(2,2))
plot(iris.som, what="obs", type="color", variable=1, print.title=TRUE, 
     main="Sepal length")
plot(iris.som, what="obs", type="color", variable=2, print.title=TRUE, 
     main="Sepal width")
plot(iris.som, what="obs", type="color", variable=3, print.title=TRUE, 
     main="Petal length")
plot(iris.som, what="obs", type="color", variable=4, print.title=TRUE, 
     main="Petal width")
plot(iris.som, what="prototypes", type="lines", print.title=TRUE)
plot(iris.som, what="obs", type="barplot", print.title=TRUE)
plot(iris.som, what="obs", type="radar", key.loc=c(-0.5,5), mar=c(0,10,2,0))

## ----irisObs-------------------------------------------------------------
plot(iris.som, what="obs", type="boxplot", print.title=TRUE)
rownames(iris)
plot(iris.som, what="obs", type="names", print.title=TRUE, scale=c(0.9,0.5))

## ----irisProto-----------------------------------------------------------
par(mfrow=c(2,2))
plot(iris.som, what="prototypes", type="3d", variable=1, main="Sepal length")
plot(iris.som, what="prototypes", type="3d", variable=2, main="Sepal width")
plot(iris.som, what="prototypes", type="3d", variable=3, main="Petal length")
plot(iris.som, what="prototypes", type="3d", variable=4, main="Petal width")

## ----irisDistProto-------------------------------------------------------
plot(iris.som, what="prototypes", type="poly.dist")
plot(iris.som, what="prototypes", type="umatrix")
plot(iris.som, what="prototypes", type="smooth.dist")
plot(iris.som, what="prototypes", type="mds")
plot(iris.som, what="prototypes", type="grid.dist")

## ----irisAdd1------------------------------------------------------------
class(iris$Species)
levels(iris$Species)
plot(iris.som, what="add", type="pie", variable=iris$Species)

## ----irisAdd2------------------------------------------------------------
plot(iris.som, what="add", type="color", variable=iris$Sepal.Length)

## ----irisMatCont, echo=FALSE---------------------------------------------
my.cont.mat <- matrix(data=c(rep(c(rep(1,50), rep(0,150)), 2), rep(1,50)), 
                      nrow=150, ncol=3)
colnames(my.cont.mat) <- levels(iris$Species)

## ----irisAdd4------------------------------------------------------------
# my.cont.mat is the contingency matrix corresponding to the variable 
# iris$Species - overview of the 5 first lines:
my.cont.mat[1:5,]
plot(iris.som, what="add", type="words", variable=my.cont.mat)

## ----irisAdd5------------------------------------------------------------
plot(iris.som, what="add", type="names", variable=rownames(iris),
     scale=c(0.9,0.5))

## ----irisAdd5bis---------------------------------------------------------
plot(iris.som, what="add", type="names", variable=iris$Species)

## ----irisQuality---------------------------------------------------------
quality(iris.som)

## ----irisSC--------------------------------------------------------------
plot(superClass(iris.som))

## ----irisSC3-------------------------------------------------------------
my.sc <- superClass(iris.som, k=3)
summary(my.sc)
plot(my.sc, plot.var=FALSE)

## ----irisSCplot, fig.width=6, fig.height=4-------------------------------
plot(my.sc, type="grid", plot.legend=TRUE)

## ----irisSCplot3d--------------------------------------------------------
plot(my.sc, type="dendro3d")

## ----irisSCplot2, fig.width=6, fig.height=4------------------------------
plot(my.sc, type="hitmap", plot.legend=TRUE)

## ----irisSCplot2B--------------------------------------------------------
plot(my.sc, type="lines", print.title=TRUE)
plot(my.sc, type="barplot", print.title=TRUE)
plot(my.sc, type="boxplot", print.title=TRUE)

## ----irisSCplot2C, fig.height=4, fig.width=6-----------------------------
plot(my.sc, type="mds", plot.legend=TRUE, cex =2)

## ----irisSCplot3---------------------------------------------------------
plot(my.sc, type="color")
plot(my.sc, type="poly.dist")
plot(my.sc, type="pie", variable=iris$Species)
plot(my.sc, type="radar", key.loc=c(-0.5,5), mar=c(0,10,2,0))

## ----irisSCplot4---------------------------------------------------------
plot(my.sc, type="color", add.type=TRUE, variable=iris$Sepal.Length)

