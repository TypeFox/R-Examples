if(!exists("Sepal.length")) {
  data(iris)
  attach(iris)
  }

require(lattice)
trellis.device("pdf", color=FALSE, file = "Examples/irisBwplot.pdf", width = 5,
               height = 4)

u = cut(Petal.Length, breaks = 3)
bwplot( ~ Sepal.Length | u * Species)
dev.off()
