# train the algorithm
## details on the iris data set are provided with 'help(iris)'
## default type: numeric
## default dimensions will be calculated by the algorithm (see 'help(trainSOM)')
## default number of iterations will also be calculated by the algorithm
## default number of intermediate savings is 0 (computational time increases when some are recorded)
iris.som <- trainSOM(x.data=iris[,1:4])

# plot the prototype values for the first variable
## either in color (2 equivalent ways to do so)
plot(iris.som, what="prototypes", type="color", variable="Sepal.Length")
plot(iris.som, what="prototypes", type="color", variable=1)
## or in 3 dimensions
plot(iris.som, what="prototypes", type="3d", variable=1)

# see the distances between prototypes
## distance between the vertices of the polygon and the grid borders 
## are proportional to the distance between the neighbouring prototypes
plot(iris.som, what="prototypes", type="poly.dist", print.title=TRUE)

# The distribution of observations on the map can be seen 
## either with a graphic
plot(iris.som, what="obs", type="hitmap")
## or with a table
table(iris.som$clustering)
## or more precisely by printing the names of the rows in their assigned neuron
### WARNING: this graphic may produce some warnings when all row names can not fit on the plot
### they can be solved by using the 'scale' argument (see 'help(wordcloud)')
plot(iris.som, what="obs", type="names")

# observation values for each neuron can be plotted with a radar plot
## the argument 'key.loc' is used to add a legend and specify its location (see 'help(stars)')
plot(iris.som, what="prototypes", type="radar", key.loc=c(-1,2), mar=c(0,10,2,0))

# the flower species distribution can be plotted by passing the 'Species' variable as an additional variable
plot(iris.som, what="add", type="pie", variable=iris$Species)

# perform a hierarchical clustering
## with 3 super clusters
iris.sc <- superClass(iris.som, k=3)
summary(iris.sc)
## plot its dendrogram
plot(iris.sc)
## identify the super clusters on the observation values on the radar plot
plot(iris.sc, type="radar", key.loc=c(-1,2), mar=c(0,10,2,0))
## plot the flower species distribution with the super cluster labels
plot(iris.sc, type="pie", variable=iris$Species)

# compute the projection quality indicators
quality(iris.som)