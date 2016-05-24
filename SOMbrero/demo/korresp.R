data(presidentielles2002)

# train the korresp algorithm
## details on the presidentielles2002 data set are provided with 'help(presidentielles2002)'
## default dimensions will be calculated by the algorithm (see 'help(trainSOM)')
presi.som <- trainSOM(x.data=presidentielles2002, type="korresp", scaling="chi2", dimension=c(8,8))

# prototypes overview can either be plotted for row variables
plot(presi.som, what="prototypes", type="lines", view="r")
# or for columns variables
plot(presi.som, what="prototypes", type="lines", view="c")

# the distances between prototypes can be displayed 
## either with smooth colors
plot(presi.som, what="prototypes", type="smooth.dist")
## or with a Multi Dimensional Scaling
plot(presi.som, what="prototypes", type="mds")

# observation distribution overviews are provided
## either with a plot
plot(presi.som, what="obs", type="hitmap")
## or with a table
table(presi.som$clustering)
## or more precisely by printing the variables names on the map
### NOTE: in the korresp SOM, both row and columns variables are considered
### WARNING: this graphic may produce some warnings when all row names can not fit on the plot
### this problem can be solved by enlarging your graphic device or by using the 'scale' argument (see 'help(wordcloud)')
plot(presi.som, what="obs", type="names", scale=c(0.6,0.7))

# hierarchical clustering is also available
## compute 3 super clusters
presi.sc <- superClass(presi.som, k=3)
## plot its 3 dimensional dendrogram
plot(presi.sc, type="dendro3d")
## view the prototypes positions in the Multi Dimensional Scaling plot
plot(presi.sc, type="mds", plot.legend=TRUE)

# two projection quality indicators are still available
quality(presi.som)
