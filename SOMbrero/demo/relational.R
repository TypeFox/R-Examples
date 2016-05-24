data(lesmis)

# plot the lesmis graph
plot(lesmis, vertex.size=0)

# train the relational algortihm
## details on the lesmis data set are provided with 'help(lesmis)'
## default dimensions will be calculated by the algorithm (see 'help(trainSOM)')
lesmis.som <- trainSOM(x.data=dissim.lesmis, type="relational")

# overview of the prototypes values
plot(lesmis.som, what="prototypes", type="barplot")

# overview of the observation distribution 
## either using the row names of the input data set as an additional variable
plot(lesmis.som, what="obs", type="names", scale=c(0.8,0.9))
## or with a table
table(lesmis.som$clustering)

# perform a hierarchical clustering with 4 super clusters
lesmis.sc <- superClass(lesmis.som, k=4)
## identify the super clusters on the map with colors
plot(lesmis.sc, type="grid")
## identify the super clusters on the prototype values barplot plot
plot(lesmis.sc, type="barplot", print.title=TRUE)

# projection quality indicators provide
## either the topographic error
quality(lesmis.som, quality.type="topographic")
## or the quantization error
quality(lesmis.som, quality.type="quantization")