#
# load the library
#
library(network)
#
# attach the sociomatrix for the Florentine marriage data
# This is not yet a graph object. 
#
data(flo)
#
# print out the sociomatrix for the Florentine marriage data
#
flo
#
# Create a network object out of the adjacency matrix and print it out
#
nflo <- network(flo,directed=FALSE)
nflo
#
# print out the sociomatrix for the Florentine marriage data
#
print(nflo,matrix.type="adjacency") 
#
# plot the Florentine marriage data
#
plot(nflo)
#
# create a vector indicating the Medici family and add it as a covariate to the
# graph object.
#
nflo <- set.vertex.attribute(nflo,"medici",c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0))
nflo
#
# create a vector indicating the Medici family for the graph
#
medici <- rep("",nrow(flo))
names(medici) <- dimnames(flo)[[1]]
medici[names(medici)=="Medici"] <- "Medici"
#
# plot the marriage data, highlighting the Medici family
#
plot(nflo,vertex.col=1+get.vertex.attribute(nflo,"medici"))

# plot the emon St. Helens network, with edge widths proportional
# to 'Frequency', and edges labeled by their id
data(emon)
par(mar=c(0,0,0,0))
plot(emon[[5]],edge.label=TRUE,edge.label.cex=0.6,
     edge.col='gray',edge.lwd=(emon[[5]]%e%'Frequency')*2)
