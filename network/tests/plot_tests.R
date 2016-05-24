# various tests for network plotting functions
# mostly recent functionality added by skyebend
require(network)
require(testthat)
# -----  test edge labels ------
ymat<-matrix(c(0,1,2,3, 0,0,0,0, 1,0,0,0, 0,0,0,0),ncol=4)
ynet<-network(ymat,ignore.eval=FALSE,names.eval='weight')
# don't do anything if no value given
plot(ynet,edge.label.col='blue',edge.label.cex='weight')
# use edge ids is if edge.label=TRUE
plot(ynet,edge.label=TRUE)

plot(ynet,edge.label='weight',edge.label.col='blue',edge.label.cex='weight')

# labels for curved edges
plot(ynet,edge.label='weight',edge.label.col='blue',edge.label.cex='weight',usecurve=TRUE)
plot(ynet,edge.label='weight',edge.label.col='blue',edge.label.cex='weight',usecurve=TRUE,edge.curve=0.5)

data(emon)
par(mar=c(0,0,0,0))
plot(emon[[5]],edge.label=TRUE,edge.label.cex=0.6,edge.col='gray',edge.lwd=(emon[[5]]%e%'Frequency')*2)

# test for labeling network with no edges #521
plot(network.initialize(1),edge.label=TRUE)

# test color stuff

col.list<-c('red','#800000','#80000505',NA)
# test is.color for vector NA processing bug #491
if(!all(is.color(col.list)[1:3] & is.na(is.color(col.list)[4]))){
  stop('is.color did not correctly recognize colors and NA values in a character vector')
}
   
col.list<-list('red','#800000','#80000505',NA)
# test is.color for list NA processing bug #491
if(!all(is.color(col.list)[1:3] & is.na(is.color(col.list)[4]))){
  stop('is.color did not correctly recognize colors and NA values in a list')
}

# ------------ as.color --------

expect_equal(as.color(c('a','b','c')),1:3)  # character
expect_equal(as.color(1:3),1:3)  # numeric
expect_equal(as.color(as.factor(c('a','b','c'))),1:3)  # factor
expect_equal(as.color(c('red','green','blue')),c('red','green','blue'))  # color name
expect_equal(as.color(c(1,0.5,0)),c("#FFFFFF", "#808080", "#000000"))# real valued  (gray)
# transparency/ opacity
expect_equal(as.color(c('red','green','blue'),0.5),c("#FF000080", "#00FF0080", "#0000FF80"))
expect_equal(as.color(1:3,0.5),c("#00000080", "#FF000080", "#00CD0080"))
expect_error(as.color(c('red','green','blue'),1.5),regexp = 'opacity parameter must be a numeric value in the range 0 to 1')


# ----- plot fixes ----

plot(network.initialize(5),vertex.lwd=c(1,2,3,5,10))

# test for expansion of label attribute name bug #785
# this should produce a plot with vertices labeled A to E, instead
# used to plot single vertex is labeled with "Label'
test<-network.initialize(5)
set.vertex.attribute(test,'Label',LETTERS[1:5])
plot(test,label='Label')

# replicates non-matching label name
plot(test,label='A')
plot(test,label=1)

# should error if all values are missing
#set.vertex.attribute(test,'bad',NA,v=1:3)
#plot(test,label='bad')

# tests for #673 plot.network.default gives error when rendering labels if two connected vertices have the same position
test<-network.initialize(2)
test[1,2]<-1
plot(test,coord=cbind(c(1,1),c(1,1)),jitter=FALSE,displaylabels=TRUE)

test<-network.initialize(3)
test[1,2:3]<-1
plot(test,coord=cbind(c(1,1,2),c(1,1,2)),jitter=FALSE,displaylabels=TRUE)

# tests for polygon sizes/sides
plot(network.initialize(7),vertex.sides=c(50,4,3,2,1,0,NA),vertex.cex=40,coord=matrix(0,ncol=7,nrow=7),jitter=F,vertex.col='#CCCCCC00',vertex.border =c('red','green','blue','orange'))
plot(network.initialize(7),vertex.sides=c(50,4,3,2,1,0,NA),vertex.cex=0)
plot(network.initialize(7),vertex.sides=c(50,4,3,2,1,0,NA),vertex.cex=NA)

