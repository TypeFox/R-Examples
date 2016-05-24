createAndPlot <-
function(data, samples=1, level="genus", colors, divisions, main, 
sub, showTipLabel=TRUE, showNodeLabel=FALSE, displayLegend=TRUE, onePerPage=FALSE){

if(missing(data))
stop("A valid data set is required.")

trees <- createTrees(data, samples, level)

if(onePerPage){
par(layout(1))
}else{ #4 per page
par(layout(matrix(c(1,3,2,4), 2, 2)))
}
plotTree(trees, colors, divisions, main, sub, 
showTipLabel, showNodeLabel, displayLegend)

par(layout(1)) #reset layout back
}
