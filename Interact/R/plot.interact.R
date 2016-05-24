plot.interact <-
function(x, numInteractions = nrow(x$interaction.ordered), ...){
  plot(1:numInteractions, (x$interaction.ordered)$qval[1:numInteractions], ylab = "FDR estimate", xlab = "Number of Interactions called significant", type = 'b',...)
}
