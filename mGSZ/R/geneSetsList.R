geneSetsList <-
function(data){
  data <- GSA.read.gmt(data)
  geneSets <- data$genesets
  names(geneSets) <- data$geneset.names
  return(geneSets)
}
