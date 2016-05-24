process.inter.cont.variable <-
function (number.of.rows, parsed.xml, variable.index, variable.name) {
  aux <- list()
  aux[[1]] <- rep("$I", number.of.rows)
  
  after.evaluator <- function(node, element.to.retrieve) {
    if (length(node["val_interv"]) == 0)
      return (NA)
    else
      return (as.numeric(xmlValue(xmlElementsByTagName(node[[1]], element.to.retrieve)[[1]])))
  }
  
  nodes <- getNodeSet(parsed.xml, paste0("/assofile/indiv_mat/ligmat/valmat[", variable.index, "]"))
  aux[[2]] <- sapply(nodes, after.evaluator, element.to.retrieve = "pmin")
  aux[[3]] <- sapply(nodes, after.evaluator, element.to.retrieve = "pmax")
  
  aux <- data.frame(aux)
  colnames(aux) <- c("$I", variable.name, variable.name)
  return (aux)
}
