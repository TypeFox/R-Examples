process.nominal.variable <-
function (labels, number.of.rows, parsed.xml, variable.index, variable.name) {
  aux <- list()
  aux[[1]] <- rep("$S", number.of.rows)
  if (labels)
    categories <- xpathSApply(parsed.xml, paste0("/assofile/variables/stvar[", variable.index, "]/nominal/nominal-desc/list-nom/label"), xmlValue)
  else
    categories <- xpathSApply(parsed.xml, paste0("/assofile/variables/stvar[", variable.index, "]/nominal/nominal-desc/list-nom/name"), xmlValue)
  
  aux[[2]] <- rep(length(categories), number.of.rows)
  nodes <- getNodeSet(parsed.xml, paste0("/assofile/indiv_mat/ligmat/valmat[", variable.index, "]"))
  
  after.evaluator <- function(node) {
    if (length(node["val_nomina"]) == 0)
      return (rep(NA, length(categories)))
    else {
      category <- as.numeric(xmlValue(node))
      return (append(rep(0, length(categories) - 1), 1, category - 1))
    }
  }
  node.categories <- t(xmlSApply(nodes, after.evaluator))
  
  aux <- data.frame(c(aux, as.data.frame(node.categories)))
  colnames(aux) <- c("$S", variable.name, categories)
  return (aux)
}
