process.mult.nominal.variable <-
function (labels, number.of.rows, parsed.xml, variable.index, variable.name) {
  aux <- list()
  aux[[1]] <- rep("$S", number.of.rows)
  if (labels)
    categories <- xpathSApply(parsed.xml, paste0("/assofile/variables/stvar[", variable.index, "]/mult_nominal/nominal-desc/list-nom/label"), xmlValue)
  else
    categories <- xpathSApply(parsed.xml, paste0("/assofile/variables/stvar[", variable.index, "]/mult_nominal/nominal-desc/list-nom/name"), xmlValue)
  
  aux[[2]] <- rep(length(categories), number.of.rows)
  nodes <- getNodeSet(parsed.xml, paste0("/assofile/indiv_mat/ligmat/valmat[", variable.index, "]"))
  
  after.evaluator <- function(node) {
    if (length(node["val_modal"]) == 0)
      return (NA)
    else {
      present.mods <- as.numeric(xmlSApply(node, xmlValue))
      modals.vector <- rep(0, length(categories) - length(present.mods))
      for (present.mod in present.mods) {
        modals.vector <- append(modals.vector, 1, present.mod - 1)
      }
      return (modals.vector)
    }
  }
  node.categories <- t(xmlSApply(nodes, after.evaluator)) 
  
  aux <- data.frame(c(aux, as.data.frame(node.categories)))
  colnames(aux) <- c("$S", variable.name, categories)
  return (aux)
}
