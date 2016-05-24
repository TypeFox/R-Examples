process.continue.variable <-
function (number.of.rows, parsed.xml, variable.index, variable.name) {
  aux <- list()
  aux[[1]] <- rep("$C", number.of.rows)
  
  after.evaluator <- function(node) {
    if (length(node["val_conti"]) == 0)
      return (NA)
    else
      return (as.numeric(xmlValue(node[[1]])))
  }
  
  aux[[2]] <- xpathSApply(parsed.xml, paste0("/assofile/indiv_mat/ligmat/valmat[", variable.index, "]"), after.evaluator)
  
  aux <- data.frame(aux)
  colnames(aux) <- c("$C", variable.name)
  return (aux)
}
