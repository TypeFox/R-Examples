process.mult.nominal.modif.variable <-
function (labels, number.of.rows, parsed.xml, variable.index, variable.name) {
  aux <- list()
  aux[[1]] <- rep("$H", number.of.rows)
  if (labels)
    categories <- xpathSApply(parsed.xml, paste0("/assofile/variables/stvar[", variable.index, "]/mult_nominal_Modif/nominal-desc/list-nom/label"), xmlValue)
  else
    categories <- xpathSApply(parsed.xml, paste0("/assofile/variables/stvar[", variable.index, "]/mult_nominal_Modif/nominal-desc/list-nom/name"), xmlValue)
  aux[[2]] <- rep(length(categories), number.of.rows)
  
  nodes <- getNodeSet(parsed.xml, paste0("/assofile/indiv_mat/ligmat/valmat[", variable.index, "]"))
  
  get.distributions <- function(node) {
    if (length(node["val_list_modal"]) == 0)
      return (rep(NA, length(categories)))
    else {  
      moda.nodes <- as.numeric(sapply(xmlSApply(node, function(x) x["no_moda"]), xmlValue))
      frequencies <- as.numeric(sapply(xmlSApply(node, function(x) x["frequency"]), xmlValue))
      missing.categories.indexes <- setdiff(1:length(categories), moda.nodes)
      for (missing.cat.index in missing.categories.indexes) {
        frequencies <- append(frequencies, 0, after = missing.cat.index - 1)
      }                  
      return (frequencies)
    }
  }
  
  all.frequencies <- t(round(sapply(nodes, get.distributions), 3))
  aux <- data.frame(c(aux, as.data.frame(all.frequencies)))
  
  colnames(aux) <- c("$H", variable.name, categories)
  return (aux)
}
