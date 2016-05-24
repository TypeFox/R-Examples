GrammarIsTerminal <- function(x) {
  
  if (class(x) != "GEPhenotype") {
    stop("Invalid Phenotype Class")
  }

  return (x$type == "T")
}
