.nodeDataVector <- function(g, att) {
  ##  data(sysdata, envir=environment())
  if (class(g)[1] != "graphNEL")
    stop("'g' has to be a 'graphNEL' object")
  stopifnot(.validateGraph(g))
  
  convert <- identity
  if (att == "atom->symbol") {
    att <- "atom"
    convert <- function(value) {
      if (is.character(value)) value
      else as.character(.chemicalElements()[value, "symbol"])
    }
  }
  else if (att == "atom->number") {
    att <- "atom"
    convert <- function(value) {
      if (is.numeric(value)) value
      ##else chemicalElements[chemicalElements$symbol == value, "number"]
      else .chemicalElements()[.chemicalElements()$symbol == value, "number"]
    }
  }

  if (is.null(nodeDataDefaults(g)[[att]]))
    stop(paste("node attribute", att, "not set"))

  raw <- nodeData(g)
  sapply(names(raw), function(v) convert(raw[[v]][[att]]))
}
