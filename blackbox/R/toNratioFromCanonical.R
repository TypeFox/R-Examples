toNratioFromCanonical <- function(parvec) { ##utility to convert parvec c(twoNmu, ..., twoNancmu) to c(Nratio, ..., twoNancmu)
  if(class(parvec)=="data.frame" || class(parvec)=="list") parvec <- unlist(parvec) ## keep names
  return(c(Nratio=parvec[["twoNmu"]]/parvec[["twoNancmu"]]))
}
