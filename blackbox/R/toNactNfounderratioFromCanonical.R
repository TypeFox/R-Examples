toNactNfounderratioFromCanonical <- function(parvec) { ##utility to convert parvec c(twoNmu, ..., twoNfoundermu, twoNancmu) to c(NactNfounderratio, ..., twoNfoundermu, twoNancmu)
  if(class(parvec)=="data.frame" || class(parvec)=="list") parvec <- unlist(parvec) ## keep names
  return(c(NactNfounderratio=parvec[["twoNmu"]]/parvec[["twoNfoundermu"]]))
}
