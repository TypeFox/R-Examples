toNfounderNancratioFromCanonical <- function(parvec) { ##utility to convert parvec c(twoNmu, ..., twoNfoundermu, twoNancmu) to c(NfounderNancratio, ..., twoNfoundermu, twoNancmu)
  if(class(parvec)=="data.frame" || class(parvec)=="list") parvec <- unlist(parvec) ## keep names
  return(c(NfounderNancratio=parvec[["twoNfoundermu"]]/parvec[["twoNancmu"]]))
}
