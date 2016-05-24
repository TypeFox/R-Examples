tolatt2Ns2 <- function(parvec) { ##util. to convert parvec : c(2Nmu, 2Nm, g) to c(2Nmu, 2Ds2[lattice])
  #if(class(parvec)=="data.frame") {parvec <- as.numeric(parvec) #as.vector(dataframe) is numeric
  #} else if(class(parvec)=="list") parvec <- as.numeric(unlist(parvec))
  if(class(parvec)=="data.frame" || class(parvec)=="list") parvec <- unlist(parvec)
  ##tmp <- condaxialS2fromg(parvec["g"])*parvec["twoNm"]
  ##return(as.numeric(c(parvec["twoNmu"], tmp)))
  return(c(twoNmu=parvec[["twoNmu"]],
           latt2Ns2 = condaxialS2fromg(parvec[["g"]], D2bool=("2D" %in% blackbox.getOption("DemographicModel"))) * parvec[["twoNm"]] ))
}
