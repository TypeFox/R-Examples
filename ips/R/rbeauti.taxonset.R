rbeauti.taxonset <- function(x, id){
  
  taxon <- lapply(x$taxon, function(x) xmlNode("taxon", 
                                               attrs = c(id = x,
                                                         spec = "Taxon")))
  taxonset <- xmlNode("taxonset",
                      attrs = c(id = x$id,
                                spec = "TaxonSet"),
                      .children = taxon)
  
  LogNormal <- xmlNode("LogNormal", attrs = c(id = paste("LogNormalDistributionModel.", x$id, sep = ""),
                                              name = "distr"),
                       .children = list(xmlNode("parameter", 1,
                                                attrs = c(estimate = "false",
                                                          id = paste("RealParameter.M", x$id, sep = ""),
                                                          name = "M")),
                                        xmlNode("parameter", 1.25,
                                                attrs = c(estimate = "false",
                                                          id = paste("RealParameter.S", x$id, sep = ""),
                                                          lower = "0.0",
                                                          name = "S",
                                                          upper = "5.0"))
                                        ))
  
  xmlNode("distribution", taxonset, LogNormal,
          attrs = c(id = paste(x$id, ".prior", sep = ""),
                    spec = "beast.math.distributions.MRCAPrior",
                    tree = paste("@Tree.t:", id, sep = "")))
}