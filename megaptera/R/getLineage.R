# extract lineage from XML
# called by: ncbiTaxonomy
# package: megaptera
# author: Christoph Heibl (at gmx.net)
# last update: 2014-10-28


getLineage <- function(xml, spec, kingdom){
  # this xPath expression avoids the 'double kingdom error' CH 2014-10-23
  xPath <- paste("/ncbi:TaxaSet/Taxon[ScientificName='", spec, 
                 "']/LineageEx[Taxon[Rank='kingdom']/ScientificName='", 
                 kingdom,"']/Taxon/", sep = "")
  out<- data.frame(name = xpathSApply(xml, paste(xPath,  "ScientificName", sep = ""), xmlValue, 
                                      namespaces =  xmlNamespaceDefinitions(xmlRoot(xml), simplify = TRUE)), 
                   rank = xpathSApply(xml, paste(xPath,  "Rank", sep = ""), xmlValue, 
                                      namespaces =  xmlNamespaceDefinitions(xmlRoot(xml), simplify = TRUE)),
                   stringsAsFactors = FALSE)
  if ( nrow(out) == 0 ){
    warning("kingdom-ambiguity-issue in ", spec)
    return(NULL)
  } else {
    return(rbind(out, c(spec, "species")))
  }
}

# c("ncbi" = "http://www.w3.org/XML/1998/namespace")
# matchNamespaces(xml, c(xml = "http://www.omegahat.org"))