SODAS.to.RSDA <-
function(XMLPath, labels = T) {
  parsed.xml <- xmlInternalTreeParse(XMLPath)
  
  containsNode <- getNodeSet(parsed.xml, "/assofile/contains")
  if (length(containsNode) == 0)
    stop("No 'contains' tag is present in the XML file")
  containsNode <- containsNode[[1]]
  if (xmlGetAttr(containsNode, "INDIVIDUALS") != "YES" || xmlGetAttr(containsNode, "VARIABLES") != "YES" || xmlGetAttr(containsNode, "RECTANGLE_MATRIX") != "YES")
    stop("Insufficient data in XML file")
  
  if (labels) {
    sym.obj.names <- xpathSApply(parsed.xml, "/assofile/individus/stindiv/label", xmlValue)
    variables.names <- xpathSApply(parsed.xml, "/assofile/variables/stvar/ident/label", xmlValue)
  }
  else {
    sym.obj.names <- xpathSApply(parsed.xml, "/assofile/individus/stindiv/name", xmlValue)
    variables.names <- xpathSApply(parsed.xml, "/assofile/variables/stvar/ident/name", xmlValue)
  }
  
  variables.types <- xpathSApply(parsed.xml, "/assofile/variables/stvar/*[2]", xmlName)
  result <- data.frame(row.names = sym.obj.names)
  number.of.rows <- nrow(result)
  
  for (i in 1:length(variables.types)) {
    
    cat(paste0("Processing variable ", i, ": ", variables.names[[i]],"\n"))
    
    switch (variables.types[[i]],
            'inter-cont' = {
              result <- cbind(result, process.inter.cont.variable(number.of.rows, parsed.xml, i, variables.names[[i]]))
            },
            'continue' = {
              result <- cbind(result, process.continue.variable(number.of.rows, parsed.xml, i, variables.names[[i]]))
            },
            'nominal' = {
              result <- cbind(result, process.nominal.variable(labels, number.of.rows, parsed.xml, i, variables.names[[i]]))
            },
            'mult_nominal' = {
              result <- cbind(result, process.mult.nominal.variable(labels, number.of.rows, parsed.xml, i, variables.names[[i]]))
            },
            'mult_nominal_Modif' = {
              type.modif <- xpathSApply(parsed.xml, paste0("/assofile/variables/stvar[", i, "]/mult_nominal_Modif/type_modif"), xmlValue)
              if (type.modif != "proba")
                cat(paste0("Unsupported type.modif in mult_nominal_Modif variable: "), type.modif, "\n")
              else
                result <- cbind(result, process.mult.nominal.modif.variable(labels, number.of.rows, parsed.xml, i, variables.names[[i]]))
            },
            cat(paste0("Variable type not supported:"), variables.types[[i]], "\n"))
  }
  
  return (newSobject(result))
}
