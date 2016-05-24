parse.graphml.internal <-
function(file, nodes, use.names) {
  doc <- xmlParse(file, useInternalNodes = TRUE)
  top <- xmlRoot(doc)
  graph <- top[["graph"]]
  if (!use.names) {
    if (xmlSize(which(xmlSApply(graph, xmlName) == "node")) != length(nodes)) stop("Incorrect number of node names")
  }
  ns <- c(ns = "http://graphml.graphdrawing.org/xmlns")
  if (use.names) {
    node.data <- getNodeSet(doc, "//ns:data[contains(@key,'d6')]/*", ns)
    for (i in 1:length(node.data)) nodes[i] <- xmlValue(node.data[[i]]["NodeLabel"]$NodeLabel[1]$text)
  }
  all <- getNodeSet(doc, "//*[position() > 1]", ns)
  keep <- getNodeSet(doc, "//ns:edge | //ns:node | //ns:graph | //ns:key[@attr.name = 'description'] | //ns:data[contains(@key,'d9')]", ns)
  remove <- getNodeSet(doc, "//ns:key[@for='port'] | //ns:key[@for='graphml'] | //ns:data[@key='d4'] | //ns:data[@key='d6'] 
 | //ns:data[@key='d7'] | //ns:data[@key='d8'] | //ns:data[@key='d10']" , ns)
  removeNodes(union(setdiff(all, keep), remove))
  temp <- tempfile(fileext = ".graphml")
  temp.xml <- saveXML(doc, file = temp)
  free(doc)
  igrph <- read.graph(temp.xml, format = "graphml")
  igrph <- set.vertex.attribute(igrph, "name", value = nodes)
  return(igrph)
}
