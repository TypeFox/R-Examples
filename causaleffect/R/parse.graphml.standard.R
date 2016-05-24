parse.graphml.standard <-
function(file, nodes, use.names) {
  doc <- xmlParse(file, useInternalNodes = TRUE)
  top <- xmlRoot(doc)
  graph <- top[["graph"]]
  if (!use.names) {
    if (xmlSize(which(xmlSApply(graph, xmlName) == "node")) != length(nodes)) stop("Incorrect number of node names")
  }
  ns <- c(ns = "http://graphml.graphdrawing.org/xmlns")
  edges <- which(xmlSApply(graph, xmlName) == "edge")
  remove.id <- 0
  removals <- list()
  for (i in 1:length(edges)) {
    current.edge <- graph[[edges[i]]]
    edge.attr <- xmlAttrs(current.edge)
    datas <- which(xmlSApply(current.edge, xmlName) == "data")
    for (j in 1:length(datas)) {
      current.data <- current.edge[[datas[j]]]
      if (xmlAttrs(current.data) == "d10") {
        arc <- current.data[[1]]
        arw.attr <- xmlAttrs(arc[["Arrows"]])
        src <- edge.attr[["source"]]
        trgt <- edge.attr[["target"]]
        two.arrows <- !identical(arw.attr[["source"]], "none") & !identical(arw.attr[["target"]], "none")
        no.arrows <- identical(arw.attr[["source"]], "none") & identical(arw.attr[["target"]], "none")
        if (two.arrows | no.arrows) {
          remove.id <- remove.id + 1
          removals[[remove.id]] <- current.edge
          e1 <- newXMLNode("edge", parent = doc, attrs = c(id = "e", source = src, target = trgt), 
                  newXMLNode("data", attrs = c(key = "d9"), cdata = TRUE, "U" ))
          e2 <- newXMLNode("edge", parent = doc, attrs = c(id = "e", source = trgt, target = src), 
                  newXMLNode("data", attrs = c(key = "d9"), cdata = TRUE, "U" ))
          addChildren(graph, kids = list(e1, e2))
        }
      }     
    }
  }
  removeNodes(removals)
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
