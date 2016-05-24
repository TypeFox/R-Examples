SavePathway2File<- function(network, layout=layout.random, name="network", file="Graph")
{   
   Split<-unlist(strsplit(as.character(network[,2]),split="|",fixed=T))
   network<-graph.data.frame(as.data.frame(cbind(first=Split[(1:(length(Split)/2))*2-1],sencond=Split[(1:(length(Split)/2))*2],network[,c(4,6)]),stringsAsFactors=FALSE),directed=F)
   .XGMML.destription <- function(name)
{
  # create top node
  # top part of xml
  top <- xmlNode("graph", attrs = c(label=name, "xmlns:dc"="http://purl.org/dc/elements/1.1/", "xmlns:xlink"="http://www.w3.org/1999/xlink", "xmlns:rdf"="http://www.w3.org/1999/02/22-rdf-syntax-ns#", "xmlns:cy"="http://www.cytoscape.org", xmlns="http://www.cs.rpi.edu/XGMML"))
  top <- append.xmlNode(top, xmlNode("att", attrs=c(name="documentVersion", value="1.1")))
  
  d <- xmlNode("rdf:Description", attrs=c("rdf:about"="http://www.cytoscape.org/"))
  d <- append.xmlNode(d,  xmlNode("dc:type", "Protein-Protein Interaction"))
  d <- append.xmlNode(d,  xmlNode("dc:description", "N/A"))
  d <- append.xmlNode(d,  xmlNode("dc:identifier", "N/A"))
  d <- append.xmlNode(d,  xmlNode("dc:date", Sys.time()))
  d <- append.xmlNode(d,  xmlNode("dc:title", name))
  d <- append.xmlNode(d,  xmlNode("dc:format", "BioNet-Cytoscape-XGMML"))
  
  c <- xmlNode("att", attrs=c(name="networkMetadata"), xmlNode("rdf:RDF", d))
  top <- append.xmlNode(top, c)
  return(top)
}
    .XGMML.nodes <- function(network)
{
  # create node-nodes
  c.node <- rep("node", length(V(network)))
  nodes <- lapply(c.node, xmlNode)
  
  # create node attributes
  a<-layout(network)[,1]
  b<-layout(network)[,2]
  V(network)$layout.x<-500*a
  V(network)$layout.y<-500*b
  attrib <- list.vertex.attributes(network)
  node.attribs <- matrix(data=NA, nrow=2, ncol=length(V(network)))
  node.attribs[1,] =paste("att type=", "\"string\"", " name=", "\"node.label\""," value=", "\"", get.vertex.attribute(network, attrib[1]), "\"", sep="")
  node.attribs[2,] =paste("graphics type=", "\"rectangle\"", " x=", "\"", get.vertex.attribute(network, attrib[2]), "\"", " y=", "\"", get.vertex.attribute(network, attrib[3]), "\"",  " fill=", "\"#C1FFC1\"", " node.border.paint=", "\"dimgray\"", " node.size=",  "\"8\"", sep="")
  node.attribs <- matrix(lapply(node.attribs, xmlNode), nrow = 2, ncol = length(V(network)))
  if(is.null(V(network)$name))
  {
    V(network)$name <- as.character(V(network))
  }
  node.label <- V(network)$name
  node.id <- as.vector(V(network))
  
  # append node attributes
  for(i in 1:length(V(network)))
  {
    nodes[[i]] <- addAttributes(nodes[[i]], label = node.label[i], id=node.id[i])
    nodes[[i]] <- append.xmlNode(nodes[[i]], node.attribs[,i])
  }
  
  return(nodes)
}

# internal method for the addition of edges to XGMML
.XGMML.edges <- function(network)
{
  # create edge-nodes
  c.edge <- rep("edge", length(E(network)))
  edges <- lapply(c.edge, xmlNode)
  
  edgelist.names <- get.edgelist(network, names=TRUE)
  edgelist.names <- paste(edgelist.names[,1], edgelist.names[,2], sep=" (pp) ")
  edgelist.ids <- get.edgelist(network, names=FALSE)
   
  # create edge attributes
  
  E(network)$color<-ifelse(E(network)$CORE_ENRICHMENT=='YES', 'red', 'dimgray');
  attrib <- list.edge.attributes(network)
  edge.attribs <- matrix(data=NA, nrow=1, ncol=length(E(network)))
  if(is.directed(network)){direct=6}else{direct=0}
  
  edge.attribs[1,] =paste("graphics width=", "\"1\"", " fill=", "\"", get.edge.attribute(network, attrib[3]), "\"", " cy:sourceArrow=", "\"0\"", " cy:targetArrow=", "\"",direct,"\""," cy:sourceArrowColor=","\"#666666\"", " cy:targetArrowColor=","\"#666666\""," cy:edgeLabelFont=","\"SanSerif-0-10\"", " cy:edgeLineType=", "\"SOLID\"", " cy:curved=","\"STRAIGHT_LINES\"",sep="")
  edge.attribs <- matrix(lapply(edge.attribs, xmlNode), nrow=1, ncol=length(E(network)))
  
  # append edge attributes
  for(i in 1:length(E(network)))
  {
    edges[[i]] <- addAttributes(edges[[i]], label=edgelist.names[i], source=edgelist.ids[i,1], target=edgelist.ids[i,2])
    edges[[i]] <- append.xmlNode(edges[[i]], c(edge.attribs[,i]))
  }
  return(edges)
}

   
	addNode <- XML::addNode
    if(is(network, "graphNEL"))
    {
      network <- igraph.from.graphNEL(network)  
    }
    top <- .XGMML.destription(name=name)
    nodes <- .XGMML.nodes(network=network)
    top <- append.xmlNode(top, nodes) 
    edges <- .XGMML.edges(network=network)
    top <- append.xmlNode(top, edges) 
    print("...writing to file")
    saveXML(top, file=paste(file, ".xgmml", sep=""), encoding="UTF-8")
	addNode <- graph::addNode
  }
  