## ----eval=FALSE----------------------------------------------------------
#  library(aop)
#  steatosis_aop <- convert_cytoscape_to_aop("path_to_file\\aop.cyjs")

## ------------------------------------------------------------------------
library(aop)
library(graph)
steatosis_json_file <- system.file("extdata", "steatosis_aop_json.cyjs", package = "aop")
steatosis_aop <- convert_cytoscape_to_aop(steatosis_json_file)

## ------------------------------------------------------------------------
steatosis_graph <- convert_aop_to_graph(steatosis_aop)
plot(steatosis_graph)

## ------------------------------------------------------------------------
getAOPNodeName(steatosis_aop, "387")

## ------------------------------------------------------------------------
aop_backdoor(steatosis_graph, "390", "388")

## ------------------------------------------------------------------------
getAOPNodeName(steatosis_aop, "389")
getAOPNodeName(steatosis_aop, "392")

