test_aop_backdoor <- function(){
  library(graph)
  library(rjson)
  steatosis_json_file <- system.file("extdata", "steatosis_aop_json.cyjs", package = "aop")
  steatosis_aop <- convert_cytoscape_to_aop(steatosis_json_file)
  steatosis_aop_graph <- convert_aop_to_graph(steatosis_aop)
  checkEquals("389", aop_backdoor(steatosis_aop_graph, "391", "388")[1])
  checkEquals("392", aop_backdoor(steatosis_aop_graph, "391", "388")[2])
}
