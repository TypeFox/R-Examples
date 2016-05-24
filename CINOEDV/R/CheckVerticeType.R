CheckVerticeType <-
function(Vertice){
  # cheak a node, or a vertice whether is a real vertice or a virtual one.
  #
  # input
  #    Vertice: a vertice in the network, namely, the graph
  # output
  #    Whether: Yes or No a real vertice
  #
  # Junliang Shang
  # 3.31/2014
  
  Whether <- TRUE
  Num <- nchar(Vertice)
  for (i in 1:Num){
    if (length(grep(":",substr(Vertice,i,i))==1)){
      Whether <- FALSE
      break
    }
  }
  
  list(Whether=Whether)
}
