removeHole <-
function(hole.poly, nodeFillingOutput) {

  #
  #  This function finds all nodes contained in hole.poly
  #  and deletes them.  It is up to the user to insure that
  #  this polygon is really contained in the larger polynomial
  #  and that it is non-intersecting.  This function can be
  #  repeated for each hole.
  #
  require(splancs)
  require(spdep)
  hole.poly <- as.matrix(hole.poly)
  if(class(nodeFillingOutput)!="nodeFillingOutput"){
    stop("Should be the output from the function nodeFilling")}
  nodes <- nodeFillingOutput$nodes
  poly <- nodeFillingOutput$poly
  nodes <- nodes[!inout(pts = nodes, poly = hole.poly),]
  nodeFillingOutput$nodes <- nodes
  return(nodeFillingOutput)
 
 }
