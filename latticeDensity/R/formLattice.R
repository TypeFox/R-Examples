formLattice <-
function(nodeFillingOutput){
   require(spdep)
   if(class(nodeFillingOutput)!="nodeFillingOutput"){
       stop("Should be the output from the functions nodeFilling or removeHole")}
   nodes <- nodeFillingOutput$nodes
   poly <- nodeFillingOutput$poly
   node.spacing <- nodeFillingOutput$node.spacing
   latt <- dnearneigh(nodes,node.spacing*0.5,node.spacing*1.5)
   formLatticeOutput <- list(EW.locs = nodeFillingOutput$EW.locs,
                              NS.locs = nodeFillingOutput$NS.locs,
                              nodes = nodes,
                              poly = poly,
                              latt = latt,
                              hole.list = nodeFillingOutput$hole.list)
   class(formLatticeOutput) <- "formLatticeOutput"
   return(formLatticeOutput)
}



