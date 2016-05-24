nodeFilling <- function(poly, node.spacing,hole.list = NULL){
   require(splancs)
   require(spdep)
   poly <- as.matrix(poly)
   if(!is.null(hole.list)){
     number.holes <- length(hole.list)
     for(k in 1:number.holes){
       hole.list[[k]] <- as.matrix(hole.list[[k]])}
     print("This function does not check to see if the holes")
     print("are nonintersecting, or whether they are contained")
     print("inside the boundary")
     flush.console()
   }
   node.spacing <- as.numeric(node.spacing)
   # 
   width.EW <- max(poly[,1]) - min(poly[,1])
   width.NS <- max(poly[,2]) - min(poly[,2])
   #
   #  Fill the region with a grid of nodes
   #
   EW.locs <- seq(min(poly[,1]),max(poly[,1]),by=node.spacing)
   NS.locs <- seq(min(poly[,2]),max(poly[,2]),by=node.spacing) 
   bound.array <- expand.grid(EW.locs,NS.locs)
   names(bound.array) <- c("x","y")
   nodes <- bound.array[inout(pts=bound.array,poly=poly),]
   nodes <- as.points(as.matrix(nodes))
   #
   #
   #
   output <- list(EW.locs = EW.locs,
                 NS.locs = NS.locs,
                 nodes = nodes,
                 poly = poly,
                 node.spacing = node.spacing,
                 hole.list = hole.list)
   class(output) <- "nodeFillingOutput"
   if(!is.null(hole.list)){for (k in 1:number.holes){
     output <- removeHole(hole.list[[k]], output)}}
   return(output)
}


