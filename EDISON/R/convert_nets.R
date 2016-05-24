#' Convert internal representation of networks.
#' 
#' Converts from representing the network as a list of target nodes to
#' representing it as a list of segments.
#' 
#' 
#' @param Ball Input network: List of target nodes, where each element is a
#' NumSegs by NumNodes matrix giving the parents for the target node in each
#' segment.
#' @param Eall Changepoints: List of target nodes, where each element contains
#' a vector of changepoints.
#' @return List with elements: \item{B_nets}{List of segments, where each
#' element contains a matrix of size NumNodes by NumNodes, representing the
#' network for that segment.} \item{segs}{Vector containing the global segment
#' boundaries.}
#' @author Frank Dondelinger
#' @export convert_nets
convert_nets <-
function(Ball, Eall) {
  
  # Find 'global' segments (set of all target-specific segments)
  segs = c()
  
  for(target in 1:length(Ball)) {
    segs = c(segs, Eall[[target]]);
  }
  
  segs = sort(unique(segs));
  
  # Initialise segment list
  B_nets = list()

  mapping = matrix(0, length(Eall), length(segs) - 1)
  
  # Build up segment list
  local.seg = rep(0, length(Ball))
  empty.matrix = matrix(0, dim(Ball[[1]])[2], length(Ball))
  
  for(segment in 1:(length(segs) - 1)) {
      
    seg.net = empty.matrix
      
    for(target in 1:length(Ball)) {
      local.seg.temp = local.seg[target]

      if(segs[segment] == Eall[[target]][local.seg.temp+1]) {
        local.seg.temp = local.seg.temp + 1
        local.seg[target] = local.seg.temp
      } 
        
      seg.net[,target] = Ball[[target]][local.seg.temp,] 
      mapping[target, segment] = local.seg.temp
    }
      
    B_nets[[segment]] = seg.net
  }

  return(list(B_nets=B_nets, seg=segs, mapping=mapping))
}

