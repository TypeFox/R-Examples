`reinforcement_tm` <- function(net){
  # Ensure that the network conforms to the tnet standard
  if(is.null(attributes(net)$tnet)) 
    net <- as.tnet(net, type="binary two-mode tnet")
  if(attributes(net)$tnet!="binary two-mode tnet")
    stop("Network not loaded properly")
  # 1-paths i1+p1
  paths <- net
  dimnames(paths)[[2]] <- c("i1","p1")
  # 2-paths i1+p1+i2
  dimnames(net)[[2]] <- c("i2","p1")
  paths <- merge(paths, net, sort=FALSE)
  paths <- paths[paths[,"i1"] != paths[,"i2"],]
  # 3-paths i1+p1+i2+p2
  dimnames(net)[[2]] <- c("i2","p2")
  paths <- merge(paths, net, sort=FALSE)
  paths <- paths[paths[,"p1"] != paths[,"p2"],]
  denominator <- nrow(paths)
  # Find which 3-paths are part of 4-cycles
  paths <- paths[order(paths[,"i1"], paths[,"p1"], paths[,"p2"]),c("i1","p1","i2","p2")]
  dimnames(net)[[2]] <- c("i1","p2")
  paths <- merge(paths, net, sort=FALSE)
  numerator <- nrow(paths)
  rm(paths)
  # Fraction
  return(numerator/denominator)
}
