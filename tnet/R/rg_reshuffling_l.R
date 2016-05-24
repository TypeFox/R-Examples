`rg_reshuffling_l` <-
function(net,keep.i=FALSE,keep.j=FALSE,seed=NULL){
  if(keep.i & keep.j)
    stop("If you keep both the i and j column, you are not randomising anything!")
  if(is.null(attributes(net)$tnet))
    net <- as.tnet(net, type="longitudinal tnet")
  if(attributes(net)$tnet!="longitudinal tnet")
    stop("Network not loaded properly")
  #Take out time column
  timecol <- net[,"t"]
  net <- as.matrix(net[,c("i","j","w")])
  #Check that there are no weaking ties
  if(sum(net[net[,"i"]!=net[,"j"],"w"] == rep(1, sum(net[,"i"]!=net[,"j"])))!=sum(net[,"i"]!=net[,"j"]))
    stop('This function cannot deal with negative ties');
  #If seed is set, set it
  if(!is.null(seed))
    set.seed(as.integer(seed))
  #Randomise ties
  N <- max(c(net[,1],net[,2]))
  active.nodes <- rep(FALSE, N)
  if(!keep.i)
    rdm.function <- "At <- which(active.nodes); At <- At[!At %in% j]; net[t,1] <- sample(At, size=1)"
  if(!keep.j)
    rdm.function <- "At <- which(active.nodes); At <- At[!At %in% i]; net[t,2] <- sample(At, size=1)"
  if(!keep.i & !keep.j)
    rdm.function <- "net[t,1:2] <- sample(which(active.nodes), size=2)"
  for(t in 1:nrow(net)) {
    tie <- as.vector(net[t,])
    i <- tie[1]
    j <- tie[2]
    if(i != j) {
      eval(parse(text=rdm.function))
    } else {
      w <- as.logical(tie[3])
      active.nodes[i] <- w
    }
  }  
  net <- data.frame(t=timecol, net, stringsAsFactors=FALSE)
  dimnames(net)[[2]] <- c("t","i","j","w");
  return(net)
}