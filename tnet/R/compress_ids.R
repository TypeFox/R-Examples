`compress_ids` <-
function(net, type=NULL){
  # Ensure type is correct
  if (is.null(attributes(net)$tnet)) {
    if(is.null(type)) {
      net <- as.tnet(net)
    } else {
      net <- as.tnet(net, type=type)
    }
  }
  type <- attributes(net)$tnet
  if(substr(type, nchar(type)-4, 100)!= " tnet") 
    stop("Network not loaded properly")

  # Define output element; [[1]] is network; [[2/3]] is node list
  out <- list()
  
  # If one-mode network
  if(attributes(net)$tnet == "weighted one-mode tnet") {
    net <- net[order(net[,"i"], net[,"j"]),]
    # Create node translation table
    out[[2]] <- unique(c(net[,"i"],net[,"j"]))
    out[[2]] <- out[[2]][order(out[[2]])]
    out[[2]] <- data.frame(old=out[[2]], new=1:length(out[[2]]))
    # Check if compression is needed
    if(sum(out[[2]][,1]==out[[2]][,2])!=nrow(out[[2]])) {
      # Update i column with new ids
      tmp <- out[[2]]
      dimnames(tmp)[[2]] <- c("i", "new")
      net <- merge(net, tmp)
      net <- net[,c("new","j","w")]
      dimnames(net)[[2]] <- c("i", "j", "w")
      # Update j column with new ids
      tmp <- out[[2]]
      dimnames(tmp)[[2]] <- c("j", "new")
      net <- merge(net, tmp)
      net <- net[,c("i","new","w")]
      dimnames(net)[[2]] <- c("i", "j", "w")
      row.names(net) <- NULL
    }
    # Put updated network in out[[1]]
    out[[1]] <- net
  
  # If binary two-mode network
  } else if(attributes(net)$tnet == "binary two-mode tnet") {
    net <- net[order(net[,"i"], net[,"p"]),]
    # Create node translation tables [[2]] and [[3]]    
    out[[2]] <- unique(net[,"i"])
    out[[2]] <- out[[2]][order(out[[2]])]
    out[[2]] <- data.frame(old=out[[2]], new=1:length(out[[2]]))
    out[[3]] <- unique(net[,"p"])
    out[[3]] <- out[[3]][order(out[[3]])]
    out[[3]] <- data.frame(old=out[[3]], new=1:length(out[[3]]))
    # Check if compression is needed
    if(sum(out[[2]][,1]==out[[2]][,2])!=nrow(out[[2]]) | sum(out[[3]][,1]==out[[3]][,2])!=nrow(out[[3]])) {
      # Update i column with new ids
      tmp <- out[[2]]
      dimnames(tmp)[[2]] <- c("i", "new")
      net <- merge(net, tmp)
      net <- net[,c("new","p")]
      dimnames(net)[[2]] <- c("i", "p")
      # Update p column with new ids
      tmp <- out[[3]]
      dimnames(tmp)[[2]] <- c("p", "new")
      net <- merge(net, tmp)
      net <- net[,c("i","new")]
      dimnames(net)[[2]] <- c("i", "p")
      row.names(net) <- NULL
    }
    # Put updated network in out[[1]]
    out[[1]] <- net
    
    # If weighted two-mode network
  }  else if(attributes(net)$tnet == "weighted two-mode tnet") {
    net <- net[order(net[,"i"], net[,"p"]),]
    # Create node translation tables [[2]] and [[3]]
    out[[2]] <- unique(net[,"i"])
    out[[2]] <- out[[2]][order(out[[2]])]
    out[[2]] <- data.frame(old=out[[2]], new=1:length(out[[2]]))
    out[[3]] <- unique(net[,"p"])
    out[[3]] <- out[[3]][order(out[[3]])]
    out[[3]] <- data.frame(old=out[[3]], new=1:length(out[[3]]))
    # Check if compression is needed
    if(sum(out[[2]][,1]==out[[2]][,2])!=nrow(out[[2]]) | sum(out[[3]][,1]==out[[3]][,2])!=nrow(out[[3]])) {
      # Update i column with new ids
      tmp <- out[[2]]
      dimnames(tmp)[[2]] <- c("i", "new")
      net <- merge(net, tmp)
      net <- net[,c("new","p","w")]
      dimnames(net)[[2]] <- c("i", "p", "w")
      # Update p column with new ids
      tmp <- out[[3]]
      dimnames(tmp)[[2]] <- c("p", "new")
      net <- merge(net, tmp)
      net <- net[,c("i","new","w")]
      dimnames(net)[[2]] <- c("i", "p", "w")
      row.names(net) <- NULL
    }
    # Put updated network in out[[1]]
    out[[1]] <- net
     
    # If longitudinal network
  }  else if(attributes(net)$tnet == "longitudinal tnet") {
    net <- net[order(net[,"i"], net[,"j"]),]
    # Create node translation table
    out[[2]] <- unique(c(net[,"i"],net[,"j"]))
    out[[2]] <- out[[2]][order(out[[2]])]
    out[[2]] <- data.frame(old=out[[2]], new=1:length(out[[2]]))
    # Check if compression is needed
    if(sum(out[[2]][,1]==out[[2]][,2])!=nrow(out[[2]])) {
      # Update i column with new ids
      tmp <- out[[2]]
      dimnames(tmp)[[2]] <- c("i", "new")
      net <- merge(net, tmp)
      net <- net[,c("t","new","j","w")]
      dimnames(net)[[2]] <- c("t","i", "j", "w")
      # Update j column with new ids
      tmp <- out[[2]]
      dimnames(tmp)[[2]] <- c("j", "new")
      net <- merge(net, tmp)
      net <- net[,c("t","i","new","w")]
      dimnames(net)[[2]] <- c("t","i", "j", "w")
      # Order network again
      net <- net[order(net[,"t"]),]
      row.names(net) <- NULL
    }
    # Put updated network in out[[1]]
    out[[1]] <- net
  } else {
    stop("type parameter misspecified")
  }
  return(out)
}
