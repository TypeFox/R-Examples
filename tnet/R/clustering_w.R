`clustering_w` <- 
function (net, measure = "am") {
  if (is.null(attributes(net)$tnet)) 
    net <- as.tnet(net, type = "weighted one-mode tnet")
  if (attributes(net)$tnet != "weighted one-mode tnet") 
    stop("Network not loaded properly")
  if (length(measure) < 1) 
    stop("You must specify a measure")
  for (i in 1:length(measure))
    if(!(measure[i] %in% c("bi", "am", "gm", "ma", "mi")))
      stop("You must specify a correct measure")

  # Sort edgelist, and name columns
  net <- net[order(net[, "i"], net[, "j"]), ]
  dimnames(net)[[2]] <- c("i", "j", "wij")
    
  # Get second part of 2-paths
  tmp <- net
  dimnames(tmp)[[2]] <- c("j", "k", "wjk")
  tmp <- merge(net, tmp)

  # Exclude reciprocated ties
  index <- tmp[, "i"] != tmp[, "k"]
  tmp <- tmp[index, ]
  tmp <- tmp[, c("i", "k", "wij", "wjk")]

  # Find closed ties
  dimnames(net)[[2]] <- c("i", "k", "wik")
  net[, "wik"] <- TRUE
  tmp <- merge(tmp, net, all.x = TRUE)
  tmp[is.na(tmp[, "wik"]), "wik"] <- FALSE

  # Calculate coefficients
  index <- as.logical(tmp[, "wik"])
  tmp <- tmp[, c("wij", "wjk")]
  output <- rep(0, length(measure))
  for (i in 1:length(measure)) {
    if (measure[i] == "am") {
      tvalues <- rowSums(tmp)
    } else if (measure[i] == "gm") {
      tvalues <- sqrt(tmp[, "wij"] * tmp[, "wjk"])
    } else if (measure[i] == "ma") {
      tvalues <- pmax.int(tmp[, "wij"], tmp[, "wjk"])
    } else if (measure[i] == "mi") {
      tvalues <- pmin.int(tmp[, "wij"], tmp[, "wjk"])
    } else if (measure[i] == "bi") {
      tvalues <- rep(1, nrow(tmp))
    } else {
      stop("measure incorrectly specified")
    }
    output[i] <- sum(tvalues[index])
    output[i] <- output[i]/sum(tvalues)
  }

  # Name and return
  names(output) <- measure
  return(output)
}
