permindex <- function(data, block=NULL, group=NULL, nsim=4999) {

  if (any(!is.null(block), !is.null(group))) data <- data[do.call(order, data[, c(block, group), drop=FALSE]),]
  rownames(data) <- NULL
  n <- nrow(data)
  if (all(is.null(block), is.null(group))) {
  	permIndex <- replicate(nsim, sample.int(n))
  }else if (all(!is.null(block), is.null(group))) {
    Bsplit <- split(seq_len(n), as.character(data[, block]))    # split 1:n seq into nB blocks
    permIndex <- replicate(nsim, unlist(lapply(Bsplit, function (x) x[sample.int(length(x))]), recursive = FALSE, use.names = FALSE))
  }else if (all(is.null(block), !is.null(group))){   # data=data; group="group"; nsim=5
    permIndex <- replicate(nsim, {Gsplit <- split(seq_len(n), as.character(data[, group]))   # split 1:n seq into separate groups
      nGsplit <- lapply(Gsplit, function (x) x[sample.int(length(x))])   # within each group permute the whole seq
      unlist(nGsplit[sample.int(length(nGsplit))], recursive = FALSE, use.names = FALSE)})          # permute the order of the whole groups and repeat process nsim times
  }else{
    ind <- seq_len(n)
    names(ind) <- as.character(data[, group])
    Bsplit <- split(ind, as.character(data[, block]))    # split 1:n seq into blocks

    permIndex <- replicate(nsim, unlist(lapply(Bsplit, function (x) {    # within each block
        Gsplit <- split(x, names(x))  # split each block seq into groups
        nGsplit <- lapply(Gsplit, function (x) x[sample.int(length(x))])   # within each group permute the whole seq
        unlist(nGsplit[sample.int(length(nGsplit))], recursive = FALSE, use.names = FALSE)          # permute the order of groups
      }), recursive = FALSE, use.names = FALSE))
  }

  colnIndex <- apply(permIndex, 2, function(x) paste(x, collapse=""))
  colnames(permIndex) <- colnIndex
  permIndex <- permIndex[, unique(colnIndex)]
  colnames(permIndex) <- NULL
  return(permIndex)
}
