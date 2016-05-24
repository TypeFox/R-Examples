netconnection <- function(treat1, treat2,
                          studlab,
                          data = NULL, subset = NULL,
                          title = "",
                          warn = FALSE) {
  
  if (is.null(data)) data <- sys.frame(sys.parent())
  ##
  ## Catch treat1, treat2, studlab from data:
  ##
  mf <- match.call()
  mf$data <- mf$subset <- mf$title <- mf$warn <- NULL
  mf[[1]] <- as.name("data.frame")
  mf <- eval(mf, data)
  ##
  ## Catch subset (possibly) from data:
  ##
  mf2 <- match.call()
  mf2$treat1 <- mf2$treat2 <- NULL
  mf2$studlab <- NULL
  mf2$data <- mf2$title <- mf2$warn <- NULL
  mf2[[1]] <- as.name("data.frame")
  ##
  mf2 <- eval(mf2, data)
  ##
  if (!is.null(mf2$subset))
    if ((is.logical(mf2$subset) & (sum(mf2$subset) > length(mf$treat1))) ||
        (length(mf2$subset) > length(mf$treat1)))
      stop("Length of subset is larger than number of comparisons.")
    else
      mf <- mf[mf2$subset,]
  
  
  treat1 <- mf$treat1
  treat2 <- mf$treat2
  ##
  if (is.factor(treat1))
    treat1 <- as.character(treat1)
  if (is.factor(treat2))
    treat2 <- as.character(treat2)
  ##
  if (length(mf$studlab) != 0)
    studlab <- as.character(mf$studlab)
  else{
    if (warn)
      warning("No information given for argument 'studlab'. Assuming that comparisons are from independent studies.")
    studlab <- seq(along = treat1)
  }
  
  
  if (any(treat1 == treat2))
    stop("Treatments must be different (arguments 'treat1' and 'treat2').")
  
  
  ##
  ## Check for correct number of comparisons
  ##
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)
      abs(x - round(x)) < tol
  tabnarms <- table(studlab)
  sel.narms <- !is.wholenumber((1 + sqrt(8 * tabnarms + 1)) / 2)
  ##
  if (sum(sel.narms) == 1)
    stop(paste("Study '", names(tabnarms)[sel.narms],
               "' has a wrong number of comparisons.",
               "\n  Please provide data for all treatment comparisons (two-arm: 1; three-arm: 3; four-arm: 6, ...).",
               sep = ""))
  if (sum(sel.narms) > 1)
    stop(paste("The following studies have a wrong number of comparisons: ",
               paste(paste("'", names(tabnarms)[sel.narms], "'", sep = ""),
                     collapse = ", "),
               "\n  Please provide data for all treatment comparisons (two-arm: 1; three-arm: 3; four-arm: 6, ...).",
               sep = ""))
  
  
  treats <- as.factor(c(as.character(treat1), as.character(treat2)))
  n <- length(unique(treats)) # Number of treatments
  m <- length(treat1)         # Number of comparisons
  Edges <- data.frame(treat1 = treats[1:m],
                      treat2 = treats[(m + 1):(2 * m)])
  B <- matrix(0, nrow = m, ncol = n) # Edge-vertex incidence matrix
  ##
  for (e in 1:m) {
    B[e, Edges[e,1]] <-  1
    B[e, Edges[e,2]] <- -1
  }
  ##
  rownames(B) <- studlab
  colnames(B) <- levels(treats)
  ##
  L.mult <- t(B) %*% B             # Laplacian matrix with multiplicity
  A <- diag(diag(L.mult)) - L.mult # Adjacency matrix
  D <- netdistance(A)              # Distance matrix
  L <- diag(rowSums(A)) - A        # Laplacian matrix without multiplicity
  ##
  n.subsets <- as.integer(table(round(eigen(L)$values,10)==0)[2])
  ##
  if (n.subsets > 1) {
    ##
    ## Block diagonal matrix in case of sub-networks
    ##
    maxdist <- dim(D)[1]
    ##
    D2 <- D
    D2[is.infinite(D2)] <- maxdist
    order.D <- hclust(dist(D2))$order
    ##
    D <- D[order.D, order.D]
    A <- A[order.D, order.D]
    L <- L[order.D, order.D]
  }
  ##
  res <- list(treat1 = treat1,
              treat2 = treat2,
              studlab = studlab,
              ##
              k = length(unique(studlab)),
              m = m,
              n = n,
              n.subnets = n.subsets,
              ##
              D.matrix = D,
              A.matrix = A,
              L.matrix = L,
              ##
              title = title,
              ##
              warn = warn,
              call = match.call(),
              version = packageDescription("netmeta")$Version
              )

  class(res) <- "netconnection"

  res
}
