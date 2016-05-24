selectdesign = function (design, MT, limit, starts=read.table(file.choose(new=FALSE))) 
{
  if (design == "CRD") {
    N <- c(rep("A", MT/2), rep("B", MT/2))
    design <- sample(N, MT, replace = FALSE)
    return(design)
  }
  if (design == "RBD") {
    ab <- c("A", "B")
    design <- numeric()
    repeat {
      design <- c(design, sample(ab, 2, replace = FALSE))
      if (length(design) == MT) 
        break
    }
    return(design)
  }
  if (design == "ATD") {
    N <- c(rep(0, MT/2), rep(1, MT/2))
    repeat {
      design <- sample(N, MT, replace = FALSE)
      check <- numeric()
      for (itr in 1:(MT - limit)) {
        check2 <- 0
        for (it in itr:(itr + limit)) {
          check2 <- check2 + design[it]
        }
        check <- cbind(check, check2)
      }
      if (sum(check == (limit + 1) | check == 0) == 0) {
        for (it in 1:(length(design))) {
          if (design[it] == 0) {
            design[it] <- "A"
          }
          else {
            design[it] <- "B"
          }
        }
        break
      }
    }
    return(design)
  }
  if (design == "AB") {
    quantity <- choose(MT - 2 * limit + 1, 1)
    design <- matrix("A", 1, MT)
    index.b <- (limit + 1):(MT - (limit - 1))
    design[1, index.b[sample(1:quantity, 1)]:MT] <- "B"
    return(design)
  }
  if (design == "ABA") {
    quantity <- choose(MT - 3 * limit + 2, 2)
    design <- matrix("A", 1, MT)
    selection <- sample(1:quantity, 1)
    index1 <- 1:(MT - 3 * limit + 1)
    index2 <- rev(index1)
    index.b.1 <- numeric()
    for (it in 1:length(index1)) {
      index.b.1 <- c(index.b.1, rep(index1[it], index2[it]))
    }
    index.b.2 <- numeric()
    for (itr in index1) {
      for (it in itr:(MT - 3 * limit + 1)) {
        index.b.2 <- c(index.b.2, 2 * limit - 1 + it)
      }
    }
    design[1, (limit + index.b.1[selection[1]]):(index.b.2[selection[1]])] <- "B"
    return(design)
  }
  if (design == "ABAB") {
    quantity <- choose(MT - 4 * limit + 3, 3)
    design <- matrix("A", 1, MT)
    selection <- sample(1:quantity, 1)
    index1 <- 1:(MT - 4 * limit + 1)
    index2 <- rev(cumsum(index1))
    index.b1.1 <- numeric()
    for (it in 1:length(index1)) {
      index.b1.1 <- c(index.b1.1, (rep((limit + index1[it]), 
                                       (index2[it]))))
    }
    index.b1.2 <- numeric()
    for (itr in index1) {
      for (it in (itr - 1):(MT - 4 * limit)) {
        index.b1.2 <- c(index.b1.2, rep((2 * limit + 
                                           it), (MT - 4 * limit + 1 - it)))
      }
    }
    design[1, (index.b1.1[selection[1]]:index.b1.2[selection[1]])] <- "B"
    indexb2 <- numeric()
    for (it in 1:length(index1)) {
      indexb2 <- c(indexb2, index1[it:length(index1)])
    }
    index.b2 <- numeric()
    for (it in 1:length(indexb2)) {
      index.b2 <- c(index.b2, indexb2[it]:length(index1))
    }
    design[1, (4 * limit - limit + index.b2[selection[1]]):MT] <- "B"
    return(design)
  }
  if (design == "MBD") {
    points <- starts
    N <- nrow(points)
    startpoints <- starts
    limits <- list()
    for (it in 1:N) {
      limits[[it]] <- startpoints[it,]
    }
    number <- numeric(N)
    for (it in 1:N) {
      number[it] <- length(limits[[it]][[1]])
    }
    startpt <- numeric(N)
    for (it in 1:N) {
      if (number[it] != 1) {
        startpt[it] <- sample(limits[[it]][[1]], 1)
      }
      else {
        startpt[it] <- limits[[it]][[1]]
      }
    }
    design <- sample(startpt, replace = FALSE)
    return(design)
  }
}
