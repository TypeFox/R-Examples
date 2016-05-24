safiModel <- function(s.d, y) {
  if (class(s.d) != "safidesign") 
    stop("object class is not safidesign")
  VP <- s.d$DoE
  if (length(y) != dim(VP)[1]) 
    stop("number of runs  in DoE and y differ")
  if (s.d$method != "SB" && s.d$method != "other") 
    stop("method not supported")
  split.points <- s.d$split.points
  sequential.bifurcation.design <- s.d$method == "SB"
  mirrored <- s.d$mirrored.runs.included
  d.f <- s.d$d.f
  n <- dim(VP)[2]
  # mirrored <- dim(VP)[1] > n+1
  if (sequential.bifurcation.design) {
    coef <- numeric(n)
    for (i in 1:n) {
      r2 <- rep(-1, n)
      if (i - 1 > 0) 
        r2[1:(i - 1)] <- 1
      r1 <- r2
      r1[i] <- 1
      # find rows
      row1 <- which(apply(mapply(data.frame(VP), r1, FUN = "=="), MARGIN = 1, FUN = all))
      row2 <- which(apply(mapply(data.frame(VP), r2, FUN = "=="), MARGIN = 1, FUN = all))
      if (mirrored) {
        r3 <- -1 * r2
        r4 <- -1 * r1
        row3 <- which(apply(mapply(data.frame(VP), r3, FUN = "=="), MARGIN = 1, FUN = all))
        row4 <- which(apply(mapply(data.frame(VP), r4, FUN = "=="), MARGIN = 1, FUN = all))
        coef[i] <- (y[row1] - y[row2] + y[row3] - y[row4])/4
      } else {
        coef[i] <- (y[row1] - y[row2])/2
      }
    }
  } else {
    coef <- lm(y ~ 0 + VP)$coefficients
    names(coef) <- NULL
  }
  len <- cumsum(sapply(split.points, length) - 1)
  norm <- sapply(split.points, diff, simplify = FALSE)
  # list of coefficient vectors corresponding to the vectors in split.points
  normalized.coefficients <- list(coef[1:len[1]]/norm[[1]])
  if (d.f > 1) {
    for (i in 2:d.f) {
      normalized.coefficients[[i]] <- coef[(len[i - 1] + 1):(len[i])]/norm[[i]]
    }
  }
  names(normalized.coefficients) <- paste("normalized.coefficients.", 1:length(split.points), sep = "")
  s.m <- list(normalized.coefficients = normalized.coefficients, DoE = s.d$DoE, split.points = s.d$split.points, 
              d.f = s.d$d.f, variable.names = s.d$variable.names, mirrored.runs.included = s.d$mirrored.runs.included, 
              method = s.d$method)
  class(s.m) <- "safimodel"
  return(s.m)
} 