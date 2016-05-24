createSafiDesign <- function(method = "SB", mirrored.runs.included = FALSE, d.f = 1, variable.names = NULL) {
  if (!is.null(variable.names) && length(variable.names) != d.f) 
    stop("number of variable names and number of variables differ")
  if (method == "SB") {
    DoE <- matrix(0, ncol = d.f, nrow = 2)
    DoE[1, ] <- 1
    DoE[2, ] <- -1
    if (d.f > 1) {
      runs <- t(sapply((1 + 1):d.f, function(x) {
        z <- rep(1, d.f)
        z[x:d.f] <- -1
        return(z)
      }))
      DoE <- rbind(DoE, runs)
      if (mirrored.runs.included) {
        DoE <- rbind(DoE, (-1) * runs)
      }
    }
  } else if (method == "other") 
    DoE <- matrix(ncol = d.f, nrow = 0) else stop("method not supported")
  if (is.null(variable.names)) {
    variable.names <- paste("x", 1:d.f, sep = "")
  }
  colnames(DoE) <- variable.names
  split.points <- replicate(d.f, c(0, 1), simplify = FALSE)
  names(split.points) <- paste("split.points.", 1:d.f, sep = "")
  object <- list(DoE = DoE, split.points = split.points, d.f = d.f, variable.names = variable.names, 
                 mirrored.runs.included = mirrored.runs.included, method = method)
  class(object) <- "safidesign"
  return(object)
} 