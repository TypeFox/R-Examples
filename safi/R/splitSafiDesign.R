splitSafiDesign <- function(s.d, new.split.points) {
  if (class(s.d) != "safidesign") 
    stop("object class is not safidesign")
  d.f <- s.d$d.f
  if (d.f != length(new.split.points)) 
    stop("new.split.points must be a list of length d.f")
  split.points <- s.d$split.points
  for (i in 1:d.f) {
    if (any(duplicated(c(split.points[[i]], new.split.points[[i]])))) 
      stop("new.split.points must contain the new points only")
  }
  method <- s.d$method
  split.points.new <- list()
  if (method != "SB" && method != "other") 
    stop("method not supported")
  for (i in 1:(s.d$d.f)) {
    split.points.new[[i]] <- sort(c(split.points[[i]], new.split.points[[i]]))
  }
  DoE <- s.d$DoE
  mirrored <- s.d$mirrored.runs.included
  d.f <- s.d$d.f
  # preperations
  anzahl <- NULL  # number of columns
  matches <- NULL  # from here new rows get a minus
  for (i in 1:d.f) {
    anzahl <- c(anzahl, diff(match(split.points[[i]], split.points.new[[i]])))
    matches <- c(matches, match(split.points.new[[i]], split.points[[i]]))
    matches <- matches[-length(matches)]
  }
  # add new columns
  DoE.new <- as.matrix(data.frame(mapply(function(x, a) kronecker(x, t(rep(1, a))), data.frame(DoE), 
                                         anzahl, SIMPLIFY = FALSE)))
  if (method == "SB") {
    # help function for new rows
    line <- function(i, M, mirrored) {
      M <- rbind(M, 1)
      M[dim(M)[1], i:(dim(M)[2])] <- -1
      if (mirrored) {
        M <- rbind(M, -1)
        M[dim(M)[1], i:(dim(M)[2])] <- 1
      }
      return(M)
    }
    # add new rows
    index <- which(is.na(matches))
    for (i in 1:length(index)) {
      DoE.new <- line(index[i], DoE.new, mirrored)
    }
  }
  names(split.points.new) <- paste("g", 1:d.f, sep = "")
  s.d$DoE <- DoE.new
  s.d$split.points <- split.points.new
  return(s.d)
} 