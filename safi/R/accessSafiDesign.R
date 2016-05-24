accessSafiDesign <- function(s.d, n.timepoints) {
  if (class(s.d) != "safidesign") 
    stop("object class is not safidesign")
  if (length(n.timepoints) == 1) 
    n.timepoints <- rep(n.timepoints, s.d$d.f)
  if (length(n.timepoints) != s.d$d.f) {
    stop("n.timepoints must be vector of length 1 or d.f")
  }
  split.points <- s.d$split.points
  DoE <- s.d$DoE
  DoE.new <- vector("list", length(n.timepoints))
  n <- 0
  
  for (i in 1:length(n.timepoints)) {
    anzahl <- diff(round(split.points[[i]] * n.timepoints[i]))
    DoE.new[[i]] <- as.matrix(data.frame(mapply(function(x, a) kronecker(x, t(rep(1, a))), data.frame(DoE[, 
                                                                                                          (n + 1):(n + length(anzahl))]), anzahl, SIMPLIFY = FALSE)))
    
    colnames(DoE.new[[i]]) <- NULL
    n <- n + length(anzahl)
  }
  names(DoE.new) <- paste("x", 1:length(split.points), sep = "")
  return(DoE.new)
} 