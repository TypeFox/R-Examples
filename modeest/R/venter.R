# Author: P. Poncet

venter <-
function(x,                  # sample (the data)
         bw = NULL,          # 'lms' parametrization : fraction of the observations to be considered for the shortest interval
         k,                  # length of the intervals
         iter = 1,           # number of iterations
         type = 1,           # -Inf, 1, 2, 3, 4, 5, 6, "-Inf", "1", "2", "3", "dalenius", "4", "shorth", "5", "ekblom", "6", "hsm"
         tie.action = "mean",
         tie.limit = 0.05)
{
################################################################################################
# Dalenius' / Venter's / LMS / Shorth mode estimator (Andrews et al., 1972)
# The estimate is the mean of the shortest interval among intervals containing 'k+1' observations.
# LMS = Least median of squares (Rousseeuw and Leroy, 1987)
################################################################################################
  
  ny <- length(x)
    
  ## Initialization
  type <- match.arg(tolower(as.character(type)), c("-inf", "1", "2", "3", "dalenius", "4", "shorth", "5", "ekblom", "6", "hsm"))
  if (type == "3") type <- "dalenius"
  if (type == "4") type <- "shorth"
  if (type == "-Inf" | type == "5") type <- "ekblom"
  if (type == "6") type <- "hsm"
  
  if (type == "hsm") return(hsm(x = x, bw = bw, k = k, tie.action = tie.action, tie.limit = tie.limit))
  
  if (missing(k) & !is.null(bw)) {
    if (bw <= 0 | bw > 1) stop("argument 'bw' must belong to (0, 1]")
    k <- ceiling(bw*ny) - 1
  } else if (missing(k) & is.null(bw)) {
    if (type == "ekblom") {
      k <- 1
    } else {
      k <- ceiling(ny/2) - 1
    }
  }
    
  if (k < 0 | k >= ny) stop("argument 'k' must belong to [0, length('x'))") 
  
  y <- sort(x)

  inf <- y[1:(ny-k)]
  sup <- y[(k+1):ny]
  diffs <- sup - inf
  i <- which(diffs==min(diffs))
  
  ## Ties?
  if (length(i) > 1) i <- .deal.ties(ny, i, tie.action, tie.limit)
  
  ## Output
  M <- switch(type,
              "1" = (y[i] + y[i+k])/2,
              "2" = y[i+floor((k+1)/2)],
              "dalenius" = median(y[i:(i+k)]),
              "shorth" = mean(y[i:(i+k)]),
              "ekblom" = ifelse(y[i+2]-y[i+1]>y[i]-y[i-1], y[i], y[i+1]))
  
  if (iter > 1) {
    M <- Recall(x=y[i:(i+k)], bw=(k+1)/ny, iter = iter-1, type=type, tie.action=tie.action, tie.limit=tie.limit)
  }
    
  #attr(M, "inf") <- y[i]
  #attr(M, "sup") <- y[i+k]
  return(M)
}

#dalenius <- venter
#Venter <- venter
#lms <- venter
shorth <- function(x, ...) venter(x, type = "shorth", ...)
