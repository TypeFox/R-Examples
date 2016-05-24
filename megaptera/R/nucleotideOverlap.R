nucleotideOverlap <- function (x) {
  
  info <- as.raw(c(136, 40, 72, 24, 48, 96, 192))
  notinfo <- as.raw(c(240, 02, 04))
  xx <- matrix(nrow = nrow(x), ncol = ncol(x), 
               dimnames = list(rownames(x), seq_along(x[1, ])))
  xx[arrayInd(which(x %in% info), .dim = dim(xx))] <- 1
  xx[arrayInd(which(x %in% notinfo), .dim = dim(xx))] <- 0
  
  nms <- rownames(x)
  id <- seq_along(nms)
  no <- matrix(nrow = nrow(x), ncol = nrow(x), 
               dimnames = list(nms, nms))
  
  ## doble loop very time-consuming!
#   system.time(
    for ( i in id ){
      for ( j in tail(id, -i) ) {
        n <- xx[i, ] + xx[j, ]
        no[i, j] <- no[j, i] <- length(which(n == 2))
      }  
    }
#     )
  no
}