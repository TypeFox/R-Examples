#' New Revealed Comparative Advantage
#' 
#' @name nrca
#' @param x A decomposed Inter-Country Input Output table as created by decompr
#' @export
#' @examples 
#' # load the decompr package
#' library(decompr)
#'            
#' # load the example data set
#' data(leather)
#' 
#' # perform Leontief decomposition
#' l <- decomp(inter,
#'             final,
#'             countries,
#'             industries,
#'             out,
#'             method = "leontief",
#'             post = "exports"    )
#' 
#' # load gvc package
#' library(gvc)
#' 
#' # perform New Revealed Comparative Advantage
#' nrca(l)
#' 


nrca <- function ( x ) {
  
 if ( attr(x, "long") == TRUE ) {
    
    # extract attributes
    k      <- attr(x, "k")
    i      <- attr(x, "i")
    G <- length(k)
    N <- length(i)
    
    x <- matrix(x[,5], nrow=G*N, byrow=TRUE)
    
    # remove anything but exports to self
    f <- rowSums(diagonals::fatdiag(diagonals::fatdiag(x, steps=G ), steps=G) )
    t <- rowSums(diagonals::fatdiag(diagonals::fatdiag(x, steps=G ), steps=G) )
    q <- sum(rowSums(diagonals::fatdiag(diagonals::fatdiag(x, steps=G ), steps=G) ))
    
    # sum across rows (source industry)
    # divide by own exports to self
    for (j in 1:G) {
      s <- seq( ((j-1)*N + 1), j*N )
      f[s] <- f[s] / sum(f[s])
    }
    
    # 
    for (i in 1:N) {
      p <- (seq(1:G)*N) - N + i
      t[p] <- sum(t[p])
    }
    
    Eij <- f
    Eit <- 1
    Enj <- t
    Ent <- q
    
  } else {
    
    # extract attributes
    k      <- attr(x, "k")
    i      <- attr(x, "i")
    G <- length(k)
    N <- length(i)
    
    # remove anything but exports to self
    f <- rowSums(diagonals::fatdiag(diagonals::fatdiag(x, steps=G ), steps=G) )
    t <- rowSums(diagonals::fatdiag(diagonals::fatdiag(x, steps=G ), steps=G) )
    q <- sum(rowSums(diagonals::fatdiag(diagonals::fatdiag(x, steps=G ), steps=G) ))
    
    # sum across rows (source industry)
    # divide by own exports to self
    for (j in 1:G) {
      s <- seq( ((j-1)*N + 1), j*N )
      f[s] <- f[s] / sum(f[s])
    }
    
    # 
    for (i in 1:N) {
      p <- (seq(1:G)*N) - N + i
      t[p] <- sum(t[p])
    }

    Eij <- f
    Eit <- 1
    Enj <- t
    Ent <- q
    
  }
  
  
  # return Bela Balassa (1965) ratio
  return( (Eij/Eit)/(Enj/Ent) )
  
}
