#' Derive the optimal lambda and c value
#'
#' Derives the optimal lambda and c-value for a given configuration.
#' 
#' @param alpha The desired confidence coefficient.
#' @param n The number of replications per population.
#' @param p The number of populations considered. This must be present if
#' min.loc is equal to 'zero'.
#' @param k The number of populations selected.
#' @param var.known A logical flag indicating if the variance of the
#' observations is known exactly. It is TRUE by default.
#' @param eps The grid size that is to be set up.
#'
#' @export
#'
#' @details This function will return the optimal lambda and c-value to be used,
#' using a grid search.
#'
#' There are essentially 2 different cases to consider. They correspond to the
#' cases when the variance is known or unknown.
#'
#' @return The function returns a list with two components, lambda and c.val
#' that are optimal.

optimalLambdaC <- function(alpha=0.05, n, p, k=1, var.known=TRUE, eps=0.1) {
  # generate the lambda and c-value sequences
  lam.seq <- seq(from=0, to=1, by=eps)
  c.val.seq <- seq(from=0, to=6, by=eps)
  lam.c.grid <- expand.grid(lam=lam.seq, c.val=c.val.seq)
  
  if(var.known) {
    # Known variance section
    eq.0 <- ((pnorm(lam.c.grid[,2]))^(p-k+1) - (pnorm(-lam.c.grid[,1]*lam.c.grid[,2]))^(p-k+1)) * 
      (pnorm(lam.c.grid[,2]) - pnorm(-lam.c.grid[,1]*lam.c.grid[,2]))^(k-1)

    eq.inf <- (pnorm(lam.c.grid[,2]) - pnorm(-lam.c.grid[,1]*lam.c.grid[,2]))^k

  } else {
    # Unknown variance section

    if(k==1) {
      eq.0 <- mapply(integrate2a, lambda=lam.c.grid$lam,
        c.val=lam.c.grid$c.val, n=n, p=p)
   
     eq.inf <- mapply(integrate2b, lambda=lam.c.grid$lam,
       c.val=lam.c.grid$c.val, n=n, p=p)

    } else { 
     eq.0 <- mapply(integrate4a, lambda=lam.c.grid$lam,
       c.val=lam.c.grid$c.val, n=n, p=p, k=k)

     eq.inf <- mapply(integrate4b, lambda=lam.c.grid$lam,
       c.val=lam.c.grid$c.val, n=n, p=p, k=k)
    }

  }

  x <- which(eq.0 > 1-alpha)
  y <- which(eq.inf > 1-alpha)
  z <- intersect(x,y)

  lam.c.grid.sub <- lam.c.grid[z,]
  interval.lengths <- (1+lam.c.grid.sub[,1])*lam.c.grid.sub[,2]
  shortest.length.id <- which.min(interval.lengths)
  return(list(lambda=lam.c.grid.sub[shortest.length.id,1], 
    c.val=lam.c.grid.sub[shortest.length.id,2]))
}
