#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2004)                              ####
####                                                 ####
#### FILE:       traceplot2.R                        ####
####                                                 ####
#### FUNCTIONS:  traceplot2                          ####
#########################################################

### ======================================
### traceplot2
### ======================================
traceplot2 <- function(x, chains, bty = "n", main, xlab, ...)
{
  if (attr(x, "class") != "mcmc") stop("This function handles only objects of class mcmc.")
  
  if (missing(chains)){
    if (is.null(attr(x, "dim"))) chains <- 1
    else                         chains <- 1:attr(x, "dim")[2]
  }
  else{
    if (is.null(attr(x, "dim")))
      chains <- 1
    else{
      ch <- chains %in% 1:attr(x, "dim")[2]
      chains <- chains[ch]
    }      
  }

  iters <- attr(x, "mcpar")
  iters <- seq(iters[1], iters[2], by = iters[3])

  if (missing(xlab)) xlab <- "Iteration"
  
  for (i in 1:length(chains)){
    name <- attr(x, "dimnames")[[2]][chains[i]]
    if (is.null(name)) name <- paste("var ", i, sep = "")
    if (missing(main)) mmain <- paste("Trace of ", name, sep = "")
    else               mmain <- main
    if (is.null(attr(x, "dim"))) sample <- x
    else                         sample <- x[,chains[i]]
    
    plot(iters, sample, type = "l", lty = 1, xlab = xlab, bty = bty, ...)
    title(main = mmain)    
  }

  return(invisible(x))    
}  
