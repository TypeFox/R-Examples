set.hlevel.stuff <-
function(x, outfile, control)
{
  for(k in 1L:length(x)) {
    if(is.null(x[[k]]$hlevel))
      x[[k]]$hlevel <- 1L
    x[[k]]$hlevel <- x[[k]]$hlevel + 1L
    x[[k]]$first <- FALSE
    x[[k]]$outfile <- outfile
    x[[k]]$iterations <- control$iterations
    x[[k]]$burnin <- control$burnin
    x[[k]]$step <- control$step
    x[[k]]$level1 <- control$level1
    x[[k]]$level2 <- control$level2
    if(!is.null(x[[k]]$h.random))
      x[[k]]$h.random <- set.hlevel.stuff(x[[k]]$h.random, outfile, control)
  }

  return(x)
}

