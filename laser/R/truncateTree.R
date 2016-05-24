`truncateTree` <-
function(x, omit.time = NULL, omit.nodes = NULL, batch = FALSE)
{
  #takes set of branching times and truncates the set at the cut-point entered
  #by user.  Parameter "omit.time" should correspond to the amount of time
  # (or genetic distance) before the present where you would like set truncated.
  #E.g., if X = (100, 80, 50, 20, 10, 5),
  # truncatetree(X, omit.time = 20) would return a vector of (80, 60, 30, 0).

  #'omit.nodes' indicates that you would like to omit the final 'omit.nodes' branching times
  #e.g., omit.nodes = 1 omits the final branching time.   E.g., X = (100, 80, 50, 20, 10, 5)
  # then truncatetree(X, omit.nodes = 1) returns x = (95, 75, 45, 15, 5).
  if (!is.numeric(x)) stop("object x not of class 'numeric'")
  if (batch == FALSE)
  {
    x <- rev(sort(x))

    if (!is.null(omit.time))
    {
      xt <- x[x >= omit.time]
      xt <- xt - omit.time
    }
    if (!is.null(omit.nodes) && omit.nodes != 0)
    {
      xt <- x[1:(length(x) - omit.nodes)]
      xt <- xt - x[(length(x) - omit.nodes + 1)]
    }
  }
  #returns an error if user inputs 'omit.nodes = 0'.  Why?
  if (batch == TRUE)
  {
    xt <- matrix(NA, nrow = nrow(x), ncol = (ncol(x) - omit.nodes))
    for (i in 1:nrow(x))
    {
      x[i, ] <- rev(sort(x[i,]))

      if (!is.null(omit.time))
      {
        stop("Can only omit nodes in batch mode\nCannot use omit.time option\n")
      }
      if (!is.null(omit.nodes) && omit.nodes != 0)
      {
        xt[i,] <- x[i, 1:(length(x[i,]) - omit.nodes)]
        xt[i,] <- xt[i,] - x[i, (length(x[i, ]) - omit.nodes + 1)]
      }

     }
  }
  return(xt)

}

