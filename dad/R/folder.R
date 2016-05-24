# Prendra en argument :
#   - soit 2 (ou +) data frames
#   - soit liste de data frames
#
# On écrira une fonction df2folder() qui prendra en argument un data.frame
# et le nom de la variable de grouupe.

folder <-
function(x1, x2 = NULL, ..., cols.select = "intersect") 
{
  # x1  :        data.frame or list of data frames; the data.
  #              - If x1 is a data.frame, it contains the first data set of the
  #                folder which will be built, and x2 must be provided.
  #              - If x1 is a list of data frames, x2 must not be provided.
  # x2  :        data.frame; the second data set.
  # cols.select : string. It can be:
  #              - "intersect": only columns (column names) that are common
  #                             to every data frames are kept.
  #                             We will have: "attr(fold, same.cols) = TRUE".
  #              - "union":     all the columns (and column names) of every
  #                             data frames are kept; the lines are completed
  #                             by NA if necessary.
  #                             We will have: "attr(fold, same.cols) = TRUE".
  # ... :        other data.frames (optional)
  
  # Checking the class or value of each argument
  
  if (!is.list(x1))
    stop("x1 must be a data frame or a list of data frames.")
  
  if ((!is.data.frame(x1)) & (!is.null(x2)))
    warning("If x1 is not a data frame, x2 is omitted.")
  
  if ((is.data.frame(x1)) & (is.null(x2)))
    stop("You cannot build a 'folder' with only one data frame.\nIf x1 is a data frame, x2 must be provided.")
  
  if ((!is.data.frame(x1)) & (length(x1) == 1))
    stop("You cannot build a 'folder' with only one data frame.\nIf x1 is a list, it must contain at least two data frames.")
  
  if ((is.data.frame(x1)) & (!is.data.frame(x2)))
    stop("If x1 is a data frame, x2 must also be a data frame.")
  
  if (!cols.select %in% c("intersect", "union"))
    stop("cols.select: wrong value. It must be either 'intersect' or 'union'")
  
  if (!is.data.frame(x1)) {
    # x1 is a list of data frames
    class.arg <- "l"
  } else {
    # x1 and x2 are data frames
    if (is.data.frame(x2)) {
      class.arg <- "d2"
    }
    dots <- list(...)
    if (length(dots) > 0) {
      if (prod(unlist(lapply(dots, is.data.frame)))) {
        # All arguments in "..." are data frames
        class.arg <- "d3"
      } else {
        warning(paste("Argument(s)", names(dots)[unlist(lapply(dots, is.data.frame))], "is/are no data frame(s).\n", names(dots), "arguments will not be used."))
      }
    }
  }
  
  switch(class.arg,
    d2 = {
      # Creation of the list
      fold <- list(x1, x2)
      
      # Number of datasets in the folder
      ndata <- 2
  
      # Names of the elements of the folder
      namesfold <- c(as.character(match.call(expand.dots = FALSE))[2],
                     as.character(match.call(expand.dots = FALSE))[3])
      names(fold) <- namesfold
    },
    d3 = {
      # Creation of the list
      fold <- list(x1, x2)
      if (!prod(unlist(lapply(dots, is.data.frame))))
        stop("All arguments in '...' must be data frames")
      fold <- c(fold, dots)
      
      # Number of data frames in the folder
      ndata <- length(fold)
      
      # Names of the elements of the folder
      namesfold <- c(as.character(match.call(expand.dots = FALSE))[2],
                     as.character(match.call(expand.dots = FALSE))[3])
      namesfold <- c(namesfold, as.character(match.call(expand.dots = FALSE)$...))
      names(fold) <- namesfold
    },
    l = {
      # Check if all elements of x1 are data frames
      if (!prod(unlist(lapply(x1, is.data.frame))))
        stop("All elements of x1 must be data frames.")
      
      # Creation of the list
      fold <- x1
      
      # Number of datasets in the folder
      ndata <- length(fold)
    }
  )
  
  # If (cols.select == "intersect"): only the columns whose names are common
  # to each element of 'fold' will be kept.
  if (cols.select == "intersect") {
    cnames.l <- lapply(fold, names)
    cnames <- cnames.l[[1]]
    for (n in 2:ndata)
      cnames <- intersect(cnames, cnames.l[[n]])
    for (n in 1:ndata) {
      fold[[n]] <- fold[[n]][cnames]
    }
    same.cols <- TRUE
  }
  
  # If (cols.select == "union"): add columns to each element of 'fold',
  # so that they all have the same columns and colnames.
  if (cols.select == "union") {
    cnames <- unique(unlist(lapply(fold, names)))
    adjcnames <- as.data.frame(t(vector(length = length(cnames))))
    colnames(adjcnames) <- cnames
    for (n in 1:ndata) {
      foldn <- merge(fold[[n]], adjcnames, all = TRUE, sort = FALSE)[1:nrow(fold[[n]]), ]
      fold[[n]] <- foldn[cnames]
    }
    same.cols <- TRUE
  }
  
  # Object which will be returned
  class(fold) <- "folder"
  attr(fold, "same.cols") <- same.cols
  attr(fold, "same.rows") <- FALSE
  
  return(fold)
}
