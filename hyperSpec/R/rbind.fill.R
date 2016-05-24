### These functions were thought to move eventually into Hadley Wickham's plyr. This is not going to
### happen, soon, however.

##' Bind matrices by row, and fill missing columns with NA
##'
##' The matrices are bound together using their column names or the column indices (in that order of
##' precedence.) Numeric columns may be converted to character beforehand, e.g. using format.
##' If a matrix doesn't have colnames, the column number is used.
##'
##' @author C. Beleites
##' @seealso   \code{\link[base]{rbind}}, \code{\link[base]{cbind}}, \code{\link[plyr]{rbind.fill}}
##' @export
##' @keywords manip
##' @rdname rbind-fill
##' @examples 
##'  A <- matrix (1:4, 2)
##'  B <- matrix (6:11, 2)
##'  A
##'  B
##'  rbind.fill.matrix (A, B)
##' 
##'  colnames (A) <- c (3, 1)
##'  A
##'  rbind.fill.matrix (A, B)
##' 
##' @return a matrix
rbind.fill.matrix <- function (...){
  matrices <- list (...)
  
  ## check the arguments
  if (! all  (sapply (matrices, is.matrix)))
    stop ("Input ", which (! sapply (matrices, is.matrix)), "is no matrix.")
  
  ## if the matrices have column names, use them 
  lcols <- lapply (matrices, .cols)
  cols  <- unique (unlist (lcols))

  ## preallocate the new spectra matrix
  pos <- sapply (matrices, nrow)
  result <- matrix (NA, nrow = sum (pos), ncol = length (cols))

  ## make an index vector for the row positions
  pos <- c (0, cumsum (pos))

  ## fill in the new matrix 
  for (i in seq_along (matrices)){
    icols <- match (lcols[[i]], cols)
    result [(pos [i] + 1) : pos [i + 1], icols] <- matrices [[i]]
  }

  attr (result, "cols") <- cols
  colnames (result) <- cols
  
  result
}

.cols <- function (x){
  cln <- colnames (x)
  if (is.null (cln)) cln <- seq_len (ncol (x))

  cln
}


##' Quick data frame.
##' Experimental version of \code{\link{as.data.frame}} that converts a
##' list to a data frame, but doesn't do any checks to make sure it's a
##' valid format. Much faster.
##'
##' @param list list to convert to data frame
##' @keywords internal
quickdf <- function(list) {
  if (is.matrix (list [[1]]))
    n <- nrow (list [[1]])
  else
    n <- length (list [[1]])
  
  structure(list,
    class = "data.frame",
    row.names = seq_len(n))
}

##' Bind matrices by row, and fill missing columns with NA
##'
##' The matrices are bound together using their column names or the column indices (in that order of
##' precedence.) Numeric columns may be converted to character beforehand, e.g. using format.  If a
##' matrix doesn't have colnames, the column number is used (via \code{\link[base]{make.names}(unique
##' = TRUE)}).
##'
##' Note that this means that a column with name \code{"X1"} is merged with the first column of a
##' matrix without name and so on.
##'
##' Vectors are converted to 1-column matrices prior to rbind.
##'
##' Matrices of factors are not supported. (They are anyways quite inconvenient.) You may convert
##' them first to either numeric or character matrices. If a character matrix is merged with a
##' numeric, the result will be character.
##'
##' Row names are ignored.
##'
##' The return matrix will always have column names.
##'
##' @aliases rbind.fill
##' @rdname rbind-fill
##' @author C. Beleites
##' @seealso   \code{\link[base]{rbind}}, \code{\link[base]{cbind}}, \code{\link[plyr]{rbind.fill}}
##' @rdname rbind-fill
##' @keywords manip
##' @examples 
##'  A <- matrix (1:4, 2)
##'  B <- matrix (6:11, 2)
##'  A
##'  B
##'  rbind.fill.matrix (A, B)
##' 
##'  colnames (A) <- c (3, 1)
##'  A
##'  rbind.fill.matrix (A, B)
##'
##'  rbind.fill.matrix (A, 99)
##' 
##' @return a matrix
rbind.fill.matrix <- function (...){
  matrices <- list (...)
  
  ## check the arguments
  tmp <- unlist (lapply (matrices, is.factor))
  if (any  (tmp))
    stop ("Input ", paste (which (tmp), collapse = ", "),
          " is a factor and needs to be converted first to either numeric or character.")

  tmp <- ! unlist (lapply (matrices, is.matrix))
  matrices [tmp] <- lapply (matrices [tmp], as.matrix)
  
  ## if the matrices have column names, use them 
  lcols <- lapply (matrices, .cols)
  cols  <- unique (unlist (lcols))

  ## the new row positions
  pos <- unlist (lapply (matrices, nrow)) # Hadley, for me nrow is about twice as fast as
                                          # .row_names_info (for matrices), the other way round for
                                          # data.frame

  ## preallocate the new spectra matrix
  result <- matrix (NA, nrow = sum (pos), ncol = length (cols))

  ## make an index vector for the row positions
  pos <- c (0, cumsum (pos))

  ## fill in the new matrix 
  for (i in seq_along (matrices)){
    icols <- match (lcols[[i]], cols)
    result [(pos [i] + 1) : pos [i + 1], icols] <- matrices [[i]]
  }

  colnames (result) <- cols
  
  result
}

.cols <- function (x){
  cln <- colnames (x)
  if (is.null (cln))
    cln <- make.names (seq_len (ncol (x)), unique = TRUE)
  cln
}

##' Combine objects by row, filling in missing columns.
##' \code{rbind}s a list of data frames filling missing columns with NA.
##'  
##' This is an enhancement to \code{\link{rbind}} which adds in columns
##' that are not present in all inputs, accepts a list of data frames, and 
##' operates substantially faster
##'  
##' @param ... data frames/matrices to row bind together
##' @keywords manip
##' @rdname rbind-fill
##' @examples
##' #' rbind.fill(mtcars[c("mpg", "wt")], mtcars[c("wt", "cyl")])
rbind.fill <- function(...) {
  dfs <- list(...)
  if (length(dfs) == 0) return(list())
  if (is.list(dfs[[1]]) && !is.data.frame(dfs[[1]])) {
    dfs <- dfs[[1]]
  }
  dfs <- dfs [!sapply (dfs, is.null)]  # compact(dfs) -> dependency plyr.
  
  if (length(dfs) == 1) return(dfs[[1]])
  
  # About 6 times faster than using nrow
  rows <- unlist(lapply(dfs, .row_names_info, 2L))
  nrows <- sum(rows)
  
  # Build up output template -------------------------------------------------
  vars <- unique(unlist(lapply(dfs, base::names)))   # ~ 125,000/s
  output <- rep(list(rep(NA, nrows)), length(vars))  # ~ 70,000,000/s
  names(output) <- vars
  
  seen <- rep(FALSE, length(output))
  names(seen) <- vars

  ## find which cols contain matrices 
  matrixcols <- unique (unlist (lapply (dfs, function (x)
                                        names (x) [sapply (x, is.matrix)])
                                ))
  seen [matrixcols] <- TRUE             # class<- will fail if the matrix is not protected by I
                                        # because 2 dims are needed
    
  for(df in dfs) {    
    if (all(seen)) break  # Quit as soon as all done

    matching <- intersect(names(df), vars[!seen])
    for(var in matching) {
      value <- df[[var]]
      if (is.factor(value)) {
        output[[var]] <- factor(output[[var]])
      } else {
        class(output[[var]]) <- class(value)
      }
    }
    seen[matching] <- TRUE
  }
  
  # Set up factors
  factors <- names(output)[unlist(lapply(output, is.factor))]
  for(var in factors) {
    all <- unique(lapply(dfs, function(df) levels(df[[var]])))  # is that unique needed?
    levels(output[[var]]) <- unique(unlist(all))
  }

  # care about matrices
  # the trick is to supply a n by 0 matrix for input without column of that name
  for (var in matrixcols){
    df <- lapply (dfs, .get.or.make.matrix, var)
    output [[var]] <- I (do.call (rbind.fill.matrix, df))
  }
  
  # Compute start and end positions for each data frame
  pos <- matrix(cumsum(rbind(1, rows - 1)), ncol = 2, byrow = TRUE)
  
  for(i in seq_along(rows)) { 
    rng <- pos[i, 1]:pos[i, 2]
    df <- dfs[[i]]
    
    for(var in setdiff (names (df), matrixcols)) {
      output[[var]][rng] <- df[[var]]
    }
  }  
  
  quickdf(output)
}

.get.or.make.matrix <- function (df, var){
  tmp <- df [[var]]
  if (is.null (tmp))
    tmp <- I (matrix (integer (), nrow = nrow (df)))
  tmp
}


