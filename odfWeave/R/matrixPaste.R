matrixPaste <- function(..., sep = rep(" ", length(list(...)) - 1))
{
   theDots <- list(...)
   if(any(unlist(lapply(theDots, function(x) !is.character(x)))))
      stop("all matrices must be character")

   numRows <-  unlist(lapply(theDots, nrow))
   numCols <-  unlist(lapply(theDots, ncol))

   if(length(unique(numRows)) > 1 | length(unique(numCols)) > 1)
      stop("all matrices must have the same dim")

   for(i in seq(along = theDots)) out <- if(i == 1) theDots[[i]] else paste(out, theDots[[i]], sep = sep[i - 1])
   matrix(out, nrow = numRows[1])
}

#mat1 <- matrix(letters[1:6], nrow = 2)
#mat2 <- matrix(LETTERS[1:6], nrow = 2)
#mat3 <- matrix(paste(1:6), nrow = 2)
#
#matrixPaste(mat1, mat2)

