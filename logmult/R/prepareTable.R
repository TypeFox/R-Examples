prepareTable <- function(tab, twoWay=TRUE, rowsup=NULL, colsup=NULL) {
  if(!twoWay & length(dim(tab)) != 3)
      stop("only three dimensional tables are supported")
  else if(!length(dim(tab)) %in% 2:3)
      stop("only two and three dimensional tables are supported")


  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback for both variable and level names
  if(is.null(rownames(tab)))
      rownames(tab) <- LETTERS[seq.int(nrow(tab))]

  if(is.null(colnames(tab)))
      colnames(tab) <- LETTERS[seq.int(ncol(tab))]

  if(length(dim(tab)) > 2 && is.null(dimnames(tab)[[3]]))
      dimnames(tab)[[3]] <- LETTERS[seq.int(dim(tab)[3])]

  if(length(names(dimnames(tab))) > 0)
      names(dimnames(tab)) <- make.names(names(dimnames(tab)), unique=TRUE)
  else
      names(dimnames(tab)) <- paste0("Var", seq.int(length(dim(tab))))


  # gnm doesn't include coefficients for row/columns that are empty or with NA name
  # so get rid of them too
  if(length(dim(tab)) == 2)
      tab <- as.table(tab[!is.na(rownames(tab)) & rowSums(tab, na.rm=TRUE) > 0,
                         !is.na(colnames(tab)) & colSums(tab, na.rm=TRUE) > 0])
  else if(length(dim(tab)) == 3)
      tab <- as.table(tab[!is.na(rownames(tab)) & apply(tab, 1, sum, na.rm=TRUE) > 0,
                          !is.na(colnames(tab)) & apply(tab, 2, sum, na.rm=TRUE) > 0,
                          !is.na(dimnames(tab)[[3]]) & apply(tab, 3, sum, na.rm=TRUE) > 0])

  if(!is.null(rowsup) && !is.matrix(rowsup))
      stop("'rowsup' must be a matrix")

  if(!is.null(colsup) && !is.matrix(colsup))
      stop("'colsup' must be a matrix")

  if(!is.null(rowsup) && ncol(rowsup) != ncol(tab))
      stop("'rowsup' must have one column for each column in model data")

  if(!is.null(colsup) && nrow(colsup) != nrow(tab))
      stop("'colsup' must have one row for each row in model data")

  tab
}
