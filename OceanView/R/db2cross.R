## =============================================================================
## From a table in "database format" (3-columns) to a "crosstab"
## =============================================================================

# Creates factors for values that are at a distance of at most df
createfactor <- function(var, df = 0) {
  if (df == 0 | is.na(df)) {
    var.unique <- sort(unique(var))
    return(list(var.input = var.unique, factor = 1 : length(var.unique), var.output = var.unique))
  } else {  
    var.unique <- sort(unique(var))
    var.diff   <- c(df+1, diff(var.unique)) # df+1 to make sure first value is selected
    var.min    <- var.unique[which(var.diff >= df)]
    var.factor <- sapply(X = var.unique, 
      FUN = function(x) which.min(abs(x-var.min)))
    var.mean   <- sapply(X = 1:length(var.min), 
      FUN = function(x) mean(var.unique[which(var.factor == x)]))
    return(list(var.input = var.unique, factor = var.factor, var.output = var.mean))
  }
}

# Creates factors for values that are closest to out
closestfactor <- function(var, out) {
    var.unique <- sort(unique(var))
    var.factor <- sapply(X = var.unique, 
      FUN = function(x) which.min(abs(x-out)))
    return(list(var.input = var.unique, factor = var.factor, var.output = out))
}

                                               
db2cross <- function(input, row = 1, col = 2, value = 3, subset = NULL, 
  df.row = NA, df.col = NA, out.row = NA, out.col = NA, 
  full.out = FALSE) {

  if (is.character(value))
    value <- which(colnames(input) == value)
  if (is.character(row))
    row <- which(colnames(input) == row)
  if (is.character(col))
    col <- which(colnames(input) == col)

  # quick and dirty
  if (!missing(subset)){
    e <- substitute(subset)
    r <- eval(e, as.data.frame(input), parent.frame())
    if (!is.logical(r))
      stop("'subset' must evaluate to logical")
    isub <- r & !is.na(r)
    input <- input[isub, ]
  }  

# Row can point to more than one column?
  if (is.character(input[, col]))
    input[, col] <- as.factor(input[, col])
    
  if (is.character(input[, row]))
    input[, row] <- as.factor(input[, row])

  if (is.character(input[, value]))
    stop(" cannot expand input; value should point to a numeric column") 

  notna <- which(!is.na(input[ , value]))
  input <- input[notna, ]
  # Check if input has only 3 columns
  IN <- cbind(as.double(input[, col]), 
              as.double(input[, row]), 
              as.double(input[, value]))
  dim.input <- dim(IN)
  if (ncol(IN) != 3)
    stop ("'input' should have three columns")

  if (! is.na(df.row))
    if (df.row < 0) stop("df.row should be a positive value or 0")

  if (! is.na(df.col))
    if (df.col < 0) stop("df.col should be a positive value or 0")

# which fortran function to use depends on df.row and df.col
  use2 <- ifelse (!is.na(df.row) | !is.na(df.col) | 
                  !all(is.na(out.row)) | !all(is.na(out.col)) , TRUE, FALSE)
  
  if (length(out.row) == 1 & is.na(out.row[1]))
    rowfac <- createfactor(input[, row], df.row)
  else
    rowfac <- closestfactor(input[, row], out.row)  

  if (length(out.col) == 1 & is.na(out.col[1]))
    colfac <- createfactor(input[, col], df.col)
  else
    colfac <- closestfactor(input[, col], out.col)

# values or names of rows and columns in crosstable
  rows <- rowfac$var.output
  cols <- colfac$var.output

  nr <- as.integer(length(rows))
  nc <- as.integer(length(cols))
  
  NAnum <- min(IN, na.rm = TRUE) - 100
  
  IN[is.na(IN)] <- NAnum
  

  if (! use2) {
  
    out <- .Fortran("crosstab", t(IN), as.integer(nrow(input)),
             as.integer(1), as.integer(2), as.integer(3),
             as.double(cols), as.double(rows), nr = nr, nc = nc, 
             cross = matrix(nrow = nr, ncol = nc, data = as.double(0)),
             count = matrix(nrow = nr, ncol = nc, data = as.integer(0)),
             NAnum = as.double(NAnum), PACKAGE = "OceanView") 
  } else {        # Need also index of each unique rowvalue/columnvalue to used row and column
    indrow <- rowfac$factor
    indcol <- colfac$factor
    
    nrow <- as.integer(length(indrow))
    ncol <- as.integer(length(indcol))
    out <- .Fortran("crosstab2", t(IN), as.integer(nrow(input)),
             as.integer(1), as.integer(2), as.integer(3),
             as.double(colfac$var.input), as.double(rowfac$var.input), 
             nr = nr, nc = nc,
             nrow = nrow, ncol = ncol, indrow = as.integer(indrow),
             indcol = as.integer(indcol),  
             cross = matrix(nrow = nr, ncol = nc, data = as.double(0)),
             count = matrix(nrow = nr, ncol = nc, data = as.integer(0)),
             NAnum = as.double(NAnum), PACKAGE = "OceanView") 
    
  }
  
  z <- out$cross
  z[out$count == 0] <- NA
  
  colnames(z) <- cols
  rownames(z) <- rows            
  out <- list(x = rows, y = cols, z = z)             
  if (full.out)
    out$map <- list(x = rowfac, y = colfac)  
    
  return(out)  
}
