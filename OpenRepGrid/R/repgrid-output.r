# trim column of a matrix to equal width
# x: matrix
trim_column_width <- function(x, just="l"){
  lengths <- apply(nchar(x), 2, max)
  for (j in seq_len(ncol(x))) {
    x[ ,j] <- format(x[ ,j], just = just, width = lengths[j])
  }
  x
}
#trim_column_width(x, just="r")


collapse_matrix <- function(x, collapse="", sep=" "){
  x <- apply(x, 1, paste, 
             collapse = collapse, sep = sep)  # collapse whole matrix
  do.call(rbind, as.list(x))                  # bind whole matrix together
}


# x:  matrix
matrix_to_single_char_matrix <- function(x, collapse="", sep=" "){
  x <- collapse_matrix(x, collapse = collapse, sep = sep)
  x <- sapply(x, strsplit, split = "")            # split to single chars
  names(x) <- NULL                                # for cleaner output
  do.call(rbind, x)                               # single chars matrix
}


# TODO: first line indent wrong
matrix_to_console <- function(x, sep=""){
  #cat(" ")  ???
  dummy <- apply(x, 1, function(x) 
                  cat(c(x, "\n"), collapse = "", sep = sep))
}

widths_matrix_columns <- function(x){
  apply(x, 2, function(y) max(nchar(y)))
}


random_df <- function(nrow=ncol, ncol=nrow, wrow=6, wcol=10){
  x <- data.frame(replicate(ncol, sample(1:5, nrow, replace=TRUE)))
  rownames(x) <- replicate(nrow, randomSentence(wrow))
  colnames(x) <- replicate(ncol, randomSentence(wcol))
  x
}


# @param  x     a single char matrix
# @param  left  number of empoty columns added on left side (default 0)
# @param  right number of empty columns added at right side (default 0)
# @return matrix
# @keywords internal
add_empty_cols <- function(x, left=0, right=0){
  x <- cbind(matrix(" ", nrow=nrow(x), ncol=left), x)
  x <- cbind(x, matrix(" ", nrow=nrow(x), ncol=right))
  x
}


# Binds two single character matrices of different size horizontally
#
# Two matrices in atomci format are binded horizontally at a specified
# position. The matrices need to be in single char format, i.e. one character per cell
# only. If the dimensions are different, the margins of the matrices are filled up with
# empty cells.
#
# @param  um        upper matrix (must be single char matrix)
# @param  lm        lower matrix (must be single char matrix)
# @param  anchors   two integers specifying at which columns matrices are aligned
# @return matrix
#
# @author            Mark Heckmann
# @keywords internal
# @examples \dontrun{
#   um <- matrix("u", ncol=10, nrow=5)
#   lm <- matrix("l", ncol=8, nrow=3)
#   bind_matrices_horizontally(um, lm, anchors=c(3,1))
# }
bind_matrices_horizontally <- function(um, lm, anchors=c(1,1)){
  diff.left <- diff(anchors)                  # add columns on left side
  if (diff.left <= 0) {             
    um.ncols.empty.left <- 0
    lm.ncols.empty.left <- abs(diff.left)
  } else {
    um.ncols.empty.left <- abs(diff.left)
    lm.ncols.empty.left <- 0
  }
  um <- add_empty_cols(um, left=um.ncols.empty.left)
  lm <- add_empty_cols(lm, left=lm.ncols.empty.left)
  
  diff.right <- diff(c(ncol(um), ncol(lm)))   # add columns on right side
  if (diff.right <= 0) {
    um.ncols.empty.right <- 0
    lm.ncols.empty.right <- abs(diff.right)
  } else {
    um.ncols.empty.right <- abs(diff.right)
    lm.ncols.empty.right <- 0
  }
  um <- add_empty_cols(um, right=um.ncols.empty.right)
  lm <- add_empty_cols(lm, right=lm.ncols.empty.right)
  
  rbind(um, lm)
}



# break at any point possible
break_output <- function(mat, ncolkeep=14, keeprows=TRUE){
   availchar <- options()$width          # get console size (problematic update)
   #print(availchar)
   #if (availchar < ncolkeep)             # set FALSE to avoid endless recursion
   #  keeprows <- FALSE
   if (ncol(mat) >= availchar) {
     mat.tmp <- mat[ , 1:(availchar - 1)]
     out.tmp <- collapse_matrix(mat.tmp, collapse = "")  # collapse rows
     matrix_to_console(out.tmp)                          # print first part to console
     cat("\n")           # empty line after print out to seperate prints
    # if (keeprows) {     # rownames after each pagebreak?
   #    mat.residual <- mat[ , c(1:(ncolkeep), availchar:ncol(mat))] 
    # } else {
       mat.residual <- mat[ , c(availchar:ncol(mat)), drop=FALSE]
     #}
     Recall(mat.residual)     # recursive output call            
   } else {
     out <- collapse_matrix(mat, collapse="")   # collapse rows
     matrix_to_console(out)                     # print to console
   }
}
   
trim_string <- function(vec, trim=NA){
  if (!is.na(trim))
    vec <- substr(vec, 1, trim)
  vec
}

make_sep_mat_atomic <- function(sep, nr){
  sep.atomic <- strsplit(sep, "")[[1]]
  matrix(sep.atomic, nrow = nr, 
                     ncol=nchar(sep), byrow=TRUE)
}




df_out <- function(df,                # data frame
                   left=NA,           # rows left
                   right=NA,          # rows right
                   showopt=1,         # options where to place left and right matrix
                                      # 0=none, 1 = left and right, 2= both left, 3=both right
                   just.rows="r",     # justification of row names
                   just.main="l",     # justification of body
                   max.char.rows=200, # max no of chars of row names to be printed
                   sep=" ",           # seperator symbol between columns
                   sep2="   ",        # seperator between row names and first column
                   equal=FALSE,       # equal width for columns (max column width)
                   prefix="",         # optional prefix before printed column name 
                                      # (e.g. "+---"). characters
                   keeprows=T,        # whether to show rows after each pagebreak
                   colstart="l",
                   margin=1,          # right margin for linebreak
                   trim=c(NA,NA),     # maximum number of character for r/c entry.
                   cut=c(NA, NA),     # maximal number of chars left and right of main matrix
                   id=c(T,T),         # id numbers at beginning/end of each row/column
                   hatform=FALSE)     # column names in hat form
{       
  # sanity checks
  if (length(trim) == 1)    # if only one parameter given, extend to the other
    trim <- recycle(trim, 2)
  if (length(cut) == 1)
    cut <- recycle(cut, 2)
  if (length(id) == 1)
    id <- recycle(id, 2)    
  if (!identical(left, NA) & !identical(right, NA)){
    if (length(left) != length(right))
      stop("left and right must have the same length")
    if (length(left) != nrow(df) | length(right) != nrow(df))
      stop("left and/or right must equal number of rows in df")
  }

  # main matrix mat.m 
  make_mat_main <- function(df) {
    mat.m <- sapply(df, as.character) # convert to character for type security
    rownames(mat.m) <- NULL           # unnecessary 
    colnames(mat.m) <- NULL           # unnecessary 
    mat.m <- as.matrix(mat.m)         # convert to matrix,
    if (nrow(df) == 1)                # re-transpose in single row case
      mat.m <- t(mat.m)
    nchar.column <- widths_matrix_columns(mat.m)   # no of chars per column
    if (equal) {                        # equal or dynamic column width
      mat.m <- format(mat.m, justify = just.main, width = max(nchar.column)) 
    } else {  
      mat.m <- trim_column_width(mat.m, just = just.main)
    }
    mat.m
  }
  
  # vec       vector of strings to be made as column matrix
  # idside    side at which id is attached (1=start, 2=end)
  # trim      number of chars to trim strings to
  # just      justification of text (l, c, r)
  make_mat_leftright <- function(vec, id=TRUE, idside=1, trim=NA, just="r"){
    if (!is.na(trim))               # trim rownames
      left <- substr(vec, 1, trim)
    if (id){                        # add id number to each row
      ids <- paste("(", seq_along(vec), ")", sep="")
      if (idside == 1)              # ids at start of string (for right side constructs)
        vec <- paste(ids, vec)
      else  vec <- paste(vec, ids)  # ids at end of string (for left side constructs)
    }
    vec <- format(vec, justify=just)   # justify rownames
    as.matrix(vec)
  }
  
  # make left and right matrices
  mat.left <- matrix("", nrow=nrow(df), ncol=0)     # default void matrix to start from
  mat.right <- matrix("", nrow=nrow(df), ncol=0)    # default void matrix to start from

  if (!identical(left, NA))                                 # trimming occures in all cases if prompted
    left <- trim_string(left, trim=trim[1])
  if (!identical(right, NA))                          
    right <- trim_string(right, trim=trim[1])  
  leftright <- paste(left, right, sep=" - ")        # join left and right strings

  # decision where and how to put left and right vectors
  if (showopt == 1) {              # #1 left to left, right to right
    if (!identical(left, NA))
      mat.left <- make_mat_leftright(left, id=id[1], idside=2, just="r")
    if (!identical(right, NA))
      mat.right <- make_mat_leftright(right, id=id[1], idside=1, just="l")
  } else if (showopt == 2) {       # #2 left and right on left side
    if (!identical(left, NA) & !identical(right, NA)){
      mat.left <- make_mat_leftright(leftright, id=id[1], idside=2, just="r") 
    } else if (identical(left, NA) & !identical(right, NA)) {
      mat.left <- make_mat_leftright(right, id=id[1], idside=2, just="r")
    } else if (!identical(left, NA) & identical(right, NA)) {
      mat.left <- make_mat_leftright(left, id=id[1], idside=2, just="r")
    }
  } else if (showopt == 3) {       # #3 left and right on right side
    if (!identical(left, NA) & !identical(right, NA)) {
      mat.right <- make_mat_leftright(leftright, id=id[1], idside=1, just="l")
    } else if (identical(left, NA) & !identical(right, NA)) {
      mat.right <- make_mat_leftright(right, id=id[1], idside=1, just="l")
    } else if (!identical(left, NA) & identical(right, NA)) {
      mat.right <- make_mat_leftright(left, id=id[1], idside=1, just="l")
    }      
  }  # #0 left and right unused, mat.left and mat.right remain void

  mat.m <- make_mat_main(df)
  mat.m.atomic <- matrix_to_single_char_matrix(mat.m, collapse=sep)

  mat.left.atomic <- matrix_to_single_char_matrix(mat.left)
  mat.right.atomic <- matrix_to_single_char_matrix(mat.right)

  widths.columns <- widths_matrix_columns(mat.m)   # vector column widths
  widths.sep1 <- nchar(sep)
  widths.sep2 <- nchar(sep2)
  
  # where to place colnames in matrix upper
  columns.start.r <- cumsum(widths.columns + widths.sep1) - widths.sep1
  columns.start.l <- columns.start.r - widths.columns + 1
  columns.start.cl <- columns.start.l + floor((widths.columns + 1) / 2)
  columns.start.cr <- columns.start.l + ceiling((widths.columns + 1) / 2)

  # select column start vector
  if (colstart == "r")
    columns.start <- columns.start.r else
  if (colstart == "cl")
    columns.start <- columns.start.cl else
  if (colstart == "cr")
    columns.start <- columns.start.cr else 
    columns.start <- columns.start.l
    
  # maximal rows of mat.u is length of column name plus starting position (plus prefix)
  names.columns <- colnames(df)                         # extract colnames
  if (!is.na(trim[2]))                                  # trim colnames
    names.columns <- substr(names.columns, 1, trim[2])

  ### hat = FALSE   (upper matrix u in descending form)
  if (!hatform){
    if (id[2]) {                                    # add id number to each col
       ids <- paste(seq_along(names.columns), "-", sep=" ")
       names.columns <- paste(ids, names.columns)
    }

    names.columns <- paste(prefix, names.columns, sep="") # add prefix (default "")
    ncol.mat.columns <- max(columns.start + 
                            nchar(names.columns) - 1)     # min no columns mat.u
    nrow.mat.columns <- length(names.columns) + 1
    mat.u.atomic <- matrix(" ", nrow=nrow.mat.columns,    # empty matrix
                           ncol=ncol.mat.columns)              

    # fill matrix upper
    names.atomic.list <- strsplit(names.columns, "")
    lengths.colnames <- nchar(names.columns)
    for (j in seq_along(columns.start)){  # vertical lines ("|") at column starts
      mat.u.atomic[(j + 1):nrow(mat.u.atomic), columns.start[j]] <- "|"
      mat.u.atomic[j, columns.start[j]:(columns.start[j] + 
                   lengths.colnames[j] -1)] <- names.atomic.list[[j]]
    }
    extra.cols.left <- 0                              # to suit results of hat=TRUE part
  }

  ### hat = TRUE  (upper matrix u in hat form)
  if (hatform){
    ncol <- length(names.columns)                     # no of columns
    midcol <- ceiling((ncol + 1) / 2)                 # determine middle column
    index.cols.left <- 1:(midcol - 1)                 # index of left columns
    index.cols.right <- midcol:ncol                   # index of right columns
    colnames.left <- names.columns[index.cols.left]   # left hat side
    colnames.right <- names.columns[index.cols.right] # right hat side 

    if (id[2]) {                                      # add id number to each col
      ids.left <- seq_along(names.columns)[index.cols.left]
      ids.right <- seq_along(names.columns)[index.cols.right]
      colnames.left <- paste(colnames.left, ids.left, sep=" - ")
      colnames.right <- paste(ids.right, colnames.right, sep=" - ")
    }  
    
    # add prefix to both sides (default "")
    colnames.left <- paste(colnames.left, strReverse(prefix), sep="")  # left side has revesred prefix 
    colnames.right <- paste(prefix, colnames.right, sep="")
    colnames.leftright <- c(colnames.left, colnames.right)
    lengths.colnames <- nchar(colnames.leftright)
  
    minpos <- min(columns.start[index.cols.left] - nchar(colnames.left))    # min pos to left
    maxpos <- max(columns.start[index.cols.right] + nchar(colnames.right))  # max pos to right

    if (minpos < 0 ) {
      extra.cols.left <- abs(minpos) 
    } else {
      extra.cols.left <- 0
    }
    ncol.mat.upper <- extra.cols.left + maxpos                                    # ncol of upper matrix
    nrow.mat.upper <- max(c(length(colnames.left), length(colnames.right))) + 1   # nrow of upper matrix
    mat.u.atomic <- matrix(" ", nrow=nrow.mat.upper,    # empty upper matrix to get filled
                                ncol=ncol.mat.upper)
                                  
    names.atomic.list.left <- strsplit(colnames.left, "")
    names.atomic.list.right <- strsplit(colnames.right, "")
    names.atomic.list.leftright <- c(names.atomic.list.left,
                                     names.atomic.list.right)

    # fill matrix u and build vertical lines for left and right side
    bottom.row <- nrow(mat.u.atomic)
    nc <- length(columns.start)
    columns.start.offsetted <- extra.cols.left + columns.start
    for (j in seq_along(columns.start)) {  # vertical lines ("|") at column starts
      if (j < ceiling((nc + 1)/ 2)) {
        mat.u.atomic[(bottom.row - j + 1):bottom.row, 
                      columns.start.offsetted[j]] <- "|"
        mat.u.atomic[(bottom.row - j), 
                     (columns.start.offsetted[j] - lengths.colnames[j] + 1):
                      columns.start.offsetted[j]] <- 
                      names.atomic.list.leftright[[j]]
      } else { 
        mat.u.atomic[(bottom.row - (nc - j) - 1):bottom.row, 
                      columns.start.offsetted[j]] <- "|"   
        mat.u.atomic[(bottom.row - (nc - j) - 1), columns.start.offsetted[j]:
                     (columns.start.offsetted[j] + lengths.colnames[j] - 1)] <- 
                      names.atomic.list.leftright[[j]]
      }
    }   # TODO: right side one row too much, maybe erase
  }

  # same part for both types
  mat.sep2.atomic <- make_sep_mat_atomic(sep2, nr=nrow(df))     # matrix to separate left and main, or main and right
  mat.lm.atomic <- cbind( mat.left.atomic, mat.sep2.atomic, mat.m.atomic, # lower matrix lm
                              mat.sep2.atomic, mat.right.atomic)
                                                    
  # join upper and lower matrix
  anchor.um <- extra.cols.left + 1
  anchor.lm <- ncol(mat.left.atomic) + ncol(mat.sep2.atomic) + 1
  mat.out.atomic <- bind_matrices_horizontally(mat.u.atomic, mat.lm.atomic,
                                               anchors=c(anchor.um, anchor.lm))
                                                                                     
  # cut output at sides if prompted
  diff.left <- diff(c(anchor.um, anchor.lm))
  if (diff.left <= 0) {             
    lm.empty.cols.left <- abs(diff.left)
  } else {
    lm.empty.cols.left <- 0
  }  
  start.main.at <- lm.empty.cols.left + ncol(cbind(mat.left.atomic, mat.sep2.atomic))
  end.main.at <- start.main.at + ncol(mat.m.atomic)
  
  if (!is.na(cut[1]) | !is.na(cut[2])){   
    if (is.na(cut[1])){
      end.left <- 1        
    } else {
      end.left <- trim_val(start.main.at - cut[1], minmax=c(1, 200))        
    }
    if (is.na(cut[2])){
      end.right <- ncol(mat.out.atomic)        
    } else {
      end.right <- trim_val(end.main.at  + cut[2],
                            minmax=c(1, ncol(mat.out.atomic)))        
    }
    mat.out.atomic <- mat.out.atomic[ , end.left:end.right]           
  }
  break_output(mat.out.atomic)
  invisible(NULL)
}

#df <- random_df(10, 25, wcol=4)
#left <- randomSentences(10, 5)
#right <- randomSentences(10, 5)
#df_out(df, left, right, h=T, cut=25, id=T, show=1)





###############################################################################
# repgrid show method

# @usage \S4method{show}{repgrid}(object)

# show method for repgrid class
# org <- list()
# org$show$cut <- 30
# org$show$showopt <- 1
# org$show$verbose <- TRUE

# method depends on the definition of the 'repgrid' object
# hence has to come before this code in COLLATE tag in DESCRIPTION

# @aliases show,repgrid-method

# Show method for repgrid
#
# @param object a \code{repgrid} object
# @docType methods
# @usage \S4method{show}{repgrid}(object)
# @author            Mark Heckmann
# @include repgrid.r
#

#' Show method for repgrid
#' 
#' @param object A \code{repgrid} object.
#' @include repgrid.r
#' 
setMethod("show", "repgrid", function(object){   
  pars <- settings() 
  trim <- c(pars$show.trim, pars$show.trim)   #trim <- c(30,30)
  cut <- c(pars$show.cut, pars$show.cut)      #cut <- c(20,20)  
  verbose <- TRUE    # what parts to print TRUE prints all information about the grid
  showopt <- 1
  id <- c(pars$c.no, pars$e.no)   #  c(T,T)
  hatform <- T
  
  x <- object
  do.bertin <- FALSE
  # verbose output displays all grid information available
  if (verbose){
    # print meta data
    if (pars$show.meta)  showMeta(x)
    if (pars$show.scale) showScale(x)    #print scale info
    cat("\nRATINGS:\n")
  }
  
  # make data frame for left and right constructs
  constructs <- getConstructNames(x)
  
  # make data frame for data
  df.ratings <- as.data.frame(x@ratings[ , ,1, drop=FALSE])     # extract scores
  colnames(df.ratings) <- getElementNames(x)                    # name columns
  left <- constructs[ ,1]
  right <- constructs[, 2]
  df_out(df.ratings, left, right, just.main="r", hatform=hatform, id=id, 
        trim=trim, cut=cut, equal=F, showopt=showopt)
  cat("\n")
  if (do.bertin)
    bertin(x)
})

# # Show method for repgrid
# # @param repgrid object
# setMethod("show", signature= "repgrid", function(object){
#   x <- object 
#   showMeta(x)
#   showScale(x)    #print scale info
# })



# output version for repertory grids:
# parameters
# 
# conside   integer to describe side where to print constructs
#           0 no constructs, 1 left side only, 2 both sides, 3 right side only






