# create table header
.header <- function (x, rows) {
  ifelse(rows, y  <- append("", colnames(x)), y  <- colnames(x))
  paste(paste(y, collapse=' & '), '\\\\')
}

# add toprule and midrule to header
.printhead <- function (top, header, mid) {
  .printif(header, paste(top, header, mid, sep='\n'))
}

# sprintf if
.printif <- function(x, str, ...){ if (!is.null(x)) sprintf(str, x, ...) }

# if booktabs else - helper function
.book <- function(opt, booktabs, true, false){
  getOption(opt, ifelse(booktabs, true, false))
}

# create table body
.body <- function (x, rows, rowsep) {
  if(rows) x  <- cbind(rownames(x), x)
  sep  <- sprintf("%s\n", rowsep)
  paste(apply(x, 1, paste, collapse=' & '), '\\\\', collapse=sep)
}

# set column definitions
.coldef <- function(x, colsep){
  paste(ifelse(sapply(x, is.numeric), 'r', 'l'), collapse=colsep)
}

# paste table
.pt <- function(x){
  cat(paste(x, collapse='\n'), '\n')
}

# triple option
.op <- function (op1, op2, def) {
  getOption(op1, getOption(op2, def))
}
