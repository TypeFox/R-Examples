bugs.log <- function (file)
{
  # Extracts the summary statistics from log.txt written by OpenBUGS.
  # Jouni Kerman 2007-01-30
  # does essentially the same thing as bugs.log() but
  # - won't crash if DIC is not there
  # - makes fewer assumptions about the structure of the matrix

  if(!file.exists(file))
    stop("Log file", file, "does not exist")
  log.txt <- readLines(file, warn=FALSE)
  extract <- function (m, line.match, skip=0, empty.left.col=TRUE) {
    start <- (skip + which(m == line.match)[1])
    if(is.na(start)) return(NA)
    if(length(start) < 1) return(NA)
    mx <- strsplit(m[-(1:start)], "\t")
    n.cols <- length(mx[[1]])
    if(n.cols < 1) return(NA)
    mxlen <- sapply(mx, length)
    end <- which(mxlen!=n.cols)[1] - 1
    if(!is.na(end)){  ### if output block does not end mx
        mx <- mx[1:end]
    }
    cm <- matrix(unlist(mx), ncol=n.cols, byrow=TRUE) # character format
    if(empty.left.col) cm <- cm[,-1]    # empty column
    col.names <- cm[1, -1]              # first column is just "node"
    row.names <- cm[,1][-1]             # first row is just ""
    col.names <- gsub("[[:space:]]+", "", col.names) # get rid of spaces
    cm <- cm[-1,-1] # delete dimname row/columns
    m <- matrix(as.numeric(cm), nrow=nrow(cm)) # convert to numeric
    dimnames(m) <- list(row.names, col.names)
    return(m)
  }
  stats <- extract(log.txt, "Summary statistics")
  DIC <- extract(log.txt, "Deviance information", skip=0, empty.left.col=FALSE)
  list(stats=stats, DIC=DIC)
}
