# takes the raw climate data as argument and checks for formatting:
# a) 4 column style as originally supported
# b) dendroclim style (12 months in one row)
# c) a list of b

climdispatch <- function(x) {
  ## is it a list?
  if (any(class(x) == "list")) {
    ## handle list case
    n <- length(x)
    minyrs <- maxyrs <- numeric(n)
    for (i in 1:n) {                    # loop through list
                                        # members, and get min and max
                                        # years for later reformatting.
      y <- x[[i]]                       # shortcut for current list
                                        # member
      if (dim(y)[2] == 13) {            # explanation see non-list
                                        # case.
        if (!any(y[,1] == seq(y[1,1], y[dim(y)[1],1], 1))) {
          stop("One member in the list supplied as climate data is not properly formatted. see '?dcc' for possible data formats.")
        } else {
          minyrs[i] <- min(y[,1])
          maxyrs[i] <- max(y[,1])
        }
      }
    }
    yrs <- max(minyrs):min(maxyrs)
    nyrs <- length(yrs)
    output.matrix <- matrix(NA, ncol = n + 2, nrow = nyrs*12)
    output.matrix[,1] <- rep(yrs, each = 12)
    output.matrix[,2] <- rep(1:12, nyrs)
    for (i in 1:n) {                    # loop through list again, and
                                        # put everything in place in
                                        # the new output matrix
      y <- x[[i]]                       # shortcut for current list
                                        # member
      for (j in 1:nyrs) {
        if (any(y[,1] == yrs[j])) {     # check, if currently selected
                                        # year is present in current
                                        # list member
          output.matrix[which(output.matrix[,1] == yrs[j]), 2+i] <-
            unlist(y[which(y[,1] == yrs[j]), 2:13]) # write elements
                                        # of specific line into i+2 th
                                        # row of output.matrix
        }                               # don't need else-case, as
                                        # matrix is already filled
                                        # with NA
      }
    }
  } else {
    ## handle non-list case
    if (dim(x)[2] == 13) {              # should have 12 months
                                        # columns and one year column
      if (!any(x[,1] == seq(x[1,1], x[dim(x)[1],1], 1))) { # check
                                        # if the first column is
                                        # perfect sequence of integer
                                        # years. if expression
                                        # evaluates to FALSE, then
                                        # this is the case, if TRUE:
                                        # stop.
        stop("unknown format of climate data")
      } else {                          # this is most probably a
                                        # dendroclim-formatted set of
                                        # climate data
        yrs <- unique(x[,1])
        nyrs <- length(yrs)
        output.matrix <- matrix(NA, ncol = 3, nrow = nyrs*12)
        output.matrix[,1] <- rep(yrs, each = 12)
        output.matrix[,2] <- rep(1:12, nyrs)
        for (i in 1:nyrs) {             # loop through years and write
                                        # respective rows in
                                        # respective columns in
                                        # output.matrix
          output.matrix[which(output.matrix[,1] == yrs[i]), 3] <-
            unlist(x[which(x[,1] == yrs[i]), 2:13])
        }
      }
    } else {                            # could still be the
                                        # originally intended format
                                        # of data.
      if (!any(x[,1] == rep(x[1,1]:x[dim(x)[1],1], each = 12))) {
                                        # check if the first column is
                                        # a perfect sequence of
                                        # integer years, each repeated
                                        # 12 times. if expression
                                        # evaluates to FALSE, then
                                        # this is the case, else stop.
        stop("unknown format of climate data")
      } else {
        if (!(any(x[,2] == rep(1:12, length(unique(x[,1])))))) {
                                        # check if the second column
                                        # is perfect sequence of 1:12
                                        # as often as there are
                                        # individual years in column
                                        # 1. if expression evaluates
                                        # to FALSE, then this is the
                                        # case, else stop.
          stop("unknown format of climate data")
        } else {
          output.matrix <- x            # pass data directly on to
                                        # pmat
        }
      }
    } 
  }
  output.matrix
}
