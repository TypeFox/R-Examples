# $Id: matchcols.R 625 2005-06-09 14:20:30Z nj7w $
# select the columns which match/don't match a set of include/omit patterns.

matchcols <- function(object, with, without, method=c("and","or"), ...)
  {
    method <- match.arg(method)
    cols <- colnames(object)

    # include columns matching 'with' pattern(s)
    if(method=="and")
      for(i in 1:length(with))
        {
          if(length(cols)>0)
            cols <- grep(with[i], cols, value=TRUE, ...)
        }
    else
      if(!missing(with))
        if(length(cols)>0)
          cols <- sapply( with, grep, x=cols, value=TRUE, ...)

    # exclude columns matching 'without' pattern(s)
    if(!missing(without))
      for(i in 1:length(without))
        if(length(cols)>0)
          {
            omit <- grep(without[i], cols, ...)
            if(length(omit)>0)
              cols <- cols[-omit]
          }

    cols
  }
