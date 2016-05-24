###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
na.omit.data.frame.mvmeta <-
function (object, ...) {
#
################################################################################
# FUNCTION TO HANDLE MISSING
#
  n <- length(object)
  omit <- FALSE
  omit2 <- TRUE
  vars <- seq_len(n)
  if(!is.null(y <- model.response(object))) vars <- vars[-1]
#
  # FIRST IN PREDICTORS (AS USUAL, TRUE IF AT LEAST ONE IS MISSING)
  for(j in vars) {
    x <- object[[j]]
    if (!is.atomic(x)) next
    x <- is.na(x)
    d <- dim(x)
    if (is.null(d)||length(d)!=2L) 
      omit <- omit | x
    else for(ii in 1L:d[2L]) omit <- omit|x[,ii]
  }
#
  # THEN IN RESPONSE (TO TRUE ONLY IF ALL ARE MISSING)
  if(!is.null(y)) {
    vars <- vars[-1]
    if (!is.atomic(y)) next
    y <- is.na(y)
    d <- dim(y)
    if (is.null(d)||length(d)!=2L) 
      omit2 <- omit2 & y
    else for(ii in 1L:d[2L]) omit2 <- omit2&y[,ii]
  }
#
  # EXCLUDE THE MISSING AND SET THE na.action ATTRIBUTE 
  omit <- omit|omit2
  xx <- object[!omit,,drop=FALSE]
  if (any(omit>0L)) {
    temp <- seq(omit)[omit]
    names(temp) <- attr(object,"row.names")[omit]
    attr(temp,"class") <- "omit"
    attr(xx,"na.action") <- temp
  }
#
  xx
}
