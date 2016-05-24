append.df2folderh <- function(fh, df, key, after = FALSE)
# fh   : folderh with tow or more data frames.
# df   : data frame to be added to the folderh fh.
# key  : name of the column of df and the first (if after = FALSE)
#        or last (if after = TRUE) element of fh containing the key.
#        - If after = FALSE: df[, key] and fh[[1]][, key] must be
#          factors with the same levels, and each one of these
#          levels must occur exactly once in fh[[1]][, key].
#        - If after = TRUE: df[, key] and fh[[length(fh)]][, key]
#          must be factors with the same levels, and each one of
#          these levels must occur exactly once in
#          fh[[length(fh)]][, key].
# after: logical.
#        - If FALSE (default), df will be added before the first
#          data frame of fh.
#        - If TRUE, df will be added after the last data frame
#          of fh.
{

  name.fh <- as.character(match.call())[2]
  name.df <- as.character(match.call())[3]

  i.add <- 1 - after + length(fh)*after

  if (!(key %in% colnames(df)))
    stop("There is no ", key, " column in ", name.df, " data frame.")

  if (!(key %in% colnames(fh[[i.add]])))
    stop("There is no ", key, " column in ", paste(names(fh)[i.add], collapse = ", "), " data frame of ", name.fh, ".")

  if (!after) {
    fh.ret <- c(list(df), fh)
    keys <- c(key, attr(fh, "keys"))
    names(fh.ret)[1] <- name.df
  } else {
    fh.ret <- c(fh, list(df))
    keys <- c(attr(fh, "keys"), key)
    names(fh.ret)[i.add+1] <- name.df
  }

  class(fh.ret) <- "folderh"
  attr(fh.ret, "keys") <- keys

  return(invisible(fh.ret))
}
