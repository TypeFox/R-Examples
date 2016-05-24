folderh <- function(df1, key1, df2=NULL, ..., na.rm = TRUE) {
  # df1 :  list of data frames (2 or more) or a data frame with at least
  #        2 column.
  # df2 :  if df1 is a data frame, df2 is a data frame. It contains
  #        the data.
  #        If df1 is a list, df2 must be NULL.
  # keys  : string vector.
  #        - If df1 and df2 are data frames, keys is a single
  #          character string, and it is the name of the column of g
  #          and x which contains the keys.
  #        - If df1 is a list and df2 = NULL, keys is a vector of
  #          strings, and length(keys) = length(df1)-1
  #          and the keys[k] is the name of a column of df1[[k]] and
  #          also the name of a column of df1[[k+1]]. It is the
  #          keys: the name of the variable which links each line of
  #          df1[[k]] to some lines of df1[[k+1]].
  #        - If df1 and df2 are data frames, and ... contains one or
  #          more data frames, then a list of the data frames:
  #          c(list(df1, df2), as.list(...)) is built, and we do the
  #          same as when df1 is a list of data frames.
  # ... :  if df1 and df2 are data frames, it can be one or more
  #        keys and data frames, ordered so: key2, df3, key3, df4...
  #
  # The object returned will be of class: "folderh".
  #
  # An objet of class "folderh" is a list of 2 or more data frames:
  # - [[1]]: data frame of at least 2 columns. One of them, whose
  #          name is the value of "key[1]" argument, contains the
  #          groups (df1 argument).
  # - [[2]]: data frame of (p+1) columns. The name of one column is
  #          the value of "key2" argument, and this column contains
  #          the groups. The other columns contain the data
  #          (df argument).
  # - and so on.
  #
  # An object of class "folderh" has one attribute named "keys"
  # - If there are only two data frames df1 and df2: key1.
  # - If there are more than two data frames: c(key1, key2, ...).

  name.1 <- as.character(match.call())[2]
  name.2 <- as.character(match.call())[4]
  dots <- list(...)
  name.others <- as.character(match.call())[-(1:4)]
  name.dots <- name.others[name.others != "na.rm"]

  name.df <- c(name.1, name.2)

  if (!is.data.frame(df1))
    stop(name.1, " is not a data frame.")

  if (!is.data.frame(df2))
    stop(name.2, " is not a data frame.")

  X <- list(df1, df2)
  keys <- key1

  if (length(dots) > 0) {
    name.X.supp <- name.dots[seq(2, length(dots), by = 2)]
    X.supp <- dots[seq(2, length(dots), by = 2)]

    name.keys.supp <- name.dots[seq(1, length(dots), by = 2)]
    keys.supp <- unlist(dots[seq(1, length(dots), by = 2)])

    name.df <- c(name.df, name.X.supp)

    is.df <- unlist(lapply(X.supp, is.data.frame))
    if (!prod(is.df))
      warning(paste("Argument(s)", name.X.supp[!is.df],
            "is/are no data frame(s).\n", name.dots,
            "arguments will not be used."))

    is.atom <- unlist(lapply(keys.supp, is.atomic))
    if (!prod(is.atom))
      warning(paste("Argument(s)", name.keys.supp[!is.df],
            "is/are no atomic.\n", name.dots,
            "arguments will not be used."))

    if (prod(is.df) & prod(is.df))
      X <- c(X, X.supp)
      names(X) <- c(name.1, name.2, name.X.supp)
      keys <- c(keys, keys.supp)
  }

  if (length(X) > 2){
    args.fh <- c(list(df2), dots, na.rm = na.rm)
    names(args.fh)[1:3] <- c("df1", "key1", "df2")
    fh.temp <- do.call(folderh, args = args.fh)
    foldh <- append.df2folderh(fh.temp, X[[1]], keys[1], after = FALSE)
  } else {
    # Checking the values of the arguments

    df1 <- X[[1]]
    df2 <- X[[2]]

    # Are df2 and df1 data frames ?
    if (!is.data.frame(df2))
      stop(paste(name.2, "is not a data frame."))
    if (!is.data.frame(df1))
      stop(paste(name.1, "is not a data frame."))

    keyg <- keys[1]

    # Is keyg the name of a column of df1 data frame?
    if (! keyg %in% colnames(df1))
      stop(paste("There is no", keyg, "variable in", name.1))

    # Is keyg the name of a column of df2 data frame?
    if (! keyg %in% colnames(df2))
      stop(paste("There is no", keyg, "variable in", name.2))

    # Checking the number of columns of df1 and df2
    if (ncol(df1) < 2)
      stop(paste(name.1, "must be a data frame with at least two column."))
    if (ncol(df2) < 2)
      stop(paste(name.2, "must be a data frame with at least two columns."))

    # Check if each group in keyg occurs only once
    if (max(table(df1[, keyg])) > 1)
      stop(paste("A level of", name.1, "[,", keyg, "] cannot occur more than once."))

    # If nag.rm = TRUE: suppression of the lines for which df1[, keyg] is NA
    if (na.rm) {
      df1.ret <- df1[!is.na(df1[, keyg]), ]
    } else {
      df1.ret <- df1
    }

    # If nax.rm = TRUE: suppression of the lines for which df2[, keyg] is NA
    if (na.rm) {
      df2.ret <- df2[!is.na(df2[, keyg]), ]
    } else {
      df2.ret <- df2
    }

    # Creation of the folder
    foldh <- list(df1 = df1.ret, df2 = df2.ret)

    names(foldh) <- c(name.1,name.2)
    class(foldh) <- "folderh"
    attr(foldh, "keys") <- keyg
  }

  names(foldh) <- name.df

  return(foldh)
}
