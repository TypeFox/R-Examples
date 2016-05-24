# FIXME: not used anywhere?
matchDataFrameSubset = function(df, ss, factors.as.chars = TRUE) {
  checkArg(df, c("list", "data.frame"))
  checkArg(ss, c("list", "data.frame"))
  if (!isProperlyNamed(df))
    stop("'df' is not proberbly named")
  if (!isProperlyNamed(ss))
    stop("'ss' is not proberbly named")
  if (any(names(ss) %nin% names(df)))
    stop("Names of 'ss' not found in 'df'")
  if (is.list(df))
    df = as.data.frame(df, stringsAsFactors = FALSE)
  if (is.list(ss))
    ss = as.data.frame(ss, stringsAsFactors = FALSE)

  df = subset(df, select = names(ss))

  if (factors.as.chars) {
    df = convertDataFrameCols(df, factors.as.char = TRUE)
    ss = convertDataFrameCols(ss, factors.as.char = TRUE)
  }

  conv = function(x) rawToChar(serialize(x, connection = NULL, ascii = TRUE))
  match(rowSapply(ss, conv, use.names = FALSE), rowSapply(df, conv, use.names = FALSE))
}
