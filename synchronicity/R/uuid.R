# UUID

#' @export
uuid=function()
{
  uuid_len = 24
  if (Sys.info()['sysname'] == "Darwin") uuid_len = 15
  paste(sample(c(letters, LETTERS), uuid_len, replace=TRUE), collapse="")
}

