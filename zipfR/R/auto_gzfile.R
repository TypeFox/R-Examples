##
## internal: automatically read/write compressed files based on filename extension (.gz)
##
## This function returns either <filename> itself for an uncompressed file, or an unopened
## gzfile() connection if <filename> ends in ".gz" (case-insensitive).
## It is safe to use file=auto.gzfile(filename) in read.table() and write.table().

auto.gzfile <- function (filename) {
  stopifnot(length(filename) == 1, is.character(filename))
  if (length(grep("\\.gz$", filename, ignore.case=TRUE)) > 0) {
    gzfile(filename)  # return connection to gzip-compressed file, which can be used for reading or writing
  }
  else {
    filename  # return filename (read.table and write.table will automatically open this file)
  }
}
