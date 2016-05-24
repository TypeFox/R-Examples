.rawToString <- function(raw, ...) {
  # This approach drops all '\0', in order to avoid warnings
  # in rawToChar().  Note, it does not truncate the string after
  # the first '\0'.  However, such strings should never occur in
  # the first place.
  raw <- raw[raw != as.raw(0)]
  rawToChar(raw)
} # .rawToString()

.readRaw <- function(con, n=1L, ...) {
  readBin(con, what="raw", n=n)
}

.readByte <- function(con, n=1L, ...) {
  readBin(con, what="integer", n=n, size=1L, signed=FALSE)
}

.readUShort <- function(con, n=1L, ...) {
  readBin(con, what="integer", n=n, size=2L, signed=FALSE, endian="little")
}

.readInt <- function(con, n=1L, ...) {
  readBin(con, what="integer", n=n, size=4L, signed=TRUE, endian="little")
}

.readFloat <- function(con, n=1L, ...) {
  readBin(con=con, what="double", n=n, size=4L, signed=TRUE, endian="little")
}

.readString <- function(con, n, ...) {
  .rawToString(.readRaw(con, n=n))
}
