message("TESTING: ASCII...")

library("R.oo")

## Display ASCII table
print(R.oo::ASCII)

idxs <- 1:255
str(idxs)

chars <- intToChar(idxs)
print(chars)
stopifnot(length(chars) == length(idxs))

idxs2 <- charToInt(chars)
str(idxs2)
stopifnot(length(idxs2) == length(chars))

stopifnot(identical(idxs2, idxs))


message("TESTING: ASCII...DONE")
