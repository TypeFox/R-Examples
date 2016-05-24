is.prefix <-
function(prefixes, sq) {
  any(sapply(prefixes, function(x) identical(sq[1:length(x)], x)))
}