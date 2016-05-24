create.pos <-
function(chr, start, end, pos) {

  if(is.null(chr)) {
    pos = list(cbind(start, end))
    names(pos) = chr
  } # if(is.null(chr))

  return(pos)
}