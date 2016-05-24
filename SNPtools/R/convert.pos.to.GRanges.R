convert.pos.to.GRanges <-
function(chr, start, end) {

  return(GRanges(seqnames = chr, ranges = IRanges(start = start, end = end)))

}
