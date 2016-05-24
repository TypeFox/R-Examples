convert.variants.to.GRanges <-
function(variants) {

  if(is.null(variants)) {
    return(NULL)
  } # if(is.null(variants))

  type = attr(variants, "type")

  gr = NULL
  if(type %in% c("snp", "indel")) {
    gr = GRanges(seqnames = variants$CHR, ranges = IRanges(start = variants$POS,
           end = variants$POS))
    metadata(gr) = variants[,-1:-2]
    attr(gr, "type") = type
  } else if(type == "sv") {
    gr = GRanges(seqnames = variants$CHR, ranges = IRanges(start = variants$START,
         end = variants$END))
    metadata(gr) = variants[,-1:-2]
    attr(gr, "type") = type
  } # else

  return(gr)

}