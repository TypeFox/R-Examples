# Given a set of genes from the get.MGI.features() function, convert them
# to GRanges.
convert.genes.to.GRanges <- 
function(mgi) {

  if(is.null(mgi)) {
    stop("convert.genes.to.GRanges: mgi cannot be null.")
  } # if(is.null(mgi))

  gr = GRanges(seqnames = mgi$seqid, ranges = IRanges(start = mgi$start,
         end = mgi$stop))
  metadata(gr) = list(Name = mgi$Name)

  return(gr)
}