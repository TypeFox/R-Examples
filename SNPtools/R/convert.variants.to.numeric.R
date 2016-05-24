convert.variants.to.numeric <-
function(variants) {

  type = attr(variants, "type")

  # Strip off the header before making the conversion.
  hdr = variants[,1:5]
  variants = as.matrix(variants[,-1:-5])

  # If there are confidence values, strip them out.
  variants = strip.quality.columns(variants)
  col.names = colnames(variants)

  # Save the NA cells.
  na.cells = which(is.na(variants))

  # See which SNPs are equal to the reference allele.
  # NOTE: if the SNPs is not bimorphic, this will make it bimorphic.
  if(type == "snp") {
    variants = variants == paste(hdr$REF, hdr$REF, sep = "")
  } else if(type == "indel") {
    variants = variants == paste(hdr$REF, hdr$REF, sep = "/")
  } else if(type == "sv") {
    variants = variants == "0"
  } # else if(type == "sv")

  # Make the SNPs a numeric matrix.
  variants = matrix(as.numeric(variants), nrow(variants), ncol(variants),
             dimnames = dimnames(variants))

  # Replace the NA cells with NAs.
  variants[na.cells] = NA

  variants = data.frame(hdr, variants, stringsAsFactors = F)
  colnames(variants) = c(colnames(hdr), col.names)
  attr(variants, "type") = type

  return(variants)

}
