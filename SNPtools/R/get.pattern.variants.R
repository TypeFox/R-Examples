get.pattern.variants <-
function(variants, strain.subset = NULL) {

  if(!is.null(strain.subset)) {

    type = attr(variants, "type")

    # Verify that all of the strains in the subset are in the colnames of variants.
    if(!all(strain.subset %in% colnames(variants))) {
      stop("all of the subset strains are not in colnames(variants).")
    } # if(!all(strain.subset %in% colnames(variants)))

    # Separate the header from the variants.
    hdr = variants[,1:5]
    variants = variants[,-1:-5]

    # If there are confidence values, strip them out.
    variants = strip.quality.columns(variants)

    # Make subset indices.
    ss = which(colnames(variants)  %in% strain.subset)
    mm = which(!colnames(variants) %in% strain.subset)

    # First, get all of the variants where the members of each subset have the
    # same allele.
    ss.allele = apply(variants[,ss, drop = F], 1, unique)
    mm.allele = apply(variants[,mm, drop = F], 1, unique)
    ss.allele = sapply(ss.allele, function(a) { a[!is.na(a)] })
    mm.allele = sapply(mm.allele, function(a) { a[!is.na(a)] })

    # Keep only the variants where each subset has only one allele.
    keep = which((sapply(ss.allele, length) == 1) &
                 (sapply(mm.allele, length) == 1))
    hdr  = hdr[keep,]
    variants = variants[keep,]
    ss.allele = unlist(ss.allele[keep])
    mm.allele = unlist(mm.allele[keep])

    # Keep only the variants that differ between groups.
    keep = ss.allele != mm.allele
    hdr  = hdr[keep,]
    variants = variants[keep,]

    variants = data.frame(hdr, variants, stringsAsFactors = F)
    attr(variants, "type") = type

  } # if(!is.null(strain.subset))

  return(variants)

}
