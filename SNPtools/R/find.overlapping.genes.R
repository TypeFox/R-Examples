# Find the overlap between genes and variants.
# Arguments: variants: data.frame of variante as returned by get.imputed.variants() with
#                  ID, CHROM and POS in the first three columns.
#            mgi.file: character, name of the Tabix indexed MGI gene feature file.
#            type: character, one of "gene" or "exon".  Currently only accepts
#                  "gene" and looks for SNPs within whole gene boundaries.
find.overlapping.genes = function(variants,
   mgi.file = "http://cgd.jax.org/tools/SNPtools/MGI/MGI.20130305.sorted.txt.gz",
   type = c("gene", "exon")) {

  gene.type = match.arg(type)
  if(gene.type == "exon") {
    stop("'exon' is not yet implemented.")
  } # if(gene.type == "exon")

  var.type = attr(variants, "type")

  retval = NULL
  if(!is.null(variants)) {

    chr = unique(variants$CHROM)

    if(var.type %in% c("snp", "indel")) {
      rng = range(variants$POS)
    } else if(var.type == "sv") {
      rng = range(c(variants$START, variants$END))
    } # else if(var.type == "sv")

    mgi = get.mgi.features(file = mgi.file, chr = chr, start = rng[1], 
                           end = rng[2], source = "MGI", type = type)
    m = convert.genes.to.GRanges(mgi)
    s = convert.variants.to.GRanges(variants)
    mtch = findOverlaps(s, m)
    
    retval = mgi[unique(subjectHits(mtch)),]

  } # if(is.null(variants))

  retval
} # find.overlapping.genes()

