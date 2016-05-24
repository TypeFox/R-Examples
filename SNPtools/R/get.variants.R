get.variants <-
function(file = "http://cgd.jax.org/tools/SNPtools/Build38/sanger.snps.NCBI38.txt.gz",
         chr, start, end, type = c("snp", "indel", "sv"), strains, polymorphic = TRUE,
         quality) {

  retval = NULL

  type = match.arg(type)

  if(missing(quality)) {
    quality = 0
  } # if(missing(quality))

  # Get the available strains.
  avail.strains = get.strains(file)
  if(missing(strains)) {
    strains = avail.strains
  } # if(is.null(strains))

  check.args(chr, start, end, strains, avail.strains)

  # Convert the start and end to bp, if neccessary.
  start = convert.pos.to.bp(start)
  end   = convert.pos.to.bp(end)

  # Create a genomic range to retrieve. 
  gr = convert.pos.to.GRanges(chr, start, end)
 
  # Open the connection to the SNP file.
  con = TabixFile(file)
  open(con)

  # Get the column names from the last row of the header info.
  hdr = headerTabix(con)
  hdr = strsplit(hdr$header, split = "\t")[[length(hdr$header)]]
  hdr = sub("^#", "", hdr)

  # Retrieve the data.
  variants = scanTabix(con, param = gr)

  # Close the connection.
  close(con)

  # Split the columns up (tab-delimited).
  variants = lapply(variants, strsplit, split = "\t")  

  # Convert the data into a matrix.  We may have gotten back no SNPs, so
  # only unpack the regions where we have SNPs.
  num.variants = sapply(variants, length)
  if(any(num.variants == 0)) {
    warning(paste("Some regions returned no variants:", paste(
            names(num.variants)[num.variants == 0], collapse = ",")))
  } # if(any(num.variants == 0))

  rng = which(num.variants > 0)

  # If we have some SNPs to process, parse them and assign alleles.
  if(length(rng) > 0) {
 
    # Select the strain columns.
    strain.columns = which(hdr %in% strains)
    if(type %in% c("snp", "indel")) {
      strain.columns = sort(c(strain.columns, strain.columns + 1))
    } else if(type == "sv") {
      strain.columns = sort(strain.columns)
    } # else if(type == "sv")

    keep = c(1:5, strain.columns)
    hdr = hdr[keep]
    for(i in rng) {
      variants[[i]] = lapply(variants[[i]], function(a) { a[keep] })
    } # for(i)

    if(type == "snp") {
      # Set the "N" alleles to NA.
      variants[[i]][variants[[i]] == "N"] = NA
    } # if(type == "snp")

    # I use two loops to reduce memory as soon as possible above.
    for(i in rng) {
      variants[[i]] = matrix(unlist(variants[[i]]), nrow = length(variants[[i]]),
                  ncol = length(variants[[i]][[1]]), byrow = T)
      variants[[i]] = data.frame(variants[[i]], stringsAsFactors = F)
      colnames(variants[[i]]) = hdr

      if(type %in% c("snp", "indel")) {
        variants[[i]]$POS = as.numeric(variants[[i]]$POS)

        # Convert the quality columns to numeric.
        qual.cols = grep("quality", colnames(variants[[i]]))
        variants[[i]][,qual.cols] = apply(variants[[i]][,qual.cols], 2, as.numeric)

        # Subset by quality score, if requested.
        if(!is.null(quality)) {
          # Get the quality columns.
          qual.cols = match(strains, colnames(variants[[i]]))
          qual.cols = qual.cols + 1

          keep = apply(apply(variants[[i]][,qual.cols], 1, ">=", quality), 2, all)
          variants[[i]] = variants[[i]][keep,]
        } # if(!is.null(quality))

      } else if(type == "sv") {
        variants[[i]]$START = as.numeric(variants[[i]]$START)
        variants[[i]]$END   = as.numeric(variants[[i]]$END)
      } # else if(type == "sv")

      # Remove non-polymorphic SNPs.
      if(nrow(variants[[i]]) > 0) {
        if(polymorphic) {
          strain.cols = match(strains, colnames(variants[[i]]))
          if(type %in% c("snp", "indel")) {
            unique.alleles = apply(variants[[i]][,strain.cols], 1, unique)
            unique.alleles = lapply(unique.alleles, function(a) { a[!is.na(a)] })
            keep = which(sapply(unique.alleles, length) != 1)
            variants[[i]] = variants[[i]][keep,]
          } else if(type == "sv") {
            # Reference alleles are all "0".
            keep = which(rowSums(variants[[i]][,strain.cols] == "0") < length(strain.cols))
            variants[[i]] = variants[[i]][keep,]
          } # else if(type == "sv")
        } # if(polymmorph)
      } # if(!is.null(variants))
    } # for(i)
  } # if(length(rng) > 0)

  # If there is only one requested region, return the data.frame, not a list.
  if(length(variants) == 1) {
    variants = variants[[1]]
  } # if(length(variants) == 1)

  # Set the class and type.
  class(variants) = c(class(variants), "variant")
  attr(variants, "type") = type

  return(variants)

} # get.variants()
