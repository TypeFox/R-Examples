get.mgi.features = function(file = "http://cgd.jax.org/tools/SNPtools/MGI/MGI.20130305.sorted.txt.gz",
  chr = NULL, start = NULL, end = NULL,
  source = c("all", "MGI", "VEGA", "ENSEMBL", "Blat", "NCBI_Gene"),
  type = c("all", "gene", "pseudogenic_transcript", "pseudogenic_exon", "pseudogene",
           "match", "match-part", "transcript", "exon", "mRNA", "five_prime_UTR",
           "start_codon", "CDS", "stop_codon", "three_prime_UTR",
           "pseudogenic_mRNA", "pseudogenic_start_codon",
           "pseudogenic_CDS", "pseudogenic_stop_codon",
           "pseudogenic_five_prime_UTR", "pseudogenic_three_prime_UTR",
           "sequence_feature")) {

  # Error checking.
  check.args(chr, start, end, "", "")

  source = match.arg(source, several.ok = T)
  type   = match.arg(type, several.ok = T)

  # Convert the start and end to bp, if neccessary.
  start = convert.pos.to.bp(start)
  end   = convert.pos.to.bp(end)

  # Create a genomic range to retrieve. 
  gr = convert.pos.to.GRanges(chr, start, end)

  # Open the connection to the locally stored Sanger SNP VCF file.
  con = TabixFile(file)
  open(con)

  # Get the column names from the last row of the header info.
  hdr = headerTabix(con)

  # Retrieve the data.
  feat = scanTabix(con, param = gr)

  # Close the connection.
  close(con)

  # Split the columns up (tab-delimited).
  feat = lapply(feat, strsplit, split = "\t")

  # Convert the data into data.frames.  We may have gotten back no features,
  # so only unpack the regions where we have features.
  num.feat = sapply(feat, length)
  if(any(num.feat == 0)) {
    warning(paste("Some regions returned no features:", paste(
            names(num.feat)[num.feat == 0], collapse = ",")))
  } # if(any(num.feat == 0))

  rng = which(num.feat > 0)
 
  # If we have some features to process, parse them and keep the type
  # requested by the user.
  if(length(rng) > 0) {
    # Convert the features to a data.frame.
    for(i in rng) {
      feat[[i]]= matrix(unlist(feat[[i]]), nrow = length(feat[[i]]),
                        ncol = length(feat[[i]][[1]]), byrow = T,
                        dimnames = list(NULL, c("seqid", "source", "type",
                        "start", "stop", "score", "strand", "phase",
                        "ID", "Name", "Parent", "Dbxref", "mgiName",
                        "bioType")))
      feat[[i]] = data.frame(feat[[i]], stringsAsFactors = F)
      feat[[i]]$start = as.numeric(as.character(feat[[i]]$start))
      feat[[i]]$stop  = as.numeric(as.character(feat[[i]]$stop))
    } # for(i)

    # Keep only the requested features.
    if(source[1] != "all") {
      for(i in rng) {
        feat[[i]] = feat[[i]][feat[[i]]$source %in% source,]
      } # for(i in rng)
    } # if(source != "all")

    if(type[1] != "all") {
      for(i in rng) {
        feat[[i]] = feat[[i]][feat[[i]]$type %in% type,]
      } # for(i in rng)
    } # if(type != "all")
  } # if(length(rng) > 0)

  if(length(feat) == 1) {
    feat = feat[[1]]
  } # if(length(feat) == 1)

  return(feat)
}
