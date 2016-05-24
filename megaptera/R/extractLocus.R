# EXTRACT ONE LOCUS FROM AN ANNOTATED GENOME
# package: megaptera
# called by: stepAP
# author: Christoph Heibl
# last change 2014-08-28

extractLocus <- function(xml, gi, locus){
  
  GBgi <- paste("gi", gi, sep = "|")
  GBlocus <- sql.wrap(locus@aliases, regex = FALSE, BOOL = "or", 
                      term = "GBQualifier_value")
  xpath <- paste("//GBSeq[GBSeq_other-seqids/GBSeqid='", GBgi, "']",
                 "//GBFeature[GBFeature_key='gene']", # or GBFeature_key='tRNA'
                 "//GBQualifier[", GBlocus, 
                 "]/../../GBFeature_intervals/GBInterval/",
                 "GBInterval_", c("from", "to"), 
                 sep = "")
  fromto <- data.frame(
    from = as.numeric(xpathSApply(xml, xpath[1], xmlValue)),
    to = as.numeric(xpathSApply(xml, xpath[2], xmlValue)))
  
  ## detect some caveats
  ## -------------------
  if (  nrow(fromto) == 0 ){
    loci <- paste("//GBSeq[GBSeq_other-seqids/GBSeqid='", GBgi, "']",
                 "//GBQualifier[GBQualifier_name='gene']/GBQualifier_value", 
                 sep = "")
    loci <- sort(unique(xpathSApply(xml, loci, xmlValue)))
    stop("failed to extract locus='", locus@aliases[1], "' from gi=", gi, 
         "\navailable loci: ", paste(loci, collapse = ", "), sep = "")
  }
  
  ## extract locus; convert to complement, if necessary
  ## --------------------------------------------------
  dna <- paste("//GBSeq[GBSeq_other-seqids/GBSeqid='", GBgi, "']",
               "//GBSeq_sequence", 
               sep = "")
  dna <- xpathSApply(xml, dna, xmlValue)
  
  cutSequence <- function(seq, pos){
    seq <- substr(seq, min(pos), max(pos))
    if ( pos[1] > pos[2] ) {
      seq <- unlist(strsplit(seq, ""))
      seq <- comp(seq)
      seq <- paste(seq, collapse = "")
    }
    seq
  }
  dna <- apply(fromto, 1, cutSequence, seq = dna)
  paste(dna, collapse = "")
}