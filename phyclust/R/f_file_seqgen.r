### Read seq-gen output.
### Return a list contains
###   nseq: Number of sequence.
###   seqlen: Sequence length.
###   seqname: Sequence's name.
###   org: Original sequence in nid. array[nseq, nsite]
read.seqgen <- function(text, byrow = TRUE, code.type = .code.type[1]){
  phylip <- list(code.type = code.type[1], info = NULL, nseq = NULL,
                 seqlen = NULL, seqname = NULL, org.code = NULL, org = NULL,
                 byrow = byrow)

### Read header.
  phylip$info <- text[1]
  tmp <- unstrsplit(phylip$info, " ")
  tmp <- tmp[tmp != ""]
  phylip$nseq <- as.numeric(tmp[1])
  phylip$seqlen <- as.numeric(tmp[2])

### Read data and transfer to nid.
  org <- do.call("rbind", lapply(text[2:length(text)], unstrsplit, ""))
  org <- as.matrix(org, nrow = phylip$nseq)
  org.names <- apply(as.matrix(org[, 1:10], nrow = phylip$nseq), 1,
                     paste, collapse = "")
  org.names <- gsub("  *", "", org.names)
  org <- org[, -(1:10)]

### Split the data by reading blocks and rejoin them by sequences.
  phylip$seqname <- org.names
  if(code.type[1] == "NUCLEOTIDE"){
    phylip$org.code <- matrix(org, nrow = phylip$nseq, ncol = phylip$seqlen)
    phylip$org <- code2nid(org)
  } else if(code.type[1] == "SNP"){
    org <- code2snp(org)
    phylip$org.code <- matrix(org, nrow = phylip$nseq, ncol = phylip$seqlen)
    phylip$org <- snp2sid(org)
  } else{
    stop("The code.type is not found.")
  }

  if(!byrow){
    phylip$org.code <- t(phylip$org.code)
    phylip$org <- t(phylip$org)
  }
  phylip$byrow <- byrow

  class(phylip) <- "seq.data"
  phylip
} # End of read.seqgen().


print.seq.data <- function(x, ...){
  seq.data <- x
  cat("code.type: ", seq.data$code.type,
      ", n.seq: ", seq.data$nseq,
      ", seq.len: ", seq.data$seqlen,
      ", byrow: ", seq.data$byrow,
      ".\n", sep = "")
}

