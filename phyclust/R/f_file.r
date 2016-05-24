### Functions to read and write files.

read.fasta <- function(filename, byrow = TRUE, code.type = .code.type[1],
    aligned = TRUE, sep = ""){
  if(is.null(code.type)){
    code.type <- "UNKNOWN"
  }
  if(code.type[1] == "CODON" && sep == ""){
    sep <- ","
  }

  ret <- read.fasta.format(filename, byrow = byrow, aligned = aligned,
                           sep = sep)

  ### transfer data.
  if(code.type[1] == "NUCLEOTIDE"){
    ret$org <- code2nid(ret$org.code)
  } else if(code.type[1] == "SNP"){
    ret$org <- snp2sid(ret$org.code)
  } else if(code.type[1] == "CODON"){
    ret$org <- as.numeric(ret$org.code)
  } else if(code.type[1] == "AMINO_ACID"){
    ret$org <- acode2aid(ret$org.code)
  }

  if(any(code.type[1] %in% .code.type)){
    ret$code.type <- code.type[1]
  }

  ret
} # End of read.fasta().

read.phylip <- function(filename, byrow = TRUE, code.type = .code.type[1],
    sep = ""){
  if(is.null(code.type)){
    code.type <- "UNKNOWN"
  }
  if(code.type[1] == "CODON" && sep == ""){
    sep <- ","
  }

  ret <- read.phylip.format(filename, byrow = byrow, sep = sep)

  ### transfer data.
  if(code.type[1] == "NUCLEOTIDE"){
    ret$org <- code2nid(ret$org.code)
  } else if(code.type[1] == "SNP"){
    ret$org <- snp2sid(ret$org.code)
  } else if(code.type[1] == "CODON"){
    ret$org <- as.numeric(ret$org.code)
  } else if(code.type[1] == "AMINO_ACID"){
    ret$org <- acode2aid(ret$org.code)
  }

  if(any(code.type[1] %in% .code.type)){
    ret$code.type <- code.type[1]
  }

  ret
} # End of read.phylip().

write.fasta <- function(seqdata, filename, classid = NULL,
    seqname = NULL, width.line = 60, lower.case = FALSE,
    code.type = .code.type[1], sep = ""){
  ### matrix for aligned, list for unaligned.
  if(!(is.vector(seqdata) || is.matrix(seqdata) || is.list(seqdata))){
    stop("The seqdata should be a vector, matirx or list.")
  }
  if(is.null(code.type)){
    code.type <- "UNKNOWN"
  }

  ### transfer data.
  if(code.type[1] == "NUCLEOTIDE"){
    seqdata <- nid2code(seqdata, lower.case = lower.case)
  } else if(code.type[1] == "SNP"){
    seqdata <- sid2snp(seqdata)
  } else if(code.type[1] == "CODON" && sep == ""){
    sep <- ","
  } else if(code.type[1] == "AMINO_ACID"){
    seqdata <- aid2acode(seqdata, lower.case = lower.case)
  }

  ret <- write.fasta.format(seqdata, filename, classid = classid,
                            seqname = seqname, width.line = width.line,
                            sep = sep)
} # End of write.fasta().

write.phylip <- function(seqdata, filename, classid = NULL,
    seqname = NULL, width.seqname = 10, width.line = 60, lower.case = FALSE,
    code.type = .code.type[1], sep = ""){
#    code.type = .code.type[1], sep = "", add.space = 10){
  if(!(is.vector(seqdata) || is.matrix(seqdata)) || is.list(seqdata)){
    stop("The seqdata should be a vector or matrix. (aligned)")
  }
  if(is.null(code.type)){
    code.type <- "UNKNOWN"
  }

  ### transfer data.
  if(code.type[1] == "NUCLEOTIDE"){
    seqdata <- nid2code(seqdata, lower.case = lower.case)
  } else if(code.type[1] == "SNP"){
    seqdata <- sid2snp(seqdata)
  } else if(code.type[1] == "CODON" && sep == ""){
    sep <- ","
  } else if(code.type[1] == "AMINO_ACID"){
    seqdata <- aid2acode(seqdata, lower.case = lower.case)
  }

  ret <- write.phylip.format(seqdata, filename, classid = classid,
                             seqname = seqname, width.seqname = width.seqname,
                             width.line = width.line, sep = sep)
#                             add.space = add.space)
} # End of write.phylip().

