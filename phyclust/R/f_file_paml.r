### Write PAML file (one of PHYLIP format).

write.paml <- function(seqdata, filename, classid = NULL,
    seqname = NULL, width.seqname = 10, width.line = 60, lower.case = FALSE,
    code.type = .code.type[1], sep = ""){
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

  ret <- write.paml.format(seqdata, filename, classid = classid,
                           seqname = seqname, width.seqname = width.seqname,
                           width.line = width.line, sep = sep)
} # End of write.paml().

write.paml.format <- function(seqdata, filename, classid = NULL,
    seqname = NULL, width.seqname = 10, width.line = 60, sep = ""){
  n.seq <- nrow(seqdata)
  tl.seq <- ncol(seqdata)

  if(is.null(seqname)){
    seqname <- as.character(1:n.seq)
  }

  if(!is.null(classid)){
    seqname <- cbind(seqname, rep("-", n.seq) , classid)
    seqname <- apply(seqname, 1, paste, collapse = "")
  }

  tmp.seqname <- strsplit(seqname, "")
  for(i in 1:n.seq){
    tl <- length(tmp.seqname[[i]])
    if(tl > width.seqname){
      tmp.seqname[[i]] <- paste(tmp.seqname[[i]][1:width.seqname],
                                collapse = "")
    } else{
      tmp.seqname[[i]] <- paste(c(tmp.seqname[[i]],
                                  rep(" ", width.seqname - tl)), collapse = "")
    }
  }
  tmp.seqname <- do.call("c", tmp.seqname)
  tmp.space <- rep(paste(rep(" ", width.seqname), collapse = ""),
                   length(tmp.seqname))

  head <- paste(n.seq, tl.seq, collapse = " ")
  write(head, file = filename)

  tl.show <- ceiling(tl.seq / width.line)
  for(i.seq in 1:n.seq){
    write(tmp.seqname[i.seq], file = filename, append = TRUE)

    for(i in 1:tl.show){
      show.range <- (1:width.line) + (i - 1) * width.line
      if(show.range[width.line] > tl.seq){
        show.range <- show.range[1:which(show.range == tl.seq)]
      }

      tmp.seqdata <- paste(seqdata[i.seq, show.range], collapse = sep)
      write(tmp.seqdata, file = filename, append = TRUE)
    }
  }
} # End of write.paml.format().

