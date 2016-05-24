### Read FASTA file.
### Return a list contains
###   nseq: Number of sequence.
###   seqlen: Sequence length.
###   seqname: Sequence's name.
###   org: Original sequence in nid. array[nsite, nseq]

read.fasta.format <- function(filename, byrow = TRUE, aligned = TRUE,
    sep = ""){
  ret <- list(code.type = "UNKNOWN", info = NULL, nseq = NULL,
              seqlen = NULL, seqname = NULL, org.code = NULL, org = NULL,
              byrow = byrow, aligned = aligned)

  tmp <- readLines(filename)
  tmp.org <- NULL
  seqname <- NULL

  ### parse data.
  i <- 1
  nseq <- 0
  tmp.i <- NULL
  repeat{
    if(regexpr("^>", tmp[i]) == 1){
      if(!is.null(tmp.i)){
#        ret$org.code[[nseq]] <- unlist(strsplit(tmp.org, sep))
        ret$org.code[[nseq]] <- unlist(strsplit(tmp[tmp.i], sep))
        tmp.i <- NULL
      }

      if(regexpr("\\|", tmp[i]) == 1){
        tmp.seqname <- unlist(strsplit(tmp[i], "\\|"))[1]
        seqname <- c(seqname, gsub(">(.*)", "\\1", tmp.seqname))
      } else{
        tmp.seqname <- unlist(strsplit(tmp[i], " "))
        tmp.id <- which(nchar(tmp.seqname) > 1)[1]
        if(tmp.id == 1){
          seqname <- c(seqname, gsub(">(.*)", "\\1", tmp.seqname[tmp.id]))
        } else{
          seqname <- c(seqname, tmp.seqname[tmp.id])
        }
      }

      nseq <- nseq + 1
    } else if(tmp[i] %in% c("", " ")){
    } else{
       tmp.i <- c(tmp.i, i)
#      tmp.org <- paste(tmp.org, tmp[i], sep = "")
    }

    i <- i + 1
    if(i > length(tmp)){
#      ret$org.code[[nseq]] <- unlist(strsplit(tmp.org, sep))
      ret$org.code[[nseq]] <- unlist(strsplit(tmp[tmp.i], sep))
      break
    }
  }

  ### check if the inputs are aligned sequences.
  flag.aligned <- aligned
  tl.seq <- do.call("c", lapply(ret$org.code, length))
  if(nseq > 1){
    if(any(tl.seq != tl.seq[1])){
      flag.aligned <- FALSE
      ret$aligned <- FALSE
    }
    tl.seq <- max(tl.seq)
  }

  ### prepare for return.
  ret$nseq <- nseq
  ret$seqlen <- tl.seq
  ret$seqname <- seqname
  if(flag.aligned){
    if(byrow){
      ret$org.code <- do.call("rbind", ret$org.code)
    } else{
      ret$org.code <- do.call("cbind", ret$org.code)
    }
  }
  ret$org <- ret$org.code

  class(ret) <- "seq.data"
  ret
} # End of read.fasta.format().


write.fasta.format <- function(seqdata, filename, classid = NULL,
    seqname = NULL, width.line = 60, sep = ""){
  if(is.vector(seqdata) && (! is.list(seqdata))){
    seqdata <- matrix(seqdata, nrow = 1)
  }

  ### matrix for aligned, list for unaligned.
  if(is.matrix(seqdata)){
    n.seq <- nrow(seqdata)
  } else if(is.list(seqdata)){
    n.seq <- length(seqdata)
  } else{
    stop("The seqdata should be a vector, matirx or list.")
  }

  if(is.null(seqname)){
    seqname <- as.character(1:n.seq)
  }

  if(!is.null(classid)){
    seqname <- cbind(seqname, rep("-", n.seq) , classid)
    seqname <- apply(seqname, 1, paste, collapse = "")
  }

  if(is.matrix(seqdata)){
    tl.seq <- ncol(seqdata)
    tl.show <- ceiling(tl.seq / width.line)
    show.range <- list()
    for(i in 1:tl.show){
      show.range[[i]] <- (1:width.line) + (i - 1) * width.line
      if(show.range[[i]][width.line] > tl.seq){
        show.range[[i]] <- show.range[[i]][1:which(show.range[[i]] == tl.seq)]
      }
    }

    ret <- NULL
    for(i in 1:n.seq){
      ret <- c(ret, paste(">", seqname[i], sep = ""))
      for(j in 1:tl.show){
        ret <- c(ret, paste(seqdata[i, show.range[[j]]], collapse = sep))
      }
    }
  } else{
    ret <- NULL
    for(i in 1:n.seq){
      tl.seq <- length(seqdata[[i]])
      tl.show <- ceiling(tl.seq / width.line)
      show.range <- list()
      for(j in 1:tl.show){
        show.range[[j]] <- (1:width.line) + (j - 1) * width.line
        if(show.range[[j]][width.line] > tl.seq){
          show.range[[j]] <- show.range[[j]][1:which(show.range[[j]] == tl.seq)]
        }
      }

      ret <- c(ret, paste(">", seqname[i], sep = ""))
      for(j in 1:tl.show){
        ret <- c(ret, paste(seqdata[[i]][show.range[[j]]], collapse = sep))
      }
    }
  }

  write.table(ret, file = filename, quote = FALSE, sep = "",
              row.names = FALSE, col.names = FALSE)
} # End of write.fasta.format().

