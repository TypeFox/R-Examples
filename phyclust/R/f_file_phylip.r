### Read PHYLIP file.
### Return a list contains
###   nseq: Number of sequence.
###   seqlen: Sequence length.
###   seqname: Sequence's name.
###   org: Original sequence in nid. array[nseq, nsite]

read.phylip.format <- function(filename, byrow = TRUE, sep = ""){
  ret <- list(code.type = "UNKNOWN", info = NULL, nseq = NULL,
              seqlen = NULL, seqname = NULL, org.code = NULL, org = NULL,
              byrow = byrow, aligned = TRUE)

### Read header.
  ret$info <- readLines(filename, n = 1)
  tmp <- unstrsplit(ret$info, " ")
  tmp <- tmp[tmp != ""]
  ret$nseq <- as.numeric(tmp[1])
  ret$seqlen <- as.numeric(tmp[2])

### Read data and transfer to nid.
  op.org <- options("stringsAsFactors")
  options(stringsAsFactors = FALSE)
  tmp <- read.table(filename, sep = "", quote = "", skip = 1, fill = TRUE)
  options(op.org)

### Split the data by reading blocks and rejoin them by sequences.
  ret$org.code <- split(tmp, gl(nrow(tmp) / ret$nseq, ret$nseq))
  ret$org.code <- do.call("cbind", ret$org.code)
  ret$org.code <- as.matrix(ret$org.code, nrow = ret$nseq)
  ret$seqname <- ret$org.code[, 1]
  ret$org.code <- apply(as.matrix(ret$org.code[, -1], nrow = ret$nseq), 1,
                        function(y){ unstrsplit(y, sep) })
  if(byrow){
    ret$org.code <- matrix(ret$org.code, nrow = ret$nseq,
                           ncol = ret$seqlen, byrow = byrow)
  } else{
    ret$org.code <- matrix(ret$org.code, ncol = ret$nseq,
                           nrow = ret$seqlen, byrow = byrow)
  }
  ret$org <- ret$org.code

  class(ret) <- "seq.data"
  ret
} # End of read.phylip.format().


write.phylip.format <- function(seqdata, filename, classid = NULL,
    seqname = NULL, width.seqname = 10, width.line = 60, sep = ""){
#    add.space = 10){
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

#   my.paste <- function(x, collapse, add.space){
#     tl.x <- length(x)
#     x.new <- NULL
#     tl.segment <- ceiling(tl.x / add.space)
#     for(i in 1:tl.segment){
#       get.range <- (i - 1) * add.space + 1:add.space
#       x.new <- c(x.new,
#                  paste(x[get.range[get.range <= tl.x]], collapse = collapse))
#     }
#     paste(x.new, collapse = " ")
#   }

  tl.show <- ceiling(tl.seq / width.line)
  for(i in 1:tl.show){
    show.range <- (1:width.line) + (i - 1) * width.line
    if(show.range[width.line] > tl.seq){
      show.range <- show.range[1:which(show.range == tl.seq)]
    }

#     tmp.seqdata <- apply(matrix(seqdata[, show.range], nrow = n.seq), 1,
#                          my.paste, collapse = sep, add.space = add.space)
    tmp.seqdata <- apply(matrix(seqdata[, show.range], nrow = n.seq), 1,
                         paste, collapse = sep)

    if(i == 1){
      tmp.seqdata <- cbind(tmp.seqname, tmp.seqdata)
    } else{
      tmp.seqdata <- cbind(tmp.space, tmp.seqdata)
    }
    tmp.seqdata <- apply(tmp.seqdata, 1, paste, collapse = "")

    write.table(tmp.seqdata, file = filename, append = TRUE, col.names = FALSE,
                row.names = FALSE, sep = "", quote = FALSE)
    write("", file = filename, append = TRUE)
  }
} # End of write.phylip.format().

