### For scuo format for Drew's scuo codes.

gen.scuo <- function(seq.string, aa.names = .CF.GV$amino.acid,
    split.S = TRUE, drop.X = TRUE, drop.MW = TRUE){
  if(split.S){
    if("S" %in% aa.names){ 
      if(! "Z" %in% aa.names){
        aa.names <- c(aa.names, "Z")
      }
    } else{
      split.S <- FALSE
    }
  } else{
    if(all(c("S", "Z") %in% aa.names)){
        split.S <- TRUE
    }
  }

  if(drop.X){
    aa.names <- aa.names[aa.names != "X"]
  }

  if(drop.MW){
    aa.names <- aa.names[!(aa.names %in% c("M", "W"))]
  }

  aa.names <- sort(aa.names)
  ret <- build.scuo(seq.string, aa.names, split.S = split.S)
  ret
} # End of gen.scuo().


build.scuo <- function(seq.string, aa.names, split.S = FALSE){
  names.seq <- names(seq.string)

  if(split.S){
    synonymous.codon <- .CF.GV$synonymous.codon.split
  } else{
    synonymous.codon <- .CF.GV$synonymous.codon
  }

  ret.id <- NULL
  ret.counts <- NULL
  for(i.aa in 1:length(aa.names)){
    tmp.id <- cbind(aa.names[i.aa], names.seq)

    scodon <- synonymous.codon[[aa.names[i.aa]]]
    tmp <- lapply(1:length(seq.string),
             function(i.gene){
               tmp.ret <- rep(0L, length(scodon))
               for(i.codon in 1:length(scodon)){
                 tmp.ret[i.codon] <- sum(seq.string[[i.gene]] == scodon[i.codon])
               }
               as.integer(tmp.ret)
             })
    tmp <- do.call("rbind", tmp)

    ### The maximum is 6 possible codons.
    add.col <- 6 - ncol(tmp)
    if(add.col > 0){
      tmp <- cbind(tmp, matrix(NA, nrow = nrow(tmp), ncol = add.col))
    }

    ret.counts <- rbind(ret.counts, tmp)
    ret.id <- rbind(ret.id, tmp.id)
  }

  rownames(ret.id) <- NULL
  rownames(ret.counts) <- NULL
  ret.id <- as.data.frame(ret.id, stringsAsFactors = FALSE)
  ret.counts <- as.data.frame(ret.counts, stringsAsFactors = FALSE)

  ### amio acid, gene name, and 6 possible codons.
  ret <- cbind(ret.id, ret.counts)
  colnames(ret) <- c("AA", "ORF", paste("C", 1:6, sep = ""))

  ret
} # End of build.scuo().
