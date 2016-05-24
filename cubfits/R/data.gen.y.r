### For wallace format.

gen.y <- function(seq.string, aa.names = .CF.GV$amino.acid,
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
  ret <- build.y(seq.string, aa.names, split.S = split.S)

  ret <- rearrange.y(ret)
  ret
} # End of gen.y().


build.y <- function(seq.string, aa.names, split.S = TRUE){
  names.seq <- names(seq.string)

  if(split.S){
    synonymous.codon <- .CF.GV$synonymous.codon.split
  } else{
    synonymous.codon <- .CF.GV$synonymous.codon
  }

  ret <- list()
  for(i.aa in 1:length(aa.names)){
    scodon <- synonymous.codon[[aa.names[i.aa]]]
    tmp <- lapply(1:length(seq.string),
             function(i.gene){
               tmp.ret <- rep(0L, length(scodon))
               for(i.codon in 1:length(scodon)){
                 tmp.ret[i.codon] <- sum(seq.string[[i.gene]] == scodon[i.codon])
               }
               as.integer(tmp.ret)
             })
    ret[[i.aa]] <- do.call("rbind", tmp)
    colnames(ret[[i.aa]]) <- scodon
    rownames(ret[[i.aa]]) <- names.seq
  }
  names(ret) <- aa.names

  ret
} # End of build.y().

