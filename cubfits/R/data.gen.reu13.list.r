### For wallace format.

gen.reu13.list <- function(seq.string, aa.names = .CF.GV$amino.acid,
    split.S = TRUE, drop.X = TRUE, drop.MW = TRUE, drop.1st.codon = TRUE){
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

  if(drop.1st.codon){
    seq.string <- lapply(seq.string, function(x){ x[-1] })
  }

  aa.names <- sort(aa.names)
  ret <- build.reu13.list(seq.string, aa.names, split.S = split.S)

  ret
} # End of gen.reu13.list().


build.reu13.list <- function(seq.string, aa.names, split.S = TRUE){
  names.seq <- names(seq.string)

  if(split.S){
    synonymous.codon <- .CF.GV$synonymous.codon.split
  } else{
    synonymous.codon <- .CF.GV$synonymous.codon
  }

  ret <- list()
  for(i.gene in 1:length(seq.string)){
    ret[[i.gene]] <- list()
    for(i.aa in 1:length(aa.names)){
      scodon <- synonymous.codon[[aa.names[i.aa]]]
      ret[[i.gene]][[i.aa]] <- list()
      for(i.codon in 1:length(scodon)){
        ret[[i.gene]][[i.aa]][[i.codon]] <-
          which(seq.string[[i.gene]] %in% scodon[i.codon])
      }
      names(ret[[i.gene]][[i.aa]]) <- scodon
    }
    names(ret[[i.gene]]) <- aa.names
  }
  names(ret) <- names.seq

  ret
} # End of build.reu13.list().

