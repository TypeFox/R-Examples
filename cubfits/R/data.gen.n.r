### For wallace format.

gen.n <- function(seq.string, aa.names = .CF.GV$amino.acid,
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
  ret <- build.n(seq.string, aa.names, split.S = split.S)

  ret <- rearrange.n(ret)
  ret
} # End of gen.n().


build.n <- function(seq.string, aa.names, split.S = TRUE){
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
               sum(seq.string[[i.gene]] %in% scodon)
             })
    ret[[i.aa]] <- as.integer(do.call("c", tmp))
    names(ret[[i.aa]]) <- names.seq
  }
  names(ret) <- aa.names

  ret
} # End of build.n().

