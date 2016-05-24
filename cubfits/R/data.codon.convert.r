### This converts between lower cases and uppler cases.
codon.low2up <- function(x){
  n.set <- c("a", "c", "g", "t")
  N.set <- c("A", "C", "G", "T")
  for(i.1 in 1:4){
    for(i.2 in 1:4){
      for(i.3 in 1:4){
        codon <- paste(n.set[i.1], n.set[i.2], n.set[i.3], sep = "")
        CODON <- paste(N.set[i.1], N.set[i.2], N.set[i.3], sep = "")
        x[x == codon] <- CODON
      }
    }
  }
  x
} # End of codon.low2up().

codon.up2low <- function(x){
  n.set <- c("a", "c", "g", "t")
  N.set <- c("A", "C", "G", "T")
  for(i.1 in 1:4){
    for(i.2 in 1:4){
      for(i.3 in 1:4){
        codon <- paste(n.set[i.1], n.set[i.2], n.set[i.3], sep = "")
        CODON <- paste(N.set[i.1], N.set[i.2], N.set[i.3], sep = "")
        x[x == CODON] <- codon
      }
    }
  }
  x
} # End of codon.up2low().


### Extension version of c2s.
codon2string <- function(i.seq.data){
  ### Check if i.seq.data has length of correct codon triplets.
  tl.i.seq.data <- length(i.seq.data)
  if(tl.i.seq.data %% 3 != 0 || tl.i.seq.data < 3){
    stop(paste("length of", seqinr::getName.SeqFastadna(i.seq.data),
               "is not a multiple of 3."))
  }

  # ret <- lapply(0:(tl.i.seq.data / 3 - 1),
  #               function(i) c2s(i.seq.data[i * 3 + (1:3)]))
  dim(i.seq.data) <- c(3, tl.i.seq.data / 3)
  ret <- lapply(1:(tl.i.seq.data / 3),
                function(i){ paste(i.seq.data[, i], collapse = "") })
  ret <- do.call("c", ret)

  ret
} # End of codon2string().



### This converts DNA lower case to uppler case.
dna.low2up <- function(x){
  n.set <- c("a", "c", "g", "t")
  N.set <- c("A", "C", "G", "T")
  for(i in 1:4){
    x[x == n.set[i]] <- N.set[i]
  }
  x
} # End of dna.low2up().

dna.up2low <- function(x){
  n.set <- c("a", "c", "g", "t")
  N.set <- c("A", "C", "G", "T")
  for(i in 1:4){
    x[x == N.set[i]] <- n.set[i]
  }
  x
} # End of dna.low2up().
