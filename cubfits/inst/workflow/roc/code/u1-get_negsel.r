### Convert delta t to negative selection.

get.negsel <- function(b.PM, id.intercept, id.slop, aa.names, b.label,
    b.ci.PM = NULL){
  ### Delta.t
  b.negsel.PM <- b.PM[id.slop]
  if(!is.null(b.ci.PM)){
    b.negsel.ci.PM <- matrix(b.ci.PM[id.slop,], ncol = 2)
  }

  ### log.mu
  b.logmu.PM <- b.PM[id.intercept]
  if(!is.null(b.ci.PM)){
    b.logmu.ci.PM <- matrix(b.ci.PM[id.intercept,], ncol = 2)
  }

  b.negsel.label <- b.label
  b.logmu.label <- b.label
  for(i.aa in aa.names){
    id.aa <- grep(paste("^", i.aa, "\\.", sep = ""), b.label)

    if(any(b.negsel.PM[id.aa] > 0)){
      ### Delta.t
      max.deltat <- max(b.negsel.PM[id.aa])
      id.max <- which.max(b.negsel.PM[id.aa])
      b.negsel.PM[id.aa] <- b.negsel.PM[id.aa] - max.deltat
      b.negsel.PM[id.aa][id.max] <- -max.deltat

      ### log.mu
      max.logmu <- b.logmu.PM[id.aa][id.max]
      b.logmu.PM[id.aa] <- b.logmu.PM[id.aa] - max.logmu
      b.logmu.PM[id.aa][id.max] <- -max.logmu

      ### Change C.I.
      if(!is.null(b.ci.PM)){
        if(length(id.aa) == 1){
          ### Delta.t
          b.negsel.ci.PM[id.aa,] <- -b.negsel.ci.PM[id.aa,]

          ### log.mu
          b.logmu.ci.PM[id.aa,] <- -b.logmu.ci.PM[id.aa,]
        } else{
          ### Delta.t
          tmp.ci <- matrix(b.negsel.ci.PM[id.aa,], ncol = 2)[id.max,]
          b.negsel.ci.PM[id.aa,] <- b.negsel.ci.PM[id.aa,] - max.deltat
          b.negsel.ci.PM[id.aa,][id.max,] <- -tmp.ci

          ### log.mu
          tmp.ci <- matrix(b.logmu.ci.PM[id.aa,], ncol = 2)[id.max,]
          b.logmu.ci.PM[id.aa,] <- b.logmu.ci.PM[id.aa,] - max.logmu
          b.logmu.ci.PM[id.aa,][id.max,] <- -tmp.ci
        }
      }

      ### Replace codon name.
      tmp <- .CF.GV$synonymous.codon[[i.aa]]
      if(sum(id.aa) != length(tmp)){
        tmp <- .CF.GV$synonymous.codon.split[[i.aa]]
      }
      ### Delta.t
      b.negsel.label[id.aa][id.max] <-
        paste(i.aa, ".", tmp[length(tmp)], sep = "")
      ### log.mu
      b.logmu.label[id.aa][id.max] <-
        paste(i.aa, ".", tmp[length(tmp)], sep = "")
    }
  }

  ### Return.
  if(!is.null(b.ci.PM)){
    ret <- list(b.negsel.PM = b.negsel.PM, b.negsel.ci.PM = b.negsel.ci.PM,
                b.negsel.label = b.negsel.label,
                b.logmu.PM = b.logmu.PM, b.logmu.ci.PM = b.logmu.ci.PM,
                b.logmu.label = b.logmu.label)
  } else{
    ret <- list(b.negsel.PM = b.negsel.PM,
                b.negsel.label = b.negsel.label,
                b.logmu.PM = b.logmu.PM,
                b.logmu.label = b.logmu.label)
  }

  ret
} # End of get.negsel().
