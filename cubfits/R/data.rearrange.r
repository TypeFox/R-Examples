### Some rearrangement utilities mainly to sort and rearrange data.

### Reording REU13's data frame such that in ORF order.
rearrange.reu13.df <- function(reu13.df){
  ret <- reu13.df
  for(i.aa in 1:length(ret)){ # i'th amino acid.
    ret[[i.aa]] <- ret[[i.aa]][order(as.character(ret[[i.aa]]$ORF)),]
    ret[[i.aa]]$ORF <- as.character(ret[[i.aa]]$ORF)
    ret[[i.aa]]$phi <- as.double(ret[[i.aa]]$phi)
    # ret[[i.aa]]$Amino.Acid <- as.character(ret[[i.aa]]$Amino.Acid)
    ret[[i.aa]]$Pos <- as.double(as.character(ret[[i.aa]]$Pos))
    ret[[i.aa]]$Codon <- as.character(ret[[i.aa]]$Codon)

    ### Add one more column for performance.
    if(is.null(ret[[i.aa]]$Codon.id)){
      ret[[i.aa]]$Codon.id <- vector(mode = "integer",
                                     length = nrow(ret[[i.aa]]))
      tmp <- sort(unique(ret[[i.aa]]$Codon))
      for(i in 1:length(tmp)){ # i'th synonymous codon.
        ret[[i.aa]]$Codon.id[ret[[i.aa]]$Codon == tmp[i]] <- as.integer(i - 1)
      }
    }
  }
  ret
} # End of rearrange.reu13.df().

### Reording y data frame such that in ORF order.
rearrange.y <- function(y){
  ret <- y
  for(i.aa in 1:length(ret)){ # i'th amino acid
    ret[[i.aa]] <- ret[[i.aa]][order(rownames(ret[[i.aa]])),]
  }
  ret
} # End of rearrange.y().

### Reording n data frame such that in ORF order.
rearrange.n <- function(n){
  ret <- n
  for(i.aa in 1:length(ret)){ # i'th amino acid
    tmp <- names(ret[[i.aa]])
    order.tmp <- order(tmp)
    ret[[i.aa]] <- as.integer(ret[[i.aa]][order.tmp])
    names(ret[[i.aa]]) <- tmp[order.tmp]
  }
  ret
} # End of rearrange.n().

### Reording x vector such that in ORF order.
rearrange.phi.Obs <- function(phi.Obs){
  tmp <- names(phi.Obs)
  order.tmp <- order(tmp)
  ret <- phi.Obs
  ret <- ret[order.tmp]
  names(ret) <- tmp[order.tmp]
  ret
} # End of rearrange.phi.Obs().

