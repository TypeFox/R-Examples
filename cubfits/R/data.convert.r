### Some format converting functions, mainly from data.frame to list.

### scodon is for synonymous codon.
convert.reu13.df.to.list <- function(reu13.df){
  names.aa <- names(reu13.df)

  tmp <- NULL
  scodon <- list()
  for(i.aa in 1:length(reu13.df)){
    tmp <- c(tmp, unique(as.character(reu13.df[[i.aa]]$ORF)))
    scodon[[i.aa]] <- sort(as.character(unique(reu13.df[[i.aa]]$Codon)))
  }
  names.ORF <- sort(unique(tmp))

  ret <- list()
  for(i.gene in 1:length(names.ORF)){ # i'th gene
    ret[[i.gene]] <- list()
    for(i.aa in 1:length(reu13.df)){ # i'th amino acid
      df <- reu13.df[[i.aa]][reu13.df[[i.aa]]$ORF == names.ORF[i.gene],
                             c("Pos", "Codon")]
      ret[[i.gene]][[i.aa]] <- list()
      for(i.scodon in 1:length(scodon[[i.aa]])){ # i'th synonymous codon
        ret[[i.gene]][[i.aa]][[i.scodon]] <-
          df$Pos[df$Codon == scodon[[i.aa]][i.scodon]]
      }
      names(ret[[i.gene]][[i.aa]]) <- scodon[[i.aa]]
      # ret[[i.gene]][[i.aa]] <- split(df$Pos, df$Codon)
    }
    names(ret[[i.gene]]) <- names.aa
  }

  names(ret) <- names.ORF
  ret
} # End of convert.reu13.df.to.list().

convert.y.to.list <- function(y){
  names.aa <- names(y)
  names.ORF <- sort(unique(rownames(y[[1]])))

  ret <- list()
  for(i.gene in 1:length(names.ORF)){ # i'th gene
    ret[[i.gene]] <- list()
    for(i.aa in 1:length(y)){ # i'th amino acid
      ret[[i.gene]][[i.aa]] <- y[[i.aa]][rownames(y[[i.aa]]) ==
                                         names.ORF[i.gene],]
    }
    names(ret[[i.gene]]) <- names.aa
  }
  names(ret) <- names.ORF
  ret
} # End of convert.y.to.list().

convert.y.to.scuo <- function(y){
  names.ORF <- rownames(y[[1]])
  names.aa <- names(y)

  ret.id <- NULL
  ret.counts <- NULL
  for(i.aa in 1:length(y)){ # i'th amino acid
    tmp.id <- cbind(names.aa[i.aa], names.ORF)
    tmp.counts <- y[[i.aa]]

    ### The maximum is 6 possible codons.
    add.col <- 6 - ncol(tmp.counts)
    if(add.col > 0){
      tmp.counts <- cbind(tmp.counts,
                          matrix(NA, nrow = nrow(tmp.counts), ncol = add.col))
    }

    ret.id <- rbind(ret.id, tmp.id)
    ret.counts <- rbind(ret.counts, tmp.counts)
  }

  rownames(ret.id) <- NULL
  rownames(ret.counts) <- NULL
  ret.id <- as.data.frame(ret.id, stringsAsFactors = FALSE)
  ret.counts <- as.data.frame(ret.counts, stringsAsFactors = FALSE)

  ### amio acid, gene name, and 6 possible codons.
  ret <- cbind(ret.id, ret.counts)
  colnames(ret) <- c("AA", "ORF", paste("C", 1:6, sep = ""))
  ret
} # End of convert.y.to.scuo().

convert.n.to.list <- function(n){
  names.aa <- names(n)
  names.ORF <- sort(unique(names(n[[1]])))

  ret <- list()
  for(i.gene in 1:length(names.ORF)){ # i'th gene
    ret[[i.gene]] <- list()
    for(i.aa in 1:length(n)){ # i'th amino acid
      ret[[i.gene]][[i.aa]] <- n[[i.aa]][names(n[[i.aa]]) == names.ORF[i.gene]]
    }
    names(ret[[i.gene]]) <- names.aa
  }
  names(ret) <- names.ORF
  ret
} # End of convert.n.to.list().

convert.seq.data.to.string <- function(seq.data){
  ret <- lapply(seq.data, codon2string)
  names(ret) <- names(seq.data)
  ret
} # End of convert.seq.data.to.string().

convert.b.to.bVec <- function(b){
  do.call("c", lapply(b, function(x) x$coefficients))
} # End of convert.b.to.bVec().

convert.bVec.to.b <- function(bVec, aa.names, model = .CF.CT$model[1]){
  n.coef <- get.my.ncoef(model)
  coef.names <- get.my.coefnames(model)

  ### Check b and aa.names.
  aa.names <- aa.names[!(aa.names %in% c("M", "W", "X"))] # single and stop codons
  if("Z" %in% aa.names){
    synonymous.codon <- .CF.GV$synonymous.codon.split[aa.names]
  } else{
    synonymous.codon <- .CF.GV$synonymous.codon[aa.names]
  }
  if(length(bVec) !=
     sum(do.call("c", lapply(synonymous.codon, length)) - 1) * n.coef){
    stop("length(bVec) is not corresponding to aa.names.")
  }

  ### Compute mutation and elong
  b.Init <- list()
  id <- 0
  for(aa in aa.names){
    n.synonymous.codon <- length(synonymous.codon[[aa]]) - 1
    tl.params <- n.coef * n.synonymous.codon

    b.Init[[aa]]$coefficients <- bVec[id + (1:tl.params)]
    b.Init[[aa]]$coef.mat <- matrix(bVec[id + (1:tl.params)],
                                   nrow = n.coef, byrow = TRUE)
    colnames(b.Init[[aa]]$coef.mat) <- synonymous.codon[[aa]][1:n.synonymous.codon]
    rownames(b.Init[[aa]]$coef.mat) <- coef.names

    id <- id + tl.params
  }

  b.Init
} # End of convert.bVec.to.b().

