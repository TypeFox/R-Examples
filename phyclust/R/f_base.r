### Object Name (Type): Definaiton (Examples)
###   codeseq (string): A nucleotide sequence. (ACGT...)
###   nidseq (vector[L]): A nucleotide id sequence. (1234....)
###                       The length should be L = 2 + 3 * n + 3 for n >= 1.


### Transfer an nucleotide sequence to a code id sequence.
### Eg. ACGT-... => 01234...
code2nid <- function(codeseq){
  UseMethod("code2nid")
} # End of code2nid().

code2nid.default <- function(codeseq){
  code.char <- as.character(.nucleotide$code)
  code.l.char <- as.character(.nucleotide$code.l)
  code.nid <- .nucleotide$nid

  nidseq <- codeseq
  for(i in 1:nrow(.nucleotide)){
    nidseq[codeseq %in% c(code.char[i], code.l.char[i])] <- code.nid[i]
  }
  nidseq[! codeseq %in% c(code.char, code.l.char)] <- .missing.code$mid
  nidseq <- as.numeric(nidseq)

  if(is.matrix(codeseq)){
    dim(nidseq) <- dim(codeseq)
  }
  nidseq
} # End of code2nid.default().

code2nid.list <- function(codeseq){
  lapply(codeseq, code2nid.default)
} # End of code2nid.list().


### Transfer a nucleotide id sequence to an nucleotide sequence.
### Eg. 01234... => ACGT-...
nid2code <- function(nidseq, lower.case = TRUE){
  UseMethod("nid2code")
} # End of nid2code().

nid2code.default <- function(nidseq, lower.case = TRUE){
  if(lower.case){
    code.char <- as.character(.nucleotide$code.l)
  } else{
    code.char <- as.character(.nucleotide$code)
  }
  code.nid <- .nucleotide$nid

  codeseq <- nidseq
  for(i in 1:nrow(.nucleotide)){
    codeseq[nidseq == code.nid[i]] <- code.char[i]
  }
  codeseq[!nidseq %in% code.nid] <- .missing.code$code

  if(is.matrix(nidseq)){
    dim(codeseq) <- dim(nidseq)
  }
  codeseq
} # End of nid2code.default().

nid2code.list <- function(nidseq, lower.case = TRUE){
  n.seq <- length(nidseq)
  n.lower.case <- length(lower.case)

  if(n.lower.case == 1){
    lower.case <- rep(lower.case, n.seq)
  } else if(n.lower.case != n.seq){
    stop("The length of lower.case is not correct.")
  }

  lapply(1:n.seq,
    function(i){ nid2code.default(nidseq[[i]], lower.case = lower.case[i]) })
} # End of nid2code.list().


### Transfer an SNP sequence to a SNP id sequence.
### Eg. 1212-... => 01012...
snp2sid <- function(snpseq){
  UseMethod("snp2sid")
} # End of snp2sid().

snp2sid.default <- function(snpseq){
  code.char <- as.character(.snp$code)
  code.sid <- .snp$sid

  sidseq <- snpseq
  for(i in 1:nrow(.snp)){
    sidseq[snpseq == code.char[i]] <- code.sid[i]
  }
  sidseq[! snpseq %in% code.char] <- .missing.code$mid
  sidseq <- as.numeric(sidseq)

  if(is.matrix(snpseq)){
    dim(sidseq) <- dim(snpseq)
  }
  sidseq
} # End of snp2sid.default().

snp2sid.list <- function(snpseq){
  lapply(snpseq, snp2sid.default)
} # End of snp2sid.list().


### Transfer a SNP id sequence to an SNP sequence.
### Eg. 01012... => 1212-..
sid2snp <- function(sidseq){
  UseMethod("sid2snp")
} # End of sid2snp().

sid2snp.default <- function(sidseq){
  code.char <- as.character(.snp$code)
  code.sid <- .snp$sid

  snpseq <- sidseq
  for(i in 1:nrow(.snp)){
    snpseq[sidseq == code.sid[i]] <- code.char[i]
  }
  snpseq[!sidseq %in% code.sid] <- .missing.code$code

  if(is.matrix(sidseq)){
    dim(snpseq) <- dim(sidseq)
  }
  snpseq
} # End of sid2snp.default().

sid2snp.list <- function(sidseq){
  lapply(sidseq, sid2snp.default)
} # End of sid2snp.list().


### Transfer a nucleotide seqeunce to a SNP sequence.
### Eg. AGCT-... => 1122-...
code2snp <- function(codeseq){
  UseMethod("code2snp")
} # End of code2snp().

code2snp.default <- function(codeseq){
  snpseq <- codeseq
  snpseq[codeseq %in% c("A", "a", "G", "g")] <- "1"
  snpseq[codeseq %in% c("C", "c", "T", "t")] <- "2"
  # snpseq[codeseq %in% c("-")] <- "-"
  snpseq[!(snpseq %in% c("1", "2", "-"))] <- "."

  if(is.matrix(codeseq)){
    dim(snpseq) <- dim(codeseq)
  }
  snpseq
} # End of code2snp.default().

code2snp.list <- function(codeseq){
  lapply(codeseq, code2snp.default)
} # End of code2snp.list().


### Transfer a SNP seqeunce to a nucleotide sequence.
### Eg. 1122-... => AACC-...
snp2code <- function(snpseq, half = TRUE){
  UseMethod("snp2code")
} # End of snp2code().

snp2code.default <- function(snpseq, half = TRUE){
  codeseq <- snpseq

  if(half){
    id <- snpseq == "1"
    tl <- sum(id)
    codeseq[id] <- rep(c("A", "G"), ceiling(tl / 2))[1:tl]
    id <- snpseq == "2"
    tl <- sum(id)
    codeseq[id] <- rep(c("C", "T"), ceiling(tl / 2))[1:tl]
  } else{
    codeseq[snpseq == "1"] <- "A"
    codeseq[snpseq == "2"] <- "C"
  }

  # codeseq[codeseq == "-"] <- "-"
  codeseq[!(snpseq %in% c("1", "2", "-"))] <- "."

  if(is.matrix(snpseq)){
    dim(codeseq) <- dim(snpseq)
  }
  codeseq
} # End of code2snp.default().

snp2code.list <- function(snpseq, half = TRUE){
  n.seq <- length(snpseq)
  n.half <- length(half)

  if(n.half == 1){
    half <- rep(half, n.seq)
  } else if(n.half != n.seq){
    stop("The length of half is not correct.")
  }

  lapply(1:n.seq, function(i){ snp2code.default(snpseq[[i]], half = half[i]) })
} # End of code2snp.list().


### Transfer nid and sid.
nid2sid <- function(nidseq){
  snp2sid(code2snp(nid2code(nidseq)))
} # End of nid2sid()

sid2nid <- function(sidseq, half = TRUE){
  code2nid(snp2code(sid2snp(sidseq), half = half))
} # End of nid2sid()


### Transfer an nucleotide id sequence to an amino acid id sequence.
### Eg. 103223.. => 3, 14, ...
nid2aid <- function(nidseq, start = 1, end = NULL, drop.gap = FALSE,
    byrow = TRUE){
  UseMethod("nid2aid")
} # End of nid2aid().

nid2aid.default <- function(nidseq, start = 1, end = NULL, drop.gap = FALSE,
    byrow = TRUE){
  code.aid.gap <- .amino.acid$aid[.amino.acid$code == "-"]
  code.nid.gap <- .nucleotide$nid[.nucleotide$code == "-"]
  code.aid <- .genetic.code$aid
  code.nid <- .genetic.code$nid

  nid2aid.a.seq <- function(a.nidseq){
    an.aidseq <- rep(code.aid.gap, length(a.nidseq))
    for(i in 1:nrow(.genetic.code)){
      an.aidseq[a.nidseq == code.nid[i]] <- code.aid[i]
    }
    an.aidseq
  }

  if(is.matrix(nidseq) && !byrow){
    nidseq <- t(nidseq)
  }
  if(is.vector(nidseq)){
    nidseq <- matrix(nidseq, nrow = 1)
  }

  n.seq <- nrow(nidseq)
  if(is.null(end)){
    end <- ncol(nidseq)
  }
  nidseq <- matrix(nidseq[, start:end], nrow = n.seq)
  end <- ncol(nidseq)

  if(! drop.gap){
    end.new <- end - end %% 3
    if(end.new < 3){
      stop("The sequence is too short.")
    }
    nidseq <- matrix(t(nidseq[, 1:end.new]), nrow = 3)

    nidseq <- apply(nidseq, 2, paste, collapse = "")
    aidseq <- matrix(nid2aid.a.seq(nidseq), nrow = n.seq, byrow = TRUE)

    if(n.seq == 1){
      aidseq <- as.vector(aidseq)
    } else{
      if(!byrow){
        aidseq <- t(aidseq)
      }
    }
  } else{
    aidseq <- list()
    for(i in 1:n.seq){
      a.nidseq <- nidseq[i,]
      a.nidseq <- a.nidseq[a.nidseq != code.nid.gap]
      end.new <- end - end %% 3
      if(end.new < 3){
        stop("The sequence is too short.")
      }
      a.nidseq <- matrix(a.nidseq[1:end.new], nrow = 3)

      a.nidseq <- apply(a.nidseq, 2, paste, collapse = "")
      aidseq[[i]] <- nid2aid.a.seq(a.nidseq)
    }
  }

  aidseq
} # End of nid2aid.default().

nid2aid.list <- function(nidseq, start = 1, end = NULL, drop.gap = FALSE,
    byrow = TRUE){
  n.seq <- length(nidseq)
  n.start <- length(start)
  n.end <- length(end)
  n.drop.gap <- length(drop.gap)

  if(n.start == 1){
    start <- rep(start, n.seq)
  } else if(n.start != n.seq){
    stop("The length of start is not correct.")
  }

  if(n.end == 0){
    end <- do.call("c", lapply(nidseq, length))
  } else if(n.end == 1){
    end <- rep(end, n.seq)
  } else if(n.end != n.seq){
    stop("The length of end is not correct.")
  }

  if(n.drop.gap == 1){
    drop.gap <- rep(drop.gap, n.seq)
  } else if(length(drop.gap) != n.seq){
    stop("The length of drop.gap is not correct.")
  }

  lapply(1:n.seq,
    function(i){ nid2aid.default(nidseq[[i]], start = start[i], end = end[i],
                                 drop.gap = drop.gap[i]) })
} # End of nid2aid.list().


### Transfer an nucleotide id sequence to a codon id sequence.
### Eg. 103223.. => 19, 43, ...
nid2cid <- function(nidseq, start = 1, end = NULL, drop.gap = FALSE,
    byrow = TRUE){
  UseMethod("nid2cid")
} # End of nid2cid().

nid2cid.default <- function(nidseq, start = 1, end = NULL, drop.gap = FALSE,
    byrow = TRUE){
  code.nid.gap <- .nucleotide$nid[.nucleotide$code == "-"]
  code.cid.gap <- .codon$cid[.codon$code == "---"]
  code.nid <- .genetic.code$nid
  code.cid <- .genetic.code$cid

  nid2cid.a.seq <- function(a.nidseq){
    a.cidseq <- rep(code.cid.gap, length(a.nidseq))
    for(i in 1:nrow(.genetic.code)){
      a.cidseq[a.nidseq == code.nid[i]] <- code.cid[i]
    }
    a.cidseq
  }

  if(is.matrix(nidseq) && !byrow){
    nidseq <- t(nidseq)
  }
  if(is.vector(nidseq)){
    nidseq <- matrix(nidseq, nrow = 1)
  }

  n.seq <- nrow(nidseq)
  if(is.null(end)){
    end <- ncol(nidseq)
  }
  nidseq <- matrix(nidseq[, start:end], nrow = n.seq)
  end <- ncol(nidseq)

  if(! drop.gap){
    end.new <- end - end %% 3
    if(end.new < 3){
      stop("The sequence is too short.")
    }
    nidseq <- matrix(t(nidseq[, 1:end.new]), nrow = 3)

    nidseq <- apply(nidseq, 2, paste, collapse = "")
    cidseq <- matrix(nid2cid.a.seq(nidseq), nrow = n.seq, byrow = TRUE)

    if(n.seq == 1){
      cidseq <- as.vector(cidseq)
    } else{
      if(!byrow){
        cidseq <- t(cidseq)
      }
    }
  } else{
    cidseq <- list()
    for(i in 1:n.seq){
      a.nidseq <- nidseq[i,]
      a.nidseq <- a.nidseq[a.nidseq != code.nid.gap]
      end.new <- end - end %% 3
      if(end.new < 3){
        stop("The sequence is too short.")
      }
      a.nidseq <- matrix(a.nidseq[1:end.new], nrow = 3)

      a.nidseq <- apply(a.nidseq, 2, paste, collapse = "")
      cidseq[[i]] <- nid2cid.a.seq(a.nidseq)
    }
  }

  cidseq
} # End of nid2cid.default().

nid2cid.list <- function(nidseq, start = 1, end = NULL, drop.gap = FALSE,
    byrow = TRUE){
  n.start <- length(start)
  n.end <- length(end)
  n.seq <- length(nidseq)

  if(n.start == 1){
    start <- rep(start, n.seq)
  } else if(n.start != n.seq){
    stop("The length of start is not correct.")
  }

  if(n.end == 0){
    end <- do.call("c", lapply(nidseq, length))
  } else if(n.end == 1){
    end <- rep(end, n.seq)
  } else if(n.end != n.seq){
    stop("The length of end is not correct.")
  }

  if(length(drop.gap) == 1){
    drop.gap <- rep(drop.gap, n.seq)
  } else if(length(drop.gap) != n.seq){
    stop("The length of drop.gap is not correct.")
  }

  lapply(1:n.seq,
    function(i){ nid2cid.default(nidseq[[i]], start = start[i], end = end[i],
                                 drop.gap = drop.gap[i]) }) 
} # End of nid2cid.list().


### Transfer an codon id sequence to an amino acid id sequence.
### Eg. 19, 43, ... => 3, 14, ...
cid2aid <- function(cidseq){
  UseMethod("cid2aid")
} # End of cid2aid().

cid2aid.default <- function(cidseq){
  code.cid <- .genetic.code$cid
  code.aid <- .genetic.code$aid
  code.aid.gap <- .amino.acid$aid[.amino.acid$code == "-"]

  aidseq <- cidseq
  for(i in 1:nrow(.genetic.code)){
    aidseq[cidseq == code.cid[i]] <- code.aid[i]
  }
  aidseq[!cidseq %in% code.cid] <- code.aid.gap

  if(is.matrix(cidseq)){
    dim(aidseq) <- dim(cidseq)
  }
  aidseq
} # End of cid2aid.default().

cid2aid.list <- function(cidseq){
  lapply(cidseq, cid2aid.default)
} # End of cid2aid.list().


### Transfer an amino acid id sequence to an amino acid code sequence.
### Eg. 3, 14, ... => DP...
aid2acode <- function(aidseq, lower.case = FALSE){
  UseMethod("aid2acode")
} # End of aid2acode().

aid2acode.default <- function(aidseq, lower.case = FALSE){
  if(lower.case){
    code.char <- as.character(.amino.acid$code.l)
  } else{
    code.char <- as.character(.amino.acid$code)
  }
  code.aid <- .amino.acid$aid

  acodeseq <- aidseq
  for(i in 1:nrow(.amino.acid)){
    acodeseq[aidseq == code.aid[i]] <- code.char[i]
  }
  acodeseq[!aidseq %in% code.aid] <- "-"

  if(is.matrix(aidseq)){
    dim(acodeseq) <- dim(aidseq)
  }
  acodeseq
} # End of aid2acode.default().

aid2acode.list <- function(aidseq, lower.case = FALSE){
  n.seq <- length(aidseq)
  n.lower.case <- length(lower.case)

  if(n.lower.case == 1){
    lower.case <- rep(lower.case, n.seq)
  } else if(n.lower.case != n.seq){
    stop("The length of lower.case is not correct.")
  }

  lapply(aidseq, aid2acode.default)
} # End of aid2acode.list().


### Transfer an amino acid code sequence to an amino acid id sequence.
### EG. DP... => 3, 14, ...
acode2aid <- function(acodeseq){
  UseMethod("acode2aid")
} # End of acode2aid().

acode2aid.default <- function(acodeseq){
  code.char <- as.character(.amino.acid$code)
  code.l.char <- as.character(.amino.acid$code.l)
  code.aid <- .amino.acid$aid

  aidseq <- acodeseq
  for(i in 1:nrow(.amino.acid)){
    aidseq[acodeseq %in% c(code.char[i], code.l.char[i])] <- code.aid[i]
  }
  aidseq[! acodeseq %in% c(code.char, code.l.char)] <-
    code.aid[code.char == "-"]
  aidseq <- as.numeric(aidseq)

  if(is.matrix(acodeseq)){
    dim(aidseq) <- dim(acodeseq)
  }
  aidseq
} # End of acode2aid.default().

acode2aid.list <- function(acodeseq){
  lapply(acodeseq, acode2aid.default)
} # End of acode2aid.list().

