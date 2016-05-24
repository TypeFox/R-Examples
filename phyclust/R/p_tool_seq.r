### Generate sequences by the HKY85 model.
gen.seq.HKY <- function(rooted.tree, pi, kappa, L, rate.scale = 1,
    anc.seq = NULL){
  if(!is.rooted(rooted.tree)){
    stop("A rooted tree is required.")
  }
  if(length(pi) != 4 || sum(pi) != 1 || any(pi <= 0) || any(pi >= 1)){
    stop("The pi is not correct.")
  }
  if(length(kappa) != 1 || kappa <= 0){
    stop("The kappa is not correct.")
  }
  if(rate.scale <= 0){
    stop("The rate.scale is not correct.")
  }
  if(!is.null(anc.seq)){
    if(any(!anc.seq %in% .nucleotide$nid[1:4]) || length(anc.seq) != L){
      stop("The anc.seq is not correct.")
    }
  }

  ttips <- rooted.tree$Nnode + 1

  ### seqgen uses the order of A, C, G, T for pi.
  ### seqgen uses ts/tv = (kappa*(freqA*freqG + freqC*freqT))/(freqR*freqY)
  ###              freqR = freqA + freqG and freqY = freqC + freqT
  ts.tv <- kappa * (pi[1] * pi[2] + pi[3] * pi[4]) /
                   ((pi[1] + pi[2]) * (pi[3] + pi[4]))

  if(is.null(anc.seq)){
    opts <- paste("-mHKY",
                  " -t", ts.tv,
                  " -f", paste(pi[c(1, 3, 2, 4)], collapse = ","),
                  " -l", L,
                  " -s", rate.scale,
                  " -u", ttips + 1,
                  " -q",
                  sep = "")
    ret <- seqgen(opts = opts, rooted.tree = rooted.tree)
  } else{
    L <- length(anc.seq)
    mu <- paste(nid2code(anc.seq, lower.case = FALSE), collapse = "")
    seqname <- paste("Ancestor  ", collapse = "")
    input <- c(paste(" 1", length(anc.seq), sep = " "),
               paste(seqname, mu, sep = ""),
               1,
               write.tree(rooted.tree, digits = 12))
    opts <- paste("-mHKY",
                  " -t", ts.tv,
                  " -f", paste(pi[c(1, 3, 2, 4)], collapse = ","),
                  " -l", L,
                  " -s", rate.scale,
                  " -u", ttips + 1,
                  " -k1",
                  " -q",
                  sep = "")
    ret <- seqgen(opts, input = input)
  }

  ret
} # End of gen.seq.HKY().


### Generate sequences by the SNP model.
gen.seq.SNP <- function(rooted.tree, pi, L, rate.scale = 1,
    anc.seq = NULL){
  if(!is.rooted(rooted.tree)){
    stop("A rooted tree is required")
  }
  if(length(pi) != 2 || sum(pi) != 1 || any(pi <= 0) || any(pi >= 1)){
    stop("The pi is not correct.")
  }
  if(rate.scale <= 0){
    stop("The rate.scale is not correct.")
  }
  if(!is.null(anc.seq)){
    if(any(!anc.seq %in% .snp$sid[1:2]) || length(anc.seq) != L){
      stop("The anc.seq is not correct.")
    }
  }

  ttips <- rooted.tree$Nnode + 1

  ### Input pi is for "1" and "2".
  ### transfer to "1" = A + G, "2" = C + T
  pi <- c(pi[1], pi[1], pi[2], pi[2]) / 2

  ### seqgen uses the order of A, C, G, T for pi.
  if(is.null(anc.seq)){
    opts <- paste("-mHKY",
                  " -f", paste(pi[c(1, 3, 2, 4)], collapse = ","),
                  " -l", L,
                  " -s", rate.scale,
                  " -u", ttips + 1,
                  " -q",
                  sep = "")
    ret <- seqgen(opts = opts, rooted.tree = rooted.tree)
  } else{
    L <- length(anc.seq)

    ### Replace half of "1" as "A" and the other half as "G"
    ### Replace half of "2" as "C" and the other half as "T"
    ### To avoid bad state of ancestor states in seqgen.
    mu <- paste(snp2code(sid2snp(anc.seq)), collapse = "")
    seqname <- paste("Ancestor  ", collapse = "")
    input <- c(paste(" 1", length(anc.seq), sep = " "),
               paste(seqname, mu, sep = ""),
               1,
               write.tree(rooted.tree, digits = 12))
    opts <- paste("-mHKY",
                  " -f", paste(pi[c(1, 3, 2, 4)], collapse = ","),
                  " -l", L,
                  " -s", rate.scale,
                  " -u", ttips + 1,
                  " -k1",
                  " -q",
                  sep = "")
    ret <- seqgen(opts, input = input)
  }

  ### Replace "A" and "G" by "1", "C" and "T" by "2".
  nseq <- length(ret) - 1
  org <- do.call("rbind", lapply(ret[2:length(ret)], unstrsplit, ""))
  org <- as.matrix(org, nrow = nseq)
  org[, 10 + (1:L)] <- code2snp(org[, 10 + (1:L)])
  org <- apply(org, 1, paste, sep = "", collapse = "")
  ret <- c(ret[1], org)
  attr(ret, "class") <- "seqgen"

  ret
} # End of gen.seq.SNP().

