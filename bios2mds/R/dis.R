dis <- function (seq1, seq2, sub.mat.id = "PAM250", gap = NULL) {

  if (!exists("sub.mat"))
    data("sub.mat", package = "bios2mds", verbose= FALSE)

  if (length(seq1) !=  length(seq2))
    stop("seq1 and seq2 are not of the same length")

 if (!is.element(sub.mat.id, names(sub.mat)))
    stop("sub.mat does not contain sub.mat.id")
        
  if (!is.null(gap) && (length(gap) != 2))
    stop("gap is not of length 2")

  x <- sub.mat[[sub.mat.id]]
  
  #turn !aa and ! gap into NA
  #aa.strict = TRUE because of sub.mat
  seq1[!as.logical(is.aa(seq1, TRUE) + is.gap(seq1))] <- NA
  seq2[!as.logical(is.aa(seq2, TRUE) + is.gap(seq2))] <- NA

  if (!is.null(gap)) {
    #turn gap into "-"
    seq1[is.gap(seq1)] <- "-"
    seq2[is.gap(seq2)] <- "-"

    #add gap score to x
    aa <- c(aa, "-")
    x <- rbind (x, "-" = rep(gap[1], nrow(x)))
    x <- cbind (x, "-" = c(rep(gap[1], ncol(x)), gap[2]))
  }
  else {
   seq1[is.gap(seq1)] <- NA
   seq2[is.gap(seq2)] <- NA
  }

  #compute Sij (score per site)
  Sij <- function(i, j) {
    ij <- na.omit(data.frame(i = i, j = j))
    return(sum(apply(ij, 1, function(i) x[i["i"], i["j"]])))
  }

  S <- Sij(seq1, seq2)

  #compute Tij (average upper limit of the score per site)
  Tij <- function(i, j) {
    ij <- na.omit(data.frame(i = i, j = j))
    return(0.5 * sum(apply(ij, 1, function(i) {x[i["i"], i["i"]] + x[i["j"], i["j"]]})))
  }

  T <- Tij(seq1, seq2)

  #compute Srandij (score per site expected from random sequences)
  Srandij <- function(i, j) {
    pair.rand <- expand.grid(list(seq_len(length(i)), seq_len(length(j))))
    return(sum(mapply(function(k, l) i[k] * j[l] * x[names(i[k]), names(j[l])], pair.rand[, 1], pair.rand[, 2])))
  }

  #compute amino acid composition for each alignment
  #count not based on aa to avoid 0 occurrences
  cases <- complete.cases(cbind(seq1, seq2))
  count1 <- summary(as.factor(seq1[cases]))
  count2 <- summary(as.factor(seq2[cases]))
  lkij <- sum(cases)
  Srand <- Srandij(count1, count2)/lkij

  #compute Vij (normalized score per site)
  if (T != Srand)
    V <- (S - Srand)/(T - Srand)
  else
    V <- NA

  #V can be negative
  if (V <= 0)
    V <- 0.001

  #dis = sim - 1
  dis <- round(1 - V, 3)
 
  return (dis)
}
