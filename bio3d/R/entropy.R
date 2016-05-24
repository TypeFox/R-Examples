"entropy" <-
function(alignment) {
  # Calculate the Shannon entropy score for each position
  # in an alignment

  if(is.list(alignment)) alignment=alignment$ali

  aa <- c("V","I","L","M",  "F","W","Y",  "S","T",
          "N","Q",  "H","K","R",  "D","E",
          "A","G",  "P",  "C",  "-","X")
  
  composition <- table(alignment)
  unk <- composition[!names( composition  ) %in% aa]
  if(length(unk) > 0) {
    warning(paste("non standard residue code:",names(unk),"mapped to X\n  "))
    for(i in 1:length(unk))
      alignment[alignment==names(unk[i])]="X"
  }
  
  len <- ncol(alignment)
  freq.22 <- matrix(0, nrow = 22, ncol = ncol(alignment),
                    dimnames = list(aa,seq(1:len)))

  freq.10 <- matrix(0, nrow = 10, ncol = ncol(alignment),
                    dimnames = list(c(1:10),c(1:len)))

  for (i in 1:len) { 
    freq.22[names(summary((as.factor(toupper(alignment[,i]))))), i] <-
      (summary(as.factor(toupper(alignment[,i])))/length(alignment[,i]))

    freq.10[1,i]  <- sum(freq.22[1:4,i])   # Hydrophobic, Aliphatic
    freq.10[2,i]  <- sum(freq.22[5:7,i])   # Aromatic
    freq.10[3,i]  <- sum(freq.22[8:9,i])   # Ser/Thr
    freq.10[4,i]  <- sum(freq.22[10:11,i]) # Polar
    freq.10[5,i]  <- sum(freq.22[12:14,i]) # Positive
    freq.10[6,i]  <- sum(freq.22[15:16,i]) # Negative
    freq.10[7,i]  <- sum(freq.22[17:18,i]) # Tiny
    freq.10[8,i]  <- sum(freq.22[19,i])    # Proline
    freq.10[9,i]  <- sum(freq.22[20,i])    # Cysteine
    freq.10[10,i] <- sum(freq.22[21:22,i]) # Gaps

  }

  entropy.22 <- vector(length = len)
  entropy.10 <- entropy.22

  for (i in 1:len) {
  #   entropy_i = sum[i] (P(X_i)log2(P(X_i)))
    entropy.22[i]  <- -1*sum(freq.22[freq.22[, i] != 0, i] *
                               log2(freq.22[freq.22[, i] != 0, i]))

    entropy.10[i] <- -1*sum(freq.10[freq.10[, i] != 0, i] *
                            log2(freq.10[freq.10[, i] != 0, i]))  
  }

  out <- list(H=entropy.22, H.10=entropy.10,
              H.norm=(1-(entropy.22/max(entropy.22))),
              H.10.norm=(1-(entropy.10/max(entropy.10))),
              freq=freq.22)
}

