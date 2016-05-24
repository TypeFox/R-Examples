"consensus" <-
function(alignment, cutoff=0.6) {
  #  Determine the consensus sequence for a given alignment

  if(is.list(alignment)) alignment=alignment$ali

  aa <- c("V","I","L","M",  "F","W","Y",  "S","T",
          "N","Q",  "H","K","R",  "D","E",
          "A","G",  "P","C",  "-","X")
  
  composition <- table(alignment)
  unk <- composition[!names( composition  ) %in% aa]
  if(length(unk) > 0) {
    warning(paste("\nnon standard residue code:",names(unk),"maped to X"))
    for(i in 1:length(unk))
      alignment[alignment==names(unk[i])]="X"
  }


  len <- ncol(alignment)

  freq <- matrix(0, nrow = 22, ncol = ncol(alignment),
                 dimnames = list(aa,seq(1:len)))

  for (i in 1:len) { 
    freq[names(summary((as.factor(toupper(alignment[,i]))))), i] <-
      (summary(as.factor(toupper(alignment[,i])))/length(alignment[,i]))
  }

  cons.freq <- apply(freq[1:20,], 2, max)
  cons.tmp  <- aa[apply(freq[1:20,], 2, which.max)]
  cons.tmp[cons.freq <= cutoff] = "-"

  return(list(seq=cons.tmp,
              freq=freq,
              seq.freq=cons.freq,
              cutoff=cutoff))
}

