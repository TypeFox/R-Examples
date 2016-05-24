m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"))
ncol.msa(m)
ncol.msa(m, names.msa(m))
