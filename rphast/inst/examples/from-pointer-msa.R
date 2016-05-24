m <- msa(seqs=c("A--ACGTAT", "AG-AGGTAA", "AGGAGGTAG"),
         names=c("human", "mouse", "rat"), pointer.only=TRUE)
m
m <- from.pointer.msa(m)
m
