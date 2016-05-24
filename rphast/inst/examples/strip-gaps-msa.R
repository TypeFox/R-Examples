m <- msa(seqs=c("A--ACGTAT-", "AG-AGGTAA-", "AGGAGGTA--"),
         names=c("human", "mouse", "rat"))
print(strip.gaps.msa(m, c("human", "mouse")), print.seq=TRUE)
print(strip.gaps.msa(m, strip.mode="any.gaps"), print.seq=TRUE)
print(strip.gaps.msa(m, strip.mode="all.gaps"), print.seq=TRUE)
print(m, print.seq=TRUE)
#' NOTE if msa stored as pointer, original object is changed
m <- as.pointer.msa(m)
temp <- strip.gaps.msa(m, "any.gaps")
print(m, print.seq=TRUE)
