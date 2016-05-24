# Here is an MSA object stored in the default mode
m1 <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
          names=c("human", "mouse", "rat"))
m2 <- m1
# NOTE seqs would not be directly accessible if stored by reference
m2$seqs[3] <- "AAAAAA"
print(m1)
print(m1, print.seq=TRUE)
print(m2, print.seq=TRUE)
