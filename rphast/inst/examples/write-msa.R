m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
         names=c("human", "mouse", "rat"))
write.msa(m, "foo.ss")
write.msa(m, "foo.fa", pretty.print=TRUE)
write.msa(m, NULL, format="PHYLIP", pretty.print=TRUE)

#clean up
unlink("foo.ss")
unlink("foo.fa")
