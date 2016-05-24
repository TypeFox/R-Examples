# read in an MSA stored in R
m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
         names=c("human", "mouse", "rat"))
print(m)
print(m, format="FASTA")
print(m, format="PHYLIP", pretty.print=TRUE)
#'
# read in an MSA stored by reference in C
m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
         names=c("human", "mouse", "rat"),
         pointer.only=TRUE)
print(m)
