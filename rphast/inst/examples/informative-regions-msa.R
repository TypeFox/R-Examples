require("rphast")
m <- msa(seqs=c("A--ACGTAT-", "AG-AGGTAA-", "AGGAGGTA--"),
         names=c("human", "mouse", "rat"))
informative.regions.msa(m, 1, refseq=NULL)
informative.regions.msa(m, 3, refseq=NULL)
informative.regions.msa(m, 3, refseq="mouse", spec=c("mouse", "rat"))
