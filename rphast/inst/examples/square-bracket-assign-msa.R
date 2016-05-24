m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
         names=c("human", "mouse", "rat"))
m[1:2,4:6] <- "G"
m[1,] <- "A"
m[,4:5] <- "-"
m["rat",] <- "ABCDEF"
