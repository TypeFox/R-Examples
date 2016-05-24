m <- msa(seqs=c("ACGTAT", "AGGTAA", "AGGTAG"),
         names=c("human", "mouse", "rat"))
print(m[c("rat", "rat", "human"), ], print.seq=TRUE)
print(m[c(3,3,1),], print.seq=TRUE)
print(m[c(TRUE, FALSE, TRUE),], print.seq=TRUE)
print(m[TRUE,], print.seq=TRUE)
print("[.msa"(m, "mouse",c(1,6,3,5)), print.seq=TRUE)
