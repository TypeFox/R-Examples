m <- msa(seqs=c("ACGT---AT", "AGGTAGTAA", "AGGAAGTAG"),
         names=c("human", "mouse", "rat"))
print(sub.msa(m, c("human", "rat"), start.col=3, end.col=6),
      print.seq=TRUE)
print(sub.msa(m, c("mouse"), keep=FALSE, refseq="human",
              start.col=3, end.col=4),
      print.seq=TRUE)
