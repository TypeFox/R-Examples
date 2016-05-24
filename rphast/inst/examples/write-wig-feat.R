f <- feat(seqname=c("chr1", "chr1", "chr1", "chr1", "chr2"),
          start=c(1, 11, 21, 100, 1),
          end=c(3, 13, 23, 102, 3),
          score=runif(5))
write.wig.feat(f)
