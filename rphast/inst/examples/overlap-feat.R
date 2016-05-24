feat1 <- feat(seqname=c(rep("chr1", 3), rep("chr2", 2)),
              start=c(1, 5, 100, 10, 20),
              end=c(7, 10, 105, 15, 30))
feat2 <- feat(seqname=c("chr1","chr2"),
              start=c(1,1),
              end=c(5,10))

overlap.feat(feat1, feat1)
overlap.feat(feat1, feat2, min.percent=0.25)
overlap.feat(feat1, feat2, min.percent=0.25, overlapping=FALSE)
overlap.feat(feat1, feat2, get.fragments=TRUE)
overlap.feat(feat1, feat2, get.fragments=TRUE)
rm(feat1, feat2)
