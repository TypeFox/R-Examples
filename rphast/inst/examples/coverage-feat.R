feat1 <- feat(seqname=c(rep("chr1", 3), rep("chr2", 2)),
              start=c(1, 5, 100, 10, 20),
              end=c(7, 10, 105, 15, 30))
feat2 <- feat(seqname=c("chr1","chr2"),
              start=c(1,1), end=c(5,10))
coverage.feat(feat1, feat2, or=FALSE)
coverage.feat(feat1, feat2, or=TRUE)
coverage.feat(feat1, feat2, get.feats=TRUE, or=TRUE)
coverage.feat(feat1, feat2, or=TRUE)
rm(feat1, feat2)
