require("rphast")
align <- msa(seqs=c("A-GTAT", "-GGTAA", "AG--AG"),
             names=c("human", "mouse", "rat"))
feats <- feat(seqname=c("MSA", "human", "human", "mouse", "mouse", "rat"),
              start=c(1, 2, 3, 1, 3, 3),
              end=c(6, 4, 4, 4, 5, 3))
convert.coords.feat(feats, align)  # convert everything to human coords
convert.coords.feat(feats, align, to="MSA")  # convert to alignment coords

# here, there is no position 6 in human alignment so feature is removed
convert.coords.feat(feat(seqname="human", start=6, end=6),
                    align, to="MSA")

# here, feature goes beyond end of MSA so it is truncated:
convert.coords.feat(feat(seqname="rat", start=2, end=100),
                    align, to="MSA")


# note that if the "to" species has gaps at the endpoints, they will
# be truncated:
align <- msa(seqs=c("A-GT-T", "ACGTGT"), names=c("human", "mouse"))
convert.coords.feat(feat(seqname="mouse", start=2, end=5),
                    align, to="human")
