require("rtfbs")
seqs <- ms(seqs=c("AAAA", "ACACACAC", "CGCCG", "ACGTACGTACGT", "CGGGGGGGGGG"),
           paste("fake", 1:5, sep=""))
gcContent.ms(seqs)
