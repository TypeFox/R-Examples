## data() sets directory
pitpropC <- as.matrix(read.table("pitpropC.tab", header=TRUE))
rownames(pitpropC) <- colnames(pitpropC)
