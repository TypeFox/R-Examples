trees <- c("((hg18, panTro2), mm9);",
           "((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):
                0.020666,canFam2:0.193569);")
prune.tree(trees, c("panTro2", "mm9"), all.but=TRUE)
prune.tree(trees, "hg18", all.but=FALSE)
