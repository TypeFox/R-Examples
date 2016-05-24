trees <- c("((hg18:1.0, panTro2:2.0):3.0, mm9:4.0);",
           "((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):
                0.020666,canFam2:0.193569);")
rescale.tree(trees, 0.5)
rescale.tree(trees, c(0.5, 2.0))
trees <- name.ancestors(trees)
rescale.tree(trees, 0.5, c("hg18-panTro2", "hg18-mm9"))
