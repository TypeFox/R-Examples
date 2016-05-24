trees <- c("((hg18, panTro2), mm9);",
           "((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385)
            :0.020666,canFam2:0.193569);")
trees <- name.ancestors(trees)
subtree(trees, c("hg18-panTro2", "mm9-rn4"))
