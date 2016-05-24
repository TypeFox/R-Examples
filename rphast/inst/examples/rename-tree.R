trees <- c("((hg18:1.0, panTro2:2.0):3.0, mm9:4.0);",
           "((hg18:0.142679,(mm9:0.083220,rn4:0.090564):0.269385):
                0.020666,canFam2:0.193569);")
rename.tree(trees,
            old.names=c("hg18", "panTro2", "mm9", "rn4", "canFam2"),
            new.names=c("human", "chimp", "mouse", "rat", "dog"))
