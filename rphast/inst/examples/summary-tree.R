tree <- "(((hg18:0.01, panTro2:0.01)hg18-panTro2:0.07,
              (mm9:0.083220,rn4:0.090564)mm9-rn4:
             0.269385)hg18-rn4:0.020666,canFam2:0.193569);"
summary.tree(tree)
summary.tree(label.subtree(tree, "mm9-rn4", "rodent", include.leading=TRUE))
