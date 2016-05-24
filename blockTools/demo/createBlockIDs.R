
data(x100)
## Block:
out <- block(x100, groups = "g", n.tr = 2, id.vars = c("id"), block.vars
             = c("b1", "b2"))
## Create vector of block IDs:
createBlockIDs(out, x100, id.var = "id")
## block ID integers are unique, even with several groups