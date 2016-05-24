
data(x100)

## Block:
b <- block(x100, groups = "g", id.vars = "id", block.vars = c("b1", "b2"))
## Assign:
a <- assignment(b)
## Examine balance:
axb <- assg2xBalance(a, x100, id.var = "id", bal.vars = c("b1", "b2"))
axb
## axb is a list with 4 elements (one for each of 3 groups, plus one for 'Overall')
