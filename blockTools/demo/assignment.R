
data(x100)

## First, block
out <- block(x100, groups = "g", n.tr = 2, id.vars = c("id"), block.vars
             = c("b1", "b2"), algorithm="optGreedy", distance =
             "mahalanobis", level.two = FALSE, valid.var = "b1",
             valid.range = c(0,500), verbose = TRUE)
## Second, assign
assg <- assignment(out, seed = 123)
## assg contains 3 data frames
