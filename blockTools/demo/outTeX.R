
data(x100)

## First, block
out <- block(x100, groups = "g", n.tr = 2, id.vars = c("id"), block.vars
             = c("b1", "b2"), algorithm="optGreedy", distance =
             "mahalanobis", level.two = FALSE, valid.var = "b1",
             valid.range = c(0,500), verbose = TRUE)
## Second, assign
assg <- assignment(out, seed = 123)

## create three .tex files of blocks
outTeX(out)
## create three .tex files of assigned blocks
##   (note: overwrites blocked .tex files)
outTeX(assg)
## create three .tex files with custom file names and captions
outTeX(assg, file.names = list("file1", "file2", "file3"), 
	   captions = list("This is caption 1.", "Caption 2.", "Caption 3?"))
