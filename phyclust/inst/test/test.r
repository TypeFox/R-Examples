library(phyclust)

### Test ms() and seqgen()
set.seed(123)
(ret.ms <- ms(nsam = 5, nreps = 1, opts = "-T"))
seqgen(opts = "-mHKY -l40 -q", newick.tree = ret.ms[3])

### Test phyclust()
X <- seq.data.toy$org
X.class <- as.numeric(gsub(".*-(.)", "\\1", seq.data.toy$seqname))
EMC.2 <- .EMControl(init.procedure = "emEM")
set.seed(1234)
(ret.2 <- phyclust(X, 4, EMC = EMC.2))
RRand(ret.2$class.id, X.class)
