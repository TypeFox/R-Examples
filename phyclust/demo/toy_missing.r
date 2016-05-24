library(phyclust, quiet = TRUE)

X <- seq.data.missing$org
X.class <- as.numeric(gsub(".*-(.)", "\\1", seq.data.missing$seqname))

### A dot plot.
windows()
plotdots(X, X.class)

### A histogram plot.
windows()
plothist(X, X.class)

### A Neighbor-Joining plot.
ret <- phyclust.edist(X, edist.model = .edist.model[3])
ret.tree <- nj(ret)
windows()
plotnj(ret.tree, X.class = X.class)

### Fit a EE, JC69 model using emEM
EMC.2 <- .EMControl(init.procedure = "emEM")
set.seed(1234)
(ret.2 <- phyclust(X, 4, EMC = EMC.2))
RRand(ret.2$class.id, X.class)
