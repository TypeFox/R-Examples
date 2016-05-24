library(phyclust, quiet = TRUE)

### Examples to analysis SGA data.
data.path <- paste(.libPaths()[1], "/phyclust/data/pony524.phy", sep = "")
my.seq <- read.phylip(data.path)
X <- my.seq$org

### Modify the em controler.
EMC <- .EMC
EMC$substitution.model <- "K80"
EMC$edist.model <- "D_HAMMING"

### Run phyclust() and find.best().
set.seed(1234)
(pony.phyclust <- phyclust(X, 2, EMC = EMC))
(pony.best <- find.best(X, 1, EMC = EMC))

### Calculate evolution distance.
pony.edist <- phyclust.edist(X, edist.model = EMC$edist.model)

### Manual initialization by hclust().
EMC$init.method <- "manualMu"
manual.id <- cutree(hclust(pony.edist), k = 2)
(pony.hc <- phyclust(X, 2, EMC = EMC, manual.id = manual.id))
