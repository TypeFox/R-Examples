library(phyclust, quiet = TRUE)

### Examples to use EM functions.
data.path <- paste(.libPaths()[1], "/phyclust/data/pony524.phy", sep = "")
my.seq <- read.phylip(data.path)
X <- my.seq$org

### Directly use phyclust().
set.seed(1234)
ret <- phyclust(X, 2)

### One EM step.
ret.em <- phyclust.em.step(X, ret)

### One E- and M- step.
ret.e <- phyclust.e.step(X, ret)
ret$Z.normalized <- ret.e
ret.m <- phyclust.m.step(X, ret)
ret.e.m <- phyclust.m.step(X, K = ret$K, Tt = ret$QA$Tt,
                           Z.normalized = ret.e,
                           substitution.model = ret$substitution.model,
                           identifier = ret$QA$identifier,
                           code.type = ret$code.type)

### Check logL.
phyclust.logL(X, ret.em)
phyclust.logL(X, ret.m)
phyclust.logL(X, ret.e.m)
