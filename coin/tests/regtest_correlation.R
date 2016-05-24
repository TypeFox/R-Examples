
### Regression tests for the correlation problem, i.e.,
### testing the independence of two numeric variables
### `x' and `y' (possibly blocked)

set.seed(290875)
library("coin")
isequal <- coin:::isequal
options(useFancyQuotes = FALSE)

### generate data
dat <- data.frame(x = rnorm(100), y = rnorm(100), block = gl(10, 10))

### not really the same, T = (rank(x) - rank(y))^2 is used here
cor.test(~ x + y, data = dat, method = "spearman")$p.value
cor.test(~ x + y, data = dat, alternative = "less", method = "spearman")$p.value
cor.test(~ x + y, data = dat, alternative = "greater", method = "spearman")$p.value

### without blocks
pvalue(spearman_test(y ~ x, data = dat))
pvalue(spearman_test(x ~ y, data = dat))
pvalue(spearman_test( ~ y + x, data = dat))
pvalue(spearman_test( ~ x + y, data = dat))

pvalue(fisyat_test(y ~ x, data = dat))
pvalue(fisyat_test(x ~ y, data = dat))
pvalue(fisyat_test( ~ y + x, data = dat))
pvalue(fisyat_test( ~ x + y, data = dat))

pvalue(quadrant_test(y ~ x, data = dat))
pvalue(quadrant_test(x ~ y, data = dat))
pvalue(quadrant_test( ~ y + x, data = dat))
pvalue(quadrant_test( ~ x + y, data = dat))

pvalue(koziol_test(y ~ x, data = dat))
pvalue(koziol_test(x ~ y, data = dat))
pvalue(koziol_test( ~ y + x, data = dat))
pvalue(koziol_test( ~ x + y, data = dat))

### with blocks
pvalue(spearman_test(y ~ x | block, data = dat))
pvalue(spearman_test(x ~ y | block, data = dat))
pvalue(spearman_test( ~ y + x | block, data = dat))
pvalue(spearman_test( ~ x + y | block, data = dat))

pvalue(fisyat_test(y ~ x | block, data = dat))
pvalue(fisyat_test(x ~ y | block, data = dat))
pvalue(fisyat_test( ~ y + x | block, data = dat))
pvalue(fisyat_test( ~ x + y | block, data = dat))

pvalue(quadrant_test(y ~ x | block, data = dat))
pvalue(quadrant_test(x ~ y | block, data = dat))
pvalue(quadrant_test( ~ y + x | block, data = dat))
pvalue(quadrant_test( ~ x + y | block, data = dat))

pvalue(koziol_test(y ~ x | block, data = dat))
pvalue(koziol_test(x ~ y | block, data = dat))
pvalue(koziol_test( ~ y + x | block, data = dat))
pvalue(koziol_test( ~ x + y | block, data = dat))

### sanity checks, those should be errors
dat <- data.frame(x = gl(2, 50), y = rnorm(100), block = rnorm(100))

try(spearman_test(y ~ x, data = dat))
try(spearman_test(y ~ x | block, data = dat))

try(fisyat_test(y ~ x, data = dat))
try(fisyat_test(y ~ x | block, data = dat))

try(quadrant_test(y ~ x, data = dat))
try(quadrant_test(y ~ x | block, data = dat))

try(koziol_test(y ~ x, data = dat))
try(koziol_test(y ~ x | block, data = dat))
