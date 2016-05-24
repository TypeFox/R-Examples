###
# Generating a subset of a full model excluding collinear variables
###

library(MuMIn)
options(na.action = na.fail)

# Fit the 'global model'
fm <- lm(y ~ (X1 + X2 + X3 + X4)^2, data = Cement)

# Suppose we want to have set of models that exclude combinations of colinear
# variables, that are significantly (p < 0.05) correlated, with Pearson
# correlation coefficient larger than r = 0.5.

is.correlated <- function(i, j, data, conf.level = .95, cutoff = .5, ...) {
	if(j >= i) return(NA)
	ct <- cor.test(data[, i], data[, j], conf.level = conf.level, ...)
	ct$p.value > (1 - conf.level) || abs(ct$estimate) <= cutoff
}
# Need vectorized function to use with 'outer'
vCorrelated <- Vectorize(is.correlated, c("i", "j"))

# Create logical matrix
smat <- outer(1:4, 1:4, vCorrelated, data = Cement)
nm <- colnames(Cement[-1])
dimnames(smat) <- list(nm, nm)

### A simpler case: exclude only pairs of variables having cor. coefficient
### r > 0.5
# smat <- abs(cor(Cement[, -5])) <= .5
# smat[!lower.tri(smat)] <- NA

# Alternatively, we can use logical expression of form:
# !((V1 && V2) || (V3 && V4))
# where V1 is collinear with V2, and V3 with V4.
# Rather that doing it by hand, we can generate it from the above matrix:
i <- as.vector(smat == FALSE & !is.na(smat))
sexpr <-parse(text = paste("!(", paste("(",
    nm[col(smat)[i]], " && ",
    nm[row(smat)[i]], ")",
    sep = "", collapse = " || "), ")"))

smat
sexpr

## =============================================================================

system.time(dd2 <- dredge(fm, subset = smat))
system.time(dd1 <- dredge(fm, subset = sexpr))

# Using the argument 'subset' in a form of matrix is usually faster.
# The results are identical:
dd1
dd2

## =============================================================================