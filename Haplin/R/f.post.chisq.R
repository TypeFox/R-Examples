f.post.chisq <- function(coeff, covar, contrast.mat){



.y <- contrast.mat %*% coeff
.cov.y <- contrast.mat %*% covar %*% t(contrast.mat)
.chisq <- t(.y) %*% solve(.cov.y) %*% .y
.df <- length(.y)
.pval <- pchisq(.chisq, df = .df, lower.tail = F)

.ut <- list(chisq = .chisq, df = .df, pval = .pval, y = .y) # probably no need for y anymore


return(.ut)

}

