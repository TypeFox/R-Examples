# The fitting routine is separated out only because it is called from multiple
#  places.
# Matched logistic regression models are fit using the coxph function,
#  unmatched ones are fit using glm.
# The very first time this is called we want to return all of the information
#  created by the underlying modeling function (fit=TRUE), other times only
#  the attributable risk itself is needed.

# Input 
#  y:      vector or matrix of response
#  x:      predictor matrix
#  w     : vector of weights
#  offset: offset vector (usually all zeros)
#  match:  the indexes for a matched analysis
#  xbase:  the matrix X(original) - X(modified)
#  fit  :  return the fit object?

attribrisk.fit <- function(x, y, w, offset, match, xbase, 
                      fit=FALSE) {
    if (is.null(match)) { # use glm
        fit1 <- glm.fit(x, y, weights=w, family=binomial(), offset=offset,
                       control=glm.control())
        cases <- which(y==1)
    }
    else {
        fit1 <- coxph.fit(x, y, offset=offset, method='exact',
                          weights=w, strata=match, 
                          rownames=dimnames(x)[[1]], 
                          control=coxph.control())
        cases <- which(y[,2] ==1)
    }

    ##LAS## if fit1 has a singularity assume a zero for the undefined coeffcient
    temp <- xbase %*% ifelse(is.na(coef(fit1)), 0, coef(fit1))
    attribrisk <- 1- sum((w* exp(-temp))[cases])/sum(w[cases])

    if (any(is.na(coefficients(fit1)))) 
        warning("at least one coefficient is na when fitting the regression for estimating attributable risk")
    

    if (fit) list(attribrisk= attribrisk, fit=fit1)
    else attribrisk
}
