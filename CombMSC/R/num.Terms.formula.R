`num.Terms.formula` <-
function(object, ...)
ifelse(length(attr(terms(object), "factors")) == 0, attr(terms(object), "intercept"),
dim(attr(terms(object), "factors"))[2] + attr(terms(object), "intercept"))

