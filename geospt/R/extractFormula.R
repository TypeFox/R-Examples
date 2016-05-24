assign("extractFormula",
function(formula, data, newdata) {
	# extract y and X from data:
    m = model.frame(terms(formula), as(data, "data.frame"))
    z = model.extract(m, "response")
    if (length(z) == 0)
        stop("no response variable present in formula")
    Terms = attr(m, "terms")
    X = model.matrix(Terms, m)

	# extract x0 from newdata:
    terms.f = delete.response(terms(formula))
    mf.f = model.frame(terms.f, newdata) #, na.action = na.action)
    x0 = model.matrix(terms.f, mf.f)
	list(z = z, X = X, x0 = x0)
}
)