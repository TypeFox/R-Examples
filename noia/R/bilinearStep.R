bilinearStep <-
function (form, X, phen, marginal, interactions, ...) 
{
    n.marginal <- names(marginal)
    n.interactions <- names(interactions)
    rrr <- nls(formula = as.formula(form), data = marginal, start = interactions, 
        ...)
    interactions <- as.list(coef(rrr))
    names(interactions) <- n.interactions
    rrr <- nls(formula = as.formula(form), data = interactions, 
        start = marginal, ...)
    marginal <- as.list(coef(rrr))
    names(marginal) <- n.marginal
    ans <- list(marginal = marginal, interactions = interactions)
    return(ans)
}
