linearRegression <-
function (phen, gen = NULL, genZ = NULL, reference = "noia", 
    max.level = NULL, max.dom = NULL, fast = FALSE) 
{
    "is.diag" <- function(mat, threshold = 1e-06) {
        return(all(mat[lower.tri(mat, diag = FALSE)] < threshold) && 
            all(mat[upper.tri(mat, diag = FALSE)] < threshold))
    }
    prep <- prepareRegression(phen = phen, gen = gen, genZ = genZ, 
        reference = reference, max.level = max.level, max.dom = max.dom, 
        fast = fast)
    nn <- NULL
    if (is.null(prep$smat)) {
        nn <- effectsNamesGeneral(prep$nloc, max.level, max.dom)
    }
    else {
        nn <- colnames(prep$smat)
    }
    regression <- lm(prep$phen ~ prep$x + 0)
    if (!is.diag(t(prep$x) %*% prep$x)) {
        warning("The decomposition of genetic effects is not orthogonal.\n")
    }
    ans <- prep
    class(ans) <- "noia.linear"
    ans$E <- coef(regression)
    names(ans$E) <- nn
    ans$variances <- effectsVariances(ans)
    ans$std.err <- effectsStdErr(regression)
    ans$pvalues <- effectsPvalues(regression)
    ans$resvar <- var(residuals(regression))
    ans$regression <- regression
    return(ans)
}
