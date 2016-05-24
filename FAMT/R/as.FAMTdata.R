as.FAMTdata <-
function (expression, covariates = NULL, annotations = NULL, 
    idcovar = 1, idannot = NULL, na.action = TRUE) 
{
    expr = expression
    n = ncol(expr)
    covar = data.frame(ID = colnames(expr), x = rep(1, n))
    annot = data.frame(ID = rownames(expr))
    if (!is.null(annotations)) {
        annot = annotations
        if ((!is.null(annotations)) & (nrow(expr) != nrow(annot))) 
            stop("Numbers of rows of expression and annotations should correspond.")
        if (is.null(idannot)) {
            if (!any(is.element(colnames(annot), "ID"))) 
                stop("One of the columns in annotations should be named ID")
        }
        if (!is.null(idannot)) {
            if (!any(is.element(1:ncol(annot), idannot))) 
                stop(paste("idannot should be in 1:", ncol(annot), 
                  sep = ""))
            names(annot)[idannot] = "ID"
        }
    }
    if (!is.null(covariates)) {
        covar = covariates
        if (ncol(expr) != nrow(covar)) 
            stop("Dimensions of expression and covariates should correspond.")
        if (!any(is.element(1:ncol(covar), idcovar))) 
            stop(paste("idcovar should be in 1:", ncol(covar), 
                sep = ""))
        if (!setequal(colnames(expr), as.character(covar[, idcovar]))) 
            stop("Names should correspond in expression and covariates")
    }
    m = nrow(expr)
    n = ncol(expr)
    nbnarow = (1:m)[apply(is.na(expr), 1, sum) > 0]
    nbnacol = (1:n)[apply(is.na(expr), 2, sum) > 0]
    na.expr = vector(length = 2, "list")
    names(na.expr) = c("Rows with missing values", "Columns with missing values")
    na.expr[[1]] = nbnarow
    na.expr[[2]] = nbnacol
    print(na.expr)
    if (na.action == TRUE) {
        if (sum(is.na(expr)) > 0) {
            expr = as.data.frame(impute.knn(as.matrix(expr))$data)
            print("Missing values were imputed using nearest neighbor averaging (impute.knn {impute})")
        }
    }
    ordcovar = order(colnames(expr))
    res = list(expression = expr[, ordcovar], covariates = covar[ordcovar, 
        ], annotations = annot, idcovar = idcovar, na.expr = na.expr)
    class(res) = list("FAMTdata", "list")
    return(res)
}
