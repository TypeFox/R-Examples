ic.est <- function (x, Sigma, ui, ci = NULL, index = 1:nrow(Sigma), meq = 0, 
    tol = sqrt(.Machine$double.eps)) 
{
    ## check inputs
    if (!is.vector(x)) 
        stop("x must be a numeric vector.")
    if (!is.matrix(Sigma)) 
        stop("Sigma must be a matrix.")
    if (!(length(x) == nrow(Sigma) & nrow(Sigma) == ncol(Sigma))) 
        stop("Check length of x and dimensions of Sigma.")
    g <- length(x)
    if (any(eigen(Sigma)$values < 0)) 
        stop("Sigma must be non-negative definite.")
    if (!(is.vector(index))) 
        stop("index must be a vector.")
    if (is.vector(ui)) 
        ui <- matrix(ui, 1, length(ui))
    if (!is.matrix(ui)) 
        stop("ui must be a matrix.")
    if (!length(index) == ncol(ui)) 
        stop("mismatch between number of columns for ui and length of index")
    if (is.null(ci)) 
        ci <- rep(0, nrow(ui))
    if (!is.vector(ci)) 
        stop("ci must be a vector.")
    if (!nrow(ui) == length(ci)) 
        stop("mismatch between number of rows in ui and elements in ci")
    namen <- names(x)
    if (is.null(namen)) {
        namen <- colnames(Sigma)
        if (is.null(namen)) 
            namen = paste("m", 1:length(x), sep = "")
        names(x) <- namen
        colnames(Sigma) <- namen
    }
    ## show linear independent subset of rows of ui
    hilf <- RREF(t(ui))
    if (hilf$rank < nrow(ui)) 
        stop(paste("Matrix ui must have full row-rank (choose e.g. rows", 
            paste(hilf$pivot, collapse = " "), ")."))

    ## expand coefficient matrix, if necessary
    uiw <- ui
    if (length(index) > g | max(index) > g) 
        stop(paste("index must be vector of index positions, at most of length ", 
            g))
    if (length(index) < g) 
        unrestr <- setdiff(1:length(x), index)
    uiw <- matrix(0, nrow(ui), g)
    uiw[, index] <- ui
    if (all(uiw %*% x - ci >= 0 * ci) & meq == 0) 
        aus <- list(b.unrestr = x, b.restr = x, Sigma = Sigma, 
            ui = ui, ci = ci, restr.index = index, meq = meq, 
            iact = NULL)
    else {
        aus <- solve.QP(Dmat = solve(Sigma), dvec = solve(Sigma, 
            x), Amat = t(uiw), bvec = ci, meq = meq)
        aus <- list(b.unrestr = x, b.restr = aus$solution, Sigma = Sigma, 
            ui = ui, ci = ci, restr.index = index, meq = meq, 
            iact = aus$iact)
    }
    aus$b.restr[abs(aus$b.restr) < tol] <- 0
    names(aus$b.restr) <- namen
    class(aus) <- "orest"
    aus
}
