orlm.default <- function (model, ui, ci = NULL, index = NULL, 
    meq = 0, tol = sqrt(.Machine$double.eps), df.error = NULL, ...) 
{
    ## this program is used, if model is a covariance matrix
    ## check model
    if (!(is.matrix(model))) 
        stop("ERROR: model must be of class lm or a covariance matrix.")
    else
    if (!(nrow(model)==ncol(model)))
        stop("ERROR: If it is not a linear model, model must be a quadratic matrix.")
    else 
    if (!(all(eigen(model,TRUE,only.values=TRUE)$values>0)))
        stop("ERROR: matrix model must be positive definite.")
    g <- nrow(model)-1
    if (is.null(df.error)) 
        stop("ERROR: df.error is required, when working from a covariance matrix")
    if (!(df.error>2)) 
        stop("ERROR: df.error must be at least 2.")
    if (is.null(index)) index <- 1:g
        else {
           if (!is.vector(index)) stop("Index must be a vector of integer numbers.")
           if (min(index) < 2) 
              stop("No restrictions on intercept possible, when working from a covariance matrix.")
           index <- index - 1
        }
        ### no intercept information obtainable from covariance matrix
        ###    i.e. index refers to x1 (2) to xp (p+1)
    ## preliminary calculations
    V <- solve(model[2:(g+1),2:(g+1)])/df.error  ## works also, if sigma^2=0
    if (is.null(colnames(V))) colnames(V) <- paste("X",1:g,sep="")
    b <- solve(model[2:(g+1),2:(g+1)],model[2:(g+1),1])
    var.y <- model[1,1]
    s2 <- c(var.y - model[1,2:(g+1)]%*%b)
    b <- c(b)
    names(b) <- colnames(V)
    orig.R2 <- model[1,2:(g+1)]%*%b/var.y
    ## check inputs
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
    hilf <- RREF(t(ui))
    if (hilf$rank < nrow(ui)) 
        stop(paste("Matrix ui must have full row-rank (choose e.g. rows", 
            paste(hilf$pivot, collapse = " "), ")."))
    ## expand ui by 0 columns, if necessary 
    uiw <- ui
    if (length(index) > g | max(index) > g) 
        stop(paste("index must be vector of index positions, at most of length ", 
            g))
    uiw <- matrix(0, nrow(ui), g)
    uiw[, index] <- ui
    ## inequality restrictions only, all fulfilled
    if (all(uiw %*% b - ci >= 0 * ci) & meq == 0) {
        aus <- list(b.unrestr = b, b.restr = b, 
            R2 = orig.R2, residuals = NULL, fitted.values = NULL, 
            weights = NULL, orig.R2 = orig.R2, 
            df.error = df.error, s2 = s2, Sigma = s2*V, 
            origmodel = NULL, ui = ui, ci = ci, iact = NULL, 
            restr.index = index, meq = meq, bootout = NULL)
            }
    else {
        ## equality restrictions involved or some inequality restrictions violated
        ## calculate restricted estimate
        aus <- solve.QP(Dmat = solve(V), dvec = solve(V, b), 
            Amat = t(uiw), bvec = ci, meq = meq)
        names(aus$solution) <- names(b)
        aus$solution[abs(aus$solution) < tol] <- 0
        ## initialize output list
        aus <- list(b.restr = aus$solution, b.unrestr = b, 
            R2 = NULL, residuals = NULL, fitted.values = NULL, 
            weights = NULL, orig.R2 = orig.R2, 
            df.error = df.error, s2 = s2, Sigma = s2*V, 
            origmodel = NULL, ui = ui, ci = ci, iact = aus$iact, 
            restr.index = index, meq = meq, bootout = NULL)
        ### R2 
        aus$R2 <- model[1,2:(g+1)]%*%t(t(aus$b.restr))/var.y
    }
    class(aus) <- c("orlm", "orest")
    aus
}
