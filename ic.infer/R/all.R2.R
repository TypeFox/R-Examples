all.R2 <- function (covmat, ui, ci = NULL, index = 2:ncol(covmat), 
    meq = 0, tol = sqrt(.Machine$double.eps), ...) 
{
    ## check for admissible model
    if (!(is.matrix(covmat))) 
        stop("ERROR: covmat must be a covariance matrix.")
    if (!(ncol(covmat)==nrow(covmat))) 
        stop("ERROR: covmat must be a square matrix.")
    if (!all(eigen(covmat,TRUE,only.values=TRUE)$values>0)) 
        stop("ERROR: covmat must be positive definite.")
    ## number of regressors 
    g <- ncol(covmat)-1
    ## check for admissible restrictions
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
    ## check for feasibility
    ## and initialize storage vector for all R2s
    if (!is.numeric(try(R2s <- rep(0, 2^g)))) 
        stop("too many regressors in model, function all.R2 needs more storage than is available")
    R2s[1] <- 0
    R2s[2^g] <- orlm(covmat, ui, ci = ci, index = index, 
        meq = meq, tol = tol, df.error=10)$R2
                    ## df.error is irrelevant, but required by orlm
    ## prepare data for calculation of sub models

    ## initialize running index
    lauf <- 1
    for (j in 1:(g - 1)) {
        sets <- nchoosek(g, j)
        for (k in 1:ncol(sets)) {
        ## reduce model to these j regressors 
        ## reduce matrix ui accordingly
        ## run orlm on reduced model for R2
            lauf <- lauf + 1
            covj <- covmat[c(1, sets[, k] + 1), c(1, sets[, k] + 1)]
            if (length(which((index-1) %in% sets[, k])) > 0) {
                uij <- matrix(ui[, which((index-1) %in% sets[, 
                  k])], nrow(ui), length(which((index-1) %in% 
                  sets[, k])))
                indexj <- index[which((index-1) %in% sets[, 
                  k])] - 1
                indexj <- which(sets[, k] %in% (indexj))
                meqj <- 0
                if (meq > 0) 
                  meqj <- length(which((index[1:meq]-1) %in% 
                    sets[, k]))
                hilf <- RREF(t(uij))
                if (hilf$rank > 0) {
                  uij <- matrix(uij[hilf$pivot, ], hilf$rank, 
                    ncol(uij))
                  cij <- ci[hilf$pivot]
                  R2s[lauf] <- orlm(covj, uij, ci = cij, 
                    index = indexj+1, meq = meqj, tol = tol, df.error=10)$R2
                    ## df.error is irrelevant, but required by orlm
                }
                else R2s[lauf] <- covj[1,-1]%*%solve(covj[-1,-1],covj[-1,1])/covj[1,1]
            }
            else R2s[lauf] <- covj[1,-1]%*%solve(covj[-1,-1],covj[-1,1])/covj[1,1]
        }
    }
    R2s
}
