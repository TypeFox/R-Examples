#### Generate asymptotic covariance matrix of correlation/covariance matrix by *column major*
asyCov <- function(x, n, cor.analysis=TRUE, dropNA=FALSE, as.matrix=TRUE,
                   acov=c("individual", "unweighted", "weighted"),
                   suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {
    if (is.list(x)) {

        ## whether to use average correlation matrix
        acov <- match.arg(acov, c("individual", "unweighted", "weighted"))
        if (acov != "individual") {
            ## Replace NA with 0 before calculations
            my.x <- lapply(x, function(x) {x[is.na(x)] <- 0; x} )

            if (acov=="unweighted") {
                ## Unweighted means = sum of correlations/no. of studies
                my.x <- Reduce("+", my.x)/pattern.na(x, show.na = FALSE)
            } else {
                my.x <- mapply("*", my.x, n, SIMPLIFY = FALSE)
                ## Weighted means = Cummulative sum of r*n/sum of n
                my.x <- Reduce("+", my.x)/pattern.n(x, n)
            }

            ## Make sure that the diagonals are 1 for correlation analysis
            if (cor.analysis) diag(my.x) <- 1
            ## Repeat it to k studies
            x <- replicate(length(x), my.x, simplify = FALSE)    
        }
        ## whether to use average correlation matrix

        
        # Check if it returns a matrix or a list
        if (as.matrix) {
          # No. of variables
          if (!identical(0, var(sapply(x, function(x){dim(x)[[1]]}))))
            stop("The dimensions of matrices should be the same in order to stack them together!")
          
          if (is.null(dimnames(x[[1]]))) {
            oldNames <- paste("x", 1:dim(x[[1]])[[1]], sep = "")
          } else {
            oldNames <- dimnames(x[[1]])[[1]]
          }
          psOldNames <- outer(oldNames, oldNames, paste, sep = "")          
          if (cor.analysis) {
            psOldNames <- vechs(psOldNames)
          } else {
            psOldNames <- vech(psOldNames)
          }
          # cov/var of psOldNames
          psCovNames <- paste("cov(", outer(psOldNames, psOldNames, paste, sep = "_"), ")", sep="")
          psCovNames <- vech(matrix(psCovNames, nrow=length(psOldNames), ncol=length(psOldNames)))

          ## Fixed a bug before v.0.7-0 that uses only the first n 
          ## out.list <- lapply(x, asyCov, n = n, cor.analysis = cor.analysis, silent = silent,
          ##                    suppressWarnings = suppressWarnings, dropNA = FALSE, ...)
          out.list <- mapply(asyCov, x, n = n, cor.analysis = cor.analysis, silent = silent,
                             suppressWarnings = suppressWarnings, dropNA = FALSE, ..., SIMPLIFY=FALSE)          
          #output
          out <- t(sapply(out.list, function(x) {(vech(x))}))
          
          ## http://stackoverflow.com/questions/17772916/using-correlation-matrices-for-meta-analytic-sem
          # Fixed a BUG when 2x2 matrices of correlation
          # The output are incorrectly arranged in a row rather than in a column.
          if (dim(out)[1]==1) out <- t(out)
          ## A list(NULL, psCovNames) is required
          dimnames(out) <- list(NULL, psCovNames)
          
          out
        } else {
          #output
          ## lapply(x, asyCov, n = n, cor.analysis = cor.analysis, silent = silent,
          ##                    suppressWarnings = suppressWarnings, dropNA = dropNA, ...)
          mapply(asyCov, x, n = n, cor.analysis = cor.analysis, silent = silent,
                             suppressWarnings = suppressWarnings, dropNA = dropNA, ..., SIMPLIFY=FALSE)          
        }
        
    } else {
      
        # Assumption: check the diagonals for missing data only
        miss.index <- is.na(Diag(x))
        x.new <- x[!miss.index, !miss.index]
        if (!is.pd(x.new)) 
            stop("x is not positive definite!\n")
        p <- nrow(x.new)
        if (is.null(dimnames(x))) {
            oldNames <- paste("x", 1:nrow(x), sep = "")
            cNames <- oldNames[!miss.index]
            dimnames(x.new) <- list(cNames, cNames)
        } else {
            oldNames <- dimnames(x)[[1]]
            cNames <- dimnames(x.new)[[1]]
        }
        # create matrix of labels for ps
        psOldNames <- outer(oldNames, oldNames, paste, sep = "")
        psMatnames <- outer(cNames, cNames, paste, sep = "")
        
        if (cor.analysis) {
            psOldNames <- vechs(psOldNames)
            acovName <- vechs(psMatnames)
            S <- mxMatrix("Stand", nrow = p, ncol = p, free = TRUE, values = jitter(vechs(cov2cor(x.new))), 
                name = "S", labels = acovName)
            D <- mxMatrix("Diag", nrow = p, ncol = p, free = TRUE, values = sqrt(Diag(x.new)), 
                name = "D")
            modelName <- "Asymptotic covariance matrix of correlation matrix"
        } else {
            psOldNames <- vech(psOldNames)
            acovName <- vech(psMatnames)
            S <- mxMatrix("Symm", nrow = p, ncol = p, free = TRUE, values = jitter(vech(x.new)), 
                name = "S", labels = acovName)
            D <- mxMatrix("Iden", nrow = p, ncol = p, name = "D")
            modelName <- "Asymptotic covariance matrix of covariance matrix"
        }
        expCov <- mxAlgebra(D %&% S, name = "expCov", dimnames = list(cNames, cNames))

        cModel <- mxModel(model = modelName, mxData(x.new, "cov", numObs = n), S, 
                          D, expCov, mxFitFunctionML(),
                          mxExpectationNormal("expCov", means=NA, dimnames = cNames))

        ## Return mx model without running the analysis
        if (run==FALSE) return(cModel)
        
        # try to run it with error message as output
        mxFit <- tryCatch(mxRun(cModel, silent=silent, suppressWarnings=suppressWarnings, ...), 
                          error = function(e) e)
        if (inherits(mxFit, "error")) {
            stop(print(mxFit))
        }
        # Need to multiply 2 to the inverse of Hessian matrix
        # http://openmx.psyc.virginia.edu/thread/360
        # Fixed a bug that all elements have to be inverted before selecting some of them
        acovS <- tryCatch(2 * solve(mxFit@output$calculatedHessian)[acovName, acovName, drop=FALSE], 
                              error = function(e) e)
        if (inherits(acovS, "error")) {
            stop(print(acovS))
        }

        ## No need to do it as [, drop=FALSE] has been added
        ## # When the dimensions are 1x1, dimnames are removed. Added them explicitly
        ## dimnames(acovS) <- list(acovName, acovName)

        if (dropNA) {          
          out <- acovS
        } else {          
          # oldNames include data for NA
          p <- length(psOldNames)
          out <- matrix(NA, nrow=p, ncol=p, dimnames=list(psOldNames, psOldNames))
          out[acovName, acovName] <- acovS
        }
        return(out)
    }
    
}
