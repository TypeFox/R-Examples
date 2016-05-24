GaussianElimination <- function(A, B, tol=sqrt(.Machine$double.eps),
    verbose=FALSE){
    # A: coefficient matrix
    # B: right-hand side vector or matrix
    # tol: tolerance for checking for 0 pivot
    # verbose: if TRUE, print intermediate steps
    # If B is absent returns the reduced row-echelon form of A.
    # If B is present, reduces A to RREF carrying B along.

    ## routine by John Fox modified to output pivot column numbers
    ## by Ulrike Groemping
    pivot.num <- integer(0)
    if ((!is.matrix(A)) || (!is.numeric(A)))
        stop("argument must be a numeric matrix")
    n <- nrow(A)
    m <- ncol(A)
    if (!missing(B)){
        B <- as.matrix(B)
        if (!(nrow(B) == nrow(A)) || !is.numeric(B))
          stop("argument must be numeric and must match the number of row of A")
        A <- cbind(A, B)
        }
    i <- j <- k <- 1
    while (i <= n && j <= m){
        while (j <= m){
            currentColumn <- A[,j]
            currentColumn[1:n < i] <- 0
            # find maximum pivot in current column at or below current row
            which <- which.max(abs(currentColumn))
            pivot <- currentColumn[which]
            if (abs(pivot) <= tol) { # check for 0 pivot
                j <- j + 1
                next
                }
            if (which > i) A[c(i, which),] <- A[c(which, i),] # exchange rows
            A[i,] <- A[i,]/pivot # pivot
            pivot.num <- c(pivot.num, j)
            k <- k + 1   
            row <- A[i,]
            A <- A - outer(A[,j], row) # sweep
            A[i,] <- row # restore current row
            if (verbose) print(round(A, round(abs(log(tol,10)))))
            j <- j + 1
            break
            }
        i <- i + 1
        }
     # 0 rows to bottom
    zeros <- which(apply(matrix(A[,1:m],nrow(A),m), 1, function(x) max(abs(x)) <= tol))
    if (length(zeros) > 0){
        zeroRows <- A[zeros,]
        A <- A[-zeros,]
        A <- rbind(A, zeroRows)
        }
    rownames(A) <- NULL
    list(E=round(A, round(abs(log(tol, 10)))), pivot=pivot.num, rank=length(pivot.num))
    }

RREF <- function(X, ...) GaussianElimination(X, ...)
    # returns the reduced row-echelon form of X and the pivot column numbers
