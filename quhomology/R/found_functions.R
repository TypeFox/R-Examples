#These two functions were originally written by Prof. John Fox.
#see: http://socserv.mcmaster.ca/jfox/Courses/R-course-Berkeley/
#and slightly adapted for this project.
#He has given me permission (by email) to use his functions for this package.

GaussianElimination <- function(A, B, tol=sqrt(.Machine$double.eps),
                                 verbose=FALSE, fractions=FALSE){
  # A: coefficient matrix
  # B: right-hand side vector or matrix
  # tol: tolerance for checking for 0 pivot
  # verbose: if TRUE, print intermediate steps
  # fractions: try to express nonintegers as rational numbers
  # If B is absent returns the reduced row-echelon form of A.
  # If B is present, reduces A to RREF carrying B along.
  if ((!is.matrix(A)) || (!is.numeric(A)))
    stop("argument must be a numeric matrix")
  n <- nrow(A)
  m <- ncol(A)
  if (!missing(B)){
    B <- as.matrix(B)
    if (!(ifelse(T==all.equal(nrow(B),nrow(A)),T,F)) || !is.numeric(B))
      stop("argument must be numeric and must match the number of row of
           A")
    A <- cbind(A, B)
  }
  i <- j <- 1
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
      row <- A[i,]
      A <- A - outer(A[,j], row) # sweep
      A[i,] <- row # restore current row
      if (verbose) if (fractions) print(fractions(A))
      else print(round(A, round(abs(log(tol,10)))))
      j <- j + 1
      break
    }
    i <- i + 1
  }
  # 0 rows to bottom
  zeros <- which(apply(A[,1:m], 1, function(x) max(abs(x)) <= tol))
  if (length(zeros) > 0){
    zeroRows <- A[zeros,]
    A <- A[-zeros,]
    A <- rbind(A, zeroRows)
  }
  rownames(A) <- NULL
  if (fractions) fractions (A) else round(A, round(abs(log(tol, 10))))
  return(A)
}


rref <- function(A, tol=sqrt(.Machine$double.eps),verbose=FALSE,
                  fractions=FALSE){
  ## A: coefficient matrix
  ## tol: tolerance for checking for 0 pivot
  ## verbose: if TRUE, print intermediate steps
  ## fractions: try to express nonintegers as rational numbers
  ## Written by John Fox
  if ((!is.matrix(A)) || (!is.numeric(A)))
    stop("argument must be a numeric matrix")
  n <- nrow(A)
  m <- ncol(A)
  for (i in 1:min(c(m, n))){
    col <- A[,i]
    col[1:n < i] <- 0
    # find maximum pivot in current column at or below current row
    which <- which.max(abs(col))
    pivot <- A[which, i]
    if (abs(pivot) <= tol) next # check for 0 pivot
    if (which > i) A[c(i, which),] <- A[c(which, i),] # exchange rows
    A[i,] <- A[i,]/pivot # pivot
    row <- A[i,]
    A <- A - outer(A[,i], row) # sweep
    A[i,] <- row # restore current row
    if (verbose)
      if (fractions) print(fractions(A))
    else print(round(A,round(abs(log(tol,10)))))
  }
  for (i in 1:n)
    if (max(abs(A[i,1:m])) <= tol)
      A[c(i,n),] <- A[c(n,i),] # 0 rows to bottom
  if (fractions) fractions (A)
  else round(A, round(abs(log(tol,10))))
}