#'Test the input correlation matrix weather it is in the feasible bounds.
#'
#'The function uses the lower and upper correlation bounds from results of the
#'function \code{valid_input_cormat} to test the users' input correlation
#'matrix.
#'
#'The function uses the lower and upper correlation bounds from results of the
#'function \code{valid_input_cormat} to test the users' input correlation
#'matrix. If all the elements in the input correlation matrix are in the bounds,
#'the function will return \code{TRUE}, otherwise it will print out the
#'elements' positions in the input correlation matrix which are out of the lower
#'and upper bounds.
#'@inheritParams valid_input_cormat
#'@param cor_matrix The input correlation matrix to be checked.
#'@return If all the elements of the input correlation matrix all in the bounds
#'  the function will return \code{TRUE}, otherwise it will print messages about
#'  the out of bounds elements' positions and then give an error message to
#'  users.
#'
#'
#'@seealso \code{\link{BoundingRA}}, \code{\link{valid_input_cormat}},
#'  \code{\link{genNORTARA}}
#'@note Because of the random samples, the results of the function may be a
#'  little different each time.
#'@examples
#'\dontrun{
#'invcdfnames <- c("qt","qpois","qnorm")
#'paramslists <- list(
#'                m1 = list(df = 3),
#'                m2 = list(lambda = 5),
#'                m3 = list(mean = 0, sd = 1)
#'                  )
#'cor_matrix_correct <- matrix(c(1,0.5,-0.3,0.5,1,0.4,-0.3,0.4,1), 3)
#'cor_matrix_wrong <- matrix(c(1,0.94,-0.3,0.94,1,0.4,-0.3,0.4,1), 3)
#'check_input_cormat(invcdfnames, paramslists, cor_matrix_correct)
#'check_input_cormat(invcdfnames, paramslists, cor_matrix_wrong)
#'}
#'@export
check_input_cormat <- function(invcdfnames, paramslists, cor_matrix){

 if (ncol(cor_matrix) != length(invcdfnames)) {
    stop("The dims of input cor_matrix does not match the length of invcdfnames!")
  }
 if (isSymmetric(cor_matrix) == FALSE) {
  stop("Input cor_matrix should be symmetric!")
 }
 ndim <- ncol(cor_matrix)
 upbound_cormat <- valid_input_cormat(invcdfnames, paramslists)$max_valid_cormat
 lowbound_cormat <- valid_input_cormat(invcdfnames, paramslists)$min_valid_cormat
 outofboundcounts <- 0
 for (i in 1:(ndim - 1))
   for (j in (i + 1):ndim){
     if ( cor_matrix[i,j] < lowbound_cormat[i,j] || cor_matrix[i,j] > upbound_cormat[i,j]) {
       outofboundcounts <- outofboundcounts + 1
       cat("The input cor_matrix[",i,",",j,"] value should be greater than ",
           lowbound_cormat[i,j]," and less than ",upbound_cormat[i,j],"!\n"
            )
    }
   }
 if (outofboundcounts > 0) {
   stop ("Please correct the above values and try again!")
 }
 TRUE
 }
