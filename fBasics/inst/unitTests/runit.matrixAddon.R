
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port: 
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file 


################################################################################
# GENERATION:               DESCRIPTION:
#  matrix                    R  Creates a matrix from the given set of values
#   diag                     R  Creates a diagonal matrix or extracts diagonals
#   triang                   M  Extracs the lower tridiagonal part from a matrix
#   Triang                   M  Extracs the upper tridiagonal part from a matrix
#   pascal                   M  Creates a Pascal matrix
#   hilbert                  M  Creates a Hilbert matrix
#   colVec                   M  Creates a column vector from a data vector
#   rowVec                   M  Creates a row vector from a data vector
#  as.matrix                 R  Attempts to turn its argument into a matrix     
#  is.matrix                 R  Tests if its argument is a (strict) matrix
#  isPositiveDefinite        M  Checks if the matrix X is positive definite
#  makePositiveDefinite      M  Forces the matrix x to be positive definite
#  dimnames                  R  Retrieves or sets the dimnames of an object
#  colnames|rownames         R  Retrieves or sets the row or column names 
#  colIds|rowIds             M  ... use alternatively
#  colIds<-|rowIds<-         M  ... for assignments
# SUBSETS:                  DESCRIPTION:
#  dim                       R  Returns the dimension of a matrix object
#  ncol|nrow                 R  Counts columns|rows of a matrix object
#  length                    R  Counts elements of a matrix object
#   "["|"[["                 R  Subsets a matrix object
#   (Arith)                  R  Elementwise Arithmetic: + - * /
#   (Lops)                   R  Elementwise logical Ops: > < >= <= == !=
#  cbind|rbind               R  Augments a matrix object by columns|rows
#  na.omit                   R  Removes NA from a matrix object
# BASIC STATISTICS:         DESCRIPTION:
#  var                       R  Returns the variance matrix
#  cov                       R  Returns the covariance matrix
#  col|rowStats              B  calculates column|row statistics 
#   col|rowMeans             R  calculates column|row means
#   col|rowAvgs              B  calculates column|row averages
#   col|rowVars              B  calculates column|row variances
#   col|rowStdevs            B  calculates column|row standard deviations
#   col|rowSkewness          B  calculates column|row skewness 
#   col|rowKurtosis          B  calculates column|row kurtosis 
#   col|rowCumsums           B  calculates column|row cumulated sums 
# LINEAR ALGEBRA:           DESCRIPTION:
#  t                         R  Returns the transposed matrix
#  det                       R  Returns the determinant of a matrix
#  inv                       M  returns the inverse of a matrix, synonyme
#  chol2inv                  R  Returns the inverse of a matrix
#  norm                      M  returns the norm of a matrix
#  rk                        M  returns the rank of a matrix
#  tr                        M  returns the trace of a matrix
#  %*%                       R  Returns the product of two matrices
#  %x%                       R  Returns the Kronecker product
#  kron                      S  returns the Kronecker product
#  vec                       M  is the operator that stacks a matrix
#  vech                      M  is the operator that stacks the lower triangle
# MORE LINEAR ALGEBRA:      DESCRIPTION:
#  chol                      R  Returns the Cholesky factor matrix
#  eigen                     R  Returns eigenvalues and eigenvectors
#  svd                       R  Returns the singular value decomposition
#  kappa                     R  Returns the condition number of a matrix
#  qr                        R  Returns the QR decomposition of a matrix
#  solve                     R  Solves a system of linear equations
#  backsolve                 R  ... use when the matrix is upper triangular
#  forwardsolve              R  ... use when the matrix is lower triangular
# TIME SERIES               DESCRIPTION:
#  tslag                     R  Lagged/leading vector/matrix of selected orders 
#  .tslag1                      Internal Function used by tslag
#  pdl                       R  Regressor matrix for polynomial distributed lags  
# NOTES:                   WHERE YOU FIND THE FUCTIONS?
#                            R  Basic R Package
#                            B  Rmetrics fBasics Package
#                            M  This Rmetrics fMultivar Package
################################################################################


test.creation =
function()
{
    # Create Pascal Matrix:
    P = pascal(3)
    P
    
    # Create lower triangle matrix
    L = triang(P)
    L
    
    # Extract diagonal part
    diag(P)
    
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------

      
test.mathOps =
function()
{
    # Create Pascal Matrix:
    P = pascal(3)
    P
    
    # Add/Subtract/Multiply/Divide:  
    X = P
    
    # Multiply matrix with a constant
    3 * X
    
    # Multiply two matrices elementwise
    X * P                     
    
    # Multiplies rows/columns of a matrix by a vector
    X %*% diag(P)            
    diag(P) %*% X
    
    # Return Value:
    return()          
}


# ------------------------------------------------------------------------------


test.subsets =
function()
{          
    # Create Pascal Matrix:
    P = pascal(3)
    P
    
    # Operate on Subsets of a Matrix:
    n = 3
    i = 2
    j = 3
    D = diag(1:3)
    
    # Return the dimension of a matrix
    dim(P)                         
    
    # Get the last colum of a matrix
    P[, ncol(P)]                   
    
    # Delete a column of a matrix
    P[, -i]                      
        
    # Permute the columns of a matrix
    P[c(3, 1, 2), ]              
    
    # Augments matrix horizontally 
    cbind(P, D)
    
    # Return Value:
    return()                           
}


# ------------------------------------------------------------------------------

         
test.apply =
function()
{
    # Apply a function to all Elements of a Matrix: 
    
    # Create Pascal Matrix:
    P = pascal(3)
    P
    
    # Return square root for each element
    sqrt(P)
    
    # Exponentiate the matrix elementwise
    exp(P)
    
    # Compute the median of each column
    apply(P, 2, "median") 
    
    # Test on all elements of a matrix       
    all( P > 2 )   
    
    # test on any element in a matrix                
    any( P > 2 )
    
    # Return Value:
    return()                  
}


# ------------------------------------------------------------------------------


test.moreOperations =
function()
{       
    # More Matrix Operations:
    
    # Create Pascal Matrix:
    P = pascal(3)
    P
    
    # Create Diagonal Matrix:
    D = diag(1:3)
    
    # Return the product of two matrices
    P %*% D   
    
    # Return the Kronecker Product                     
    P %x% D                        
    
    # Return the transposed matrix
    t(P)                           
    
    # Return the inverse of a matrix
    inv(P)  
    
    # Return the norm of a matrix                      
    norm(P)                        
    
    # Return the determinante of a matrix
    det(P)                         
    
    # Return the rank of a matrix
    rk(P)                            
    
    # Return trace of a matrix
    tr(P)                          
    
    # Return the variance matrix
    var(P)     
    
    # Return the covariance matrix                   
    cov(P) 
    
    # Stack a matrix
    vec(P) 
    
    # Stack the lower triangle
    vech(P)
    
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------

 
test.linearAlgebra =
function()
{  
    # More Linear Algebra:
    
    # Create Pascal Matrix:
    P = pascal(3)
    P
    
    # Example Matrix and Vector
    X = P
    b = c(1, 2, 3)
    
    # Return the Cholesky factor matrix
    chol(X)                        
    
    # Return eigenvalues and eigenvectors
    eigen(X)                       
    
    # Return the singular value decomposition
    svd(X)                         
    
    # Return the condition number of a matrix
    kappa(X)                       
    
    # Return the QR decomposition of a matrix
    qr(X)                          
    
    # Solve a system of linear equations
    # ... use backsolve when the matrix is upper triangular
    # ... use forwardsolve when the matrix is lower triangular
    solve(X, b)  
    backsolve(Triang(X), b)
    solve(Triang(X), b)                 
    forwardsolve(triang(X), b) 
    solve(triang(X), b)
    
    # Return Value:
    return()
}  


################################################################################

