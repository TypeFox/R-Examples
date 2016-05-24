rsmvnorm <-
function(R=R,cor.matrix=cor.matrix)
{
 if(!is.numeric(R) | R<1) 
    stop("'R' must be greater than or equal to one")
 R <- as.integer(R)
 if(!is.numeric(cor.matrix))  
    stop("'cor.matrix' must be numeric")
 cor.matrix <- as.matrix(cor.matrix)
 if(!isSymmetric(cor.matrix)) 
    stop("'cor.matrix' must be a symmetric matrix") 
 if(any(diag(cor.matrix)!=1)) 
    stop("the diagonal elements of 'cor.matrix' must be equal to one")
 if(any(cor.matrix> 1) | any(cor.matrix< -1))
    stop("all the elements of 'cor.matrix' must be on [-1,1]")
 if(any(eigen(cor.matrix,symmetric=TRUE,only.values=TRUE)$values<=0))
    stop("'cor.matrix' must be a positive definite matrix") 
 p <- ncol(cor.matrix)
 ans <- matrix(rnorm(R*p),R,p)%*%chol(cor.matrix)
 ans
}
