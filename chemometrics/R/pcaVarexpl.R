"pcaVarexpl" <-
function(X,a,center=TRUE,scale=TRUE,plot=TRUE, ...)
{
# Diagnostics of the explained variance per variable in PCA
#
# INPUT: 
# X ... X data
# a ... maximum number of components to use, e.g. the number of X-variables
# center ... centring of X (FALSE or TRUE)
# scale ... scaling of X (FALSE or TRUE)
# plot ... make plot with explained variance
# ... additional plotting arguments
#
# OUTPUT:
# ExplVar ... explained variance for each variable using "a" components


    if (a < 1 || a > min(nrow(X) - 1, ncol(X)))
            stop("Invalid number of components, a")

    ###################################################################################################
    X <- as.matrix(scale(X,center=center,scale=scale))

    # PCA
      if (ncol(X)>nrow(X)){
     	e <- eigen(X%*%t(X))
	T <- e$vectors %*% diag(sqrt(e$values))
        P <- t(X) %*% T %*% diag(1/e$values)
      }
      else {
	X_svd <- svd(X)
	P <- X_svd$v
	T <- X %*% P
      }

      varexpl=1-apply((X-T[,1:a]%*%t(P[,1:a]))^2,2,sum)/apply(X^2,2,sum)

    if (plot){
	barplot(varexpl,ylab="Explained variance",ylim=c(0,1),...)
		
    }

    list(ExplVar=varexpl)
}

