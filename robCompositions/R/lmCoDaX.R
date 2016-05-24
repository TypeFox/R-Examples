#' Classical and robust regression of non-compositional response on
#' compositional predictors
#' 
#' Delivers appropriate inference for regression of y on a compositional matrix
#' X.
#' 
#' Compositional explanatory variables should not be directly used in a linear
#' regression model because any inference statistic can become misleading.
#' While various approaches for this problem were proposed, here an approach
#' based on the isometric logratio (ilr) transformation is used.
#' 
#' @aliases lmCoDaX ilrregression robilrregression
#' @param y The response which should be non-compositional
#' @param X The compositional predictors as a matrix, data.frame or numeric
#' vector
#' @param method If robust, LTS-regression is applied, while with method equals
#' \dQuote{classical}, the conventional least squares regression is applied.
#' @return An object of class \sQuote{lts} or \sQuote{lm} and two summary
#' objects.
#' @author Peter Filzmoser
#' @seealso \code{\link{lm}}
#' @references Filzmoser, P., Hron, K., Thompsonc, K. (2012) Linear regression
#' with compositional explanatory variables.  \emph{Journal of Applied
#' Statistics}, 39, 1115-1128.
#' @keywords models
#' @export
#' @examples
#' 
#' ## How the total household expenditures in EU Member
#' ## States depend on relative contributions of 
#' ## single household expenditures:
#' data(expendituresEU)
#' y <- as.numeric(apply(expendituresEU,1,sum))
#' lmCoDaX(y, expendituresEU, method="classical")
#' lmCoDaX(y, expendituresEU, method="robust")
#' 
lmCoDaX <- function(y, X, method = "robust"){

	ilrregression <- function(X,y){
	# delivers appropriate inference for regression of y on a compositional matrix X
	# PF, Aug 18, 2010
	#
	# classical regression
		d <- data.frame(y=y,X=X)
		lmcla <- lm(y~.,data=d)
		lmcla.sum <- summary(lmcla)
	# ilr regressions:
		ilr.sum <- lmcla.sum
		for (j in 1:ncol(X)){
			Zj <- isomLR(cbind(X[,j],X[,-j]))
			dj <- data.frame(y=y,Z=Zj)
			res <- lm(y~.,data=dj)
			res.sum <- summary(res)
			if (j==1){
				ilr.sum$coefficients[1:2,] <- res.sum$coefficients[1:2,]
				ilr.sum$residuals <- res.sum$residuals
				ilr.sum$sigma <- res.sum$sigma
				ilr.sum$r.squared <- res.sum$r.squared
				ilr.sum$adj.r.squared <- res.sum$adj.r.squared
				ilr.sum$fstatistic <- res.sum$fstatistic
			}
			else{
				ilr.sum$coefficients[j+1,] <- res.sum$coefficients[2,]
			}
		}
		list(lm = lmcla, lm=lmcla.sum,ilr=ilr.sum)
	}
	
	robilrregression <- function(X,y){
	# delivers appropriate inference for robust regression of y on a compositional matrix X
	# PF, Aug 18, 2010
	#
	# classical regression
		d <- data.frame(y=y,X=X)
		lmcla <- robustbase::ltsReg(y~.,data=d)
		lmcla.sum <- summary(lmcla)
	# ilr regressions:
		ilr.sum <- lmcla.sum
		for (j in 1:ncol(X)){
			Zj <- isomLR(cbind(X[,j],X[,-j]))
			dj <- data.frame(y=y,Z=Zj)
			res <- robustbase::ltsReg(y~.,data=dj)
			res.sum <- summary(res)
			if (j==1){
				ilr.sum$coefficients[1:2,] <- res.sum$coefficients[1:2,]
				ilr.sum$residuals <- res.sum$residuals
				ilr.sum$sigma <- res.sum$sigma
				ilr.sum$r.squared <- res.sum$r.squared
				ilr.sum$adj.r.squared <- res.sum$adj.r.squared
				ilr.sum$fstatistic <- res.sum$fstatistic
			}
			else{
				ilr.sum$coefficients[j+1,] <- res.sum$coefficients[2,]
			}
		}
		list(lm = lmcla, lm=lmcla.sum,ilr=ilr.sum)
	}
	
	if( method =="classical"){
		reg <- ilrregression(X,y)	
	} else if(method=="robust"){
		reg <- robilrregression(X,y)
	}
	
	return(reg)

}

