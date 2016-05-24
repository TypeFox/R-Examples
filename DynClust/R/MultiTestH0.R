
#' @name MultiTestH0
#' @title Statistical test of zero mean for dynamics
#' @description MultiTestH0 tests if each column vectors of a matrix seen as a noisy dynamic is of zero mean (\eqn{H_0}) or not. The multiple statistical test assumes known variance and is based on a multiple \eqn{\chi^2} test.
#'
#' @param proj.matrix a matrix whose colummns are tested to have zero mean or not.
#' @param data.var a numeric providing the known variance
#' @param thrs a numeric vector of thresholds specified as the \eqn{1-\alpha} quantiles of the multiple test under the null on each partition. If thrs=NULL (default) return the global p-value which is the minimum of all p-values obtained on each partition.
#' @return If thrs is provide returns a Boolean vector with length the number of columns of proj.matrix. Element j is TRUE if the null hypothesis (no difference with the null vector) is accepted for column j of proj.matrix. Otherwise, return one p-value per column.
#' @references Baraud Y., Huet S., Laurent B. \emph{Ann. Stat.} (2003) Durot C., Rozenholc Y. \emph{Methods Math. Stat.} (2006)
#' @author Tiffany Lieury, Christophe Pouzat, Yves Rozenholc
#' @export
MultiTestH0 <-
    function
### 'MultiTestH0' tests if the vector or matrix 'proj.matrix' given as an argument is significantly different from a null vector or matrix at a level alpha.
(
    proj.matrix,   
### a matrix of signals to be tested with as many columns as signals to test
    data.var,
### a numeric indicating the variance
    thrs=NULL
### a numeric vector indicating the thresholds (for each partitions) of the statistical multitest H0
### if thrs=NULL return the global p-value which is the minimum of all p-values obtained on each partitions
    ){
    MultiTestOneColumn <- function(i){
        X = proj.matrix[,i]
        test <- TRUE 
        pvalue = 1
        k <- 0
        while(test && (2^(k+1) <= length(X)+1)) {
            k <- k+1
            norm2.proj = sum(X[2^(k-1):(2^k-1)]^2)/data.var[i]
            if (is.null(thrs)) pvalue <- min(pvalue,1-pchisq(norm2.proj,df=2^k-2^(k-1)))
            else test <- (norm2.proj<=thrs[k])
        }
        if (is.null(thrs)) return(pvalue) else return(test)
    }   
    
    if (is.vector(proj.matrix)) proj.matrix=as.matrix(proj.matrix)
    
    if (length(data.var)==1) data.var=rep(data.var,dim(proj.matrix)[2])
    
                                        # test return a vector of length ncol(proj.matrix)
                                        # test[k] is TRUE if proj.matrix[,k] is statistically of zero mean
                                        # or if thrs is NULL return the p-values
    test <- sapply(1:dim(proj.matrix)[2],MultiTestOneColumn)
    
    test
}
