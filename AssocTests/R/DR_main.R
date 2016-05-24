##' Conduct the distance regression with or without the adjustment of
##' the covariates to detect the association between a distance matrix
##' and some independent variants of interest.
##'
##' The pseudo \emph{F} statistic based on the distance regression with
##' or without the adjustment of the covariates detects the
##' association between a distance matrix and some independent
##' variants of interest. A distance matrix can be transformed into a
##' similarity matrix easily.
##' @title Distance regression
##' @param simi.mat a similarity matrix among the subjects.
##' @param null.space a numeric vector to show the column numbers of
##' the null space in the \code{x.mat}.
##' @param x.mat the covariate matrix which combines the null space
##' and the matrix of interest.
##' @param permute logical. If \code{TRUE}, the Monter Carlo sampling is used
##' without replacement; otherwise, with replacement. The default is
##' \code{TRUE}.
##' @param n.MonterCarlo the number of times for the Monter Carlo
##' procedure. The default is \code{1000}.
##' @param seed if it is not \code{NULL}, set the random number generator
##' state for random number generation. The default is \code{NULL}.
##' @return A list with class "\code{htest}" containing the following components:
##' \tabular{llll}{
##' \code{statistic} \tab \tab \tab \cr
##' \tab \tab \tab the observed value of the test statistic.\cr
##' \code{p.value} \tab \tab \tab \cr
##' \tab \tab \tab the p-value for the test.\cr
##' \code{alternative} \tab \tab \tab \cr
##' \tab \tab \tab a character string describing the alternative hypothesis.\cr
##' \code{method} \tab \tab \tab \cr
##' \tab \tab \tab a character string indicating the type of test performed.\cr
##' \code{data.name} \tab \tab \tab \cr
##' \tab \tab \tab a character string giving the names of the data.
##' }
##' @author Lin Wang, Wei Zhang, and Qizhai Li.
##' @references Q Li, S Wacholder, DJ Hunter, RN Hoover, S Chanock, G
##' Thomas, and K Yu. Genetic Background Comparison Using
##' Distance-Based Regression, with Applications in Population
##' Stratification Evaluation and Adjustment. \emph{Genetic
##' Epidemiology}. 2009; 33(5): 432-441.
##' @references J Wessel and NJ Schork. Generalized Genomic
##' Distance-Based Regression Methodology for Multilocus Association
##' Analysis. \emph{American Journal of Human Genetics}. 2006; 79(5):
##' 792-806.
##' @references MA Zapala and NJ Schork. Multivariate Regression
##' Analysis of Distance Matrices for Testing Associations Between
##' Gene Expression Patterns and Related Variables. \emph{Proceedings
##' of the National Academy of Sciences of the United States of
##' America}. 2006; 103(51): 19430-19435.
##' @examples
##' data(drS.eg)
##' null.space <- 1
##' x.mat <- matrix(c(rep(1, 600), rep(0, 200)), ncol=2)
##' dr(drS.eg, null.space, x.mat, permute = TRUE, n.MonterCarlo = 50, seed = NULL)
##' @export
dr <- function(simi.mat, null.space, x.mat, permute=TRUE, n.MonterCarlo=1000, seed=NULL)
{
    if (!is.null(seed))
    {
        set.seed(seed)
    }

    x1 <- x.mat[,null.space]

    if (length(null.space)==1)
    {
        x1 <- matrix(x1, ncol=1)
    }

    x.hat <- x.mat %*% solve(t(x.mat)%*%x.mat) %*% t(x.mat)

    x1.hat <- x1 %*% solve(t(x1)%*%x1) %*% t(x1)

    n <- nrow(simi.mat)

    I.n <- diag(n)

    cent <- I.n - matrix(1,nrow=n,ncol=n)/n

    i.x  <- I.n - x.hat

    i.x1 <- I.n - x1.hat

    Q <- i.x1 %*% cent %*% simi.mat %*% cent %*% i.x1

    alter.hat <- x.hat - x1.hat

    F.obs <- sum(alter.hat*Q) / sum(i.x*Q)

    U <- 1:n

    F.star <- rep(NA, n.MonterCarlo)

    for (i in 1:n.MonterCarlo)
    {
        id.sam <- sample(U, replace=!permute)

      	Q.star <- Q[id.sam, id.sam]

	      F.star[i] <- sum(alter.hat*Q.star) / sum(i.x*Q.star)
    }

    pv <- sum(F.star >= F.obs)/n.MonterCarlo

    a <- deparse(substitute(simi.mat))
    b <- deparse(substitute(x.mat))
    
    structure( 
    list(statistic=c(F = F.obs), 
        p.value = pv, 
        alternative = "the pair-wise similarity is influenced by the variants of interest", 
        method = "Distance regression", 
        data.name = paste(a, "and", b, sep=" ")
        ), 
    .Names=c("statistic", "p.value", "alternative", "method", "data.name"), 
    class="htest"
    ) 
}
