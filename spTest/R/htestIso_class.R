#' An S4 class to hold output from a hypothesis test of isotropy. Extends the S3 class 'htest'.
#'
#' @slot p.value.finite A length-one numeric vector containing the p-value approximated by using a finite sample adjustment.
#' @slot sigma.hat A matrix containing the estimated asymptotic variance-covariance matrix.
#' @slot n.subblocks A length-one numeric vector containing the number of subblocks used in estimating the asymptotic variance-covariance.
#' @slot n.boot A length-one numeric vector containing the number of bootstrap samples used in estimating the asymptotic variance-covariance.
htestIso <- setClass("htestIso", slots = list(p.value.finite = "numeric", p.value.refl = "numeric", p.value.comp = "numeric", sigma.hat = "matrix", n.subblocks = "numeric", n.boot = "numeric"), contains = setOldClass('htest'))

#' Print hypothesis test results.
#'
#' Print the results from a nonparametric hypothesis test of isotropy/symmetry.
#'
#' @method print htestIso
#' @export
#' @keywords external
#' @param x An object of class 'htestIso' from the \code{\link{LuTest}}, \code{\link{GuanTestGrid}}, \code{\link{GuanTestUnif}}, or \code{\link{MaityTest}} functions
#' @param digits Number of significant digits in printed output.
#' @param prefix A character defining the prefix used for lines of output. Default is a tab.
#' @param ... Other arguments to print.
#'
#' @return Summary results of the hypothesis test.
#'
#' @examples
#' library(mvtnorm)
#' set.seed(1)
#' #number of rows and columns
#' nr <- 12
#' nc <- 18
#' n <- nr*nc
#' #Set up the coordinates
#' coords <- expand.grid(0:(nr-1), 0:(nc-1))
#' coords <- cbind(coords[,2], coords[,1])
#' #compute the distance between sampling locations
#' D <- as.matrix(dist(coords))
#' #Set parameter values for exponential covariance function
#' sigma.sq <- 1
#' tau.sq <- 0.0
#' phi <- 1/4
#' R <- sigma.sq * exp(-phi*D)
#' R <- R + diag(tau.sq, nrow = n, ncol = n)
#' #Simulate Gaussian spatial data
#' z <- rmvnorm(1,rep(0,n), R, method = c("chol"))
#' z <-  z-mean(z)
#' z <- t(z)
#' mydata <- cbind(coords, z)
#' mylags <-  rbind(c(1,0), c(0, 1), c(1, 1), c(-1,1))
#' myA <-  rbind(c(1, -1, 0 , 0), c(0, 0, 1, -1))
#' tr <- GuanTestGrid(mydata, delta = 1, mylags, myA, df = 2, window.dims = c(3,2), 
#' pt.est.edge = TRUE, sig.est.edge = TRUE, sig.est.finite = TRUE )
#' print.htestIso(tr) #print the summary
#' tr #can also print it using this command

## S3 method for objects of class 'htestIso'
print.htestIso = function(x, digits = getOption("digits"), prefix = "\t", ...) 
{
    cat("\n")
    cat(strwrap(x$method, prefix = prefix), sep = "\n")
    cat("\n")
    cat("data:  ", x$data.name, "\n", sep = "")
    
if(!is.null(x$p.value.refl))
{
    	cat(strwrap("Cramer-Von-Mises GoF Test:"))
    	cat("\n")
    	
    	out <- character()
        fp <- format.pval(x$p.value.refl, digits = max(1L, digits - 
            3L))
        out <- c(out, paste("p-value (reflection symmetry) = ", fp) )
        
        if (!is.null(x$p.value.comp)) 
        {
        	fp <- format.pval(x$p.value.comp, digits = max(1L, digits - 
            3L))
        	out <- c(out, paste("p-value (complete symmetry) = ", fp) )
        }    	
     
     cat( strwrap(paste(out, collapse = ", ") ) , "\n")
     cat(strwrap("alternative hypothesis: periodogram ratios do not follow F(2,2) distribution, i.e., process is not refl./comp. symmetric"))
 	
}
    
if(is.null(x$pvalue.refl))
{
    out <- character()
    if (!is.null(x$statistic)) 
        out <- c(out, paste(names(x$statistic), "=", format(signif(x$statistic, 
            max(1L, digits - 2L)))))
    if (!is.null(x$parameter)) 
        out <- c(out, paste(names(x$parameter), "=", format(signif(x$parameter, 
            max(1L, digits - 2L)))))
    if (!is.null(x$p.value)) {
        fp <- format.pval(x$p.value, digits = max(1L, digits - 
            3L))
        out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == 
            "<") fp else paste("=", fp)))
    }
 
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")

    out <- character()
    if (!is.null(x$p.value.finite)) {
        fp <- format.pval(x$p.value.finite, digits = max(1L, digits - 
            3L))
        out <- c(out, paste("p-value (finite adj.)", if (substr(fp, 1L, 1L) == 
            "<") fp else paste("=", fp)))
    }
    if (!is.null(x$n.subblocks)) {
        nblk <- x$n.subblocks
        out <- c(out, paste("number of subblocks: ", nblk, sep = ""))
    }
    if (!is.null(x$n.boot)) {
        nboot <- x$n.boot
        out <- c(out, paste("number of bootstraps: ", nboot, sep = ""))
    }
    
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")

    if (!is.null(x$alternative)) {
        cat("alternative hypothesis: ")
        if (!is.null(x$null.value)) {
            if (length(x$null.value) == 1L) {
                alt.char <- switch(x$alternative, two.sided = "not equal to", 
                  less = "less than", greater = "greater than")
                cat("true ", names(x$null.value), " is ", alt.char, 
                  " ", x$null.value, "\n", sep = "")
            }
            else {
                cat(x$alternative, "\nnull values:\n", sep = "")
                print(x$null.value, digits = digits, ...)
            }
        }
        else cat(x$alternative, "\n", sep = "")
    }
    
    if (!is.null(x$conf.int)) {
        cat(format(100 * attr(x$conf.int, "conf.level")), " percent confidence interval:\n", 
            " ", paste(format(c(x$conf.int[1L], x$conf.int[2L])), 
                collapse = " "), "\n", sep = "")
    }
    
    if (!is.null(x$estimate)) {
        cat("\n")
        cat("sample estimates: (lag value)\n")
        print(x$estimate, digits = digits, ...)
    }
    
    if (!is.null(x$sigma.hat)) {
        cat("\n")
        cat("estimated asymp. variance-covariance matrix:\n")
        print(x$sigma.hat, digits = digits, ...)
    }
    
    cat("\n")
    invisible(x)
}

}

#' @keywords internal
setMethod("print", signature = "htestIso", definition = print.htestIso)
