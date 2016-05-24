# Bootstrap test for the Generalized Pareto distribution

`gpd.test` <-
function(x,J=999)
{
dname <- deparse(substitute(x))
if( min(x) < 0 )
	stop("There are negative observations. \nAll data must be positive real numbers.")
n <- length(x)
gammap <- .aml(x,n,ceiling(.2*n))
gamman <- .combined(x)[1]
Fn <- seq(1,n,1)/n
r1 <- .r1(x,n)
r2 <- .r2(x,n,Fn)
p.value1 <- sum( .cc1(n,gamman,J) < r1)/J
p.value2 <-  sum( .cc2(n,gammap,Fn,J) < r2 )/J
p.value <- max(p.value1,p.value2)
results <- list("p.value"=p.value, method = "Bootstrap test for the generalized Pareto distribution", data.name=dname)
    pvalues <- as.matrix(c(p.value1,p.value2),ncol=1)
    colnames(pvalues) <- c("p.value")
    rownames(pvalues) <- c("H_0^-: x has a gPd with negative shape parameter ","H_0^+: x has a gPd with positive shape parameter ")
    r <- as.matrix(c(r1,r2),ncol=1)
    colnames(r) <- c("R-statistic")	
class(results) = "htest"
return(list("boot.test"=results,"p.values"=cbind(pvalues,r)))
}

