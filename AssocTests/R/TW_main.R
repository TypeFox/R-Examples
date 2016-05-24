##' Find the significant eigenvalues of a matrix.
##'
##' @title Tracy-Wisdom test
##' @param eigenvalues a numeric vector whose elements are the
##' eigenvalues of a matrix. The values should be sorted in the
##' descending order.
##' @param eigenL the number of the eigenvalues.
##' @param criticalpoint a numeric value corresponding to the
##' significance level. If the significance level is 0.05, 0.01,
##' 0.005, or 0.001, the critical point should be set to be \code{0.9793},
##' \code{2.0234}, \code{2.4224}, or \code{3.2724}, accordingly. The default is \code{2.0234}.
##' @return A list with class "\code{htest}" containing the following components:
##' \tabular{llll}{
##' \code{statistic} \tab \tab \tab \cr
##' \tab \tab \tab a vector of the Tracy-Wisdom statistics.\cr
##' \code{alternative} \tab \tab \tab \cr
##' \tab \tab \tab a character string describing the alternative hypothesis.\cr
##' \code{method} \tab \tab \tab \cr
##' \tab \tab \tab a character string indicating the type of test performed.\cr
##' \code{data.name} \tab \tab \tab \cr
##' \tab \tab \tab a character string giving the name of the data. \cr
##' \code{SigntEigenL} \tab \tab \tab \cr
##' \tab \tab \tab the number of the significant eigenvalues.
##' }
##' @author Lin Wang, Wei Zhang, and Qizhai Li.
##' @references N Patterson, AL Price, and D Reich. Population
##' Structure and Eigenanalysis. \emph{PloS Genetics}. 2006; 2(12):
##' 2074-2093.
##' @references CA Tracy and H Widom. Level-Spacing Distributions and
##' the Airy Kernel. \emph{Communications in Mathematical
##' Physics}. 1994; 159(1): 151-174.
##' @examples
##' tw(eigenvalues = c(5, 3, 1, 0), eigenL = 4, criticalpoint = 2.0234)
##' @export
tw <- function(eigenvalues, eigenL, criticalpoint=2.0234)
{
    a <- deparse(substitute(eigenvalues))
    dex <- which(eigenvalues <= 1e-8)
    eigenvalues[dex] <- 1e-8

    L1 <- rev(cumsum(rev(eigenvalues)))
    L2 <- rev(cumsum(rev(eigenvalues^2)))
    N <- eigenL:1
    S2 <- N^2*L2/(L1^2)
    v <- N*(N+2)/(S2-N) # Effective number of markers

    L <- N*eigenvalues/L1

    v.st <- sqrt(v-1)
    N.st <- sqrt(N)

    mu  <- (v.st+N.st)^2/v
    sig <- (v.st+N.st)/v * (1/v.st+1/N.st)^(1/3)

    twstat <-(L-mu)/sig

    #sink(output)
    #cat("TWstat = ", twstat, '\n')
    #sink()

    d <- which(twstat < criticalpoint)[1]

    if (length(d)==0)
    {
        d <- -100
    }else
    {
        d <- d-1
    }

    structure( 
    list(statistic = c(TW = twstat), 
        alternative = "the eigenvalue is significant", 
        method = "Tracy-Wisdom test", 
        data.name = a,
        SigntEigenL = d
        ), 
    .Names=c("statistic", "alternative", "method", "data.name", "SigntEigenL"), 
    class="htest"
    )
}
