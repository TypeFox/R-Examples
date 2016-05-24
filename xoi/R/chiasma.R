## chiasma.R

#' Estimate chiasma distribution from crossover counts
#'
#' Fit several models, with an assumption of no chromatid interference, to
#' crossover count data to obtain fitted distributions of the number of
#' chiasmata.
#'
#' @details
#' See Broman and Weber (2000) for details of the method.
#'
#' We use R's \code{\link[stats]{integrate}} function for numerical integrals,
#' \code{\link[stats]{optimize}} for optimizing the likelihood, and
#' \code{\link[stats]{uniroot}} for identifying the endpoints of the likelihood
#' support interval.
#'
#' @param xo Vector of non-negative integers; the number of crossovers in a set
#' of meiotic products.
#' @param max.chiasma Maximum number of chiasmata to allow.
#' @param n.iter Maximum number of iterations in the EM algorithm.
#' @param tol Tolerance for convergence of the EM algorithm.
#' @param verbose If TRUE, print number of interations for each of the 4 models
#' at the end.
#' @return A list with three components.
#'
#' First, a matrix containing the observed distribution of the numbers of
#' crossovers, followed by the fitted distributions under the Poisson model,
#' the truncated Poisson model (assuming an obligate chiasma), the obligate
#' chiasma model, and the freely varying model.  In all cases we assume no
#' chromatid interference.
#'
#' Second, a matrix containing the estimated distributions of the number of
#' chiasmata on the four-strand bundle for the above four models.
#'
#' Third, the estimated average number of crossovers under the Poisson and
#' truncated Poisson models.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{fitGamma}}, \code{\link[qtl]{fitstahl}},
#' \code{\link{countxo}}
#' @references Broman, K. W. and Weber, J. L. (2000) Characterization of human
#' crossover interference. \emph{Am. J. Hum. Genet.} \bold{66}, 1911--1926.
#'
#' Broman, K. W., Rowe, L. B., Churchill, G. A. and Paigen, K. (2002) Crossover
#' interference in the mouse. \emph{Genetics} \bold{160}, 1123--1131.
#'
#' Yu, K. and Feinbold, E. (2001) Estimating the frequency distribution of
#' crossovers during meiosis from recombination data.  \emph{Biometrics}
#' \bold{57}, 427--434.
#' @keywords models
#' @examples
#'
#' data(bssbsb)
#'
#' # estimated number of crossovers on chr 1
#' nxo <- countxo(bssbsb, chr=1)
#'
#' # estimate chiasma distribution
#' \dontrun{chiasma(nxo)}
#' \dontshow{chiasma(nxo, tol=0.001)}
#'
#' @useDynLib xoi
#' @export
chiasma <-
    function(xo, max.chiasma=max(xo)*2+5, n.iter=10000, tol=1e-6, verbose=FALSE)
{
    n.xo <- length(xo)

    res <- .C("chiasma",
              as.integer(xo),
              as.integer(n.xo),
              as.integer(max.chiasma),
              p.ch = as.double(rep(0,(max.chiasma+1)*4)),
              p.xo = as.double(rep(0,(max.chiasma+1)*4)),
              lambda = as.double(c(0,0)),
              as.double(rep(0,(max.chiasma+1)*n.xo+2+2*(max.chiasma+1))),
              n.iter = as.integer(c(n.iter,0,0,0,0)),
              as.double(tol),
              PACKAGE="xoi")

    if(verbose)
        cat("Done!  number of iterations = ",
            paste(as.character(res$n.iter[-1]),collapse=" "), "\n")

    xo.table <- rbind(table(factor(xo,levels=0:max.chiasma))/length(xo),
                      matrix(res$p.xo,nrow=4,byrow=TRUE))
    ch.table <- matrix(res$p.ch,nrow=4,byrow=TRUE)
    colnames(ch.table) <- 0:(ncol(ch.table)-1)

    rownames(ch.table) <- c("truncPois","Pois","oblchi","free")
    rownames(xo.table) <- c("observed", "truncPois","Pois","oblchi","free")

    out <- list(xo.table=xo.table*n.xo,ch.table=ch.table,
                lambda = res$lambda)
    attr(out, "n.iter") <-  res$n.iter[-1]
    out
}
