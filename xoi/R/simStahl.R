## simStahl.R

#' Simulate crossover locations under the Stahl model
#'
#' Simulate crossover locations under the Stahl model.
#'
#' The Stahl model is an extension to the gamma model, in which chiasmata occur
#' according to two independent mechanisms.  A proportion \eqn{p} come from a
#' mechanism exhibiting no interference, and a proportion 1-\eqn{p} come from a
#' mechanism in which chiasma locations follow a gamma model with interference
#' parameter \eqn{\nu}{nu}.
#'
#' @param n.sim Number of meiotic products to simulate.
#' @param nu The interference parameter in the gamma model.
#' @param p The proportion of chiasmata coming from the no-interference
#' mechanism.
#' @param L Chromosome length (in cM).
#' @param obligate_chiasma Require an obligate chiasma (requires \code{nu} to
#' be an integer; if nu is not an integer, it is rounded.
#' @param n.bins4start We approximate the distribution of the location of the
#' first crossover from the mechanism exhibiting interference using a even grid
#' with this many bins. (Only if \code{nu} is not an integer.)
#' @return A vector of length \code{n.sim}, each element being empty (for
#' products with no crossovers) or a vector of crossover locations, in cM.  An
#' attribute, \code{L}, contains the chromosome length in cM.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link{fitGamma}}, \code{\link[qtl]{sim.cross}}
#' @references Copenhaver, G. P., Housworth, E. A. and Stahl, F. W. (2002)
#' Crossover interference in Arabidopsis. \emph{Genetics} \bold{160},
#' 1631--1639.
#'
#' Housworth, E. A. and Stahl, F. W. (2003) Crossover interference in humans.
#' \emph{Am J Hum Genet} \bold{73}, 188--197.
#' @keywords datagen
#' @examples
#'
#' # simulations with no interference, chromosome of length 80 cM
#' xoNI <- simStahl(100, nu=1, p=0, L=80)
#'
#' # simulations under gamma model with nu=7.6
#' xogamma <- simStahl(100, nu=7.6, p=0, L=80)
#'
#' # simulations under Stahl model with nu=7.6, p=0.1
#' xostahl <- simStahl(100, nu=7.6, p=0.1, L=80)
#'
#' # simulations under chi-square model with nu=11 (m=10) and obligate chiasma
#' xo_oblchi <- simStahl(100, nu=11, p=0, L=80, obligate_chiasma=TRUE)
#'
#' # simulations under Stahl model with nu=11, p=0.1, and obligate chiasma
#' xo_oblchi_stahl <- simStahl(100, nu=11, p=0.1, L=80, obligate_chiasma=TRUE)
#' @importFrom stats qpois dpois
#' @useDynLib xoi
#' @export
simStahl <-
    function(n.sim, nu=1, p=0, L=100,
             obligate_chiasma=FALSE, n.bins4start=10000)
{
    if(nu <= 0) stop("nu should be positive.")
    if(p < 0 || p > 1) stop("p should be in [0,1].")
    if(p==1) { # if p is 1, might as well take nu == 1
        nu <- 1L
        p <- 0
    }
    if(n.sim <= 0) stop("n should be a positive integer.")
    if(L < 0) stop("L should be positive.")
    if(n.bins4start < 1000) {
        warning("n.bins4start should be large.  Using 1000.")
        n.bins4start <- 1000
    }
    if(nu %% 1 < 1e-6) { # if nu is very close to integer, just make it an integer
        nu <- as.integer(nu)
    }

    max.nxo <- qpois(1-1e-10, L)*10

    if(obligate_chiasma || is.integer(nu)) { # use integer version

        if(!is.integer(nu)) {
            nu <- as.integer(round(nu))
            warning("Simulations with obligate chiasma require that nu be an integer; ",
                    "rounding nu to ", nu)
        }

        if(obligate_chiasma)
            Lstar <- calc_Lstar(L, nu-1, p)
        else Lstar <- L

        out <- .C("R_simStahl_int",
                  as.integer(n.sim),
                  as.integer(nu-1),
                  as.double(p),
                  as.double(L),
                  as.double(Lstar),
                  nxo = as.integer(rep(0,n.sim)), # number of crossovers
                  loc = as.double(rep(0,n.sim*max.nxo)),
                  as.integer(max.nxo),
                  as.integer(obligate_chiasma),
                  PACKAGE="xoi")

        out <- lapply(as.data.frame(rbind(out$nxo, matrix(out$loc, nrow=max.nxo))),
                      function(a) {if(a[1]==0) return(numeric(0)); a[(1:a[1])+1] })
    }
    else {
        out <- .C("simStahl",
                  as.integer(n.sim),
                  as.double(nu),
                  as.double(p),
                  as.double(L/100),
                  nxo = as.integer(rep(0,n.sim)), # number of crossovers
                  loc = as.double(rep(0,n.sim*max.nxo)),
                  as.integer(max.nxo),
                  as.integer(n.bins4start),
                  PACKAGE="xoi")

        out <- lapply(as.data.frame(rbind(out$nxo, matrix(out$loc*100, nrow=max.nxo))),
                      function(a) {if(a[1]==0) return(numeric(0)); a[(1:a[1])+1] })
    }

    attr(out, "L") <- L
    names(out) <- NULL
    out
}

# calculate reduced length to give expected no. chiasmata when conditioning on >= 1
#
# L and Lstar are in cM
calc_Lstar <-
    function(L, m=0, p=0)
{
    if(L <= 50) stop("Must have L > 50")
    if(m < 0) stop("Must have m==0")
    if(!is.integer(m) && m %% 1 > 1e-8) {
        warning("m must be an non-negative integer; rounding")
        m <- round(m)
    }
    if(p < 0 || p > 1)
        stop("p must be in [0, 1]")
    if(p==1) { # if p == 1, might as well take m=0, p=0
        m <- 0
        p <- 0
    }

    func_to_zero <- function(Lstar, L, m=0, p=0) {
        if(m==0)
            denom <- 1 - exp(-Lstar/50)
        else {
            lambda1 <- Lstar/50 * (m+1) * (1-p)
            lambda2 <- Lstar/50 * p
            denom <- 1 - sum(dpois(0:m, lambda1) * (m+1 - (0:m))/(m+1)) * exp(-lambda2)
        }

        2*L - 2*Lstar / denom
    }

    stats::uniroot(func_to_zero, c(1e-8, L), L=L, m=m, p=p,
                   tol=sqrt(.Machine$double.eps))$root
}
