##=============================##
##  Bonferroni-Holm-Shaffer's  ##
##  correction for p.value     ##
##=============================##
#' Multiplicity adjustment by Bonferroni-Holm-Shaffer's rule
#' 
#' @param pValues 
#'     \code{numeric} vector of \emph{p}-vaules 
#' @return 
#'     \code{numeric} vector of corrected \emph{p}-vaules
#' @author Federico Mattiello  <federico.mattiello@@gmail.com>
#' @references 
#'     Shaffer J.P. (1986) 
#'     Modified Sequentially Rejective Multiple Test Procedures, 
#'     \emph{Journal of the American Statistical Association}, 
#'     \bold{81}, 826--831.\cr
#' @seealso \code{\link{p.adjust}}, \code{\link{p.adjust.methods}}
#' @export
#' @examples
#' set.seed(123)
#' p.raw <- runif(10, max = 0.2)
#' rbind(p.raw, p.adj = SOUP::BHS(p.raw))
#' 
BHS <- function (pValues)
{
    K <- length(pValues)
    nCmprs <- as.integer(round(.5 + sqrt(.25 + 2 * K)))
    S <- array(0, dim = c(K, nCmprs + 1))
    R <- array(NA, dim = dim(S))
    ch <- seq_len(nCmprs) * (seq_len(nCmprs) - 1)/2
    for (n in seq_len(nCmprs)[-1L])
    {
        x <- sort(unique(c(ch[1:n] + t(S[, n:1]))), decreasing = TRUE)
        y <- rep(NA, x[1])
        for (cc in 1:length(x[-1])) y[(x[cc] <= x[1]:1) & is.na(y)] <- x[cc]
        R[1:x[1], n + 1] <- y
        S[1:length(x), n + 1] <- x
    }# END:for-n
    den <- R[, nCmprs + 1]
    ord <- order(pValues)
    pBhs <- pValues[ord] * den
    pBhs[ord] <- ifelse(pBhs > 1, 1, pBhs)
    names(pBhs) <- names(pValues)
    return(pBhs)
}#=END=

##==============================================##
##  p.values correction with the "minP" method  ##
##==============================================##
#' FWE Adjustment Using Permutation and NPC
#' 
#' Multiplicity correction controlling the Family-Wise Error using the 
#' permutation \emph{p}-values and NonParametric Combination with \emph{minP} 
#' as combining function.
#' 
#' @title FWE Adjustment Using Permutation
#' @param Pmat 
#'     \code{matrix} of \emph{p-}values where comparisons are on the columns
#' @return 
#'     \code{numeric} vector of corrected p.values
#' @author Dario Basso and Federico Mattiello  <federico.mattiello@@gmail.com>
#' @references 
#'     Pesarin, F. and Salmaso, L. (2010) 
#'     \emph{Permutation Tests for Complex Data}.
#'     Wiley: United Kingdom \cr
#'     
#'     Finos, L. and Pesarin, F. and Salmaso, L. (2003) 
#'     Test combinati per il controllo della 
#'     {molteplicit\`a} mediante procedure di closed testing, 
#'     \emph{Statistica Applicata}, \bold{15}, 301--329.
#' @seealso \code{\link{p.adjust}}, \code{\link{p.adjust.methods}}
#' @export
#' @examples
#' set.seed(123)
#' P <- matrix(runif(1010), nrow = 101, ncol = 10,
#'   dimnames = list(c("p-obs", paste("p-*", 1L:100)), LETTERS[1L:10]))
#' P[1L, 1L:4] <- 1/100
#' FWEminP(P)
#' 
FWEminP <- function (Pmat)
{
    B <- dim(Pmat)[1] - 1
    nVars <- dim(Pmat)[2]
    pRes <- array(0, dim = c(nVars, 1))
    
    pObs <- Pmat[1, ]
    ord <- order(pObs, decreasing = FALSE)
    pOrd <- pObs[ord]
    PmatOrd <- Pmat[, ord]
    
    Tmat <- apply(PmatOrd, MARGIN = 1, FUN = min)
    pRes[1] <- mean(Tmat[-1] <= Tmat[1])
    
    if (nVars > 2)
    {
        for (j in seq_len(nVars)[-c(1, nVars)])
        {
            Tmat <- apply(PmatOrd[, j:nVars], MARGIN = 1, FUN = min)
            pRes[j] <- max(mean(Tmat[-1] <= Tmat[1]), pRes[j - 1])
        }# END:for
    } else {}# END: if - nVars
    
    pRes[nVars] <- max(pOrd[nVars], pRes[nVars - 1])
    pRes[ord] <- pRes
    rownames(pRes) <- colnames(Pmat)
    return(pRes)
}# END: FWEminP
##- documenting it
# # prompt(FWEminP, file = "man/FWEminP.Rd")
