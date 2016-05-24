##===========================##
##  Passing from statistics  ##
##  to permutation p.values  ##
##===========================##
#' Transforming Test Statistics to Permutation \emph{p}-Values
#' 
#' @title From Test Statistics To \emph{p}-Values
#' @param Tmat 
#'     \eqn{3}-dimensional (\eqn{2}-dimensional) \code{array} containing 
#'     the test statistics where the first horizontal slice (first row) 
#'     contains the observed value 
#' @param obsOnly 
#'     \code{logical}, if \code{FALSE} (default) the whole permutation 
#'     distribution of the computed \emph{p}-values is returned, if 
#'     \code{TRUE} only the observed ones are returned.
#' @return 
#'     if \code{obsOnly} is \code{FALSE} an \code{array} of the same dimension
#'     of the input matrix \code{Tmat}, otherwise only the first row 
#' @author Livio Finos, Aldo Solari and Federico Mattiello 
#'     <federico.mattiello@@gmail.com>
#' @seealso \code{\link[flip:permutationSpace]{permutationSpace}}, 
#'     \code{\link{NPC}}
#' @export
#' 
t2p <- function(Tmat, obsOnly = FALSE)
{
    ### dimensions check
    dimT <- dim(Tmat)
    ### number of permutations
    B <- dimT[1] #- 1
    ### observed value
    P <- apply(Tmat, MARGIN = seq_along(dimT)[-1], 
            FUN = function(z)
            {
                mean(z >= z[1])
            })
    ### p.values distribution
    if(obsOnly)
    {
        if(is.null(dim(P)))
        {
            names(P) <- dimnames(Tmat)[[2]]
        } else
        {
            dimnames(P) <- dimnames(Tmat)[-1]
        }# END:if-vector-P
    } else
    {
#        obs <- P
        Ptemp <- apply(-Tmat, MARGIN = seq_along(dimT)[-1], 
                FUN = function(x)
                {
                    rank(x, ties.method = "min")
                }
        )/B
        
#        P <- array(NA, dim = dimT)
        switch(as.character(length(dimT)),
            "2" = {
                P <- Ptemp
#                P[1L, ] <- obs
                # P[2:(B + 1), ] <- Ptemp[-1, ]
            },
            "3" = {
                P <- Ptemp
#                P[1L, , ] <- obs
#                P[2L:(B), , ] <- Ptemp
                # P[2L:(B + 1), , ] <- Ptemp[-1, , ]
            }
        )# END:switch-dim-Ptemp
        dimnames(P) <- dimnames(Tmat)
        # dimnames(P)[[1L]] <- c("p-obs", paste("p-*", 1L:B, sep = ""))        
    }# END:if-obsOnly
    
    return(P)
    
}# END: t2p
### documenting it
# # prompt(t2p, file = "man/t2p.Rd")
