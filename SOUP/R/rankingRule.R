################################################################################
### Title: rankingRule.R                                                     ###
###                                                                          ###
### Project: SOUP                                                            ###
###                                                                          ###
### Version: 0.1 - 13/lug/2013                                               ###
###                                                                          ###
### Description: Constructs the ranking with a suitable rule and making      ###
###     use of pairwise hypotheses p.values                                  ###
###                                                                          ###
### Authors: Federico Mattiello <Federico.Mattiello@UGent.be>                ###
###                                                                          ###
### Maintainer: Federico Mattiello <Federico.Mattiello@UGent.be>             ###
###                                                                          ###
### Versions:                                                                ###
###     > 0.1 - 13/lug/2013 - 15:24:00:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################
###

#' Performs the ranking
#' 
#' 
#' @param permSpace 
#'     object of the class \code{\linkS4class{PermSpace}} containing the 
#'     permutation space of the test statistic
#' @param alpha 
#'     \code{numeric} significance level to be employed, in case it is 
#'     a vector, a ranking is computed for each value of \code{alpha}
#' @param multAdjMethod 
#'     multiplicity adjustment method to be used for the pairwise hypotheses 
#'     \emph{p}-values
#' @param groupsLabs
#'     \code{character} vector containing the groups' labels
#' @return 
#'     object of the class \code{\linkS4class{RankResults}} containing 
#'     results of the ranking procedure and other information
#' @author Federico Mattiello <federico.mattiello@@gmail.com>
#' @export
#' 
rankingRule <- function(permSpace, alpha = 0.05, multAdjMethod, groupsLabs)
{
    ### multiplicity correction check
    multAdjMethod <- match.arg(multAdjMethod, 
            c("FWEminP", "BHS",  p.adjust.methods))
    
    ### check input class
    if(!is(permSpace, "PermSpace"))
    {
        stop("input object is not of class \"PermSpace\"")
    } else {}# END: if - PermSpace
    
    ### converting "p.values" permSpace into a "list" if it is not
    if(is.list(permSpace@P.H0Low))
    {
        P.H0Low <- permSpace@P.H0Low
        P.H0Gre <- permSpace@P.H0Gre
    } else
    {
        P.H0Low <- list(permSpace@P.H0Low)
        P.H0Gre <- list(permSpace@P.H0Gre)
    }# END: if - permSpace@stats
    
    ### global variables
    K <- ncol(P.H0Low[[1L]])
    G <- as.integer(round(.5 + sqrt(.25 + 2 * K)))
    
#    tmp <- dimnames(P.H0Low[[1L]])[[2L]]
#    if((missing(groupsLabs) || is.null(groupsLabs)) && !is.null(tmp))
#    {
#        groupsLabs <- tmp
#    } else
#    {
#        groupsLabs <- seq_len(G)
#    }# END: ifelse - missing groupsLab
    # CM <- .DesM(rep.int(1L, C))
    
    ### outputs
    rankRes <- vector("list", length(P.H0Low))
    
    ### looping along "stats" list
    for(i in seq_along(P.H0Low))
    {
        # T <- stats[[i]] %*% CM
        # sgn <- sign(T[1, ])
        # P <- t2p(T)
        # sgn <- sign(stats[[i]][1, , drop = FALSE] %*% CM)
        pairComp <- matrix(
            NA, nrow = G, ncol = G,
            dimnames = list(groupsLabs, groupsLabs)
        )
        ## p.values multiplicity adjustment
        
        switch(multAdjMethod,
            "BHS" = {
                p.valuesLow <- BHS(pValues = P.H0Low[[i]][1, ])
                p.valuesGre <- BHS(pValues = P.H0Gre[[i]][1, ])
            },
            "FWEminP"   = {
                p.valuesLow <- t(FWEminP(P.H0Low[[i]]))
                p.valuesGre <- t(FWEminP(P.H0Gre[[i]]))
            },
            {##  other types of multiplicity adjustments
                p.valuesLow <- p.adjust(P.H0Low[[i]][1, ], 
                        method = multAdjMethod)
                p.valuesGre <- p.adjust(P.H0Gre[[i]][1, ], 
                        method = multAdjMethod)
            }
        )#END: switch - multAdjMethod
        
        ### pairwise hypotheses matrix
        pairComp[lower.tri(pairComp)] <- p.valuesLow
        pairComp <- t(pairComp)
        pairComp[lower.tri(pairComp)] <- p.valuesGre
        
        ### rankings
        ranks <- array(NA, dim = c(length(alpha), G),
            dimnames = list(paste("alpha=", alpha, sep = ""), groupsLabs))
        
        for(aa in seq_along(alpha))
        {
            pairHypMat <- ifelse((2 * pairComp) < alpha[aa], 1L, 0L)
            r <- colSums(pairHypMat, na.rm = TRUE) - 
                    rowSums(pairHypMat, na.rm = TRUE)
            ranks[aa, ] <- rank(r, ties.method = "min")
        }# END:for-alpha
        
        ### list of "RankResults" objects
        rankRes[[i]] <- new("RankResults",
            alpha = alpha,
            ranks = ranks,
            # p.values = t(p.values),
            p.values = pairComp,
            p.adj.method = multAdjMethod
        )
    }# END:for-stats
    ### unlisting rankRes if is a one-length list
    if(length(P.H0Low) == 1L)
    {
        Res <- rankRes[[1L]]
    } else
    {
        Res <- rankRes
    }# END:ifelse-rankRes
    
    return(Res)
}#=END=
