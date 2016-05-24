###---------------------------------------------------------------###
###- MODIFICARE PER FARE I TESTS NELL'ALTRA DIREZIONE: t2p(-T2)! -###
###---------------------------------------------------------------###
#' NonParametric Combination of the test statistic matrix, mostly for 
#' internal use.
#' 
#' It takes as input the 3-ways matrix containing the raw test statistics and 
#' perform the NPC with the possibility to add sets of weights for weighting 
#' variables differently, and/or to select subsets of variables to which NPC 
#' has to be applied separately. Note that \code{weights} and \code{subsets} 
#' are placed in the ``\code{\dots}'' argument of the \pkg{SOUP} function 
#' call hence they are not documented in the \pkg{SOUP} help.
#' 
#' @title NonParametric Combination of Test Statistic
#' @param rawStats 
#'     3-ways \code{array} containing the test statistic computed on 
#'     all pairwise comparisons (both directions) in all variables, 
#'     with dimensions \eqn{B+1 \times p \times K}{B+1 x p x K}, 
#'     where \code{B} is the number of permutations, \code{p} is the number 
#'     of variables and \code{K} is the number of pairwise comparisons.
#' @param combFun 
#'     \code{character} string indicating the combining function to be used, 
#'     can take as input \emph{p}-values or the raw test statistic depending
#'     on the choice, can be one of  \code{"Fisher"}, \code{"Liptak"}, 
#'     \code{"minP"}, \code{"Tippett"}, \code{"maxT"}, \code{"sumT"}
#'     and \code{"direct"}; note that \code{"direct"} is equivalent to 
#'     \code{"sumT"} and \code{"Tippett"} is equivalent to \code{"minP"}.
#'     \code{"Liptak"} use the normal quantile function. See references for 
#'     details.
#' @param seed 
#'     \code{integer} seed for the \code{RNG} (random number generator), 
#'     taken from the input for reproducibility purposes
#' @param p.values 
#'     \code{logical}, if \code{TRUE} means that the input matrix 
#'     \code{rawStats} contains \emph{p}-values rather than raw test  
#'     statistics and so has to be treated differently, \emph{i.e.} the first  
#'     passage from test statistics to \emph{p}-values is omitted (function 
#'     \code{\link{t2p}})
#' @param tails
#'     \code{integer} vector of \eqn{\pm 1}{{+1,-1}} containing the 
#'     alternatives for response variables: \code{+1} means ``the higher the 
#'     better'', \code{-1} means ``the lower the better'' (direction of 
#'     preference), if \code{NULL} (default) all variables are considered 
#'     to be of the type ``the higher the better'
#' @param subsets 
#'     \code{integer}, \code{character} or \code{list} where each element 
#'     contains the subset of column indices that need to be treated 
#'     separately, if \code{NULL} (default) all input variables are considered
#' @param weights 
#'     \code{integer}, \code{character} or \code{list} where each element 
#'     contains the weights of the variables, one for each \code{subset}, 
#'     if \code{NULL} (default) variables are treated equally \emph{i.e.} 
#'     all have the same weight
#' @param iteratedNPC 
#'     \code{logical}, if \code{TRUE} it performs the iterated combination 
#'     procedure running function \code{\link{iterNPC}} and the choice of 
#'     \code{comb.funct} becomes thus irrelevant; otherwise just perform the 
#'     requested NPC.
#' @return 
#'     an object of class \code{\linkS4class{PermSpace}} or a \code{list} 
#'     containing \code{PermSpace} objects, in the case of multiple weights 
#'     and/or subsets.
#' @author Federico Mattiello <federico.mattiello@@gmail.com>
#' @references 
#'     Pesarin, F. and Salmaso, L. (2010) \emph{Permutation Tests for 
#'     Complex Data}. Wiley: United Kingdom \cr
#'     \cr
#'     Pesarin F. (2001) \emph{Multivariate Permutation Tests with Applications 
#'     in Biostatistics}. Wiley: New York.
#' @note 
#'     This function is mainly taken from the function 
#'     \code{\link[flip:permutationSpace]{npc}} in the package 
#'     \code{\link[flip]{flip}}.
#' @seealso 
#'     \code{\link{t2p}}
#' @keywords 
#'     array, manip
#' @export
#' 
NPC <- function(rawStats, combFun = "Fisher", seed, p.values = FALSE,
    tails = NULL, subsets = NULL, weights = NULL, 
    iteratedNPC = FALSE)
{
    ##- combining function matching
    combFun <- match.arg(
        tolower(combFun[1]),
        c("fisher", "liptak", "minp", "tippett", "maxt", "sumt", "direct")
    )
    if(combFun == "tippett")
    {
        combFun <- "minp"
    }
    if(combFun == "direct")
    {
        combFun <- "sumt"
    }
    ##- check dimensions of rawStats
    if(length(dim(rawStats)) == 2) {
        rawStats <- array(rawStats,
            dim      = c(nrow(rawStats), 1, ncol(rawStats)),
            dimnames = list(rownames(rawStats), NULL, colnames(rawStats))
        )
    }# END:checkDim-rawStats

    ##- number of groups
    K <- dim(rawStats)[3]
    C <- as.integer(round(.5 + sqrt(.25 + 2 * K)))
    ##-    number and names of variables
    nVar <- dim(rawStats)[2]
    if(is.null(dimnames(rawStats)[[2]])) {
        varNames <- as.character(seq_len(nVar))
    } else {
        varNames <- dimnames(rawStats)[[2]]
    }# END:if-null-varNames
    ##- label of Pairwise Comparisons
    labsPC <- dimnames(rawStats)[[3]]
        
    ##--------------------------##
    ##  SUBSETS & WEIGHTS       ##
    ##--------------------------##
    {# < START: SUBSETS & WEIGHTS >
        ##- subsets    
        if(missing(subsets) || is.null(subsets)) {
            one.subset <- TRUE
            subsets <- list(seq_len(nVar))
        } else {
            one.subset   <- !is.list(subsets)
            many.subsets <- !one.subset
            subsets <- as.list(subsets)
        }# END:ifelse-subsets
        ##- weights
        if(missing(weights) || is.null(weights)) {
            ## one subset, with many subsets weights must be "!missing"
            one.weight <- TRUE
            weights <- list(rep.int(1L, length(subsets[[1]])))
        } else {
            one.weight   <- !is.list(weights)
            many.weights <- !one.weight
            weights <- as.list(weights)
        }# END:ifelse-weights
        ##- subsets and weights cases
        if(one.subset) {
            ## one subsets, one OR many weights
            if(one.weight) {case <- "1"} else {case <- "2"}
        } else {
            ## many subsets, many weights
            if(many.weights) {case <- "3"}
        }# END:cases-subsets&weights
        ##- check matching lengths of subsets and weights
        wts.len <- sapply(weights, FUN = length)
        sub.len <- sapply(subsets, FUN = length)
        if(case %in% c("1", "2")) {
            if(any(wts.len != nVar)) {
                stop("lengths of \"weights\" and number of variables do not match")
            } else { }# END:if
        } else {# case == "3"
            if(length(subsets) != length(weights)) {
                stop("number of \"subsets\" and \"weights\" do not match")
            } else {
                if(any(wts.len != sub.len)) {
                    stop("lengths of \"subsets\" and \"weights\" do not match")
                } else { }# END:if
            }# END:ifelse-lengths
            if(is.character(subsets[[1L]])) {
                sub.ind <- lapply(subsets,
                    FUN = function(x, names) {seq_along(names)[names == x]}, names = varNames
                )
            } else {
                sub.ind <- subsets
            }# END:matching-subsets-names
        }# END:ifelse-case
    }# < END: SUBSETS & WEIGHTS >
    ##--------------------------##
    ##  END: SUBSETS & WEIGHTS  ##
    ##--------------------------##
    ##- setting tails
    if(p.values) {
        rawStats[, tails == -1, ] <- 1 - rawStats[, tails == -1, ]
    } else {
        rawStats[, tails == -1, ] <- -rawStats[, tails == -1, ]
    }# END:asymptotic-p.values
    
    ##================================================##
    ##== NUOVA rawStats PER NUOVO METODO DI RANKING ==##
    ##================================================##
    ##- Tpd = T of pairwise differences,
    # Tpd1 <- rawStats[, , 1:K]
    
    ##- Trp = T of ranking parameters, resulting from the
    ##  sum of comparisons 'i > h' with 'h != i'
    ##  if "p.values == TRUE" Fisher's comb. fun. else "direct" comb. fun.
    # if(p.values) {
        # Trp1 <- tensor(-log(rawStats), combM, 3, 1)
    # } else {
        # Trp1 <- tensor(rawStats, combM, 3, 1)
    # }# END:if-asymptotic-p.values
    
    
    ##- p.value transformation
    ## distinguish between iterated NPC and single NPC
    if (iteratedNPC)
    {
        if (p.values)
        {
            T1.H0Low <-  rawStats
            T1.H0Gre <- 1 - rawStats
        } else
        {
            T1.H0Low <- t2p( rawStats)
            T1.H0Gre <- t2p(-rawStats)
        }# END - if else: p.values or raw statistics
    } else
    {
        switch(combFun,
            'fisher' = {
                if(p.values) {
                    T1.H0Low <- -log(rawStats)
                    T1.H0Gre <- -log(1 - rawStats)
                } else {
                    T1.H0Low <- -log(t2p( rawStats))
                    T1.H0Gre <- -log(t2p(-rawStats))
                }# END:if-asymptotic-p.values
            },
            'liptak' = {
                if(p.values) {
                    T1.H0Low <- qnorm(1 - rawStats)
                    T1.H0Gre <- qnorm(rawStats)
                } else {
                    T1.H0Low <- qnorm(1 - t2p(rawStats))
                    T1.H0Gre <- qnorm(1 - t2p(-rawStats))
                }# END:if-asymptotic-p.values
            },
            'sumt'   = {
                if(p.values) {
                    T1.H0Low <- -log(rawStats)
                    T1.H0Gre <- -log(1 - rawStats)
                } else {
                    T1.H0Low <-  rawStats
                    T1.H0Gre <- -rawStats
                }# END:if-asymptotic-p.values
            },
            'maxt'   = {
                if(p.values) {
                    T1.H0Low <- 1 / rawStats
                    T1.H0Gre <- 1 / (1 - rawStats)
                } else {
                    T1.H0Low <-  rawStats
                    T1.H0Gre <- -rawStats
                }# END:if-asymptotic-p.values
            },
            'minp'   = {
                if(p.values) {
                    T1.H0Low <- 1 / rawStats
                    T1.H0Gre <- 1 / (1 - rawStats)
                } else {
                    T1.H0Low <- 1 / t2p( rawStats)
                    T1.H0Gre <- 1 / t2p(-rawStats)
                }# END:if-asymptotic-p.values
            }
        )# END:switch-combFun
    }# END - if else: iteratedNPC
    
    ##- permutation space of ranking statistics, i.e. weighting, subsetting
    ##  and combining of pairwise comparisons
    switch(case, 
        '1' = {
            wts <- aperm(array(weights[[1L]],
                dim = c(nVar, dim(T1.H0Low)[-2])), c(2, 1, 3)
            )
            T2.H0Low <- T1.H0Low * wts
            T2.H0Gre <- T1.H0Gre * wts
            
            #############################################################
            #############################################################
            #############################################################
            #############################################################
            
            ### distinguish between iterated or single NPC 
            if (iteratedNPC)
            {
                P.H0Low <- apply(T2.H0Low, MARGIN = 3, FUN = iterNPC,
                        plotIt = FALSE, onlyCombined = TRUE)
                P.H0Gre <- apply(T2.H0Gre, MARGIN = 3, FUN = iterNPC, 
                        plotIt = FALSE, onlyCombined = TRUE)
            } else
            {
                if(combFun %in% c("fisher", "liptak", "sumt")) {
    #                T3.H0Low <- tensor(T2.H0Low, rep.int(1, nVar), 2, 1)
    #                T3.H0Gre <- tensor(T2.H0Gre, rep.int(1, nVar), 2, 1)
                    T3.H0Low <- colSums(aperm(T2.H0Low, c(2L, 1L, 3L)))
                    T3.H0Gre <- colSums(aperm(T2.H0Gre, c(2L, 1L, 3L)))
                } else {
                    T3.H0Low <- apply(T2.H0Low, MARGIN = c(1, 3), FUN = max)
                    T3.H0Gre <- apply(T2.H0Gre, MARGIN = c(1, 3), FUN = max)
                }# END:ifelse-sumCombFun
                
                ##- permutation space of p.values
                P.H0Low <- t2p(T3.H0Low)
                P.H0Gre <- t2p(T3.H0Gre)
                colnames(T3.H0Low) <- labsPC
                colnames(T3.H0Gre) <- labsPC
            }# END: if else iteratedNPC
            
            colnames(P.H0Low) <- labsPC
            colnames(P.H0Gre) <- labsPC
        },
        '2' = {
            T3.H0Low <- T3.H0Gre <- vector("list", length(weights))
            P.H0Low <- P.H0Gre <- vector("list", length(weights))
            
            for(i in seq_along(weights)) {
                wts <- aperm(array(weights[[i]],
                    dim = c(nVar, dim(T1.H0Low)[-2])), c(2, 1, 3)
                )
                nm.wt <- paste(paste(
                        paste(varNames, "wt", sep = ":"),
                        weights[[i]], sep = "="), collapse = ", "
                )# END:names
                
                T2.H0Low <- T1.H0Low * wts
                T2.H0Gre <- T1.H0Gre * wts
                
                ### distinguish between iterated or single NPC 
                if (iteratedNPC)
                {
                    P.H0Low[[i]] <- apply(T2.H0Low, MARGIN = 3, FUN = iterNPC,
                            plotIt = FALSE, onlyCombined = TRUE)
                    P.H0Gre[[i]] <- apply(T2.H0Gre, MARGIN = 3, FUN = iterNPC, 
                            plotIt = FALSE, onlyCombined = TRUE)
                    
                    colnames(P.H0Low[[i]]) <- colnames(P.H0Gre[[i]]) <- labsPC
                    names(P.H0Low)[i] <- names(P.H0Gre)[i] <- nm.wt
                } else
                {
                    if(combFun %in% c("fisher", "liptak", "sumt")) {
    #                    T3.H0Low[[i]] <- tensor(T2.H0Low, rep.int(1, nVar), 2, 1)
    #                    T3.H0Gre[[i]] <- tensor(T2.H0Gre, rep.int(1, nVar), 2, 1)
                        T3.H0Low[[i]] <- colSums(aperm(T2.H0Low, c(2L, 1L, 3L)))
                        T3.H0Gre[[i]] <- colSums(aperm(T2.H0Gre, c(2L, 1L, 3L)))
                    } else {
                        T3.H0Low[[i]] <- apply(T2.H0Low, MARGIN = c(1, 3), FUN = max)
                        T3.H0Gre[[i]] <- apply(T2.H0Gre, MARGIN = c(1, 3), FUN = max)
                    }# END:ifelse-sumCombFun
                    colnames(T3.H0Low[[i]]) <- colnames(T3.H0Gre[[i]]) <- labsPC
                    names(T3.H0Low)[i] <- names(T3.H0Gre)[i] <- nm.wt
                }# END: if else iteratedNPC
            }# END:for - weights
            
            ##- permutation space of p.values
            ### distinguish between iterated or single NPC 
            if (!iteratedNPC)
            {
                P.H0Low <- lapply(T3.H0Low, FUN = t2p)
                P.H0Gre <- lapply(T3.H0Gre, FUN = t2p)
                names(P.H0Low) <- names(P.H0Gre) <- names(T3.H0Low)
            } else {}
            
        },
        '3' = {
            T3.H0Low <- T3.H0Gre <- vector("list", length(subsets))
            P.H0Low <- P.H0Gre <- vector("list", length(subsets))
            
            for(i in seq_along(subsets)) {
                wts <- aperm(array(weights[[i]],
                    dim = c(sub.len[i], dim(T1.H0Low)[-2])), c(2, 1, 3)
                )
                nm.wt <- paste(paste(paste(varNames[sub.ind[[i]]], "wt", sep = ":"),
                        weights[[i]], sep = "="), collapse = ", "
                )# END:names
                
                T2.H0Low <- T1.H0Low[, sub.ind[[i]], , drop = FALSE] * wts
                T2.H0Gre <- T1.H0Gre[, sub.ind[[i]], , drop = FALSE] * wts
                
                ### distinguish between iterated or single NPC 
                if (iteratedNPC)
                {
                    P.H0Low[[i]] <- apply(T2.H0Low, MARGIN = 3, FUN = iterNPC,
                            plotIt = FALSE, onlyCombined = TRUE)
                    P.H0Gre[[i]] <- apply(T2.H0Gre, MARGIN = 3, FUN = iterNPC, 
                            plotIt = FALSE, onlyCombined = TRUE)
                    
                    colnames(P.H0Low[[i]]) <- colnames(P.H0Gre[[i]]) <- labsPC
                    names(P.H0Low)[i] <- names(P.H0Gre)[i] <- nm.wt
                } else
                {
                    if(combFun %in% c("fisher", "liptak", "sumt")) {
    #                    T3.H0Low[[i]] <- tensor(T2.H0Low, rep.int(1, sub.len[i]), 2, 1)
    #                    T3.H0Gre[[i]] <- tensor(T2.H0Gre, rep.int(1, sub.len[i]), 2, 1)
                        T3.H0Low[[i]] <- colSums(aperm(T2.H0Low, c(2L, 1L, 3L)))
                        T3.H0Gre[[i]] <- colSums(aperm(T2.H0Gre, c(2L, 1L, 3L)))
                    } else {
                        T3.H0Low[[i]] <- apply(T2.H0Low, MARGIN = c(1, 3), FUN = max)
                        T3.H0Gre[[i]] <- apply(T2.H0Gre, MARGIN = c(1, 3), FUN = max)
                    }# END:ifelse-sumCombFun
                    colnames(T3.H0Low[[i]]) <- colnames(T3.H0Gre[[i]]) <- labsPC
                    names(T3.H0Low)[i] <- names(T3.H0Gre)[i] <- nm.wt
                }# END: if else iteratedNPC
            }# END:for
            
            ##- permutation space of p.values
            ### distinguish between iterated or single NPC 
            if (!iteratedNPC)
            {
                P.H0Low <- lapply(T3.H0Low, FUN = t2p)
                P.H0Gre <- lapply(T3.H0Gre, FUN = t2p)
                names(P.H0Low) <- names(P.H0Gre) <- names(T3.H0Low)
            } else {}
        }
    )# END:switch
    #>
    # browser()
    #<
    
    ### distinguish between iterated or single NPC 
    if (iteratedNPC)
    {
        T3.H0Low <- array(0, dim = c(0, 0, 0))
        T3.H0Gre <- array(0, dim = c(0, 0, 0))
        
    } else {}
    
    
    ##- outputs: object "PermSpace"
    permSpace <- new(
        Class      = "PermSpace",
        seed       = seed,
        T.H0Low    = T3.H0Low,
        T.H0Gre    = T3.H0Gre,
        P.H0Low    = P.H0Low,
        P.H0Gre    = P.H0Gre,
        rawStats   = rawStats,
        comb.funct = combFun
    )
    return(permSpace)
}#=END=

##- documenting it
# prompt(NPC, file = "C:/Users/Public/Documents/Docs_VMware/myRlibrary/soup2/man/NPC.Rd")
