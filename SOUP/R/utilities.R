##===========================================================##
##  Constructing Design    Matrix for pairwise comparisons   ##
##===========================================================##
.allPerms <- function(x)
{
    countPerm <-  function(m)
        round(exp(lfactorial(sum(m)) - sum(lfactorial(m))))
    
    xVal <- unique(x)
    xTab <- table(x)
    
    if (length(xVal) == 1)
    {
        out <- xVal
    } else
    {
        n <- sum(xTab)
        out <- matrix(0 , n, countPerm(xTab))
        where <- 0
        for (i in seq_along(xVal))
        {
            if (xTab[i] == 1)
            {
                xTabNew <- xTab[-i]
                xValNew <- xVal[-i]
            } else
            {
                xTabNew <- xTab
                xTabNew[i] <- xTabNew[i] - 1
                xValNew <- xVal
            }
            range <- where + seq_len(countPerm(xTabNew))
            where <- max(range)
            out[1L, range] <- xVal[i]
            out[2L:n, range] <- Recall(rep.int(xValNew, times = xTabNew))
        }# END: for 
    }# END: if-else finished or not
    out
}# END: .allPerms




##===========================================================##
##  Constructing Design    Matrix for pairwise comparisons   ##
##-----------------------------------------------------------##
## @author: Federico Mattiello                               ##
## @date:   17/10/2011                                       ##
## @version: 3.0                                             ##
## @notes: unbalanced case taken into account                ##
## @notes: now 4 times faster thanks to the use of "matrix"  ##
##         instead of "kronecker"                            ##
##-----------------------------------------------------------##
## Inputs:                                                   ##
## - N: number of replications of the experiment, either     ##
##      an integer vector or a "table", both of length "C"   ##
##-----------------------------------------------------------##
## Outputs:                                                  ##
## - M: matrix of {-1, 0, 1}, of dimensions (sum(N) x K)     ##
##      where K = C*(C-1)/2; if dataset has "nObs" rows and  ##
##     "p" columns then "t(dataset) %*% M = P", where P is   ## 
##      a (p x K) matrix of pairwise differences of sums     ##
##===========================================================##
#' Construct Design Matrix that pre-multiplied to the dataset gives the 
#' pairwise mean differences (wrt the groups)
#' 
#' @title Design Matrix For Pairwise Differences
#' @rdname DesM
#' @param N
#'     number of replication of the experiment for each group, either an 
#'     \code{integer} vector or a \code{table} of length \code{G} 
#'     \emph{i.e.} number of groups/treatments 
#' @return 
#'     matrix of 
#' @author Federico Mattiello <federico.mattiello@@gmail.com>
#' 
.DesM <- function(N)
{
    M <- NULL
    C <- length(N)
    sq <- c(0, cumsum(N))
    for (i in seq_len(C - 1))
    {
        tmp <- NULL
        negD <- -diag(C - i)
        for (j in (i + 1):C)
        {
            tmp <- rbind(tmp,
                matrix(negD[j - i, ], nrow = N[j], ncol = C - i, byrow = TRUE)
            )# END:tmp
        }# END:for-j
        A <- rbind(
            array(0, c(sq[i], C - i)), 
            array(1, c(N[i], C - i)),
            tmp
        )# END:A
        M <- cbind(M, A)
    }# END:for-i
    return(M)
}#=END=



###- "Dummyzing" categorical variables -###
#' Transform a \code{factor} or a factor-like vector of \code{character} into
#' a matrix of dichotomous dummy variables
#' 
#' @title From Factor To Dummy Variables
#' @rdname dummyze
#' @param x 
#'     \code{factor}, \code{character} or \code{integer} to be \emph{dummyzed}
#' @return 
#'     a matrix of \{-1, 0, 1\} of dimension \code{length(x)} \eqn{\times}{x}
#'     \code{G} where \code{G} is the number of levels (distinct values) of 
#'     the factor (vector)
#' @author Federico Mattiello <federico.mattiello@@gmail.com>
#' @export
#' 
.dummyze <- function(x)
{
    u <- sort(unique(x))
    d <- matrix(0, nrow = length(x), ncol = length(u), 
            dimnames = list(NULL, u))
    for(k in seq_along(u))
    {
        d[, k] <- as.integer(x %in% u[k])
    }
    return(d)
}#=END=
##- documenting it
# # prompt(.dummyze, file = "man/.dummyze.Rd")

##================================================##
##  Rearrange p.values from vector to a 3D array  ##
##------------------------------------------------##
## Input: "p x 2K" matrix of univariate p.values  ##
## Output: 3D array                               ##
##================================================##
## AGGIORNATA
## DEVE PRENDERE IN INGRESSO LA "P" E ANCHE FARE LA CORREZIONE PER MOLTEPLICITA"
#' Take as input the \eqn{3}-ways array of \emph{p}-values and return the 
#' \eqn{G \times G \times V}{G x G x V} matrix of observed \emph{p}-values 
#' adjusted for multiplicity; \eqn{G, V} are, respectively, 
#' the number of groups/treatments and the number of variables.
#' 
#' @title Construct Univariate \emph{p}-Values Matrix
#' @rdname makePValueMat
#' @param P 
#'     the \eqn{3}-ways \code{array} containing the \emph{p}-values for each 
#'     pairwise comparison in each variable
#' @param multAdjMethod 
#'     the \code{character} string indicating which multiplicity correction 
#'     must be used
#' @param groupsLabs
#'     \code{character} vector containing the groups' labels 
#' @return 
#'     an object of the class \code{\linkS4class{PValueMat}}
#' @author Federico Mattiello <federico.mattiello@@gmail.com>
#' 
.makePValueMat <- function(P, multAdjMethod, groupsLabs)
{
    #>
    # browser()
    #<
    ##- multiplicity correction check
    multAdjMethod <- match.arg(multAdjMethod, 
            c("FWEminP", "BHS",  p.adjust.methods))
    nVar <- dim(P)[2]
    K <- dim(P)[3]
    C <- as.integer(round(.5 + sqrt(.25 + 2 * K)))
    M.raw <- M.adj <- array(NA, 
        dim = c(C, C, nVar), 
        dimnames = list(groupsLabs, groupsLabs, dimnames(P)[[2]]))
    ##- raw p.values
    tmp <- array(NA, dim = c(C, C))
    for(j in seq_len(nVar))
    {
        tmp[lower.tri(tmp)] <- P[1, j, ]
        tmp <- t(tmp)
        tmp[lower.tri(tmp)] <- 1 - P[1, j, ]
        M.raw[, , j] <- tmp
    }# END:for
    ##- p.values multiplicity adjustment
    switch(multAdjMethod,
        # "bonferroni" = {
            # P.adj <- P[1, , ] * K
            # for(j in seq_len(nVar)) {
                # tmp[lower.tri(tmp)] <- P.adj[j, 1:K]
                # tmp <- t(tmp)
                # tmp[lower.tri(tmp)] <- P.adj[j, (K + 1):(2 * K)]
                # M.adj[, , j] <- tmp
                #- only for pretty output
                # multAdjMethod <- "Bonferroni"
            # }# END:for
        # },
        "BHS" = {
            for(j in seq_len(nVar))
            {
                P.adj <- BHS(pValues = P[1, j, ])
                tmp[lower.tri(tmp)] <- P.adj
                tmp <- t(tmp)
                tmp[lower.tri(tmp)] <- 1 - P.adj
                M.adj[, , j] <- tmp
            }# END:for
        },
        "FWEminP" = {
            for(j in seq_len(nVar))
            {
                P.adj <- FWEminP(P[, j, ])
                tmp[lower.tri(tmp)] <- P.adj
                tmp <- t(tmp)
                tmp[lower.tri(tmp)] <- 1 - P.adj
                M.adj[, , j] <- tmp
            }# END:for
        }, 
        {# one of p.adjust.methods
            if(!(multAdjMethod %in% p.adjust.methods))
            {
                multAdjMethod <- "bonferroni"
            }# END:check-p.adjust.methods
            for(j in seq_len(nVar))
            {
                P.adj <- p.adjust(P[1, j, ], method = multAdjMethod)
                tmp[lower.tri(tmp)] <- P.adj
                tmp <- t(tmp)
                tmp[lower.tri(tmp)] <- 1 - P.adj
                M.adj[, , j] <- tmp
            }# END:for
        }
    )# END:switch-multAdjMethod
    ##- outputs
    res <- new("PValueMat", 
        raw.p.values = M.raw,
        adj.p.values = M.adj,
        p.adj.method = multAdjMethod
    )
    return(res)
}#=END=



##====================================##
##  Constructs the permutation space  ##
##  of the row indexes (IDs)          ##
##====================================##
#' Constructs the permutation space of row indices depending on the type of 
#' analysis
#' 
#' @title Constructs Row Indices Permutation Space
#' @rdname makePermSpaceID
#' @param nObs
#'     \code{integer} total number of observations
#' @param analysisType 
#'     \code{character} type of the analysis to be performed
#' @param strata 
#'     \code{character} vector or \code{factor} containing the covarite used 
#'     for stratification; if \code{analysisType} is \code{"strata"} 
#'     then this parameter is provided
#' @param seed 
#'     optional \code{integer} seed for the \code{RNG}
#' @param nPerms 
#'     \code{integer} number of permutations to be performed
#' @return 
#'     either a \code{matrix} in the case of \code{"simple"} 
#'     \code{analysisType} or a \eqn{3}-ways \code{array}, containing the 
#'     permuted row indices 
#' @author Federico Mattiello <federico.mattiello@@gmail.com>
#' 
.makePermSpaceID <- function(nObs, analysisType, strata, seed, nPerms)
{
    if (!missing(seed) && !is.null(seed))
    {
        set.seed(seed)
    } else {}# END: if - set.seed
    
    ### start
    switch(analysisType,
        "simple" = {
            permIndexMat <- apply(
                array(seq_len(nObs), dim = c(nObs, nPerms)), 
                2, FUN = sample, replace = FALSE
            )
        },
        "strata" = {
            if(missing(strata))
            {
                stop("strata can not be missing ", 
                        "with analysis type \"strata.aov\"")
            } else {
                S <- length(tab.strata <- table(strata))
                strataLevs <- unique(strata)
                if(nObs != sum(tab.strata)) {
                    stop("number of observations and length ", 
                            "of \"strata\" variable differs")
                } else {}# END:if
                
                tempInd <- array(NA, c(max(tab.strata), nPerms, S))
                seqObs <- seq_len(nObs)
                for(ss in seq_len(S))
                {
                    tempInd[1:(tab.strata[[ss]]), , ss] <- 
                            seqObs[strata == strataLevs[ss]]
                }# END:for
                permIndexMat <- apply(tempInd, MARGIN = c(2, 3), 
                        FUN = .sampleNA, replace = FALSE) 
            }# END:ifelse-missing-strat
        },
        regres.aov = { }
    )# END:switch-analysisType
    
    return(permIndexMat)
}#=END=



##=================================================##
##  Modified "sample" function: mySample(3L) == 3  ##
##=================================================##
.mySample <- function(x, ...)
{
    x[sample.int(length(x), ...)]
}#=END=


##=======================================##
##  Function for residualizing response  ##
##  variables by numerical covariates    ##
##=======================================##
#' Residualises the response variables \emph{w.r.t.} the matrix of covariates
#' like in the linear models, therefore projecting on the space \eqn{I - H} 
#' where \eqn{H} is the \emph{hat} matrix of \eqn{X}
#' 
#' @rdname orthoX
#' @param Y 
#'     the \code{matrix} or \code{data.frame} of response variables
#' @param X 
#'     the \code{matrix} of covariates
#' @return 
#'     the residualised \code{Y} \code{matrix} 
#' @author Federico Mattiello <federico.mattiello@@gmail.com>
#' @export
#' 
.orthoX <- function(Y, X)
{
    ## matrix (I - H)
    # IP0 <- chol(chol2inv(qr(X)$qr))# = (X"X)^(-1)
    IP0 <- diag(nrow(X)) - X %*% qr.coef(qr(X), diag(nrow(X)))
    
    ## spectral decomposition of (I - H)
    ei <- eigen(IP0)
    ## optional removal of eigenvector associated with "below the tolerance" 
    ## eigenvalues (all eigenvalues are 0 or 1)
    # ei$vectors <- ei$vectors[, (ei$values > 1e-1)] 
    
    Y <- crossprod(ei$vectors, Y)
    # X <- t(ei$vectors) %*% X
    
    return(Y)
}# END:orthoX

##==============================================================##
##  Function that constructs matrix for direct combination      ##
##  over pairwise comparisons of the type:                      ##
##  1>2, 1>3, ..., C-1>C, 2>1, 3>1, ..., C>C-1                  ##
##  T[,j,] %*% combMat(C) = r, where T[,j,] is the matrix of    ##
##  statistics for all "2K" pairwise comparisons and "r" is the ##
##  "B x C" matrix containing sums of i>j with j != i           ##
##==============================================================##
# .pairDiffSumMat <- function(C) {
    # K <- C * (C - 1)/2
    # v <- NULL
    # m <- NULL
    # for (i in seq_len(C - 1)) {
        # v <- c(v, rep(1, C - i), rep(0, K))
        # m <- rbind(m, cbind(array(0, c(C - i, i)), diag(C - i)))
    # }#END:for
    # v <- c(v, rep(0, K))
    # r <- rbind(array(v, c(K, C)), m)
    # return(r)
# }#=END=


##===================================================##
##  sample only available data (the 3D matrix of     ##
##  indexes in the case of a categorical covariate)  ##
##===================================================##
.sampleNA <- function(x, ...)
{
    x[!is.na(x)] <- .mySample(x[!is.na(x)], ...)
    return(x)
}#=END=

##====================================================##
##  function that calculates AOV's residuals sigma^2  ##
##====================================================##
.sigma <- function(x, X.lm)
{
    # x.hat <- X.lm %*% qr.coef(qr(X.lm, LAPACK = TRUE), x)
    x.hat <- X.lm %*% qr.coef(qr(X.lm), x)
    # s2 <- sum((x - x.hat)^2)/(length(x) - ncol(X.lm))
    # sqrt(s2)
    s2 <- crossprod((x - x.hat)^2, rep.int(1, NROW(x))) /
        (NROW(X.lm) - ncol(X.lm))
    return(drop(sqrt(s2)))
}# END:sigma
# .sigma2 <- function(y, X.lm){
    # y.hat <- X.lm %*% solve(t(X.lm)%*%X.lm) %*% t(X.lm) %*% y
    # s2 <- sum((y - y.hat)^2)/(length(y) - ncol(X.lm))
    # sqrt(s2)
# }# END:sigma

##=================================================##
##  From a factor variable to an integer variable  ##
##  (only if this makes sense)                     ##
##=================================================##
.unfactor <- function(fac)
{
    if(is.factor(fac))
    {
        ord <- order(fac)
        res <- rep.int(as.integer(levels(fac)), times = table(fac))
        res[ord] <- res
        return(res)
    } else {
        return(fac)
    }# END:ifelse
}#=END=
##- documenting it
# # # prompt(.unfactor, file = "man/.unfactor.Rd")
