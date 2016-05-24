###
### TODO: try rowMeans instead of tensor 
###
#' 
#' Main function of the package, interface for every analysis. 
#' The dataset can be balanced or not for almost all possible choices of the 
#' input parameters. The function allows also for the presence of one or more 
#' continuous covariates or for stratified analysis.
#' 
#' Depending on the chosen p-values type and on the analysis type, only some 
#' options can be selected:
#' \itemize{
#' \item{}{
#'     with \code{"simple"} or \code{"regres"} analysis and 
#'     \code{"asymptotic"} \emph{p}-values, \code{"Hotelling"} and 
#'     \code{"Ttest"}; with \code{permutation} \emph{p}-values \code{"AD"}, 
#'     \code{"Hotelling"} and \code{"meanDiff"} can be selected.}
#' \item{}{
#'     With \code{"strata"} analysis and \code{"asymptotic"} \emph{p}-values, 
#'     \code{"lmCoef"} and \code{"Ttest"}; with \code{"permutation"} 
#'     \emph{p}-values \code{"AD"} and \code{"meanDiff"} can be selected.}
#' }
#' @title SOUP Main Function
#' @param Y 
#'     input \code{matrix} where each column is a response variables.
#' @param covars 
#'     it can be a \code{matrix}, a \code{data.frame} or a \code{formula}, 
#'     in the first two cases it must contains at least the labels of groups, 
#'     in the latter case it has to be a right-sided \code{formula} 
#'     (\emph{e.g.} \code{~ v1 + v2}) specifying the model to extract from 
#'     the \code{data} input.
#' @param data 
#'     optional \code{data.frame} containing covariates requested by \code{covars}, 
#'     if \code{covars} is not a formula this input is useless.
#' @param analysisType
#'     \code{character}, type of the analysis to be performed: it can be 
#'     \code{"simple"} if the only covariate is the labels of groups, 
#'     \code{"strata"} if there is also a stratifying (categorical) covariate, 
#'     \code{"regres"} if there is one or more (numerical or not) covariate(s) 
#'     besides labels of groups. In the latter case the linear effect of the 
#'     covariates is removed from the response variables are 
#'     residualised by the matrix \eqn{V^{-1/2}} obtained from 
#'     \eqn{V = I - H} (where \eqn{I} is the identity matrix and \eqn{H} is 
#'     the ``hat'' matrix of the OLS, by means of a spectral decomposition. 
#' @param p.adj.method 
#'     \code{character} string containing the type of required \emph{p}-value 
#'     adjustment
#' @param p.valuesType 
#'     \code{character} string indicating the type of \emph{p}-value to be 
#'     used, it can be \code{"permutation"} or \code{"asymptotic"} 
#' @param testStatistic 
#'     \code{character} string indicating the test statistic to be used, it 
#'     depends on both \code{analysisType} and on \code{p.valuesType} and 
#'     the alternatives are:
#'     \describe{
#'         \item{\code{AD, meanDiff}}{
#'             for all \code{analysisType} but only using \code{permutation} 
#'             \emph{p}-values}
#'         \item{\code{Ttest}}{
#'             for all \code{analysisType} but only using \code{asymptotic} 
#'             \emph{p}-values}
#'         \item{\code{Hotelling}}{
#'             with both \code{permutation} and \code{asymptotic} 
#'             emph{p}-values, with \code{"simple"} and \code{"regres"} but 
#'             not with \code{"strata"} \code{analysisType}}
#'         \item{\code{lmCoef}}{
#'             only with \code{"strata"} \code{analysisType} and with 
#'             \code{"asymptotic"} \emph{p}-values}
#'     }
#' @param combFunct 
#'     \code{character} string containing the desired combining function to be 
#'     used, choices are: 
#'     \describe{
#'     \item{\code{Fisher}}{
#'         the famous Fisher's \emph{p}-values combining function}
#'     \item{\code{Liptak}}{
#'         it uses the quantile function of the Normal distribution to combine 
#'         \emph{p}-values}
#'     \item{\code{minP, tippett}}{
#'         combine \emph{p}-values by taking the minimum across the set}
#'     \item{\code{maxT}}{
#'         combines directly the test statistics by taking the maximum across 
#'         the set}
#'     \item{\code{direct, sumT}}{
#'         combine the test statistics by summing them}
#'     \item{\code{sumT2}}{
#'         combines the test statistics by squaring and summing them}
#'     }
#'     See the references for more details about their properties.
#' @param univ.p.values 
#'     \code{logical}, if \code{TRUE} (default) \emph{p}-values are returned 
#'     for each variable separately in a 3-ways \code{array}, the chosen 
#'     multiplicity correction is performed independently for each variable 
#' @param tails 
#'     \code{integer} vector of \eqn{\pm 1}{\{+1,-1\}} containing the 
#'     alternatives for response variables: \code{+1} means ``the higher the 
#'     better'', \code{-1} means ``the lower the better'' (direction of 
#'     preference), if \code{NULL} (default) all variables are considered 
#'     to be of the type ``the higher the better''
#' @param linearInter 
#'     \code{logical}, if \code{TRUE} the presence of linear interaction is 
#'     assumed between levels of the stratifying covariate and response 
#'     variables, this affects only the \code{"lmCoef"} test statistic in the 
#'     (in the \code{"strata"} \code{analysisType}), 
#'     basically the contrasts matrix of groups is multiplied by the levels 
#'     of the stratifying factor. 
#' @param returnPermSpace 
#'     \code{logical} if \code{TRUE} (default) the whole permutation space is 
#'     returned, class \code{\linkS4class{PermSpace}}, otherwise it is an empty 
#'     instance of the class.
#' @param nPerms 
#'     \code{integer} number of permutation to be performed
#' @param alpha 
#'     \code{numeric} desired significance level, \emph{i.e.} type-I error
#' @param seed 
#'     \code{integer} seed for the Random Number Generator
#' @param iteratedNPC
#'     \code{logical}, single or iterated Non-Parametric Combination, see \
#'     code{\link{iterNPC}} for details.
#' @param ... 
#'     put here the optional \code{weights} and \code{subsets} for the 
#'     \code{\link{NPC}} function and the permutation space of rows indexes 
#'     \code{permSpaceID}. 
#'     The latter allows to exactly reproduce a previous analysis, if all 
#'     other inputs are kept equal, or to see what happens changing for 
#'     example only the \code{testStatistic}.
#' @return 
#'    an object of class \code{\linkS4class{SoupObject}}.
#' @author Federico Mattiello <federico.mattiello@@gmail.com>
#' @references 
#'     Pesarin, F. and Salmaso, L. (2010) 
#'     \emph{Permutation Tests for Complex Data}.
#'     Wiley: United Kingdom.\cr
#'     
#'     Pesarin F. (2001) 
#'     \emph{Multivariate Permutation Tests with Applications in Biostatistics} 
#'     Wiley: New York.\cr
#'     
#'     Federico Mattiello (2010) 
#'     \emph{Some resampling-based procedures for ranking of multivariate 
#'     populations}, Master's Thesis, Faculty of Statistical Sciences: Padova.
#' 
#' @export
#' @examples
#' ### 
#' ### testing SOUP
#' ### 
#' rm(list = ls()); gc(reset = TRUE)
#' 
#' require(SOUP)
#' n <- 5L         # replication of the experiment
#' G <- 4L         # number of groups
#' nVar <- 10L     # number of variables
#' shift <- 1.5    # shift to be added to group 3
#' alpha <- c(0.01, 0.05, 0.1)        # significance levels
#' 
#' ## groups factor
#' groups <- gl(G, n, labels = paste("gr", seq_len(n), sep = "_"))
#' 
#' set.seed(12345)
#' Y <- matrix(rnorm(n * G * nVar), nrow = n * G, ncol = nVar)
#' colnames(Y) <- paste("var", seq_len(nVar), sep = "_")
#' ind1 <- groups == unique(groups)[3L]
#' Y[ind1, ] <- Y[ind1, ] + shift
#' 
#' res <- SOUP(Y = Y, covars = as.matrix(groups), analysisType = "simple", 
#'         testStatistic = "meanDiff", combFunct = "Fisher",
#'         alpha = alpha, 
#'         subsets = list("first" = 1:5, "second" = 6:10), 
#'         weights = list(
#'                 "firstW" = c(.1, .2, .1, .5, .1), 
#'                 "secondW" = rep.int(1, 5)),
#'         p.valuesType = "permutation", p.adj.method = "FWEminP")
#' res
#' 
SOUP <- function(Y, covars, data = NULL, analysisType, p.adj.method, p.valuesType,
    testStatistic, combFunct, univ.p.values = TRUE, tails = NULL, 
    linearInter = FALSE, returnPermSpace = TRUE,  
    nPerms = 999L, alpha = 0.05, seed, iteratedNPC, ...)
{

    ### 
    ### Matching arguments
    ### 
    
    ### iteratedNPC, if missing is set to FALSE
    if (missing(iteratedNPC))
    {
        iteratedNPC <- FALSE
    } else {}
    
    ### analysisType
    if(missing(analysisType))
    {
        stop("\"analysisType\" is missing, must be one of ", 
                "\'simple\', \'strata\' or \'regres\'"
        )# END:stop
    } else {
        analysisType <- match.arg(analysisType, c("simple", "strata", "regres"))
    }# END:if-analysisType
    
    
    ### covars, i.e. covariate(s)
    if(missing(covars) || is.null(covars))
    {
        stop("\"covars\" is empty. It must contains at least the ",
                "vector of groups\' labels")
    } else {}# END:covars-present
    
    
    ### p.adj.method
    if(missing(p.adj.method))
    {
        ## default = bonferroni
        warning ("p.value adjustment method is missing, ",
            "using \"BHS\" as default", call. = FALSE)
        p.adj.method <- "BHS"
    } else {
        p.adj.method <- match.arg(p.adj.method, c("BHS", "FWEminP", p.adjust.methods))
    }# END:if-p.adj.method
    
    
    ### p.value.type
    if(missing(p.valuesType))
    {
        stop("\"p.valuesType\" is missing, must be either ", 
                "\"asymptotic\" or \"permutation\"")
    } else {
        p.valuesType <- match.arg(p.valuesType, c("asymptotic", "permutation"))
    }# END: missing - p.valuesType
    
    
    ### selecting the proper "testStatistic": depends on "analysisType" 
    ### and "p.valuesType"
    if(missing(testStatistic))
    {
        stop(
            "\"testStatistic\" must be one of: ",
            "\"AD\", \"lmCoef\", \"Hotelling\", \"meanDiff\" or \"Ttest\""
        )
    } else
    {
        if(p.valuesType == "asymptotic")
        {
            switch(analysisType,
                "simple" = {
                    testStatistic <- match.arg(testStatistic, 
                            c("Hotelling", "Ttest"))
                },
                "strata" = {
                    testStatistic <- match.arg(testStatistic, 
                            c("lmCoef", "Ttest"))
                },
                "regres" = {
                    testStatistic <- match.arg(testStatistic, 
                            c("Hotelling", "Ttest"))
                }
            )# END: switch - analysisType
        } else ##  permutation p.values
        {
            switch(analysisType,
                "simple" = {
                    testStatistic <- match.arg(testStatistic, 
                            c("AD", "Hotelling", "meanDiff"))
                },
                "strata" = {
                    testStatistic <- match.arg(testStatistic, 
                            c("AD", "meanDiff"))
                },
                "regres" = {
                    testStatistic <- match.arg(testStatistic, 
                            c("AD", "Hotelling", "meanDiff"))
                }
            )# END:switch-analysisType
        }# END:if-p.valuesType
    }# END:missing-testStatistic
    
    
    ### checks for Hotelling's test-statistic
    if(testStatistic == "Hotelling")
    {
        if(univ.p.values)
        {
            warning ("When using \"Hotelling\" test-statistic ",
                    "univariate p-values are not calculated.",
                call. = FALSE)# END:warning
            univ.p.values <- FALSE
        } else {}# END:if-univ.p.values
        
        if(analysisType != "simple")
        {
            warning ("When using \"Hotelling\" test-statistic no ",
                    "covariates are considered, so only \"simple\" ",
                    "analysis can be performed.", call. = FALSE)
            analysisType <- "simple"
        } else { }# END: if - simple
    }# END: check - hotelling
    
    
    ### check combining function
    if(missing(combFunct))
    {
        if (testStatistic != "Hotelling")
        {
            message ("combining function is missing, ",
                    "using \"Fisher\" as default")
            combFunct <- "fisher"
        } else {}
    } else {
        combFunct <- match.arg(tolower(combFunct),
            c("fisher", "liptak", "minp", "tippett", 
                    "maxt", "sumt", "direct", "sumt2"))
    }# END: missing - combFunct
    
    
    ### check of tails
    if(missing(tails) || is.null(tails))
    {
        tails <- rep.int(1L, NCOL(Y))
    } else
    {
        if(length(tails) != NCOL(Y))
        {
            warning ("Number of \"tails\" differs from number of ",
                    "variables. \"tails\" are all set to 1.", 
                    call. = FALSE)
            tails <- rep.int(1L, NCOL(Y))
        } else {}
    }# END: if - check tails
    
    
    ### extract arguments in '...'
    dots <- list(...)
    
    
    ### check dim of Y
    if(NCOL(Y) == 1L)
    {
        dim(Y) <- c(NROW(Y), 1)
    } else {}# END:dim-Y
    
    
    
    ### 'covars' is a formula of variables in 'data'
    if(is(covars, "formula"))
    {
        if(length(covars) > 2L)
        {
            warning("\"covars\" is a \"formula\" with response, ",
                    "it has to be a right-sided formula ",
                    "(responses are in \"Y\").")
                covars <- covars[-2]
        } else {}# END: if - formula check
        
        if(!is.null(data))
        {
            covars <- model.frame(covars, data = data)
        } else {
            stop("\"data\" argument is missing and \"covars\" is a formula, ", 
                    "not a matrix.")
        }# END: ifelse - data check
    } else {}# END: ifelse - covars is a formula
    
    
    ### "covars" is a matrix
    if(!is.data.frame(covars))
    {
        covars <- data.frame(covars)
    } else {}# END:if
    
    
    
    ### set the seed for the random number generator
    if(missing(seed))
    {
        seed <- round(1e3 * runif(1), digits = 0)
    }# END:ifelse-set.seed
    set.seed(seed)
    
    ### 
    ### END Matching
    ### 
    
    
    ### regression aov: residualise Y with a linear regression of Y on the 
    ### covariate(s); after that is like the "simple"
    if(analysisType == "regres")
    {
        Xmm <- model.matrix(model.frame(covars[, -1L, drop = FALSE]), data = data)
        Y <- .orthoX(Y, Xmm)
        analysisType <- "simple"
        covars <- covars[, 1L, drop = FALSE]
    } else {}# END: if - regres
    
    
    
    ### do the analysis
    switch(analysisType,
            
        ### simple aov: AD = Anderson-Darling, Hotelling = Hotelling's T^2,
        ##  meanDiff = differences of means
        "simple" = {
            ### check number of columns of covars
            if(NCOL(covars) > 1L)
            {
                stop("selected \"simple\" so \"covars\" must have only one column")
            }# END:if-no-covariate(s)
            
            ### ordering of data
            groups <- covars[, 1]
            
            ord <- order(groups)
            Y <- Y[ord, ]
            groups <- groups[ord]
            
            ### permutation space of row IDs
            if(is.null(dots$permSpaceID))
            {
                permIndexMat <- .makePermSpaceID(nObs = NROW(Y),
                    analysisType = analysisType, seed = seed, nPerms = nPerms)
            } else
            {   ## take input permSpace of indexes (for reproducibility)
                permIndexMat <- dots$permSpaceID
                seed <- integer()
            }# END: indexes - permSpace
            
            switch(testStatistic,
                "AD"        = {
                    T <- .simpleAD (
                        dataset = Y, groups = groups, indexMat = permIndexMat
                    )# END:AD
                },
                "Hotelling" = {
                    T <- .simpleHotelling (
                        dataset = Y, groups = groups, indexMat = permIndexMat,
                        p.valuesType = p.valuesType
                    )# END:hotelling
                },
                "meanDiff"  = {
                    T <- .simpleMeanDiff (
                        dataset = Y, groups = groups, indexMat = permIndexMat
                    )# END:meanDiff
                },
                "Ttest"  = {
                    T <- .simpleTtest(
                        dataset = Y, groups = groups, indexMat = permIndexMat
                    )# END:Ttest
                },
                {
                    stop("can not match the test statistics")
                }
            )# END:switch-simple
        },
        
        ### stratified aov: "AD" = Anderson-Darling, "lmCoef" = anova test 
        ### for the effect of the "groups" and pairwise differences between 
        ### coefficients estimated by "lm", variable by variable, 
        ### "meanDiff" = differences of mean, "Ttest" = t-test as in
        ### the case "simple"-"meanDiff"-asymptotic p.values but with 
        ### different denominator (the MSE take into account also the 
        ### stratifying variable considered as factor)
        "strata" = {
            ### check number of columns of covars
            if(NCOL(covars) > 2)
            {
                stop("selected \"strata\" so \"covars\" must have 2 columns, ",
                        "both factor variables")
            }# END:if-no-covariate(s)
            
            ### ordering of data
            # if((NCOL(covars) == 2) && is.factor(covars[, 1]) && is.factor(covars[, 2])) {
            if(NCOL(covars) == 2L)
            {
                ord    <- order(covars[, 2], covars[, 1])
                covars <- covars[ord, ]
                Y      <- Y[ord, ]
                groups <- covars[, 1]
                strata <- .unfactor(covars[, 2])
            }# END:if-1-covariate
            ### permutation space of rows IDs
            if(is.null(dots$permSpaceID)) {
                permIndexMat <- .makePermSpaceID(nObs = NROW(Y),
                    analysisType = analysisType, strata = strata, 
                    seed = seed, nPerms = nPerms
                )
            } else {# take input permSpace of indexes (for reproducibility)
                permIndexMat <- dots$permSpaceID
                seed <- integer()
            }# END:indexes-permSpace
            ### select test statistic
            switch(testStatistic,
                "AD"       = {
                    T <- .strataAD(
                        dataset = Y, groups = groups, strata = strata, 
                        indexMat = permIndexMat
                    )# END:AD
                },
                "lmCoef"    = {
                    T <- .strataLmCoef(
                        dataset = Y, groups = groups, strata = strata, 
                        indexMat = permIndexMat, linearInter = linearInter
                    )# END:lmCoef
                },
                "meanDiff" = {
                    T <- .strataMeanDiff(
                        dataset = Y, groups = groups, strata = strata, 
                        indexMat = permIndexMat, linearInter = linearInter
                    )# END:meanDiff
                },
                "Ttest"    = {
                    T <- .strataTtest(
                        dataset = Y, groups = groups, strata = strata, 
                        indexMat = permIndexMat, linearInter = linearInter
                    )# END:Ttest
                },
                {
                    stop("can not match the test statistics")
                }
            )# END:switch-testStatistic-strata
        }
    )# END:switch-analysisType
    ###-----------------------------##
    
    
    ### from raw statistics to univariate p.values
    if(p.valuesType == "asymptotic")
    {
        P <- T
    } else
    {
        P <- t2p(T)
    }# END:ifelse-asympotic-p.values

    ### create object PValueMat
    if(univ.p.values)
    {
        pValueMat <- .makePValueMat(
            P = P, multAdjMethod = p.adj.method, groupsLabs = unique(groups)
        )
    } else
    {
        pValueMat <- new("PValueMat")
    }# END:if-univ.p.values
    
    
    ### NPC and create object "PermSpace"; subsets and weights are 
    ### extracted from "..."
    if(testStatistic == "Hotelling")
    {
        T.H0Low <- T
        if(p.valuesType == "asymptotic")
        {
            T.H0Gre <- 1 - T.H0Low
        } else
        {
            T.H0Gre <- -T.H0Low
        }# END:asymptotic
        
        permSpace <- new(
            Class      = "PermSpace",
            seed       = seed,
            T.H0Low    = T.H0Low,
            T.H0Gre    = T.H0Gre,
            P.H0Low    = t2p(T.H0Low),
            P.H0Gre    = t2p(T.H0Gre),
            rawStats   = array(0, dim = c(0, 0, 0)),
            comb.funct = combFunct)
    } else
    {
        ### change behaviour: if iteratedNPC then add permSpace2 to the output list
#        permSpace <- NPC (rawStats = T, combFun = combFunct, seed = seed, 
#                p.values = (p.valuesType == "asymptotic"), 
#                tails = tails, subsets = dots$subsets, weights = dots$weights, 
#                iteratedNPC = FALSE)
        permSpace <- NPC (rawStats = T, combFun = combFunct, seed = seed, 
                p.values = (p.valuesType == "asymptotic"), 
                tails = tails, subsets = dots$subsets, weights = dots$weights, 
                iteratedNPC = iteratedNPC)
    }# END:if-permSpace-generation
    permSpace@IDs <- permIndexMat
    
    ### create object "RankResults"
    rankRes <- rankingRule(permSpace = permSpace, alpha = alpha, 
            multAdjMethod = p.adj.method, groupsLabs = unique(groups))
    
    ### conditional removal of PermSpace
    if (!returnPermSpace)
    {
        permSpace <- new(Class = "PermSpace", seed = seed, rawStats = T)
    } else {}
    
    ### create object "SoupObject"
    soupRes <- new("SoupObject",
        call        = match.call(),
        rankResults = as.list(rankRes),
        pValueMat  = pValueMat,
        permSpace   = permSpace)
    
    return(soupRes)
    
}#=END=


