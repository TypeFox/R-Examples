#' Apply flips and swaps to a record. 
#'
#' @param flips Index vector in r. r will be sign-flipped at flips
#' @param swaps nx2 matrix denoting value swaps in r.
#' @param r numerical record.
#' @return r, with flips and swaps applied
#'
#' @keywords internal
applyFix <- function(flips, swaps, r){
    if ( length(flips) > 0 )
        r[flips] <- -r[flips]
    if ( !is.null(nrow(swaps)) ) 
        apply(swaps, 1, function(sw){r[sw] <<- r[sw[2:1]]})
    return(r)
}

#' Workhorse for correctSigns
#'
#' 
#' @title workhorse for correctSigns
#'
#' @param r The numerical record to correct
#' @param A1 Equality matrix
#' @param C1 Constant vector for equalities
#' @param eps Tolerance for equality-checking
#' @param A2 Inequality matrix
#' @param C2 Constant vector for inequalities
#' @param epsvec Vector to check against. (\code{.Machine.double.eps} for \code{>} inequalities, otherwise 0.)
#' @param flip indices in \code{r}, where r may be sign-flipped
#' @param swap \eqn{n\times2} matrix with indices in \code{r}, each row indicating a possible value swap.
#' @param w weight vector of \code{length(flips)+nrow(swaps)} if \code{swapIsOneFlip==TRUE}, otherwise of \code{length(r)}
#' @param swapIsOneFlip logical. If \code{TRUE}, weights are assigned to each flip or swap action. If \code{FALSE} weights are assigned to every changed variable. 
#' @param maxActions Maximum number of \code{flips+swaps} to try.
#' @param maxCombinations the maximum number of flip/swap combinations to try. See the description in \code{\link{correctSigns}}
#'
#' @return a list containing
#'  \tabular{ll}{
#'      \code{S} \tab \code{n x length(r)} array with corrected versions of \code{r} \cr
#'      \code{weight} \tab vector of length \code{n} with total weight for each solution \cr
#'      \code{nFlip} \tab number of sign flips for every solution \cr
#'      \code{nSwap} \tab number of value swaps for every solution\cr
#' }
#'
#' @seealso \code{\link{correctSigns}}
#'
#'
getSignCorrection <- function( r, A1, C1, eps, A2, C2, epsvec, flip, swap, w, 
    swapIsOneFlip, maxActions, maxCombinations ){
    nflip <- length(flip)
    ntot <- nflip + nrow(swap) 
     
    S <- array(dim=c(0,length(r)))
    weight <- nFlip <- nSwap <- numeric(0)
    i <- 0
    for( nActions in 1:maxActions ){
        if ( choose(ntot, nActions) > maxCombinations ) break 
        I <- combn(ntot, nActions)
        for ( k in 1:ncol(I) ){
            flips <- flip[I[I[ ,k] <= nflip, k]]
            swaps <- swap[I[I[ ,k] >  nflip, k]-nflip,,drop=FALSE]
            s <- applyFix(flips, swaps, r)
            if ( all(abs(A1 %*% s - C1) < eps) && all(A2 %*% s - C2 < epsvec) ){
                i <- i + 1
                S <- rbind(S,s)
                if ( swapIsOneFlip ){
                    weight[i] <- sum(w[I[,k]])
                } else {
                    iFlip <- I[I[ ,k] <= nflip, k] 
                    iSwap <- I[I[ ,k] >  nflip, k]-nflip 
                    weight[i] <- sum(w[flip[I[iFlip,k]]]) + sum(w[swap[I[iSwap,k]]])
                }
                nFlip[i] <- length(flips)
                nSwap[i] <- nrow(swaps)
            }
        }
        if (length(S) > 0) break
    }
    return(list(S=S, weight=weight, nFlip=nFlip, nSwap=nSwap))
}

#' Correct sign errors and value interchanges in data records. 
#'
#' This algorithm tries to correct records violating linear equalities by sign flipping and/or value interchanges.
#' Linear inequalities are taken into account when judging possible solutions. If one or more inequality restriction
#' is violated, the solution is rejected. It is important to note that the \code{\link{status}} of a record has
#' the following meaning:
#' 
#' \tabular{ll}{
#' \code{valid} \tab The record obeys all equality constraints on entry. No error correction is performed. \cr
#' \code{}      \tab It may therefore still contain inequality errors.\cr
#' \code{corrected} \tab Equality errors were found, and all of them are solved without violating inequalities.\cr
#' \code{partial}\tab Does not occur\cr
#' \code{invalid} \tab The record contains equality violations which could not be solved with this algorithm\cr
#' \code{NA} \tab record could not be checked. It contained missings.
#' }
#'
#' The algorithm applies all combinations of (user-allowed) flip- and swap combinations to find a solution, and minimizes 
#' the number of actions (flips+swaps) that have to be taken to correct a record. When multiple solutions are found, the
#' solution of minimal weight is chosen. The user may provide a weight vector with weights for every flip and every swap,
#' or a named weight vector with a weight for every variable. If the weights do not single out a solution, the first one
#' found is chosen.
#'
#' If arguments \code{flip} or \code{swap} contain a variable not in \code{E}, these variables will be ignored by the algorithm.
#'
#'
#'
#' @title Correct sign errors and value interchanges in data records
#'
#' @param E An object of class \code{\link[editrules:editmatrix]{editmatrix}}
#' @param dat \code{data.frame}, the records to correct.
#' @param ... arguments to be passed to other methods.
#' @param flip A \code{character} vector of variable names who's values may be sign-flipped
#' @param swap A \code{list} of \code{character} 2-vectors of variable combinations who's values may be swapped
#' @param maxActions The maximum number of flips and swaps that may be performed
#' @param maxCombinations The number of possible flip/swap combinations in each step of the algorithm is \code{choose(n,k)}, with \code{n}
#'      the number of \code{flips+swaps}, and \code{k} the number of actions taken in that step. If \code{choose(n,k)} exceeds \code{maxCombinations},
#'      the algorithm returns a record uncorrected.
#' @param eps Tolerance to check equalities against. Use this to account for sign errors masked by rounding errors.
#' @param weight weight vector. Weights can be assigned either to actions (flips and swap) or to variables.
#'      If \code{length(weight)==length(flip)+length(swap)}, weights are assiged to actions, if \code{length(weight)==ncol(E)}, weights
#'      are assigned to variables. In the first case, the first \code{length{flip}} weights correspond to flips, the rest to swaps. 
#'      A warning is issued in the second case when the weight vector is not named. See the examples for more details.
#' @param fixate a \code{character} vector with names of variables whos values may not be changed
#' @return a \code{\link{deducorrect-object}}. The \code{status} slot has the following columns for every records in \code{dat}.
#'
#'  \tabular{ll}{
#'      \code{status}\tab a \code{\link{status}} factor, showing the status of the treated record.\cr
#'      \code{degeneracy}\tab the number of solutions found, \emph{after} applying the weight\cr
#'      \code{weight}\tab the weight of the chosen solution\cr
#'      \code{nflip}\tab the number of applied sign flips\cr
#'      \code{nswap}\tab the number of applied value interchanges\cr
#'  }
#' @example ../examples/correctSigns.R
#' @references
#' Scholtus S (2008). Algorithms for correcting some obvious
#' inconsistencies and rounding errors in business survey data. Technical
#' Report 08015, Netherlands.
#' @seealso \code{\link{deducorrect-object}}
#' @export
correctSigns <- function(E,dat, ...){
    UseMethod("correctSigns")
}

#'
#' @method correctSigns editset
#' @rdname correctSigns
#' @export
correctSigns.editset <- function(E, dat, ...){
    correctAndRevert(correctSigns.editmatrix, E, dat, ...)
}   


#' @method correctSigns editmatrix
#' @rdname correctSigns
#' @export
correctSigns.editmatrix <- function(
    E, 
    dat,
    flip = getVars(E),
    swap = list(),
    maxActions = length(flip)+length(swap),
    maxCombinations = 1e5,
    eps=sqrt(.Machine$double.eps),
    weight = rep(1,length(flip)+length(swap)),
    fixate = NA,
    ...){

    ops <- getOps(E)
    if ( !isNormalized(E) ) E <- normalize(E)
  
    vars <- getVars(E)
 
    # remove variables in flip, not occuring in E
    flip <- flip[flip %in% vars]
    # fixation variables. 
    if ( !is.na(fixate) ) flip <- setdiff(flip,fixate)
    # remove fixated variables and pairs with a variable not in E from swaplist
    if (length(swap)>0){
        for ( i in 1:length(swap) ){ 
            if (    any(swap[[i]] %in% fixate) 
                 || any(!(swap[[i]] %in% vars) )
               ) swap[[i]] <- NULL
        }
    }                   
   

    # determine if flipIsOneSwap=TRUE or not by length of *weight*
    if (length(weight) == length(flip) + length(swap)){
        swapIsOneFlip <- TRUE
    } else if (length(weight) == length(vars) ){
        swapIsOneFlip <- FALSE
        if ( all(names(weight) %in% vars) & length(names(weight))==length(vars )){
            weight <- weight[vars]
        } else {
            names(weight) <- vars
            warning(paste("Weight vector has no names. Assuming same order as getVars(E)"))
            
        }
    } else {
        stop("Length of weight vector does not correspond with number of actions or number of columns in editmatrix")
    }

    # prepare matrices and constants
    eq <- ops == "=="
    A1 <- getA(E)[eq, ,drop=FALSE]
    C1 <- getb(E)[eq]
    A2 <- getA(E)[!eq, ,drop=FALSE]
    C2 <- getb(E)[!eq]
    epsvec <- ifelse(ops[!eq]=="<=", 100*.Machine$double.eps, 0)
    
    D <- as.matrix(dat[, getVars(E)])
    # from flip and swap names to indices
    flip <- sapply(flip, function(fl) which(vars==fl))
    swap <- sapply(swap, function(sw) c(which(vars==sw[1]), which(vars==sw[2])))
    swap <- array(t(swap), dim=c(length(swap)/2, 2))
    status <- status(nrow(dat))
    wgt <- degeneracy <- nflips <- nswaps <- numeric(nrow(dat))
    corrections <- data.frame(row=numeric(),variable=factor(levels=vars),old=numeric(),new=numeric())
    
    for ( i in which(complete.cases(D)) ){
        r <- D[i, ]
        iViolated <- abs(A1 %*% r - C1)  > eps
        if ( !any(iViolated) ){
            status[i] <- "valid"
            next
        }
        adapt <- which(colSums(abs(A1[iViolated, ,drop=FALSE]))>0)
        flipable <- flip[flip %in% adapt]
        swapable <- array(dim=c(0,2))
        apply(swap, 1, function(sw) if ( all(sw %in% adapt) ) swapable <<- rbind(swapable,sw) )
        corr <- getSignCorrection(
            r  = r, 
            A1 = A1, 
            C1 = C1, 
            eps= eps, 
            A2 = A2, 
            C2 = C2, 
            epsvec = epsvec, 
            flip = flipable, 
            swap = swapable, 
            w = weight, 
            swapIsOneFlip = swapIsOneFlip, 
            maxActions = min(maxActions,length(flipable)+nrow(swapable)), 
            maxCombinations = maxCombinations
        )
        if ( nrow(corr$S) > 0 ){
            iMin <- which.min(corr$weight)
            D[i,] <- corr$S[iMin,]
            status[i] <- "corrected"
            wgt[i] <- min(corr$weight)
            degeneracy[i] <- length(corr$weight)
            nflips[i] <- corr$nFlip[iMin]
            nswaps[i] <- corr$nSwap[iMin]
            changed <- D[i,] != r
            corrections <- rbind(corrections,
                data.frame(
                    row = rep(i,sum(changed)),
                    variable = vars[changed],
                    old = r[changed],
                    new = D[i,changed,drop=TRUE]))
        } else {
            status[i] <- "invalid"
        }
    }
    dat[,vars] <- D
    rownames(corrections) <- NULL
    return(newdeducorrect(
        corrected   = dat,
        status      = data.frame(
            status  = status,
            weight  = wgt,
            degeneracy = degeneracy, 
            nflip = nflips, 
            nswap=nswaps),
            corrections = corrections,
            ...
    ))
}
