##*****************************************************************************
##      DO NOT ROXYGENISE THIS!!!
## 
##' Compute the empirical survival values and the empirical return
##' periods at the observations of an object. These are used as
##' plotting positions in several plots.
##'
##' When the object contains historical information (\code{MAX} or
##' \code{OTS}), the computation is an adaptation Hirsch and Stedinger
##' (1987) for the Marked Process (MP) context. The original method is
##' devoted to block maxima and interpolates the survival which are
##' first computed at the thresholds. For MP, the interpolation is
##' done for the \emph{inverse return periods}, and the survival values are
##' deduced from those of the inverse return periods.
##' 
##' @title Compute empirical survivals and return periods
##'
##' @param object The object containing the data, with class
##' \code{"Renouv"} or \code{"Rendata"}.*
##' 
##' @param points Option for the computation of the plotting
##' positions.  When \code{points} is set to \code{"p"}, the
##' \eqn{p}-points formula is used with the selected value of
##' \code{a}. This formula is used to compute the survival from which
##' the return period is computed. When instead \code{points} is set
##' to \code{"H"}, \emph{Nelson's formula} is used to compute the
##' return periods, the survival value still being given by the
##' \eqn{p}-points formula. When the data is heterogenous, i.e.  when
##' \code{object} contains \code{MAX} and/or \code{OTS} data, Nelson's
##' formula is used only to compute the return periods of the upper
##' slice of levels.
##' 
##' @param a Parameter used in the interpolation formula for the
##' inverse return periods as in Hirsch and Stedinger (1987).
##'
##' 
##' @return A list with the following elements
##' 
##' \item{x}{Numeric vector containing the ordered values from all
##' the available sources in the object: main sample, historical
##' periods either 'MAX' or 'OTS'.  }
##'
##'\item{group, groupNames}{ Integer and character vectors giving the
##' source of the values in \code{x} (in the same order of the
##' values. For instance, \code{group[10]} gives the group form which
##' \code{x[10]} was extracted, and the name of this group is
##' \code{groupNames[group[10]]}.}
##'
##'\item{S, T}{Numeric vectors of the same length as \code{x} and
##' containing the corresponding estimation of the survival value and
##' of the return period.}
##' 
##' \item{thesh, lambda.thresh, S.thresh, T.thresh}{Vector of
##' thresholds and the corresponding estimation for the event rate,
##' survival and return period. All the estimations are \emph{for the
##' threshold values}. The value of \code{T.thresh[i]} for a threshold
##' \code{thresh[i]} results from a simple computation: divide the sum of the
##' durations for blocks with thresholds \code{>= thresh[i]} by the number of events
##' for these blocks.}
##' 
##' @references
##' The original method for block maxima is described in  
##'  
##' Hisch R.M. and Stedinger J.R.(1887) Plotting Positions for Historical Floods
##' and their precision. \emph{Water Ressources Research}, vol. 23, N. 4 pp. 715-727.
##'
##' Millard S. and Neerchal N. (2001). \emph{Environmental Statistics with S-Plus}. CRC Press
##'
##' The adaptation for the Marked Process context is described in the
##' \emph{Renext Computing Details} document.
##' @author Yves Deville
##' @examples
##' ## use an object with class "Rendata"
##' ST1 <- SandT(object = Garonne)
##' ## basic return level plot
##' plot(ST1$T, ST1$x, col = ST1$group, log = "x")
##' ## use an object with class "Renouv"
##' fit <- Renouv(x = Garonne)
##' ST2 <- SandT(object = fit)
##' plot(ST2$T, ST2$x, col = ST2$group, log = "x")
SandT <- function(object,
                  points = c("p", "H"),
                  a = 0,
                  naive = FALSE) {
    
    if (naive) {
        if (!missing(a) || !missing(points)) {
            warning("when 'naive' is TRUE, formals 'a' and 'points' are ignored")
        }
        return(oldSandT(object))
    }
    
    points <- match.arg(points)
    
    if (class(object) == "Rendata") {
        vn  <- object$info$varName
        if (!is.null(object$OTdata)){
            newObj <- list(effDuration = object$OTinfo$effDuration,
                           x.OT = object$OTdata[ , vn])
        } else newObj <- list()
        newObj$min <- min(newObj$x.OT)
        OTS <- makeOTSdata(object)
        newObj$history.OTS <- OTS
        MAX <- makeMAXdata(object)
        newObj$history.MAX <- MAX
        object <- newObj
    } else if (class(object) == "Renouv") {
        ## clean and change object
        x.OT <- object$threshold + object$y.OT
        object <- list(effDuration = object$effDuration,
                       x.OT = x.OT,
                       min = object$threshold,
                       history.OTS = object$history.OTS,
                       history.MAX = object$history.MAX)
    } else {
        stop("'object' must have (old)class \"Rendata\" or \"Renouv\"")
    }
    
    if (!is.null(object$x.OT) && length(object$x.OT) > 0) {
        vals <- object$x.OT
        samp <- rep(1L, length(object$x.OT))
        names(samp) <- names(vals) <- rep("OT", length(object$x.OT))
        mins <- object$min
        durs <- object$effDuration
        groupNames <- "OT"
        sourceNames <- "OT"
        names(mins) <- names(durs) <- groupNames[1L]
        nBlock <- 1L
    }  else {
        vals <- numeric(0)
        samp <- numeric(0)
        mins <- numeric(0)
        durs <- numeric(0)
        groupNames <- character(0)
        sourceNames <- character(0)
        nBlock <- 0L
    }
    
    ## manage in the same way 'OTS' periods
    if (!is.null(object$history.OTS) && object$history.OTS$flag){
        r <- object$history.OTS$r
        nms <- paste("OTS.block", 1L:length(r), sep = "")
        groupNames <- c(groupNames, nms)
        sourceNames <- c(sourceNames, rep("OTS", length(r)))
        newDurs <- object$history.OTS$effDuration
        newMins <- object$history.OTS$threshold
        names(newDurs) <- names(newMins) <- nms
        durs <- c(durs, newDurs)
        mins <- c(mins, newMins)
        newSamp <- nBlock + rep(1L:length(r), times = r)
        newVals <- unlist(object$history.OTS$data)
        names(newSamp) <- names(newVals) <- rep(nms, times = r)
        samp <- c(samp, newSamp)
        vals <- c(vals, newVals)
        nBlock <- nBlock + length(r)
    }
    
    ## manage 'MAX' periods
    if (!is.null(object$history.MAX) && object$history.MAX$flag){
        r <- object$history.MAX$r
        nms <- paste("MAX.block", 1L:length(r), sep = "")
        groupNames <- c(groupNames, nms)
        sourceNames <- c(sourceNames, rep("MAX", length(r)))
        newSamp <- nBlock + rep(1L:length(r), times = r)
        newVals <- unlist(object$history.MAX$data)
        names(newSamp) <- names(newVals) <- rep(nms, times = r)
        samp <- c(samp, newSamp)
        vals <- c(vals, newVals)
        nBlock <- nBlock + length(r)
        newDurs <- object$history.MAX$effDuration
        d <- diff(sort(vals))
        smallVal<- min(d[d > 0]) * 0.5
        newMins <- unlist(lapply(object$history.MAX$data, min)) - smallVal
        names(newDurs) <- names(newMins) <- nms
        durs <- c(durs, newDurs)
        mins <- c(mins, newMins)
    }
    
    ## order the thresholds
    ind <- order(mins)
    mins <- mins[ind]
    durs <- durs[ind]
    
    ## use more algorithmic notations as usesd in the 'Renext Computind Details'
    ## document, and remove ex aequo if necessary
    ##
    ## BUG fixed here >= 2.1-6 old code here :
    ## keep <- c(TRUE, diff(mins) > 0); u <- mins[keep]; W <- cumsum(durs)[keep]
    u <- mins
    W <- cumsum(durs)
    J <- length(u)

    ## for objects such as those created with 'RenouvNoEst'
    if (J == 0) {
        return(list(x = numeric(0),
                    group = numeric(0),
                    groupNames = character(0),
                    sourceNames = character(0),
                    lambda = NA,
                    S = numeric(0),
                    T = numeric(0),
                    thresh = c(object$threshold, Inf),
                    lambda.thresh = numeric(0),
                    S.thresh = numeric(0),
                    T.thresh = numeric(0)))
    }
    
    ind <- order(vals)
    vals <- vals[ind]
    samp <- samp[ind]
    
    ##===========================================================================
    ## First pass: compute A, D, rate and inverse return period AT THE
    ## THRESHOLDS store the indices of vals falling betwwen the
    ## thresholds.
    ## ===========================================================================
    interv <- findInterval(x = vals, vec = u, rightmost.closed = TRUE)
    A <- D <- invT.thresh <- lambda.thresh <- rep(0, J + 1L)
    ind <- list()
    for (j in J:1L) {
        ind[[j]] <- (interv == j)
        A[j] <- sum(ind[[j]])
        lambda.thresh[j] <- A[j] / W[j] 
        D[j] <- D[j + 1L] + A[j] + lambda.thresh[j] * (W[J] - W[j])
        invT.thresh[j] <- D[j] / W[J] 
    }
    
    lambda <- sum(lambda.thresh[1:J])
    ## cat("lambda = ", lambda, "\n")
    
    ##===========================================================================
    ## Second pass: compute the inverse retun period at the data by
    ## interpolation, as well as the survival.  The 'Hpoints' options
    ## can be used only (at the time) for the upper slice of values.
    ## ===========================================================================
    T <- rep(NA, length(vals))
    S <- rep(NA, length(vals))
    Aprec <- FALSE
    
    for (j in J:1L) {
        if (A[j]) {
            ## reverse order
            prov <- invT.thresh[j + 1L] +
                (invT.thresh[j] - invT.thresh[j + 1L]) * ((A[j]:1L) - a) / (A[j] -2*a + 1)
            S[ind[[j]]] <- prov / lambda
            if (!Aprec && points == "H") {
                Sprov <- invT.thresh[j] / lambda
                Hprov <- -log(Sprov)
                Hprov <- Hprov + Hpoints(A[j])
                T[ind[[j]]]  <- exp(Hprov) / lambda
            } else {
                ## cat(sprintf("j = %d\n", j))
                T[ind[[j]]]  <- 1 / S[ind[[j]]] / lambda
            }
            Aprec <- TRUE
        }
    }
    
    list(x = vals,
         group = samp,
         groupNames = groupNames,
         sourceNames = sourceNames,
         lambda = lambda,
         S = S,
         T = T,
         thresh = c(u, Inf),
         lambda.thresh = lambda.thresh,
         S.thresh = invT.thresh / lambda,
         T.thresh = 1 / invT.thresh)
    
}

##*****************************************************************************
##  Private function to maintain compatibility.
##*****************************************************************************

oldSandT <- function(object) {

    Sfun <- function(r) 1.0 / (r:1L)
    
    if (!is.null(object$x.OT) && length(object$x.OT) > 0) {
        vals <- sort(object$x.OT)
        nx.OT <- length(object$x.OT)
        samp <- rep(1L, nx.OT)
        names(samp) <- names(vals) <- rep("OT", length(object$x.OT))
        mins <- min(object$x.OT)
        durs <- object$effDuration
        groupNames <- "OT"
        sourceNames <- "OT"
        names(mins) <- names(durs) <- groupNames[1L]
        nBlock <- 1L
        S <- 1 - (1L:nx.OT) / (nx.OT + 1)
        lambda <- nx.OT / durs
        T <- 1 / lambda / S
    }  else {
        stop("'naive' plotting postions can be used only ",
             "when there is 'OT' data combined with 'MAX' or 'OTS' data")
    }
    
    ## manage in the same way 'OTS' periods
    if (!is.null(object$history.OTS) && object$history.OTS$flag){
        r <- object$history.OTS$r
        nms <- paste("OTS.block", 1L:length(r), sep = "")
        groupNames <- c(groupNames, nms)
        sourceNames <- c(sourceNames, rep("OTS", length(r)))
        newDurs <- object$history.OTS$effDuration
        newMins <- object$history.OTS$threshold
        names(newDurs) <- names(newMins) <- nms
        Npred <- lambda * newDurs
        Npred[Npred < r] <- r[Npred < r]
        newSamp <- nBlock + rep(1L:length(r), times = r)
        newVals <- unlist(lapply(object$history.OTS$data, sort))
        newS <- unlist(sapply(r, Sfun))
        a <- (Npred + 1L) / Npred
        newT <- rep(newDurs * a, times = r) * newS
        names(newSamp) <- names(newVals) <- rep(nms, times = r)
        samp <- c(samp, newSamp)
        vals <- c(vals, newVals)
        T <- c(T, newT)
        nBlock <- nBlock + length(r)
    }
    
    ## manage 'MAX' periods
    if (!is.null(object$history.MAX) && object$history.MAX$flag){
        r <- object$history.MAX$r
        nms <- paste("MAX.block", 1L:length(r), sep = "")
        groupNames <- c(groupNames, nms)
        sourceNames <- c(sourceNames, rep("MAX", length(r)))
        newDurs <- object$history.MAX$effDuration
        newSamp <- nBlock + rep(1L:length(r), times = r)
        newVals <- unlist(lapply(object$history.MAX$data, sort))
        Npred <- lambda * newDurs
        Npred[Npred < r] <- r[Npred < r]
        newS <- unlist(lapply(r, Sfun))
        a <- (Npred + 1L) / Npred
        newT <- rep(newDurs * a, times = r) * newS
        names(newSamp) <- names(newVals) <- rep(nms, times = r)
        samp <- c(samp, newSamp)
        vals <- c(vals, newVals)
        T <- c(T, newT)
        nBlock <- nBlock + length(r)
    }
    
    list(x = vals,
         group = samp,
         groupNames = groupNames,
         sourceNames = sourceNames,
         lambda = lambda,
         T = T)
    
}
