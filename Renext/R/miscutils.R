
## To be used later in a summary or coef

coefftable <- function(fit) {
  s <- sqrt(diag(fit$cov))
  res <- cbind(fit$estimate, s, fit$estimate/s)
  colnames(res) <- c("est.", "sd", "t")
  res
}

##' Describe a period of time knowing its start and end using
##' a concise format.
##'
##' Depending on the length of the period, the date formats for
##' the start and end will be more or less detailled.
##' 
##' @title Describe a period of time knowing its start and end using
##' a concise format
##' 
##' @param start A vector that can be coerced into \code{POSIXct}. 
##'
##' @param end A vector that can be coerced into \code{POSIXct} and
##' with the same length as \code{start}.
##'
##' @return A character vector of the same length as \code{start}
##' and \code{end} with a period description.
##'
##' @examples
##' Renext:::formatPeriod(start = rep("1900-01-01", 3),
##'                       end = c("1900-01-10", "1901-01-10", "2010-01-10"))
##' 
formatPeriod <- function(start, end) {
    n <- length(start)
    if (length(end) != n) stop("'start' and 'end' must have the same length")
    start <- as.POSIXct(start)
    end <- as.POSIXct(end)
    dur <- difftime(end, start, units = "day") / 365
    text <- rep("", n)
    ind <- (dur > 100)
    if (any(ind)){
        text[ind] <- paste(format(start[ind], "%Y"),
                           format(end[ind], "%Y"), sep = "-")
    }
    ind <- (dur > 10) & (dur <= 100)
    if (any(ind)){
        text[ind] <- paste(format(start[ind], "%Y-%m"),
                           format(end[ind], "%Y-%m"), sep = " to ")
    }
    ind <- (dur <= 10) 
    if (any(ind)){
        text[ind] <- paste(format(start[ind], "%Y-%m-%d"),
                           format(end[ind], "%Y-%m-%d"), sep = " to ")
    }
    text
}
##******************************************************************************
##' Replace NULL names by empty names.
##'
##' @title Replace NULL names by empty names
##' 
##' @param x A vector.
##' 
##' @return A character vector with the same length as \code{x}.
##'
fillNames <- function(x) {
    if (!is.atomic(x)) stop("'x' is intended to be atomic here")
    nm <- names(x)
    if (is.null(nm)) nm <- rep("", length(x))
    nm
}

##******************************************************************************
##' Transform a 'show' indication in a versatile fashion.
##'
##' @title Transform a 'show' indication in a versatile fashion
##'
##' @param show a logical or character vector indicating which
##' elements must be shown.
##'
##' @param x an object of class \code{"Renouv"} containing the information
##' about \code{MAX} or \code{OTS} blocks.
##' 
##' @param type character. Can be \code{"MAX"} or {"OTS"}.
##' 
##' @return A logical vector with one element for each block.
##'
##' @examples
##' fit1 <- Renouv(x = rexp(100), threshold = 0, effDuration = 100,
##'                MAX.effDuration = c(100, 40),
##'                MAX.data = list("deluge" = c(1.2, 3.4, 5.1), "dryness"= c(0.1, 0.3)))
##'
##' res1a <- Renext:::transShow(show = TRUE, x = fit1)
##' res1b <- Renext:::transShow(show = "deluge", x = fit1)
##'
##' fit2 <- Renouv(x = rexp(100), threshold = 0, effDuration = 100,
##'                OTS.effDuration = c(100, 40),
##'                OTS.threshold = c(1, 0.05),
##'                OTS.data = list("deluge" = c(1.2, 3.4, 5.1), "dryness"= numeric(0)) )
##'
##' res2a <- Renext:::transShow(show = TRUE, x = fit2, type = "OTS")
##' res2b <- Renext:::transShow(show = "deluge", x = fit2, type = "OTS")
##' 
##' 
transShow <- function(show, x, type = "MAX") {

    h <- x[[paste("history", type, sep = ".")]]
    if (!is.null(show) && h$flag) {
        nB <- nlevels(h$block)
        if (is.character(show)) {
            show1 <- rep(FALSE, nB)
            for (ib in 1L:length(show)) {
                ind <- grepl(show[ib], h$blockNames)
                show1[ind] <- TRUE
            }
            show <- show1
        } else if (is.logical(show)) {
            if (length(show) == 1L)  {
                show <- rep(show, nB)
            } else if (length(show) == nB) {
                if (!is.null(names(show))) {
                    if (min(nchar(names(show))) > 0L) {
                        ind <- match(names(show), h$blockNames)
                        if (any(is.na(ind))) {
                            stop("in show[[\"", type, "\"]], invalid names: ",
                                 names(show)[is.na(ind)])
                        }
                        show <- show[ind]
                    } 
                }
            } else stop("a logical 'show' can have length 1 or ", nB)
        }
        names(show) <- paste("block", 1L:nB, sep = "")
        return(show) 
    } else {
        return (FALSE)
    }
}

##' Replace elements in a hierarchical list as \code{RLPar}
##'
##' @title Replace elements in hierarchical list
##'
##' @param from The list of change elements.
##'
##' @param to The original list which is copied and edited.
##'
##' @param trace Level of verbosity.
##' 
##' @return A hierarchical list containing all elements of \code{to}
##' with the same (hierarchical) names and their value possibly
##' modified when they appear in \code{from}
##'
##' @examples
##'
##' .par <- RLpar(mono = mono)
##' From <-  list(quant = c(col = "azure"), OT = c(pch = 1))
##' Renext:::rReplace(from = From, to = .par, trace = 1)$quant
##' 
`rReplace` <- function(from, to, trace = 0) {
  
  ulFrom <- unlist(from)
  ulTo <- unlist(to)
  nm0 <- match(names(ulFrom), names(ulTo))
    
  if (any(is.na(nm0))) {
    warning("'from' contains unused names:\n",
            sprintf("\"%s\"", names(ulFrom)[is.na(nm0)]))
  }
  
  Names0 <- names(ulFrom)[!is.na(nm0)]

  ## This can not be vectorized since nmVec is a vector with elements for
  ## the hierarchical levels 1, 2, ...
  
  for (Name0 in Names0){
    nmVec <- unlist(strsplit(Name0, split = "\\."))
    if (trace) cat("name = ", nmVec, ", value = ", from[[nmVec]], "\n")
    to[[nmVec]] <- from[[nmVec]]
  }
  
  to 
  
}
