library(shiny)
library(DVHmetrics)
library(markdown)

trimWS <- function(x, side="both")  {
    side <- match.arg(side, c("left", "right", "both"))
    pattern <- switch(side, left="^\\s+", right="\\s+$", both="^\\s+|\\s+$")
    gsub(pattern, "", x)
}

getStrIDs <- function(x, what=c("structure", "patient"), choices=FALSE) {
    UseMethod("getStrIDs")
}

getStrIDs.DVHLst <- function(x, what=c("structure", "patient"), choices=FALSE) {
    what <- match.arg(what)

    strids <- if(what == "structure") {
        unique(vapply(x, function(y) y$structure, character(1)))
    } else {
        unique(vapply(x, function(y) y$patID,     character(1)))
    }

    revList <- as.list(seq_along(sort(strids)))
    if(choices) {
        return(setNames(revList, strids))
    } else {
        return(strids)
    }
}

getStrIDs.DVHLstLst <- function(x, what=c("structure", "patient"), choices=FALSE) {
    what <- match.arg(what)

    strids <- if(what == "structure") {
        unique(unlist(lapply(x, function(y) {
            lapply(y, function(z) z$structure) })))
    } else {
        unique(unlist(lapply(x, function(y) {
            lapply(y, function(z) z$patID) })))
    }

    revList  <- as.list(seq_along(sort(strids)))
    if(choices) {
        return(setNames(revList, strids))
    } else {
        return(strids)
    }
}

constrOut <- c("patID"="1", "structure"="2", "constraint"="3", "observed"="4",
               "compliance"="5", "deltaV"="6", "deltaVpc"="7", "deltaD"="8",
               "deltaDpc"="9", "dstMin"="10", "dstMinRel"="11",
               "ptMinD"="12", "ptMinV"="13")
constrOutInv <- c("1"="patID", "2"="structure", "3"="constraint", "4"="observed",
                  "5"="compliance", "6"="deltaV", "7"="deltaVpc", "8"="deltaD",
                  "9"="deltaDpc", "10"="dstMin", "11"="dstMinRel",
                  "12"="ptMinD", "13"="ptMinV")

## Hodges-Lehman pseudo median
pseudoMed <- function(x) {
    wilcox.test(x, conf.int=TRUE)$estimate
}
