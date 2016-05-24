P4.4 <- function (ID, digits = 4, rela = FALSE, parft=FALSE, parftdf=FALSE, 
    detailed=FALSE)
{
    ## function used for calculations in the paper
    ## final updated function P4.4 for R package DoE.base will omit "newdf" 
    ##     and will rename "new" to "parft"

    ## function to calculate pattern of numbers of generalized length 4 words
    ## for projections into four factors
    if (!is.logical(rela))
        stop("rela must be TRUE or FALSE")
    if (!is.logical(parft))
        stop("parft must be TRUE or FALSE")
    if (!is.logical(parftdf))
        stop("parftdf must be TRUE or FALSE")
    ## retrieve child array or array identified by character string
          ## gsub for case where ID is character string
    IDname <- gsub("\"", "", deparse(substitute(ID)))
    if (all(IDname %in% oacat$name)) {
        if (!exists(IDname))
            ID <- eval(parse(text = paste("oa.design(", IDname,
                ")")))
        else if (is.character(ID))
            ID <- eval(parse(text = paste("oa.design(", IDname,
                ")")))
    }
    if ((rela || parft || parftdf) & !(isTRUE(all.equal(length2(ID), 0)) & isTRUE(all.equal(length3(ID),
        0))))
        stop(IDname, " is not a strength 3 array, \nP4.4 with rela, parft or parftdf TRUE is inadequate.")
    if (!(is.data.frame(ID) | is.matrix(ID)))
        stop("ID must be a data frame or a matrix")
    if (is.matrix(ID))
        ID <- as.data.frame(ID)
    if (!ncol(ID) >= 4)
        return(NULL)
    hilf <- length4(ID, J=TRUE)
    fhilf <- factor(names(hilf), levels=unique(names(hilf))) ## bug fix 2 Sep 2013
                       ## hilf was in unexpected order before, 
                       ## yielding wrong calculations for rela in designs 
                       ## with many columns due to character instead of 
                       ## numeric sorting
    hilf <- sapply(split(hilf, fhilf), function(obj) sum(obj^2))
    if (rela) {
        waehl <- nchoosek(ncol(ID), 4)
        nlevels <- sapply(ID, function(obj) length(unique(obj)))
        div <- apply(waehl, 2, function(obj) min((nlevels[obj] -
            1)))
    }
    if (parft) {
        waehl <- nchoosek(ncol(ID), 4)
        nlevels <- sapply(ID, function(obj) length(unique(obj)))
        div <- apply(waehl, 2, function(obj) 4/sum(1/(nlevels[obj] -
            1)))  ## divisor, multiplier is 1/divisor
    }
    if (parftdf) {
        waehl <- nchoosek(ncol(ID), 4)
        nlevels <- sapply(ID, function(obj) length(unique(obj)))
        div <- apply(waehl, 2, function(obj) mean(nlevels[obj] -
            1))  ## divisor, multiplier is 1/divisor
    }

    aus <- table(round(hilf, digits))
    if (rela || parft || parftdf) aus <- table(round(hilf/div, digits))
    
    ## formatting the table for output
    aus <- cbind(length4 = as.numeric(names(aus)), frequency = aus)
    if (rela)
        colnames(aus) <- c("length4.rela", "frequency")
    if (parft)
        colnames(aus) <- c("length4.parft", "frequency")
    if (parftdf)
        colnames(aus) <- c("length4.parftdf", "frequency")

    rownames(aus) <- rep("", nrow(aus))

    ## attaching attributes
    attr(aus, "A4") <- A4 <- sum(hilf)
    if (!(rela || parft || parftdf) & detailed & A4>0) 
        attr(aus, "detail") <- round(hilf, digits)
    if (rela) {
        attr(aus, "rA4") <- rA4 <- sum(hilf/div)
        if (rA4 > 0) attr(aus, "GR") <- round(4+1-sqrt(max(hilf/div)),digits)
          else attr(aus, "GR") <- ">=5"
        if (rA4 > 0 & detailed) attr(aus, "detail") <- round(hilf/div, digits)
        }
      if (parft){
        attr(aus, "sumPARFT4") <- sumPARFT4 <- sum(hilf/div)
        if (sumPARFT4 > 0 & detailed) 
            attr(aus, "detail") <- round(hilf/div, digits)
        }
      if (parftdf){
        attr(aus, "sumPARFTdf4") <- sumPARFTdf4 <- sum(hilf/div)
        if (sumPARFTdf4 > 0 & detailed)
            attr(aus, "detail") <- round(hilf/div, digits)
        }
    aus
}
