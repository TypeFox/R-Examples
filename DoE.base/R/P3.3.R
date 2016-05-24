P3.3 <- function (ID, digits = 4, rela = FALSE, parft=FALSE, parftdf=FALSE, 
    detailed=FALSE)
{
    ## function used for calculations in the paper
    ## final updated function P3.3 for R package DoE.base will omit "newdf" 
    ##     and will rename "new" to "parft"
    
    ## function to calculate pattern of numbers of generalized length 3 words
    ## for projections into three factors
    if (!is.logical(rela))
        stop("rela must be TRUE or FALSE")
    if (!is.logical(parft))
        stop("parft must be TRUE or FALSE")
    if (!is.logical(parftdf))
        stop("parft must be TRUE or FALSE")
    
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
    if ((rela || parft || parftdf) & !(isTRUE(all.equal(length2(ID), 0))))
        stop(IDname, " is not an orthogonal array, \nP3.3 with rela, parft or parftdf TRUE is inadequate.")
    if (!(is.data.frame(ID) | is.matrix(ID)))
        stop("ID must be a data frame or a matrix")
    if (is.matrix(ID))
        ID <- as.data.frame(ID)
    if (!ncol(ID) >= 3)
        return(NULL)
    hilf <- length3(ID, J = TRUE)
    fhilf <- factor(names(hilf), levels=unique(names(hilf))) ## bug fix 2 Sep 2013
                       ## hilf was in unexpected order before, 
                       ## yielding wrong calculations for rela in designs 
                       ## with many columns due to character instead of 
                       ## numeric sorting
    hilf <- sapply(split(hilf, fhilf), function(obj) sum(obj^2))
    if (rela) {
       # divisors are created from nlevels by relative numbers
        waehl <- nchoosek(ncol(ID), 3)
        nlevels <- sapply(ID, function(obj) length(unique(obj)))
        div <- apply(waehl, 2, function(obj) min((nlevels[obj] -
            1)))
    }
    if (parft) {
        waehl <- nchoosek(ncol(ID), 3)
        nlevels <- sapply(ID, function(obj) length(unique(obj)))
        div <- apply(waehl, 2, function(obj) 3/sum(1/(nlevels[obj] -
            1)))  ## divisor, multiplier is 1/divisor
    }
    if (parftdf) {
        waehl <- nchoosek(ncol(ID), 3)
        nlevels <- sapply(ID, function(obj) length(unique(obj)))
        div <- apply(waehl, 2, function(obj) mean(nlevels[obj] -
            1))  ## divisor, multiplier is 1/divisor
    }
    
    aus <- table(round(hilf, digits))
    if (rela || parft || parftdf)
        aus <- table(round(hilf/div, digits))
    
    ## formatting the table for output
    aus <- cbind(length3 = as.numeric(names(aus)), frequency = aus)
    if (rela)
        colnames(aus) <- c("length3.rela", "frequency")
    if (parft)
        colnames(aus) <- c("length3.parft", "frequency")
    if (parftdf)
        colnames(aus) <- c("length3.parftdf", "frequency")
    rownames(aus) <- rep("", nrow(aus))
    
    ## attaching attributes
    attr(aus, "A3") <- A3 <- sum(hilf)
    if (detailed & A3 > 0 & !(rela || parft || parftdf)) 
            attr(aus, "detail") <- round(hilf, digits)
    
    if (rela) {
            attr(aus, "rA3") <- rA3 <- sum(hilf/div)
            if (detailed & rA3 > 0) attr(aus, "detail") <- round(hilf/div, digits)
            if (rA3 > 0) attr(aus, "GR") <- round(3+1-sqrt(max(hilf/div)),digits)
                else attr(aus, "GR") <- ">=4"
        }
    if (parft){  
            attr(aus, "sumPARFT3") <- sumPARFT3 <- sum(hilf/div)
            if (detailed & sumPARFT3 > 0) 
                attr(aus, "detail") <- round(hilf/div, digits)
        }
    if (parftdf){  
            attr(aus, "sumPARFTdf3") <- sumPARFTdf3 <- sum(hilf/div)
            if (detailed & sumPARFTdf3 > 0) 
                attr(aus, "detail") <- round(hilf/div, digits)
        }
    aus
}
