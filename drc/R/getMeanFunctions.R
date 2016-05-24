"getMeanFunctions" <- function(noParm = NA, fname = NULL, flist = NULL, display = TRUE)
{
    if (is.null(flist))
    {
        fctList <- list(
        LL.2(), LL.3(), LL.3u(), LL.4(), LL.5(),
        W1.2(), W1.3(), W1.4(),
        W2.2(), W2.3(), W2.4(),
        BC.4(), BC.5(), 
        LL2.2(), LL2.3(), LL2.3u(), LL2.4(), LL2.5(),                        
        AR.2(), AR.3(),
        MM.2(), MM.3() 
#        baro5(), 
#        boltzmann(), 
#        CRS.4a(), CRS.4b(), CRS.4c(),  # cedergreen(alpha = 1)
#        CRS.6(), expDecay(), gompertzd(), 
#        L.3(), L.4(), L.5(),
#        richards(), 
#        UCRS.4a(), UCRS.4b(), UCRS.4c(), UCRS.5a(), UCRS.5b(), UCRS.5c(),  # ucedergreen(alpha = 1), 
        )
    } else {
        fctList <- flist
    }

#    grepFct <- function(x){grep(x, "Weibull")}
    textVec <- NULL
    lapFct <- function(x) {c(x$"name", x$"text")} 
    if (!is.null(fname))
    {
        textVec <- fname       
        sapFct <- function(x, object){grep(x, object$"name", fixed = TRUE)}
    }
#    if (!is.null(ftext))
#    {
#        textVec <- ftext
#        lapFct <- function(x) {x$"text"}        
#        sapFct <- function(x, object){grep(x, object$"text", fixed = TRUE)}
#    }
    

    displayFunction <- function(object)
    {
#        && (is.null(ftext) || (length(grep(ftext, object$"text")) > 0)) )    
        if ( ((is.na(noParm)) || (noParm == object$"noParm")) 
        && (is.null(textVec) || (any(textVec %in% object$"name"))) )
#        && (is.null(textVec) || (sum(unlist(sapply(textVec, sapFct, object = object))) > 0)) )
        {
            if (display)
            {
                cat(object$"text", "\n")
                cat(paste("(", object$"noParm", " parameters)", sep = ""), "\n")
#                cat("Equation: ", object$"equation", "\n")  # duplicate with help page ... skip
#                cat("Reference:", object$"reference", "\n")  # duplicate with help page ... skip
                cat("In 'drc': ", object$"name", "\n\n")
            }
            return(object)
        } else {
            return(NULL)
        }
    }
    lapList1 <- lapply(fctList, displayFunction)
    
    ## Removing NULL elements in the list 'lapList1'
    lapList2 <- list()
    counter <- 1
    for (i in 1:length(lapList1))
    {
        if (!is.null(lapList1[[i]]))
        {
            lapList2[[counter]] <- lapList1[[i]]
            counter <- counter + 1
        }
    }
    
    if (!is.null(textVec) || !is.na(noParm))
    {
#        invisible(lapply(fctList, displayFunction))
        invisible(lapList2)
    } else {
        invisible(lapply(lapList2, lapFct))
    }
}


#"getDatasets" <- function()
#{
#    data(package = "drc")
#}
