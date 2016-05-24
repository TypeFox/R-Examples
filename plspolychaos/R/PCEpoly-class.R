#########################################################
# Class to store the result of functions 'polyLeg' and
# 'analyticsPolyLeg'
#########################################################
setClass("PCEpoly", slots = c(.Data = "matrix", STRUC = "PCEdesign", nvx = "numeric", 
    call = "call"))

#########################################################
# print method
# "all": option to fix the display. Vector. Valid values are:
# "FALSE" : the design only
# "TRUE": the design and the lhs
# " ..." : extend=TRUE/FALSE (option of print.PCEdesign
#          TRUE to display all the monomes)
#         and all options passed as it to print

print.PCEpoly <- function(x, all = FALSE, ...) {
    print.PCEdesign(x@STRUC, all, ...)
    cat("Number of rows:", nrow(x), "\n")
    
    if (all == TRUE) {
        cat("Created by:\n")
        print(x@call, ...)
    }
    return(invisible())
}  # end print.PCEpoly

#########################################################
## show method
show.PCEpoly <- function(object) {
    print.PCEpoly(object)
    return(invisible())
}  # end show.PCEpoly


# --------------------------------------
setMethod("show", signature(object = "PCEpoly"), definition = show.PCEpoly)

#########################################################
## getNames method
getNames.PCEpoly <- function(object) {
    
    slotnames <- slotNames(object)
    
    for (a in slotnames) {
        cat(" Slot: ", a, ".", sep = "")
        cde <- paste("class(object@", a, ")", sep = "")
        cat(" Class: \"", eval(parse(text = cde)), "\".", sep = "")
        cde <- paste("dim(object@", a, ")", sep = "")
        z <- eval(parse(text = cde))
        if (!is.null(z)) {
            cat(" Dimension: ", paste("(", paste(z, collapse = ", "), ")", sep = ""), 
                ".", sep = "")
        } else {
            cde <- paste("length(object@", a, ")", sep = "")
            z <- eval(parse(text = cde))
            if (!is.null(z)) {
                cat(" Length: ", paste("(", paste(z, collapse = ", "), ")", sep = ""), 
                  ".", sep = "")
            }
        }
        
        switch(a, .Data = {
            cat(" Legendre polynomial")
        }, STRUC = {
            cat(" The polynomial structure")
        }, nvx = {
            cat(" The number of inputs")
        }, call = {
            cat(" The command which creates the object")
        }, cat("Unknown slot ", a, "."))  # fin switch
        
        cat("\n")
        
        
    }  # fin a
    return(invisible())
}  # fin getNames


# --------------------------------------
setMethod("getNames", signature(object = "PCEpoly"), definition = getNames.PCEpoly) 
