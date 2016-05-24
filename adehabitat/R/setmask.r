"setmask" <- function(x, mask)
{
    ## Verifications
    if ((!inherits(x, "asc"))&(!inherits(x, "kasc")))
        stop("x should be an object of class \"asc\" or \"kasc\"")
    if (!inherits(mask, "asc"))
        stop("mask should be of class \"asc\"")
    if (attr(x, "xll")!=attr(mask, "xll"))
        stop("Objects should have the same xll attribute")
    if (attr(x, "yll")!=attr(mask, "yll"))
        stop("Objects should have the same yll attribute")
    if (attr(x, "cellsize")!=attr(mask, "cellsize"))
        stop("Objects should have the same cellsize attribute")

    if (inherits(x, "kasc")) {

        ## If x is of class kasc
        ## Additional verifications
        if (attr(x, "nrow")!=ncol(mask))
            stop("Maps should have the same number of columns")
        if (attr(x, "ncol")!=nrow(mask))
            stop("Maps should have the same number of rows")

        ## Sets the mask thanks to managNAkasc
        x$mask0012<-as.vector(mask)
        so<-managNAkasc(x)
        so$mask0012<-NULL
        sorties<-so

    } else {

        ## If x is of class asc
        ## Additional verifications
        if (nrow(x)!=nrow(mask))
            stop("Objects should have the same attributes")
        if (ncol(x)!=ncol(mask))
            stop("Objects should have the same attributes")

        ## Sets the mask thanks to managNAkasc
        u<-as.kasc(list(x=x, mas=mask))
        so<-managNAkasc(u)
        sorties<-getkasc(so, "x")
    }
    return(sorties)
}

