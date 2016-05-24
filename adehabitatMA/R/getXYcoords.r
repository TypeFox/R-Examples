".getXYcoords" <- function(w)
{
    ## Verifications
    if ((((!inherits(w, "asc"))&
          (!inherits(w, "kasc")))&
         (!inherits(w,"sahrlocs")))&
        (!inherits(w,"mapattr")))
        stop("non convenient object")

    ## Gets the attributes
    cs<-attr(w, "cellsize")
    xll<-attr(w, "xll")
    yll<-attr(w, "yll")

    ## Computation of the number of rows and columns of the matrix
    if (inherits(w,"asc")) {
        nr<-nrow(w)
        nc<-ncol(w)
    }
    if (((inherits(w,"kasc"))|(inherits(w, "sahrlocs")))|
        (inherits(w, "mapattr"))){
        nc<-attr(w, "nrow")
        nr<-attr(w, "ncol")
    }
    ## The results
    x<-xll+c(0:(nr-1))*cs
    y<-yll+c(0:(nc-1))*cs
    return(list(x=x, y=y))
}

