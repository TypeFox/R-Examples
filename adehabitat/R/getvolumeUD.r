"getvolumeUD" <- function(x)
{
    ## Verifications
    if ((!inherits(x, "khrud"))&(!inherits(x, "kbbhrud")))
        stop("x should be an object of class \"khrud\" or \"kbbhrud\"")

    ## for each animal
    for (i in 1:length(x)) {

        ## gets the UD of the animal
        asc<-x[[i]]$UD
        cs<-attr(asc,"cellsize")

        ## computes the volume for each pixel
        ## thanks to a call to the C function calcvolume
        v<-.C("calcvolume", as.double(t(asc)), as.integer(ncol(asc)),
              as.integer(nrow(asc)), as.double(cs), PACKAGE="adehabitat")[[1]]

        ## standardize it so that the total volume is 1 over the area
        index<-1:length(v)
        vord<-v[order(v, decreasing=TRUE)]
        indord<-index[order(v, decreasing=TRUE)]
        vsu<-cumsum(vord)
        vreord<-vsu[order(indord)]*100

        ## output
        u<-matrix(vreord, ncol=ncol(asc), byrow=TRUE)
        UD <-getascattr(asc,u)
        attr(UD,"UD") <- "volume"
        x[[i]]$UD <- UD
    }
    ## OUTPUT
    class(x)<-c("khrvol", "khr")
    return(x)
}

