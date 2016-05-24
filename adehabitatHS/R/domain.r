domain <- function(x, pts, type=c("value", "potential"),
                     thresh=0.95)
{
    ## Verifications
    if (!inherits(x, "SpatialPixelsDataFrame"))
        stop("should be an object of class SpatialPixelsDataFrame")
    gridded(x) <- TRUE
    gr <- gridparameters(x)
    if (nrow(gr) > 2)
          stop("x should be defined in two dimensions")
    if ((gr[1, 2] - gr[2, 2])> get(".adeoptions", envir=.adehabitatMAEnv)$epsilon)
        stop("the cellsize should be the same in x and y directions")
    if (!inherits(pts, "SpatialPoints"))
        stop("should inherit from class \"SpatialPoints\"")

    typ <- c("numeric","factor")[as.numeric(sapply(slot(x, "data"),
                                                   is.factor))+1]
    if (!all(typ=="numeric"))
        stop("All variables in kasc should be of mode numeric")
    type<-match.arg(type)

    ## Preparation of the data to be passed to the C function "fctdomain"
    ## 1. spatial join of the points
    ptsmod <- as.matrix(join(pts, x))

    ## 2. deletes the missing values
    kascmod <- as.matrix(slot(x,"data"))
    if (any(is.na(kascmod)))
        stop("the same area should be provided for all variables")
    ptsmod <- ptsmod[!is.na(ptsmod[,1]),]

    ## 3. Computation of the range of environmental variables
    rg <- apply(kascmod, 2, function(x) diff(range(x)))

    ## Call to the C function
    toto<-.C("fctdomain", as.double(t(kascmod)), as.double(t(ptsmod)),
             as.double(rg), as.integer(nrow(ptsmod)),
             as.integer(nrow(kascmod)), as.integer(ncol(ptsmod)),
             double(nrow(kascmod)), PACKAGE="adehabitatHS")[[7]]

    ## Should the value or the potential habitat be exported in the output ?
    if (type!="value") {
        toto[toto<=thresh]<-NA
        toto[toto>thresh]<-1
    }

    toto <- data.frame(domain=toto)
    coordinates(toto) <- coordinates(x)
    gridded(toto) <- TRUE
    return(toto)
}
