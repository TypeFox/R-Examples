"domain" <-
function(kasc, pts, type=c("value", "potential"),
                 thresh=0.95)
{
    ## Verifications
    if (!inherits(kasc, "kasc"))
        stop("should be an object of class \"kasc\"")
    if (ncol(pts)!=2)
        stop("pts should have 2 columns")
    typ<-""
    for (i in 1:length(kasc)) {
        if (is.factor(kasc[[i]])) {
            typ[i] <- "factor"
        }
        else {
            typ[i] <- "numeric"
        }
    }
    if (!all(typ=="numeric"))
        stop("All variables in kasc should be of mode numeric")
    type<-match.arg(type)

    ## Preparation of the data to be passed to the C function "fctdomain"
    ## 1. spatial join of the points
    ptsmod<-as.matrix(join.kasc(pts, kasc))

    ## 2. deletes the missing values
    kasct<-kasc2df(kasc)
    kascmod<-as.matrix(kasct$tab)
    if (any(is.na(kascmod)))
        stop("the same area should be provided for all variables")

    ## 3. Computation of the range of environmental variables
    rg<-apply(kascmod, 2, function(x) range(x)[2] - range(x)[1])

    ## Call to the C function
    toto<-.C("fctdomain", as.double(t(kascmod)), as.double(t(ptsmod)),
             as.double(rg), as.integer(nrow(ptsmod)),
             as.integer(nrow(kascmod)), as.integer(ncol(ptsmod)),
             double(nrow(kascmod)), PACKAGE="adehabitat")[[7]]

    ## Transfo of the output vector into a map (equivalent to df2kasc)
    N <- nrow(kasc)
    indw <- c(1:N)
    n1 <- length(toto)
    compl <- rep(NA, N - n1)
    output <- c(toto, compl)
    indcompl <- indw[is.na(match(indw, kasct$index))]
    indtot <- c(kasct$index, indcompl)
    output <- output[sort(indtot, index.return = TRUE)$ix]
    output<-matrix(output, attr(kasc,"ncol"))

    ## Should the value or the potential habitat be exported in the output ?
    if (type!="value") {
        output[output<=thresh]<-NA
        output[output>thresh]<-1
    }

    ## Output
    attr(output, "xll") <- attr(kasc, "xll")
    attr(output, "yll") <- attr(kasc, "yll")
    attr(output, "cellsize") <- attr(kasc, "cellsize")
    attr(output, "type") <- "numeric"
    class(output)<-"asc"
    return(output)
}

