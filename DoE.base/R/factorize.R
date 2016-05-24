## relies on the generic factorize and factorize.default from package conf.design

factorize.factor <- function (x, name = deparse(substitute(x)), extension = letters,
    drop = FALSE, sep = "", ...)
{
    ## adapted from package conf.design
    ## modifications:
    ##     default drop=FALSE
    ##     modified output: not factor or data frame of factors
    ##       but integer vector (in case of drop) or integer matrix
    ##     (function design did not work,
    ##      class design is not helpful)
    
    ### uses the default method of function factorize in conf.design
    llev <- conf.design::factorize(length(levels(x)))
    if (length(llev) == 1)
        return(if (drop) as.integer(x)-1 else {
            nm <- name
            x <- matrix(as.integer(x)-1,ncol=1)
            colnames(x) <- nm
            x
        })
    D <- NULL
    for (i in llev) {
        E <- D
        D <- NULL
        for (j in 1:i) D <- rbind(D, cbind(E, j))
    }
    l <- matrix(NA, nrow=length(x), ncol=ncol(D))
    for (i in seq(along = llev)) l[,i] <- D[, i][x] - 1
    colnames(l) <- paste(name, extension[1:length(llev)], sep = sep)
    l
}
factorize.design <- function (x, extension = letters, sep = ".", long=FALSE, ...)
{
    if (!"design" %in% class(x)) stop("applicable to class design objects only")
    if (is.null(design.info(x)))
            stop("applicable to uncorrupted class design objects from the class defined in package DoE.base only")
    di <- design.info(x)
    if (!di$type %in% c("full factorial","oa")) stop("x must be an unblocked full factorial or a design created with function oa.design")

    namen <- names(factor.names(x))
    aus <- numeric(0)
    for (nm in namen) {
       aus <- cbind(aus, factorize.factor(as.factor(x[[nm]]),name=nm, extension=extension, sep=sep, ...))
       }
    if (long) aus else apply(aus,2,max)+1
}

factorize.data.frame <- function (x, extension = letters, sep = ".", long=FALSE, ...)
{
    namen <- colnames(x)
    aus <- numeric(0)
    for (nm in namen) {
       aus <- cbind(aus, factorize.factor(as.factor(x[[nm]]),name=nm, extension=extension, sep=sep, ...))
       }
    if (long) aus else apply(aus,2,max)+1
}