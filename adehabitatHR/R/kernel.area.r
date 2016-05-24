"kernel.area" <- function(x, percent = seq(20, 95, by = 5),
                          unin = c("m", "km"),
                          unout = c("ha", "km2", "m2"), standardize=FALSE)
{
    if (!inherits(x, "estUDm")&(!inherits(x, "estUD")))
        stop("x should be of class \"estUD\" or \"estUDm\"")
    unin <- match.arg(unin)
    unout <- match.arg(unout)

    if (inherits(x, "estUD")) {
        if (!slot(x, "vol"))
            x <- getvolumeUD(x, standardize=standardize)
        tmp <- as.data.frame(x)[,1]
        ## Check that all the contour are within the study area limits
        gp <- gridparameters(x)[,3]
        tmpm <- matrix(tmp, ncol=gp[2], nrow=gp[1], byrow = TRUE)
        ma <- min(c(tmpm[c(1:nrow(tmpm)), c(1,ncol(tmpm))],
                    tmpm[c(1,nrow(tmpm)), c(1:ncol(tmpm))]))
        if (any(percent>=ma))
            warning(paste("The grid is too small to allow the estimation of home-range\nfor the following value of percent: ", paste(percent[percent>=ma], collapse = ","), ". You should rerun kernelUD with a larger extent parameter", sep=""))

        area <- unlist(lapply(percent, function(o) length(tmp[tmp<=o])*(gridparameters(x)[1,2]^2)))
        if (unin == "m") {
            if (unout == "ha")
                area <- area/10000
            if (unout == "km2")
                area <- area/1e+06
        }
        if (unin == "km") {
            if (unout == "ha")
                area <- area * 100
            if (unout == "m2")
                area <- area * 1e+06
        }
        names(area) <- as.character(percent)
        return(area)
    } else {
        area <- do.call("data.frame", lapply(x, function(j) {
            kernel.area(j, percent, unin, unout)
        }))

        class(area) <- c("hrsize", "data.frame")
        return(area)

    }
}

