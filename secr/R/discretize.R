
## 2015-11-26 adapted to also handle traps objects
## 2015-12-06 added transect functionality
## 2016-01-07 cell.overlap, type

discretize <- function (object, spacing = 5, outputdetector = c('proximity','count','multi'),
                        tol = 0.001, cell.overlap = FALSE, type = c('centre','any','all'), ...) {
    ## convert capthist data from polygon detectors to point detector
    outputdetector <- match.arg(outputdetector)
    type <- match.arg(type)
    if (ms(object)) {
        CHlist <- lapply(object, discretize, spacing, outputdetector, tol, cell.overlap, type, ...)
        class(CHlist) <- class(object)
        CHlist
    }
    else {

        onepolytraps <- function(trapsCH) {
            trps <- make.mask(trapsCH, buffer = spacing*2, spacing = spacing)
            ## arrange shared grid
            mx <- min(trapsCH[,1]); dx <- mx - trunc(mx/spacing)*spacing
            my <- min(trapsCH[,2]); dy <- my - trunc(my/spacing)*spacing
            trps$x <- trps$x - dx
            trps$y <- trps$y - dy
            if (type == 'centre') {
                poly <- inflate(trapsCH, 1+tol)  ## inflate is fn in utility.r
                OK <- pointsInPolygon(trps,poly)
            }
            else {
                poly <- inflate(trapsCH, 1+tol)  ## inflate is fn in utility.r
                cell <- matrix(c(-1,-1,1,1,-1,-1,1,1,-1,-1), ncol = 2) * spacing/2
                trpsc <- rbind(sweep(trps, STATS=c(-1,-1)*spacing/2, MARGIN=2, FUN='+'),
                               sweep(trps, STATS=c(-1,1)*spacing/2,  MARGIN=2, FUN='+'),
                               sweep(trps, STATS=c(1,1)*spacing/2,   MARGIN=2, FUN='+'),
                               sweep(trps, STATS=c(1,-1)*spacing/2,  MARGIN=2, FUN='+')
                               )
                fn <- get(type)  ## 'all' and 'any' are function names
                OK <- apply(matrix(pointsInPolygon(trpsc,poly), ncol = 4), 1, fn)
            }
            trps <- subset(trps, OK)
            temptraps <- read.traps(data = trps, detector = outputdetector, spacing = spacing)
            rownames(temptraps) <- 1:nrow(temptraps)
            if (cell.overlap) {
                ## cell overlap with polygon
                spp <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(trapsCH))), ID=1)))
                cell <- matrix(c(-1,-1,1,1,-1,-1,1,1,-1,-1), ncol = 2) * spacing/2
                onecell <- function(xy) {
                    cell <- sweep(cell, STATS = xy, FUN = '+', MARGIN = 2)
                    cell <- SpatialPolygons(list(Polygons(list(Polygon(cell)), ID=1)))
                    if (!requireNamespace('rgeos', quietly = TRUE))
                        stop ("package rgeos is required for area of overlap")
                    rgeos::gArea(rgeos::gIntersection(cell,spp))
                }
                overlap <- apply(trps,1,onecell)/spacing^2
            }
            else overlap <- rep(1, nrow(trps))

            if (!is.null(usage(trapsCH))) {
                usage(temptraps) <- matrix (usage(trapsCH), byrow = TRUE,
                                            nrow = nrow(trps), ncol = ncol(object)) * overlap
            }
            else {
                usage(temptraps) <- matrix (overlap, nrow = nrow(trps), ncol = ncol(object))
            }
            if (!is.null(covariates(trapsCH))) {
                covdf <- as.data.frame(covariates(trapsCH)[rep(1,nrow(temptraps)),])
                rownames(covdf) <- rownames(temptraps)
                covariates(temptraps) <- covdf
            }
            else
                covariates(temptraps) <- data.frame(polyID=rep(NA,nrow(temptraps)))

            if (!is.null(markocc(trapsCH))) {
                markocc(temptraps) <- markocc(trapsCH)
            }
            temptraps
        }

        trps <- if (inherits(object, 'traps')) object else traps(object)

        ## first make combined traps object
        if (!(detector(trps) %in% c('polygon', 'polygonX','transect', 'transectX')))
            stop ("discretize is for polygon or transect data")

        if (detector(trps) %in% c('transect', 'transectX')) {
            ## for transects, snip and reduce do all the work
            temp <- snip(object, by = spacing, ...)
            if (inherits(object, 'capthist'))
                temp <- reduce(temp, outputdetector = outputdetector)
            return(temp)
        }
        else {

            polylevels <- levels(polyID(trps))
            polytrps <- split(trps, polylevels, bytrap = TRUE)
            tmp <- lapply (polytrps, onepolytraps)
            trps <- if (length(tmp)==1) tmp[[1]] else do.call(rbind, tmp) ## but loses numbering
            covariates(trps)$polyID <- rep(polylevels, sapply(tmp, nrow))

            if ((inherits(object, 'traps')))
                return(trps)
            else {

                ## now assemble captures dataframe
                trpnum <- nearesttrap(xy(object), trps)

                df <- data.frame(
                    trap = factor(trpnum, levels =1:nrow(trps)),
                    occ = factor(occasion(object), levels = 1:ncol(object)),
                    ID = factor(animalID(object, names = F)),
                    alive = alive(object))

                if (outputdetector %in% c('multi')) {
                    alivesign <- df$alive*2 - 1
                    tempnew <- matrix(0, nrow = nrow(object), ncol = ncol(object))
                    tempnew[cbind(df$ID, df$occ)] <- trpnum * alivesign
                    if (detector(traps(object)) %in% 'polygon')
                        warning("polygon data converted to multi-catch; ",
                                "information may be lost", call. = FALSE)
                }
                else {
                    tempnew <- table(df$ID, df$occ, df$trap)
                    alivesign <- tapply(df$alive, list(df$ID,df$occ,df$trap),all)
                    alivesign[is.na(alivesign)] <- TRUE
                    alivesign <- alivesign * 2 - 1
                    if (! (outputdetector %in% c('count'))
                        && (length(tempnew)>0)) {
                        ## convert 'proximity' to binary
                        tempnew[tempnew>0] <- 1
                        warning("count data converted to binary; ",
                                "information may be lost", call. = FALSE)
                    }
                    tempnew <- tempnew * alivesign
                }

                ## restore attributes
                rownames(tempnew) <- rownames(object)
                class(tempnew) <- 'capthist'
                session(tempnew) <- session(object)
                attr(tempnew, 'n.mash') <- attr(object, 'n.mash')
                attr(tempnew, 'centres') <- attr(object, 'centres')
                if (!is.null(covariates(object)))
                    if (nrow(covariates(object)) == nrow(tempnew))
                        covariates(tempnew) <- covariates(object)
                traps(tempnew) <- trps

                ## unmarked/nonID sightings cannot be differentiated by detector: save total
                if (!is.null(Tu(object)))
                    Tu(tempnew) <- sum(Tu(object))
                if (!is.null(Tm(object)))
                    Tm(tempnew) <- sum(Tm(object))
                return(tempnew)
            }
        }
    }
}
