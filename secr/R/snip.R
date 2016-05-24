###############################################################################
## package 'secr'
## snip.R
## 2012-12-10,11,12,13 slice transects into small discrete units
## uses xyontransect from verify.r
## transectX, transect
###############################################################################

#------------------------------------------------------------------------
# divide a transect into many equal-sized parts, each 'by' metres long
slice <- function (object, from = 0, by = 1000, length.out = NULL, keep.incomplete = TRUE) {

    ## case 1 : from = 0, keep.incomplete = TRUE
    ## case 2 : from = 0, keep.incomplete = FALSE
    ## case 3 : from > 0, keep.incomplete = TRUE
    ## case 4 : from > 0, keep.incomplete = FALSE

    sliceone <- function (df) {
        d <- sqrt(diff(df$x)^2 + diff(df$y)^2)
        d <- c(0,cumsum(d))           ## distance of vertices along old
        cumd <- d[length(d)]          ## total length

        case <- integer(0)
        if ((keep.incomplete) & (from==0))
            case <- 1
        else if ((!keep.incomplete) & (from==0))
            case <- 2
        else if ((keep.incomplete) & (from>0))
            case <- 3
        else if ((!keep.incomplete) & (from>0))
            case <- 4
        else stop ("unrecognised case")

        if (!is.null(length.out)) {
            if (length.out > 0)
                by <- (cumd - from) / length.out
        }

        ## exit early if no computation required
        if ((by > cumd) & !keep.incomplete) {
            warning ("zero-length transect", call. = FALSE)
            df <- df[1,, drop=FALSE]
            return(df)
        }

        ## new transect number corresp each old vertex
        seg <- floor((d-from)/by)+1
        if (case == 3) seg <- seg + 1
        # rownames(df) <- paste(seg, sprintf("%04d",as.numeric(rownames(df))), sep='.')
        digits <- trunc(log10(nrow(df)) + 1)
        fmt <- paste("%0", digits, "d", sep = "") 
        rownames(df) <- paste(seg, sprintf(fmt,1:nrow(df)), sep='.')
    
        ## proceed by defining intermediate points between new segments
        breaks <- seq(from,cumd,by)
        if (from == 0)
            breaks <- breaks[-1]
        if (keep.incomplete)
            to <- cumd
        else
            to <- breaks[length(breaks)]

        if (length(breaks) > 0) {

            if ((from>0) & (!keep.incomplete))
                warning ("initial ", as.character(round(from)),
                         " metres discarded", call.=FALSE)
            if ((breaks[length(breaks)] < cumd) & (!keep.incomplete))
                warning ("final ", as.character(round(cumd - breaks[length(breaks)])),
                         " metres discarded", call.=FALSE)

            ## for each break determine the preceding old vertex
            breakstartold <- sapply(breaks, function(z) which.max(d>=z)) - 1
            breakstartold[breakstartold==0] <- 1  ## first is special case

            ## interpolate break points along transect between old vertices
            oldseglength <- d[breakstartold+1] - d[breakstartold]
            dseg <- breaks - d[breakstartold]  ## distance of termini within segment
            f <- dseg / oldseglength
            dfbreaks <- data.frame(
                x = df$x[breakstartold] + f * (df$x[breakstartold + 1] - df$x[breakstartold]),
                y = df$y[breakstartold] + f * (df$y[breakstartold + 1] - df$y[breakstartold]))
            k <- nrow(dfbreaks)

            ##
            if (k > 0) {
                ## potential start and end vertices for new transects
                indices <- c(1:k, 1:k)
                tempname <- c((1:k) - 0.0001, (1:k) + 0.0001)
                ## drop start or end vertices for terminal tansects as needed
                drop <- switch (case, numeric(0), 2*k, numeric(0), c(1,2*k))
                if (length(drop)>0) {
                    indices <- indices[-drop]
                    tempname <- tempname[-drop]
                }
                dfbreaks <- dfbreaks[indices,]
                if (case %in% c(1,2,3)) tempname <- tempname + 1
                rownames(dfbreaks) <- tempname
                if (!keep.incomplete) {
                    df <- df[(d >= from) & (d <= to),]
                }
                df <- rbind(df,dfbreaks)
            }
            else if (k==0) warning ("snip of zero length")
        }
        df <- df[order(as.numeric(rownames(df))),]
        df
    }
    IDfn <- function(x) {
        if (length(x) > 1)
            paste(x[1], sprintf("%04d", as.numeric(x[2])), sep='.')
        else
            x
    }
    ############
    ## main line
    if (!detector(object) %in% c('transect','transectX'))
        stop ("requires 'transect' input")
    temp <- split(as.data.frame(object), transectID(object))
    temp <- lapply(temp, sliceone)
    temp <- do.call(rbind, temp)
    oldID <- as.numeric(do.call(rbind, strsplit(rownames(temp),'.', fixed=T))[,1])
    ID <- sapply(strsplit(rownames(temp),'.', fixed=T), IDfn)
    temp <- split(temp, ID)
    newobj <- make.transect(transectlist = temp, exclusive =
                            detector(object) == 'transectX')
    if (!is.null(usage(object))) {
## TO BE FIXED 2012-12-22
        usagematrix <- matrix (0, ndetector(newobj), ncol(usage(object)))
        usagematrix <- matrix(usagematrix, nrow = nrow(newobj))
        usage(newobj) <- usagematrix[oldID,]
    }
    if (!is.null(covariates(object))) {
        covariates(newobj) <- covariates(object)[oldID,]
    }
    newobj
}
#------------------------------------------------------------------------

## for 'object' either traps or capthist:
## modelled in part on reduce.capthist
snip <- function (object, from = 0, by = 1000, length.out = NULL, keep.incomplete = TRUE) {

    if (ms(object)) {
        ## for each component session
        temp <- lapply (object, snip, from = from, by = by, length.out = length.out,
                        keep.incomplete = keep.incomplete)
        if (inherits(object,'capthist'))
            class(temp) <- c('list', 'capthist')
        else
            class(temp) <- c('list', 'traps')
        return(temp)
    }
    else {

        if (inherits(object, 'traps'))
            return (slice(object,  from = from, by = by, length.out = length.out,
                          keep.incomplete = keep.incomplete))
        else if (inherits(object, 'capthist')) {
            newtraps <- slice(traps(object), from = from, by = by, length.out =
                              length.out, keep.incomplete = keep.incomplete)
            newtrap <- xyontransect(xy(object), newtraps)
            newtrap <- factor(newtrap, levels = 1:length(levels(polyID(newtraps))))
            old.row.names <- row.names(object)
            df <- data.frame(
                             trap = trap(object, names = F),
                             occ = occasion(object),
                             ID = animalID(object, names = F),
                             alive = alive(object),
                             x = xy(object)[,1],
                             y = xy(object)[,2],
                             newtrap = newtrap)

            if (detector(traps(object)) == 'transect') {
                newobj <- table(df$ID, df$occ, df$newtrap)
                alivesign <- tapply(df$alive, list(df$ID,df$occ,df$newtrap),all)
                alivesign[is.na(alivesign)] <- TRUE
                alivesign <- alivesign * 2 - 1
                newobj <- newobj * alivesign
            }
            else if (detector(traps(object)) == 'transectX') {
                alivesign <- df$alive*2 - 1
                alivesign[is.na(alivesign)] <- TRUE
                newobj <- matrix(0, nrow = nrow(object), ncol = ncol(object))
                newobj[cbind(df$ID, df$occ)] <- as.numeric(df$newtrap) * alivesign
            }
            else stop ("unsuitable detector type")

            class(newobj) <- 'capthist'
            traps(newobj) <- newtraps
            rownames(newobj) <- old.row.names
            if (detector(traps(object)) == 'transectX')
                xy(newobj) <- xy(object)
            else {
                detectorder <- order(df$newtrap, df$occ, df$ID)
                xy(newobj) <- df[detectorder,c('x','y'), drop = FALSE]
            }
            covariates(newobj) <- covariates(object)   # OK because all animals transfer
            if (!keep.incomplete)
                newobj <- subset(newobj)   ## drops null capture histories if transects trimmed
            return(newobj)
        }
    }
}

#---------------------------------------------------------------------------------------

