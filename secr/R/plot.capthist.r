##############################################################################
## package 'secr'
## plot.capthist.R
## 2013-11-20
## 2015-10-11 type = 'sightings'
##############################################################################

plot.capthist <- function(x, rad = 5,
   hidetraps = FALSE, tracks = FALSE,
   title = TRUE, subtitle = TRUE,
   add = FALSE,
   varycol = TRUE, icolours = NULL, randcol = FALSE,
   lab1cap = FALSE, laboffset = 4,
   ncap = FALSE,
   splitocc = NULL, col2 = 'green',
   type = c("petal", "n.per.detector", "n.per.cluster", "sightings"),
   cappar = list(cex=1.3, pch=16, col='blue'),
   trkpar = list(col='blue', lwd=1),
   labpar = list(cex=0.7, col='black'),
   ...)

    # see also version in d:\single sample with stems=F, mst=F 2009 02 22

{

    ## recursive if list of capthist
    if (ms(x)) {
        sapply (x, plot.capthist,
            rad = rad, hidetraps = hidetraps, tracks = tracks,
            title = title, subtitle = subtitle, add = add, varycol = varycol, icolours =
            icolours, randcol = randcol, lab1cap = lab1cap, laboffset =
            laboffset, ncap = ncap, splitocc = splitocc, col2 = col2,
            type = type, cappar = cappar, trkpar = trkpar, labpar = labpar, ...)
    }
    else {

        plotprox <- function (x) {
            x <- abs(x)   # do not distinguish deads for now
            # some dupl points will be over plotted - could increase rad for
            # captures after first at a trap on a given day
            dx  <- rep((cos((1:nocc) * 2 * pi / nocc) * rad), ncol(x))[x>0]
            dy  <- rep((sin((1:nocc) * 2 * pi / nocc) * rad), ncol(x))[x>0]
            if (varycol) icol <<- icol+1  ## 2009 10 02
            par(trkpar)
            if (varycol) par(col=icol)    ## override 2009 10 02
            if (tracks) lines (traps$x[x]+dx, traps$y[x]-dy)
            par(cappar)
            if (varycol) par(col=icol)    ## override 2009 10 02
            points (traps$x[x]+dx, traps$y[x]-dy)
        }
        plotpolygoncapt <- function (xy) {
            # some dupl points will be over plotted - could increase rad for
            # captures after first at a trap on a given day
            if (varycol) icol <<- icol+1  ## 2009 10 02
            par(trkpar)
            if (varycol) par(col=icol)    ## override 2009 10 02
            if (tracks) lines (xy$x, xy$y)
            par(cappar)
            if (varycol) par(col=icol)    ## override 2009 10 02
            points (xy$x, xy$y)
        }
        plotcapt <- function (x) {
            x <- abs(x)   # do not distinguish deads for now
            dx <- (cos((1:nocc) * 2 * pi / nocc) * rad)
            dy <- (sin((1:nocc) * 2 * pi / nocc) * rad)
            occ2 <- (1:nocc) %in% splitocc
            x0 <- x>0
            x1 <- x0 & !(occ2)
            x2 <- x0 & occ2
            x[!x0] <- 1   # fool into keeping length..
            usesplit <- !is.null(splitocc) & (sum(x2)>0)
            if (varycol) icol <<- icol+1
            par(trkpar)
            if (varycol) par(col=icol)   # override

            if (tracks) {
                if (usesplit) {
                    par(col=col2)
                    lines ((traps$x[x]+dx)[x0], (traps$y[x]-dy)[x0])  # all
                    par(trkpar)
                    if (varycol) par(col=icol)   # override
                    lines ((traps$x[x]+dx)[x1], (traps$y[x]-dy)[x1])  # pre-split
                }
                else
                lines ((traps$x[x]+dx)[x0], (traps$y[x]-dy)[x0])  # all
            }

            par(cappar)
            if (varycol) par(col=icol)   # override
            points ((traps$x[x]+dx)[x0], (traps$y[x]-dy)[x0])
            if (usesplit) {
              par(col=col2)
              points ((traps$x[x]+dx)[x2], (traps$y[x]-dy)[x2])
            }
        }
        labcapt <- function (n) {
            if ( detectr %in% c('proximity', 'count', 'polygonX',
                'transectX', 'signal', 'signalnoise', 'polygon',
                'transect', 'unmarked', 'presence') ) {
                warning ("labels not implemented for this detector type")
            }
            else {
                t1 <- abs(x[n,])
                o1 <- sum(cumsum(t1)==0)+1   # first occasion
                t1 <- t1[o1]                 # first trap site
                dx  <- (cos((1:nocc) * 2 * pi / nocc) * rad)[o1]
                dy  <- (sin((1:nocc) * 2 * pi / nocc) * rad)[o1]
                par(labpar)
                if (varycol) par(col=n)   # override
                laboffsety <- ifelse (length(laboffset)==2, laboffset[2], laboffset[1])
                text (traps$x[t1]+dx+laboffset[1], traps$y[t1]-dy+laboffsety, row.names(x)[n])
            }
        }
        labhead <- function (n, df) {
                par(labpar)
                if (varycol) par(col=n)   # override
                laboffsety <- ifelse (length(laboffset)==2, laboffset[2], laboffset[1])
                text (head(df[[n]],1)$x++laboffset[1], head(df[[n]],1)$y+laboffsety,
                      row.names(x)[n])
        }
        ncapt <- function (x) {
            if (detectr %in% .localstuff$detectors3D){
               temp <- t(apply (abs(x),c(2,3),sum)) # capts/trap/day)
            }
            else if (detectr %in% .localstuff$polydetectors){
               stop ("ncap does not work with polygon and similar detectors")
            }
            else {
               fx   <- factor(x, levels=0:nrow(traps))
               temp <- table (fx, col(x)) # capts/trap/day
               temp <- temp[-1,,drop=F] # drop zeros
            }
            dx  <- rep(cos((1:nocc) * 2 * pi / nocc) * rad, rep(nrow(traps),nocc))
            dy  <- rep(sin((1:nocc) * 2 * pi / nocc) * rad, rep(nrow(traps),nocc))
            par(labpar)
            par(adj=0.5)
            OK <- temp>0
            text ((traps$x[row(temp)]+dx)[OK], (traps$y[row(temp)]-dy)[OK], as.character(temp[OK]))
        }

        plotsignal <- function (df, minsignal, maxsignal,n) {
            # df is a dataframe for one animal
            # some dupl points will be over plotted - could increase rad for
            # captures after first at a trap on a given day

            ## function (signal, occasion, trap, minsignal, maxsignal, n)
            .localstuff$i <- .localstuff$i+1
            dx <- rep( (cos((.localstuff$i) * 2 * pi / n) * rad), nrow(df))
            dy <- rep( (sin((.localstuff$i) * 2 * pi / n) * rad), nrow(df))
            sq <- order(df$signal)     # plot darkest points last
            df <- df[sq,]
            df$trap <- as.character(df$trap)
            if (maxsignal>minsignal)
                greycol <- grey(0.7 * (1 - (df$signal-minsignal)/(maxsignal-minsignal)))
            else
                greycol <- grey(0.5)
            if (tracks) lines (traps$x[df$trap]+dx, traps$y[df$trap]-dy, col=greycol)
            par(cappar)
            points (traps[df$trap,'x']+dx, traps[df$trap,'y']-dy, col = greycol)
        }
        plotsightings <- function (x) {
            ## plot sightings; assumes proximity or count detectors
            Tu <- Tu(x)
            Tu0 <- Tu
            marking <- Tu
            if (is.null(Tu)) stop ("sightings type requires sighting data")
            markocc <- markocc(traps(x))
            Tu0[Tu!=0] <- NA
            Tu0[,markocc>0] <- NA
            Tu[Tu==0] <- NA
            marking[,markocc<1] <- NA
            dx  <- rep((cos((1:nocc) * 2 * pi / nocc) * rad), each = nrow(Tu))
            dy  <- rep((sin((1:nocc) * 2 * pi / nocc) * rad), each = nrow(Tu))
            par(cappar)
            dx0 <- dx; dx0[is.na(Tu0)] <- NA
            
            if (detector(traps(x)) %in% c('polygon')) {
                centres <- split(traps(x), polyID(traps(x)))
                ## assume each polygon closed, so first vertex redundant
                centres <- lapply(centres, function(xy) apply(xy[-1,,drop=FALSE], 2, mean))
                trapxy <- data.frame(do.call(rbind,centres))
                names(trapxy) <- c('x','y')
            }
            else trapxy <- traps(x)
            text (rep(trapxy$x, nocc) + dx, rep(trapxy$y, nocc) - dy, Tu, cex = 0.7)
            points (rep(trapxy$x, nocc) + dx0, rep(trapxy$y, nocc) - dy, pch = 1, cex = 0.7)
            ## marking occasions shown as dot
            dx[is.na(marking)] <- NA
            points (rep(trapxy$x, nocc) + dx, rep(trapxy$y, nocc) - dy, pch=16, cex=0.4)
        }
        
        ###########
        ## MAINLINE

        ## suggested by Mike Meredith 2013-05-24
        opal <- palette() ; on.exit(palette(opal))
        type <- match.arg(type)
        traps <- traps(x)
        detectr <- detector(traps)
        nocc <- ncol(x)
        nanimal <- nrow(x)

        if (type == 'petal')
            cappar <- replacedefaults (list(cex=1.3, pch=16, col='blue'), cappar)
        if (type == 'sightings')
            cappar <- replacedefaults (list(cex=1, pch=16, col='blue'), cappar)
        if (type %in% c('n.per.cluster','n.per.detector'))
            cappar <- replacedefaults (list(cex = 3, pch = 21), cappar)

        trkpar <- replacedefaults (list(col='blue', lwd=1), trkpar)
        labpar <- replacedefaults (list(cex=0.7, col='black'), labpar)
        initialpar <- par(cappar)

        if (!add) plot(traps, hidetr=hidetraps, ...)

        ## check added 2012-10-25
        if (nrow(x) == 0) {
            warning("no detections in capthist object")
            type <- 'null'
        }

        if (type == 'petal') {
            if (is.null(icolours)) icolours <- topo.colors((nanimal+1)*1.5)
            if (varycol) {
                if (randcol) icolours <- sample(icolours)
                test <- try (palette(icolours))  ## too many?
                if (inherits(test, 'try-error'))
                    stop ("requested too many colours; ",
                          "try with varycol = FALSE")
                icol <- 0
            }
            if ((nocc == 1) & ! (detectr %in% c('signal','signalnoise'))) rad <- 0

            if ( detectr %in% c('proximity', 'count', 'unmarked', 'presence') )
            {
                ## detector number
                w <- apply(x,1:2,function(x) (abs(x)>0) * (1:length(x)))
                w <- array(w, dim=dim(x)[c(3,1,2)])  ## 2013-11-20, in case dim dropped
                multiples <- any(apply(x,1:2, function(y) sum(abs(y)))>1)
                if (multiples & tracks)
                    warning("track for repeat detections on same occasion",
                            " joins points in arbitrary sequence")
                w <- aperm(w, c(2,3,1))
                apply( w, 1, plotprox )
            }
            else
            if ( detectr %in% .localstuff$polydetectors ) {
                ## occasions not distinguished
                lxy <- split (xy(x), animalID(x, names=FALSE))
                lapply (lxy, plotpolygoncapt)
            }
            else
            if ( detectr %in% c('signal','signalnoise') )
            {
                .localstuff$i <- 0
                temp <- data.frame(ID = animalID(x), occ = occasion(x),
                                   trap=trap(x), signal = signal(x))
                lsignal <- split(temp, animalID(x, names = FALSE))
                lapply(lsignal, plotsignal, minsignal = min(temp$signal),
                    maxsignal = max(temp$signal), n=nanimal)
            }
            else  {   ## single, multi-catch traps
                apply( x, 1, plotcapt )
            }

            if (lab1cap) {
                if ( detectr %in% .localstuff$polydetectors ) {
                    lxy <- split (xy(x), animalID(x, names = FALSE))
                    sapply(1:nanimal, labhead, df=lxy)
                }
                else sapply(1:nanimal, labcapt)
            }

            if (ncap) { ncapt(x)}

        }
        else if (type %in% c('n.per.cluster','n.per.detector')) {
            if (type == 'n.per.detector') {
                ## never yields zeros
                temp <- table(trap(x), animalID(x))>0
                nj <- apply(temp,1,sum)
                centres <- traps(x)[names(nj),]
            }
            else if (type == 'n.per.cluster') {
                nj <- cluster.counts(x)
                centres <- cluster.centres(traps)
                # hide zeros, if present
                centres <- centres[nj>0,]
                nj <- nj[nj>0]
            }
            else stop ("unrecognised type")

            if (is.null(icolours)) {
                icolours <- topo.colors(max(nj)*1.5)
            }
            palette (icolours)
            npal <- length(icolours)
            if (max(nj) < npal) {
                cols <- npal-nj
            }
            else {
                cols <- round(npal * (1-nj/max(nj)))
            }
            if (cappar$pch == 21)
                fg <- 'black'
            else
                fg <- cols
            par(cappar)
            points(centres, col = fg, bg = cols, pch = cappar$pch, cex = cappar$cex)

            if (ncap) {
                par(labpar)
                par(adj=0.5)
                vadj <- diff(par()$usr[3:4])/500  ## better centring!
                text(centres$x, centres$y + vadj, nj)
            }

            ## should export data for legend:
            ## count classes 0, 1-, 2-,...,-max
            ## and corresponding colours
            tempcol <- npal- (1:max(nj))
            output <- data.frame(
                legend = 1:max(nj),
                col = tempcol,
                colour = icolours[tempcol],
                stringsAsFactors = FALSE
            )
        }
        else if (type == 'sightings') {
            plotsightings(x)
        }
        else
            if (type != 'null') stop ("type not recognised")
        

        ####################################################
        ## Titles

        nd <- if(detectr %in% .localstuff$exclusivedetectors) 
            sum(abs(x)>0) else sum(abs(x))
        if (type == 'sightings') {
            markocc <- markocc(traps(x))
            Tu <- Tu(x)
            nd <- sum(Tu)
            nocc <- sum(markocc<1)
        }
        
        if (is.logical(title)) {
            txt <- ifelse (is.null(session(x)), paste(deparse(substitute(x)),
                       collapse=''), session(x))
            title <- ifelse(title, txt, '')
        }
        if (title != '') {
            par(col='black')
            mtext(side=3,line=1.2, text = title, cex=0.7)
        }
        if (is.logical(subtitle)) {

            subtitle <- if (subtitle) {
                if (type == 'sightings') {
                    if (any(markocc<0))
                        paste(nocc, 'sighting occasions,', nd, 'sightings of marked and unmarked animals')
                    else
                        paste(nocc, 'sighting occasions,', nd, 'sightings of unmarked animals')
                }
                else
                    paste(nocc, 'occasions,', nd, 'detections,',
                          nanimal, 'animals')
            }
                else ''
        }
        if (subtitle != '') {
            par(col='black')
            mtext(text = subtitle, side=3,line=0.2, cex=0.7)
        }
        ####################################################

        par(initialpar)   # restore
        if (type %in% c('n.per.detector','n.per.cluster'))
            invisible(output)
        else
            invisible(nd)
    }
}
############################################################################################

plotMCP <- function(x, add = FALSE, col = 'black', fill = NA, lab1cap = FALSE,
                    laboffset = 4, ncap = FALSE, ...) {
    plotone <- function (df, col, fill) {
        par(fg=col)
        polygon(df, col=fill)
    }
    mcp <- function (df) {
        df <- df[chull(df[,1], df[,2]),]
        rbind(df, df[1,])
    }
    if (ms(x)) {
        lapply(x, plotMCP, add = add, col = col, ...)
    }
    else {
        xyl <- telemetryxy(x)
        if (is.null(xyl)) {
            if (detector(traps(x)) %in% c('polygon','polygonX','telemetry'))
                xyl <- split(xy(x), animalID(x))
            else
                stop ("requires polygon or telemetry data")
        }
        if (!add)
            plot(traps(x), ...)
        fg <- par()$fg
        on.exit(par(fg=fg))
        if (missing(col)) col <- 'black'
        xymcp <- lapply(xyl, mcp)
        mapply(plotone, xymcp, col, fill)

        if (lab1cap | ncap) {
            labhead <- function(xy, nam, num, col) {
                laboffsety <- ifelse (length(laboffset)==2, laboffset[2], laboffset[1])
                if (ncap)
                    text (xy[1,1] + laboffset[1], xy[1,2] + laboffsety, num, col=col)
                else
                    text (xy[1,1] + laboffset[1], xy[1,2] + laboffsety, nam, col=col)
            }
            mapply(labhead, xymcp, names(xyl), sapply(xyl,nrow), col)
        }
        invisible(xyl)
    }
}
############################################################################################
# plotMCP(AB2004,col=1:12, gridl=F)

