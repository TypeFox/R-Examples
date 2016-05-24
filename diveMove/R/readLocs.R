
"readLocs" <- function(locations, loc.idCol, idCol, dateCol, timeCol=NULL,
                       dtformat="%m/%d/%Y %H:%M:%S", tz="GMT", classCol,
                       lonCol, latCol, alt.lonCol=NULL, alt.latCol=NULL, ...)
{
    ## Value: A data frame with ARGOS locations.
    ## --------------------------------------------------------------------
    ## Arguments: locations=quoted file name, including path, of file to
    ## read, or data.frame with data to read, or a text-mode connection,
    ## loc.idCol=column number containing the location id, idCol=column
    ## number identifying locations belonging to different groups,
    ## dateCol=column number containing dates and, optionally, times,
    ## timeCol=optional column number containing times, latCol and
    ## lonCol=latitude and longitude column numbers, respectively,
    ## alt.latCol and alt.lonCol=alternative latitude and longitude
    ## columns, respectively, classCol=ARGOS classification; ...= passed to
    ## read.csv()
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    if (inherits(locations, "connection") ||
        (is.character(locations) && file.exists(locations))) {
        srcfile.name <- ifelse(inherits(locations, "connection"),
                               basename(summary(locations)$description),
                               basename(locations))
        inLocs <- read.csv(locations, ...)
    } else {
        if (! is.data.frame(locations)) {
            stop ("'locations' must be a data.frame, path to a file, or a connection")
        } else {inLocs <- locations}
    }
    if (missing(loc.idCol)) {
        loc.id <- seq(nrow(inLocs))
    } else loc.id <- inLocs[, loc.idCol]
    if (missing(idCol)) {
        id <- rep(1, nrow(inLocs))
    } else id <- inLocs[, idCol]
    dtpasted <- paste(inLocs[, dateCol], inLocs[, timeCol])
    datetime <- as.POSIXct(strptime(dtpasted, format=dtformat), tz=tz)
    locs <- data.frame(loc.id=loc.id, id=id, time=datetime,
                       lon=inLocs[, lonCol], lat=inLocs[, latCol],
                       class=inLocs[, classCol])
    if (!is.null(alt.lonCol)) locs$alt.lon <- inLocs[, alt.lonCol]
    if (!is.null(alt.latCol)) locs$alt.lat <- inLocs[, alt.latCol]
    comment(locs) <- ifelse(exists("srcfile.name"), srcfile.name,
                            paste(deparse(match.call()), collapse=""))

    locs[order(locs[, 2], locs[, 3]), ] # sort by seal id and time
}

## TEST ZONE --------------------------------------------------------------
