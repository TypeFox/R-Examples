
getKMLcoordinates <- function(kmlfile, ignoreAltitude = FALSE) {
    if (missing(kmlfile)) 
        stop("kmlfile is missing")
    kml <- paste(readLines(kmlfile, encoding = "UTF-8"), 
        collapse = " ")
    ## ++ new code Mike Sumner 120509
    ## remove tabs first
    kml <- gsub("[[:blank:]]+", " ", kml)
    ##   

    re <- "<coordinates> *([^<]+?) *<\\/coordinates>"
    mtchs <- gregexpr(re, kml)[[1]]
    coords <- list()
    for (i in 1:(length(mtchs))) {
        kmlCoords <- unlist(strsplit(gsub(re, "\\1", substr(kml, 
# Kent Johnson bugfix (added -1) 151013
            mtchs[i], (mtchs[i] + (attr(mtchs, "match.length")[i]-1))), 
            perl = TRUE), split = " "))
        m <- t(as.matrix(sapply(kmlCoords, function(x) as.numeric(unlist(
            strsplit(x, ","))), USE.NAMES = FALSE)))
        if (!ignoreAltitude && dim(m)[2] != 3) 
            message(paste("no altitude values for KML object", i))
        coords <- append(coords, 
            ifelse(ignoreAltitude, list(m[, 1:2]), list(m)))
    }
    coords
}
