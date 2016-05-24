read.coord <- function (inputfile = "", type = "geographical")
{
    syst <- pmatch(type, c("geographical", "cartesian"))
    if (is.na(syst))
        stop("Invalid coordinate system")
    if (length(syst) != 1)
        stop("Ambiguous coordinate system")
    if(inputfile == "") inputfile <- choose.files(filters = Filters["txt",], multi = FALSE)
    n <- nchar(inputfile)
    stopifnot(length(n) != 0)
    stopifnot(n != 0)
    if (substr(inputfile, n - 3, n) != ".txt") {
        cat("The input file is not .txt \n")
        return(invisible())
    }
    h <- scan(inputfile, sep = ",", what = character(0), nlines = 1)
    if (length(h) != 3) {
        cat(c("Inadequate format of input file \n", "Tip 1: Fields must be separated by comma \n",
            "Tip 2: Include the following header: ID, Longitude, Latitude \n",
            "Tip 3: The first field corresponds to species ID. The others fields are interchangeable \n",
            "Tip 4: Put coordinates in decimal format. \n"))
        return(invisible())
    }
    h <- tolower(h)
    lo <- grep("lo", h)
    la <- grep("la", h)
    if (length(lo) != 1 | length(la) != 1) {
        cat(" No recognizable field identity. Verify that they are well written.  \n")
        return(invisible())
    }
    else if (la == lo) {
        cat(" No recognizable field identity. Verify that they are well written.  \n")
        return(invisible())
    }
    if (any(c(lo, la) == 1)) {
        cat(" Remember that IDs should be included in the first field \n")
        return(invisible())
    }
    puntos <- scan(inputfile, sep = ",", what = list(character(0),
        double(0), double(0)), skip = 1)
    numpoints <- length(puntos[[1]])
    o <- order(puntos[[1]])
    puntos[[1]] <- puntos[[1]][o]
    puntos[[2]] <- puntos[[2]][o]
    puntos[[3]] <- puntos[[3]][o]
    if (any(is.na(puntos[[lo]])) | any(is.na(puntos[[la]]))) {
        cat(" There are some lines wihtout coordinates associated to \n")
        return(invisible())
    }
    if(syst == 1) {
      if (any(abs(puntos[[lo]])/180 > 1)) {
          cat(" Some longitude values fall outside the allowable range [-180,180] \n")
          return(invisible())
      }
      if (any(abs(puntos[[la]])/90 > 1)) {
          cat(" Some latitude values fall outside the allowable range [-90,90] \n")
          return(invisible())
      }
    }
    codigo <- unique(puntos[[1]])
    k <- tapply(puntos[[1]], puntos[[1]], length)
    rslts <- list(Numpoints = numpoints, Points = data.frame(cbind(IDsp = rep(1:length(codigo), k),
        Longitud = as.double(puntos[[lo]]), Latitud = as.double(puntos[[la]]))),
        Label = codigo, Type = ifelse(syst == 1, "geographical", "cartesian"))
    class(rslts) <- "dnpoint"
    rslts
}