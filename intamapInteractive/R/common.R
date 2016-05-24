


 spdf2list2 = function (data) {
# Modified function to convert a SpatialPolygonsDataFrame into a list of polygons
# This is a new version, also checking if there are several polygons in each Polygons
# Original function taken from GeoXp
    if (class(data)[1] != "SpatialPolygonsDataFrame") 
        stop("Must be a SpatialPolygonsDataFrame object")
    poly <- data@polygons
    n <- length(poly)
    for (i in 1:n) {
      np = length(poly[[i]]@Polygons)
      for (j in 1:np) {
        if (i ==1 & j ==1) {
          X <- poly[1][[1]]@Polygons[[1]]@labpt[1]
          Y <- poly[1][[1]]@Polygons[[1]]@labpt[2]
          contours <- rbind(NA, NA, NA, poly[1][[1]]@Polygons[[1]]@coords)
        } else {
          X <- rbind(X, poly[i][[1]]@Polygons[[j]]@labpt[1])
          Y <- rbind(Y, poly[i][[1]]@Polygons[[j]]@labpt[2])
          contours = rbind(contours, NA, NA, NA, poly[i][[1]]@Polygons[[j]]@coords)
        }
      }
    }
    contours = rbind(contours, NA, NA, NA)
    return(list(X = X, Y = Y, poly = contours))
}




dSolve = function(diffs) {
  D = as.data.frame(diffs$D)
  xnam <- paste("x", 1:dim(D)[2], sep="")
  names(D) = xnam
  fmla <- as.formula(paste("Q ~ ", paste(xnam, collapse= "+"),"-1"))
  Q = diffs$Q
  V = diffs$V
  dat = cbind(D,Q=Q)
  ols = as.data.frame(summary(lm(fmla,dat))$coefficients)
  wls = as.data.frame(summary(lm(fmla,dat,weights = 1/V))$coefficients)
  names(ols) = c("ols","ols.std","ols t value", "ols.p")
  names(wls) = c("wls","wls.std","wls t value", "wls.p")
  return(cbind(ols,wls))
}

#EJP: single comment sign added to all lines in this function:
#(double comment means this was already out-commented

#zip.file.extract.dir <- function (file, tmpdir, zipname = "R.zip", unzip = getOption("unzip"),zipdir)
#{
#    if (!is.character(unzip) || length(unzip) != 1)
#        stop("'unzip' must be a single character string")
#    if (!nzchar(unzip))
#        unzip <- "internal"
#    path <- dirname(file)
#    topic <- basename(file)
#    if (file.exists(file.path(path, zipname))) {
##        if (missing(zipdir)) {
##          tmpd <- tempdir()
##        } else {
##          tmpd <- zipdir
##        }
#        tmpd = ifelse(missing(zipdir),tempdir(),zipdir)
#        if (unzip != "internal") {
#            cmd <- paste(unzip, "-oq", shQuote(file.path(path,
#                zipname)), topic, " -d ", tmpd)
#            res <- if (.Platform$OS.type == "windows")
#                system(cmd, invisible = TRUE)
#            else system(paste(cmd, "> /dev/null"))
#            if (!res)
#                file <- file.path(tmpd, topic)
#        }
#        else {
#            rc <- .Internal(int.unzip(file.path(path, zipname),
#                topic, tmpd))
#            if (rc == 0)
#                file <- file.path(tmpd, topic)
#        }
#    }
#    file
#}

commonArea = function(...) {
  intamap:::commonArea(...)
}