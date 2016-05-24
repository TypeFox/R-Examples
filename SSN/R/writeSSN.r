################################################################################
# Write SSN object
#
################################################################################

writeSSN <- function(ssn, filename, o.write = FALSE) {

  suppressWarnings(warning("writeSSN"))

  if(!file.exists(filename)) {
      dir.create(filename)
  } else {
      if(o.write == FALSE) stop("file exists and o.write = FALSE")
      if(o.write == TRUE) {
          unlink(filename, recursive = TRUE)
          dir.create(filename) }
  }

  oldwd <- getwd()
  setwd(filename)

  on.exit({

	setwd(oldwd)
  })


  ssn.tmp <- ssn
  pred.len <- length(ssn.tmp@predpoints@ID)


  ssn.tmp@path <- getwd()


  ssn.list <- list.files(ssn@path)
  ind.shp <- file_ext(ssn.list) == "shp"
  sub.list <- substr(ssn.list[ind.shp], 1, nchar(ssn.list[ind.shp])-4)

  ind.list <- vector(mode = "logical", length = length(ssn.list))

  for(j in 1:length(ssn.list)) {
    tmp <- unlist(strsplit(ssn.list[j], "[.]"))
    if((tmp[1] %in% sub.list) == TRUE) ind.list[j] <- TRUE
  }
  ssn.list <- ssn.list[!ind.list]

  for (i in 1:length(ssn.list)) {
    fn.old <- file.path(ssn@path, ssn.list[i])
    if (basename(fn.old) != "distance")  {
        #if (substr(basename(fn.old), 1, 5) != "sites") {

        fn.new <- file.path(ssn.tmp@path, ssn.list[i])

        file.copy(fn.old, fn.new, overwrite = TRUE)
        #}
    }
  }
  rm(fn.old, fn.new)

  ## Copy SITES--------------------------------------------
  # create spatial points data frame
  ## coords <- ssn.tmp@obspoints@SSNPoints[[1]]@point.coords

  ## proj4string <- ssn.tmp@proj4string
  ## data.tmp <- ssn.tmp@obspoints@SSNPoints[[1]]@point.data

  ## ind.xy <- names(data.tmp) == "coords_x1" | names(data.tmp) == "coords_x2"
  ## if (sum(ind.xy) > 0) {
  ##     data.tmp <- data.tmp[,!ind.xy]}

  ## sites.sub <- SpatialPointsDataFrame(coords = coords, data = data.tmp, proj4string = proj4string)
  sites.sub <- as.SpatialPointsDataFrame(ssn.tmp)

  writeSpatialShape(sites.sub, "sites")
  write.dbf.SSN(sites.sub@data, "sites")

  ## Copy edges --------------------------------------------
  edges.sub <- as.SpatialLinesDataFrame.SpatialStreamNetwork(ssn.tmp)
  writeSpatialShape(edges.sub, "edges")
  write.dbf.SSN(edges.sub@data, "edges")

  ## Copy prediction sites

  if (pred.len > 0) {
      pred.name.list <- ssn.tmp@predpoints@ID
      for(i in 1:pred.len) {
          pred.name <- pred.name.list[i]
          ## coords <- ssn.tmp@predpoints@SSNPoints[[i]]@point.coords

          ## proj4string <- ssn.tmp@proj4string
          ## data.tmp <- ssn.tmp@predpoints@SSNPoints[[i]]@point.data

          ## ind.xy <- names(data.tmp) == "coords_x1" | names(data.tmp) == "coords_x2"
          ## if (sum(ind.xy) > 0) {
          ##     data.tmp <- data.tmp[,!ind.xy]}

          ## preds.sub <- SpatialPointsDataFrame(coords = coords, data = data.tmp,
          ##     proj4string = proj4string)

          preds.sub <- as.SpatialPointsDataFrame(ssn.tmp, pred.name)
          writeSpatialShape(preds.sub, pred.name)
          write.dbf.SSN(preds.sub@data, pred.name)

          #rm(pred.name, coords, data.tmp, ind.xy, preds.sub)
          rm(pred.name, preds.sub)
      }
  }

  ## Import SSN

  if (pred.len ==0) {
    ssn.tmp <- importSSN(ssn.tmp@path, o.write = T)}

  if (pred.len > 0) {
      ssn.tmp <- importSSN(ssn.tmp@path, o.write = T)

      for (j in 1:pred.len) {
        ssn.tmp <- importPredpts(ssn.tmp, pred.name.list[j], obj.type = "ssn")}
  }

  setwd(oldwd)

  ##ssn.tmp
}

