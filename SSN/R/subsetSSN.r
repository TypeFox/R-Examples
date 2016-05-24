################################################################################
# Subset SSN object
#
# Observed sites are selected based on a logical expression. The results are
# exported to a new .ssn folder. There is also an option to clip edges and
# prediction sites, based on the NetworkID values contained in the observed
# sites selection set.
################################################################################

subsetSSN <- function(ssn, filename, subset, clip = FALSE) {

  dir.create(filename)
  oldwd <- getwd()
  setwd(filename)

  on.exit({
	setwd(oldwd)
  })


  ssn.tmp <- ssn
  pred.len <- length(ssn.tmp@predpoints@ID)

  #Check to see if attribute exists
  s <- deparse(substitute(subset))

  ind <- eval(substitute(subset), ssn.tmp@obspoints@SSNPoints[[1]]@point.data)
  ind.na <- is.na(ind)
  ind[ind.na] <- FALSE
  rm(ind.na)

  if (sum(ind) == 0) {
    stop("No records were selected based on logical expression")}

  ssn.tmp@obspoints@SSNPoints[[1]]@point.data<- ssn@obspoints@SSNPoints[[1]]@point.data[ind,]
  #ssn.tmp@obspoints@SSNPoints[[1]]@network.point.coords<-ssn@obspoints@SSNPoints[[1]]@network.point.coords[ind,]

  if(sum(ind)==1)
    ssn.tmp@obspoints@SSNPoints[[1]]@point.coords<- as.matrix(t(ssn@obspoints@SSNPoints[[1]]@point.coords[ind, ]))
  else ssn.tmp@obspoints@SSNPoints[[1]]@point.coords<- ssn@obspoints@SSNPoints[[1]]@point.coords[ind,]

  ssn.tmp@path <- getwd()
  #ssn.tmp@path <- filename

  rm(ind)

  if (clip == FALSE) {
    ssn.list <- list.files(ssn@path)
    for (i in 1:length(ssn.list)) {
      fn.old <- file.path(ssn@path, ssn.list[i])
      if (basename(fn.old) != "distance")  {
        if (substr(basename(fn.old), 1, 5) != "sites") {
          #fn.new <- file.path(filename, ssn.list[i])
          fn.new <- file.path(ssn.tmp@path, ssn.list[i])

          file.copy(fn.old, fn.new, overwrite = TRUE)}}
    }
    rm(fn.old, fn.new)

  } else {
  # Clip everything to the new extent
    #subset edges
    ind.edges <- eval(substitute(subset), ssn.tmp@data)
    ind.na <- is.na(ind.edges)
    ind.edges[ind.na] <- FALSE
    rm(ind.na)
    sl <- ssn.tmp@lines[ind.edges]
    proj4string <- ssn.tmp@proj4string
    edges.sl <- SpatialLines(sl, proj4string)
    edges.sub <- SpatialLinesDataFrame(edges.sl, ssn.tmp@data[ind.edges,], match.ID = FALSE)

    #edges.sub<- edges.sldf[edges.sldf$netID %in% netID.list,]
    writeSpatialShape(edges.sub, "edges")
    write.dbf.SSN(edges.sub@data, "edges", max_nchar = 30)

    ind.dup <- !duplicated(edges.sub@data$netID)
    netID.list <- edges.sub@data$netID[ind.dup]

    rm(proj4string, sl, edges.sl, edges.sub, ind.edges)

    #subset pred points
    #pred.len <- length(ssn.tmp@predpoints@ID)
    if (pred.len > 0) {
      for (i in 1:pred.len) {
        pred.name <- ssn.tmp@predpoints@ID[[i]]

        ind.preds <- eval(substitute(subset), ssn.tmp@predpoints@SSNPoints[[i]]@point.data)
        ind.na <- is.na(ind.preds)
        ind.preds[ind.na] <- FALSE
        rm(ind.na)

        if(sum(ind.preds)==1)
          coords<- as.matrix(t(ssn@predpoints@SSNPoints[[i]]@point.coords[ind.preds, ]))
        else coords<- ssn@predpoints@SSNPoints[[i]]@point.coords[ind.preds,]
        #coords <- ssn.tmp@predpoints@SSNPoints[[i]]@point.coords[ind.preds,]
        proj4string <- ssn.tmp@predpoints@SSNPoints[[i]]@proj4string
        data.tmp <- ssn.tmp@predpoints@SSNPoints[[i]]@point.data[ind.preds,]

        ind.xy <- names(data.tmp) == "coords_x1" | names(data.tmp) == "coords_x2"
        if (sum(ind.xy) > 0) {
          data.tmp <- data.tmp[,!ind.xy]}

        preds.sub <- SpatialPointsDataFrame(coords = coords, data = data.tmp, proj4string = proj4string)
        writeSpatialShape(preds.sub, pred.name)
        write.dbf.SSN(preds.sub@data, pred.name, max_nchar = 30)
        rm(coords, proj4string, data.tmp, preds.sub, ind.preds, ind.xy)
    }}

    #copy netID files
    for (i in 1:length(netID.list)) {
      fn.old <- file.path(ssn@path, paste("netID", netID.list[i], ".dat", sep = ""))
      fn.new <- file.path(ssn.tmp@path, paste("netID", netID.list[i], ".dat", sep = ""))
      file.copy(fn.old, fn.new, overwrite = TRUE)
    }
    rm(fn.new, fn.old)
  }

  # create spatial points data frame...
  coords <- ssn.tmp@obspoints@SSNPoints[[1]]@point.coords

  proj4string <- ssn.tmp@proj4string
  data.tmp <- ssn.tmp@obspoints@SSNPoints[[1]]@point.data

  ind.xy <- names(data.tmp) == "coords_x1" | names(data.tmp) == "coords_x2"
  if (sum(ind.xy) > 0) {
      data.tmp <- data.tmp[,!ind.xy]}

  sites.sub <- SpatialPointsDataFrame(coords = coords, data = data.tmp, proj4string = proj4string)
  writeSpatialShape(sites.sub, "sites")
  write.dbf.SSN(sites.sub@data, "sites", max_nchar = 30)

  #create binaryID database in new .ssn folder
  ########### This could be improved. Only import the rids for the segments in
  # the subset, instead of everything in the netID tables.
  #createBinaryID(ssn.tmp, o.write = T)


  if (pred.len ==0) {
    ssn.tmp <- importSSN(ssn.tmp@path, o.write = T)}

  if (pred.len > 0) {
      pred.name.list <- ssn.tmp@predpoints@ID
      ssn.tmp <- importSSN(ssn.tmp@path, o.write = T)

      for (j in 1:pred.len) {
        ssn.tmp <- importPredpts(ssn.tmp, pred.name.list[j], obj.type = "ssn")}
  }

  setwd(oldwd)

  ssn.tmp
}








