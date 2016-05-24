imageRename <- function(inDir,
                        outDir,
                        hasCameraFolders,
                        keepCameraSubfolders,
                        copyImages = FALSE,
                        writecsv = FALSE){

  wd0 <- getwd()
  on.exit(setwd(wd0))
  
  file.sep <- .Platform$file.sep
  
  stationCol <- "Station"
  cameraCol  <- "Camera"

  # check call for consistency
  stopifnot(is.logical(copyImages))
  stopifnot(is.logical(writecsv))
  stopifnot(is.logical(hasCameraFolders))
  
    if(isTRUE(hasCameraFolders)){
      stopifnot(hasArg(keepCameraSubfolders))
      stopifnot(is.logical(keepCameraSubfolders))
   } else {
    keepCameraSubfolders <- FALSE
   }
   if(hasArg(hasCameraFolders) & hasArg(keepCameraSubfolders)){
    if(keepCameraSubfolders == TRUE & hasCameraFolders == FALSE){stop("If hasCameraFolders is FALSE, keepCameraSubfolders must be FALSE too")}
   }
  
  
  stopifnot(length(inDir) == 1)
  if(hasArg(outDir)){
    stopifnot(length(outDir) == 1)
    if(isTRUE(all(unlist(strsplit(tolower(inDir), split = file.sep)) %in%
                  unlist(strsplit(tolower(outDir), split = file.sep))))) stop("outDir may not be identical to or a subdirectory of inDir")

    list.files.tmp <- list.files(outDir, recursive = TRUE)
    if(copyImages == TRUE & length(list.files.tmp) >= 1) stop("outDir must be empty if you wish to copy images")
  }
  if(copyImages == TRUE){
    
    if(any(c(grep("/$", inDir) == 1, grep("/$", outDir) == 1))) stop("inDir and outDir may not end with /")

  } else {
    if(isTRUE(grepl("/$", inDir))) stop("inDir may not end with /")
  }
  if(Sys.which("exiftool") == "") stop("cannot find ExifTool")
  if(isTRUE(writecsv) & hasArg(outDir) == FALSE) stop("writecsv is TRUE. Please specify outDir")




  # list of subdirectories of inDir
  dirs <- list.dirs(inDir, full.names = TRUE, recursive = FALSE)
  dirs_short <- list.dirs(inDir, full.names = FALSE , recursive = FALSE)

  # make sure none is empty
  list_n_files <- lapply(dirs, list.files, pattern = ".jpg$|.JPG$", recursive = TRUE)

  if(any(unlist(lapply(list_n_files, length)) == 0)){
    warning("at least one station directory contains no JPEGs:  ", paste(names(which(lapply(list_n_files, length) == 0)), collapse = "; "), call. = FALSE)
  }


  # function body

  copy.info.table <- data.frame()

  for(i in which(lapply(list_n_files, length) >= 1)){     # for all non-empty directories

    # check if there are any images not in camera trap subfolders
    if(hasCameraFolders == TRUE && length(list.files(dirs[i], pattern = ".jpg$|.JPG$", recursive = FALSE)) >= 1){
      stop(paste("Directory ", dirs[i], " contains images not sorted into Camera Trap subfolders. Check argument 'hasCameraFolders'"))
    }

    # run exiftool to get image date and time

    command.tmp <- paste('exiftool -q -f -t -r -Directory -FileName -EXIF:DateTimeOriginal -ext JPG "', dirs[i], '"', sep = "")
    tmp1 <-  strsplit(system(command.tmp, intern=TRUE), split = "\t")
    colnames.tmp <- c("Directory", "FileName", "DateTimeOriginal")
    metadata.tmp <- as.data.frame(matrix(unlist(lapply(tmp1, FUN = function(X){X[2]})),
                                         ncol = length(colnames.tmp),
                                         byrow = TRUE),
                                  stringsAsFactors = FALSE)

    colnames(metadata.tmp) <- colnames.tmp
    rm(colnames.tmp, tmp1)

    if(length(metadata.tmp) == 0){
      length.tmp <- length(list.files(dirs[i], pattern = ".jpg$|JPG$", ignore.case = TRUE, recursive = TRUE))
      warning(paste(dirs[i], "seems to contain no images;", " found", length.tmp, "jpgs"))    # give message if station directory contains no jpgs
    } else {

      message(paste(dirs_short[i], ":", nrow(metadata.tmp), "images"))

      if(isTRUE(hasCameraFolders)){
        filenames_by_subfolder <- lapply(as.list(list.dirs(dirs[i], full.names =TRUE, recursive = FALSE)),
                                         FUN = list.files, pattern = ".jpg$|.JPG$", recursive = TRUE, ignore.case = TRUE)
        metadata.tmp$CameraID <- rep(list.dirs(dirs[i], full.names = FALSE, recursive = FALSE),
                                     times = unlist(lapply(filenames_by_subfolder, length)))
        metadata.tmp$Station <- rep(dirs_short[i], times = nrow(metadata.tmp))
        colnames(metadata.tmp)[grep("CameraID", colnames(metadata.tmp))] <- cameraCol
        colnames(metadata.tmp)[grep("Station", colnames(metadata.tmp))] <- stationCol
      } else {
        metadata.tmp$CameraID <- "NA"                                                    # "camera" name if no camera id available (only station ids)
        metadata.tmp$Station <- rep(dirs_short[i], times = nrow(metadata.tmp))
        colnames(metadata.tmp)[grep("CameraID", colnames(metadata.tmp))] <- cameraCol
        colnames(metadata.tmp)[grep("Station", colnames(metadata.tmp))] <- stationCol
      }

      # make time readable
      metadata.tmp$DateTimeOriginal <- as.POSIXct(strptime(x = as.character(metadata.tmp$DateTimeOriginal),
                                                           format = "%Y:%m:%d %H:%M:%S", tz = "UTC"))
      metadata.tmp$DateTimeOriginal2 <- as.POSIXlt(metadata.tmp$DateTimeOriginal)

      # exclude images for which no DateTimeOriginal was found
      suppressWarnings(na.date.rows <- which(is.na(metadata.tmp$DateTimeOriginal)))
      if(length(na.date.rows) != 0){
        #metadata.tmp.na.date <- metadata.tmp[na.date.rows,]
        warning(paste("couldn't read DateTimeOriginal tag of ",
                    file.path(metadata.tmp2$Directory,  metadata.tmp2$FileName)[na.date.rows]))
        metadata.tmp <- data.frame(metadata.tmp, DateReadable = NA)
        metadata.tmp$DateReadable[-na.date.rows] <- TRUE
        metadata.tmp$DateReadable[na.date.rows] <- FALSE
      } else {
        metadata.tmp$DateReadable <- TRUE
      }
      rm(na.date.rows)

      # rearrange column order
      metadata.tmp <- metadata.tmp[,c("Directory",  "FileName", stationCol, cameraCol,
                                      "DateTimeOriginal", "DateTimeOriginal2", "DateReadable")]

      # find images taken within 1 minute of one another (to append number)
      metadata.tmp$DateTimeOriginal2$sec <- 0

      if(isTRUE(hasCameraFolders)){
        metadata.tmp.split <- split(x = metadata.tmp, f = list(metadata.tmp[,cameraCol],
                                                               as.POSIXct(metadata.tmp$DateTimeOriginal2)
        ), drop = TRUE)
      } else {
        metadata.tmp.split <- split(x = metadata.tmp, f = list(metadata.tmp[,stationCol],
                                                               as.POSIXct(metadata.tmp$DateTimeOriginal2)
        ), drop = TRUE)
      }

      metadata.tmp.split2 <- lapply(metadata.tmp.split, FUN = function(X){
        X2 <- X[with(X, order(DateTimeOriginal)), ]
        cbind(X2, minute.append = seq(from = 1, to = nrow(X2), by = 1))
      })

      metadata.tmp2 <- do.call("rbind", metadata.tmp.split2)    # reassemble
      if(any(metadata.tmp$DateReadable == FALSE)){
        metadata.tmp2 <- rbind(cbind(metadata.tmp[metadata.tmp$DateReadable == FALSE,], minute.append = NA) ,metadata.tmp2)
      }
      rm(metadata.tmp.split, metadata.tmp.split2)

      # convert time object character vector to and format for outfilename
      time.tmp <- gsub(pattern = ":", replacement = "-", metadata.tmp2$DateTimeOriginal)
      time.tmp2 <- gsub(pattern = " ", replacement = "__", time.tmp)
      metadata.tmp2$DateTime_for_filename <- time.tmp2
      rm(time.tmp, time.tmp2)

      # create outfilename
      if(isTRUE(copyImages)){
        if(keepCameraSubfolders == TRUE){
          metadata.tmp2$outDir <- file.path(outDir, metadata.tmp2[,stationCol], metadata.tmp2[,cameraCol])
        } else {
          metadata.tmp2$outDir <- file.path(outDir, metadata.tmp2[,stationCol])
        }
      }

      if(isTRUE(hasCameraFolders)){
        metadata.tmp2$filename_new <- paste(paste(metadata.tmp2[,stationCol],
                                                  metadata.tmp2[,cameraCol],
                                                  paste(metadata.tmp2$DateTime_for_filename, "(",metadata.tmp2$minute.append, ")", sep = ""),
                                                  sep = "__"),
                                            ".JPG", sep = "")
      } else {
        metadata.tmp2$filename_new <- paste(paste(metadata.tmp2[,stationCol],
                                                  paste(metadata.tmp2$DateTime_for_filename, "(",metadata.tmp2$minute.append, ")", sep = ""),
                                                  sep = "__"),
                                            ".JPG", sep = "")
      }

      # copy images
      if(isTRUE(copyImages)){
        if(keepCameraSubfolders == TRUE){
          wd_out_tmp <- file.path(outDir, dirs_short[i], unique(metadata.tmp2[,cameraCol]))
        } else {
          wd_out_tmp <- file.path(outDir, dirs_short[i])
        }
        sapply(wd_out_tmp, dir.create, recursive = TRUE)
        metadata.tmp2$CopyStatus[metadata.tmp2$DateReadable == TRUE] <- file.copy(from = apply(metadata.tmp2[metadata.tmp2$DateReadable == TRUE, c("Directory", "FileName")], MARGIN = 1, FUN = paste, collapse = file.sep),
                                                                               to = apply(metadata.tmp2[metadata.tmp2$DateReadable == TRUE, c("outDir", "filename_new")], MARGIN = 1, FUN = paste, collapse = file.sep),
                                                                               overwrite = FALSE)
        metadata.tmp2$CopyStatus[metadata.tmp2$DateReadable == FALSE] <- FALSE
        rm(wd_out_tmp)


      } else {
        metadata.tmp2$CopyStatus <- FALSE
      }
      copy.info.table <- rbind(copy.info.table, metadata.tmp2)
      rm(metadata.tmp2)
    }
  }

  rownames(copy.info.table) <- NULL
  copy.info.table <- copy.info.table[,-which(names(copy.info.table) %in%
                                               c("DateTimeOriginal2", "DateTime_for_filename", "minute.append"))]
  # save table
  if(writecsv == TRUE){
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
    setwd(outDir)
    write.csv(copy.info.table, file = paste("_renaming_table_", Sys.Date(), ".csv", sep = ""),
              row.names = FALSE)
  }
  return(copy.info.table)
}