getSpeciesImages <- function(species,
                             inDir,
                             outDir,
                             createStationSubfolders = FALSE,
                             IDfrom,
                             metadataSpeciesTag,
                             metadataHierarchyDelimitor = "|")
{

  file.sep <- .Platform$file.sep

  stopifnot(is.logical(createStationSubfolders))
  if (isTRUE(all(unlist(strsplit(tolower(inDir), split = file.sep)) %in%
                 unlist(strsplit(tolower(outDir), split = file.sep)))))
    stop("outDir may not be identical to or a subdirectory of inDir")

  if (length(list.files(file.path(outDir, species),
                        recursive = TRUE, full.names = TRUE, pattern = ".jpg$|.JPG$")) !=   0){
    stop(paste(outDir, file.sep, species, " is not empty",  sep = ""))
  }


  if(hasArg(metadataSpeciesTag)){
    if(class(metadataSpeciesTag) != "character"){stop("metadataSpeciesTag must be of class 'character'")}
    if(length(metadataSpeciesTag) != 1){stop("metadataSpeciesTag must be of length 1")}
  }


  if(class(IDfrom) != "character"){stop("IDfrom must be of class 'character'")}
  if(IDfrom %in% c("metadata", "directory") == FALSE) stop("'IDfrom' must be 'metadata' or 'directory'")

  if(IDfrom == "metadata"){    
    stopifnot(metadataHierarchyDelimitor %in% c("|", ":"))
    metadata.tagname <- "HierarchicalSubject"
  } 

multiple_tag_separator = "_&_"

  # define directories to save images to
  dirs.out0 <- file.path(outDir, species)


  copy.table <- data.frame()

  if(IDfrom == "directory"){

    dir.list <- list.dirs(inDir, recursive = TRUE, full.names = TRUE)         # recursive to also find species directories in camera subdirectories
    dir.list.short <- list.dirs(inDir, recursive = TRUE, full.names = FALSE)

    # find directories of requested species
    spec.index <- grep(species, dir.list, ignore.case = TRUE)

    if (length(spec.index) == 0)      stop(paste("found no directory named  '", species, "' in ", inDir, sep = ""))
    dir.list2 <- dir.list[spec.index]
    dir.list.short2 <- dir.list.short[spec.index]

    # find the station directories that have species records

    for (i in 1:length(dir.list2)) {
      files.tmp <- list.files(dir.list2[i], recursive = FALSE,
                              full.names = TRUE, include.dirs = FALSE, pattern = ".jpg$|.JPG$")
      files.tmp.short <- list.files(dir.list2[i], recursive = FALSE,
                                    full.names = FALSE, include.dirs = FALSE, pattern = ".jpg$|.JPG$")


      copy.table.tmp <- as.data.frame(matrix(nrow = length(files.tmp),
                                             ncol = 5))
      colnames(copy.table.tmp) <- c("DirectoryOrig", "FilenameOrig", "DirectoryCopy",  "FilenameCopy", "copy_ok")

      copy.table.tmp$DirectoryOrig <- dir.list2[i]
      copy.table.tmp$FilenameOrig <- files.tmp.short

      if(isTRUE(createStationSubfolders))  copy.table.tmp$DirectoryCopy <- file.path(dirs.out0, unlist(strsplit(dir.list.short2[i], split = file.sep))[1])
      if(!isTRUE(createStationSubfolders)) copy.table.tmp$DirectoryCopy <- dirs.out0
      copy.table.tmp$FilenameCopy <- files.tmp.short

      copy.table <- rbind(copy.table, copy.table.tmp)

      rm(copy.table.tmp)
    }

  } else {

    dirs <- list.dirs(inDir, recursive = FALSE, full.names = TRUE)     # exiftool will work through directory recursively, so set to FALSE here
    dirs_short <- list.dirs(inDir, recursive = FALSE, full.names = FALSE)


    for(i in 1:length(dirs_short)){

      command.tmp <- paste('exiftool -t -q -r -f -Directory -FileName -HierarchicalSubject -ext JPG "', dirs[i], '"', sep = "")
      colnames.tmp <- c("Directory", "FileName", "HierarchicalSubject")

      # run exiftool and make data frame
      tmp1 <- strsplit(system(command.tmp, intern=TRUE), split = "\t")

      metadata.tmp <- as.data.frame(matrix(unlist(lapply(tmp1, FUN = function(X){X[2]})),
                                           ncol = length(colnames.tmp),
                                           byrow = TRUE),
                                    stringsAsFactors = FALSE)
      colnames(metadata.tmp) <- colnames.tmp



      if(class(metadata.tmp) == "data.frame"){

        message(paste(dirs_short[i], ": ", nrow(metadata.tmp), "images"))

        # add metadata
        metadata.tmp <- addMetadataAsColumns (intable                    = metadata.tmp,
                                              metadata.tagname           = metadata.tagname,
                                              metadataHierarchyDelimitor = metadataHierarchyDelimitor,
                                              multiple_tag_separator     = multiple_tag_separator)

        # assign species ID
        metadata.tmp <- assignSpeciesID (intable                = metadata.tmp,
                                         IDfrom                 = IDfrom,
                                         metadataSpeciesTag     = metadataSpeciesTag,
                                         speciesCol             = "species",
                                         dirs_short             = dirs_short,
                                         i_tmp                  = i,
                                         multiple_tag_separator = multiple_tag_separator
        )

        # filter by species
        metadata.tmp <- metadata.tmp[metadata.tmp$species == species,]

        if(nrow(metadata.tmp) >= 1){
          copy.table.tmp <- data.frame(DirectoryOrig = metadata.tmp$Directory,
                                       FilenameOrig  = metadata.tmp$FileName,
                                       DirectoryCopy = NA,
                                       FilenameCopy  = metadata.tmp$FileName)

          if(isTRUE(createStationSubfolders))  copy.table.tmp$DirectoryCopy <- file.path(dirs.out0, dirs_short[i])
          if(!isTRUE(createStationSubfolders)) copy.table.tmp$DirectoryCopy <- dirs.out0

          copy.table <- rbind(copy.table, copy.table.tmp)
        }

      } # end if(class(metadata.tmp) == "data.frame"){
    }   # end for(i in 1:length(dirs_short)){
  }     # end else (IDfrom)

  # if species not found, exit
  if(nrow(copy.table) == 0) {stop(paste("did not find species '", species, "' in  ", inDir, sep = ""), call. = FALSE)}

  # check if there are any duplicate file names (would cause overwriting of images)
  if(any(tapply(copy.table$FilenameCopy,
                INDEX = copy.table$DirectoryCopy,
                FUN   = function(X){any(duplicated(X))}))){
    message(table(copy.table$FilenameCopy)[table(copy.table$FilenameCopy) >= 2])
    stop("Cannot copy. there are duplicate filenames", call. = FALSE)
  }


  # create directories
  sapply(unique(copy.table$DirectoryCopy), FUN = dir.create, recursive = TRUE, showWarnings = FALSE)

  # copy images
  copy.table$copy_ok <- file.copy(from = file.path(copy.table$DirectoryOrig, copy.table$FilenameOrig),
                                  to   = file.path(copy.table$DirectoryCopy, copy.table$FilenameCopy))

  # little summary line
  message(paste(toupper(species), " - copied", sum(copy.table$copy_ok),
              "out of", nrow(copy.table), "images in", inDir, "  to  ",
              file.path(outDir, species)))

  copy.table <- cbind(species = species, copy.table)
  return(copy.table)
}