#' @title Merge DAS Data Into Acoustic Detections
#' @description Fill in sighting information for acoustic detections from DAS file.
#'
#' @param acoust.file filename of a .csv file of acoustic detections exported from Access database.
#' @param das.file filename of a DAS file.
#' @param out.file desired output .csv filename.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom utils read.csv write.csv
#' @export
#' 
das.merge.acoustic <- function(acoust.file, das.file, out.file) {
  das.df <- das.read(das.file)
  sight.df <- das.df[das.df$Event == "S", ]
  rownames(sight.df) <- sight.df$Sight
  a.df <- das.df[das.df$Event == "A", ]
  rownames(a.df) <- a.df$Sight
  for(i in rownames(sight.df)) {
    sight.df[i, "Spp1"] <- a.df[i, "Spp1"]
    sight.df[i, "Spp2"] <- a.df[i, "Spp2"]
    sight.df[i, "Spp3"] <- a.df[i, "Spp3"]
  }

  acoust.df <- read.csv(acoust.file, stringsAsFactors = FALSE)
  acoust.df$first_dist_units <- gsub("[[:blank:]]", "", acoust.df$first_dist_units)
  acoust.df$detection_type <- toupper(gsub("[[:blank:]]", "", acoust.df$detection_type))
  acoust.df$first_detection <- toupper(gsub("[[:blank:]]", "", acoust.df$first_detection))

  need.das <- which(acoust.df$first_detection == "V")

  for(i in need.das) {
    vid <- as.character(acoust.df$vis_id[i])
    acoust.df[i, "date_time_start"] <- format(sight.df[vid, "Date"], "%m/%d/%Y %T")
    acoust.df[i, "latlong_LAT"] <- sight.df[vid, "Lat"]
    acoust.df[i, "latlong_LON"] <- sight.df[vid, "Long"]
    acoust.df[i, "species1_class1"] <- as.numeric(sight.df[vid, "Spp1"])
    acoust.df[i, "species2_class1"] <- as.numeric(sight.df[vid, "Spp2"])
    acoust.df[i, "species1_class2"] <- as.numeric(sight.df[vid, "Spp3"])
    acoust.df[i, "class1"] <- "V"
    acoust.df[i, "class2"] <- ifelse(is.na(sight.df[vid, "Spp3"]), NA, "V")
    acoust.df[i, "first_angle"] <- sight.df[vid, "Angle"]
    acoust.df[i, "first_dist"] <- sight.df[vid, "Dist"]
    acoust.df[i, "loc_type"] <- ""
    acoust.df[i, "first_dist_units"] <- "nmi"
    acoust.df[i, "beam_dist"] <- -999
    acoust.df[i, "det_dist"] <- -999
    acoust.df[i, "loc_method"] <- "Vis"
  }

  write.csv(acoust.df, file = out.file, na = "", row.names = FALSE)
  invisible(acoust.df)
}
