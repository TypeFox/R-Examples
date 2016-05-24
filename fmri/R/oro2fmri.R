############################################################################
## as("nifti", "fmridata")
############################################################################

## setAs("nifti", "fmridata", function(from) { oro2fmri(from) })

oro2fmri <- function(from, value=NULL, level=0.75, setmask=TRUE) {
  ## Convert nifti() S4 object to "fmridata" S3 object
  header <- vector("list", length(slotNames(from)))
  for (i in 2:length(header)) {
    header[[i]] <- slot(from, slotNames(from)[i])
  }
  names(header) <- sub("_", "", slotNames(from))
  index <- ! (slotNames(from) %in%
              c(".Data","trail","extensions","extender","reoriented"))
  header <- header[index]
  index <- names(header) %in% "datatype"
  names(header)[index][1] <- "datatype1"
  index <- names(header) %in% "dim"
  names(header)[index] <- "dimension"
  header$endian <- "little"
  header$extension <-  vector("raw", 4) # NULL
  dx <- header$dimension[2]
  dy <- header$dimension[3]
  dz <- header$dimension[4]
  dt <- header$dimension[5]
  dd <- ifelse(header$dimension[1] == 5, header$dimension[6], 1)
  if (min(abs(header$pixdim[2:4])) != 0) {
    weights <- abs(header$pixdim[2:4] / min(abs(header$pixdim[2:4])))
  } else {
    weights <- NULL
  }
  ttt <- from@.Data
  if (dd == 1) {
    dim(ttt) <- c(dx,dy,dz,dt)
  } else {
    if (dt == 1) {
      dim(ttt) <- c(dx,dy,dz,dd)
    } else {
      dim(ttt) <- c(dx,dy,dz,dt,dd)
    }
  }
  if (dd == 1) {
    mask <- array(TRUE, c(dx,dy,dz))
    if (setmask) {
      mask[ttt[,,,1] < quantile(ttt[,,,1], level, na.rm=TRUE)] <- FALSE
      dim(ttt) <- c(prod(dim(ttt)[1:3]), dim(ttt)[4])
      na <- ttt %*% rep(1, dim(ttt)[2])
      mask[is.na(na)] <- FALSE
      ttt[is.na(na), ] <- 0
      dim(mask) <- c(dx, dy, dz)
      mask <- connect.mask(mask)
    }
    z <- list(ttt = writeBin(as.numeric(ttt), raw(), 4), 
              format = "NIFTI",
              delta = header$pixdim[2:4],
              origin = c(header$qoffsetx, header$qoffsety, header$qoffsetz),
              orient = NULL, 
              dim = header$dimension[2:5],
              dim0 = header$dimension[2:5],
              roixa = 1,
              roixe = dx,
              roiya = 1,
              roiye = dy,
              roiza = 1, 
              roize = dz,
              roit = 1:dd,
              weights = weights,
              header = header,
              mask = mask)
    class(z) <- "fmridata"
  } else {
        z <- list(ttt = writeBin(as.numeric(ttt), raw(), 4),
              format = "NIFTI",
              delta = header$pixdim[2:4],
              origin = c(header$qoffsetx, header$qoffsety, header$qoffsetz),
              orient = NULL,
              dim = c(dx, dy, dz, dd),
              dim0 = c(dx, dy, dz, dd),
              roixa = 1,
              roixe = dx,
              roiya = 1,
              roiye = dy,
              roiza = 1,
              roize = dz,
              roit = 1:dd,
              weights = weights,
              header = header)
  }
  attr(z, "file") <- ""
  return(z)
}

fmri2oro <- function(from, value=NULL, verbose=FALSE, reorient=FALSE,
                     call=NULL) {
  ## Convert "fmridata" S3 object to nifti() S4 object
  require(oro.nifti)
  nim <- nifti()
  nim@"sizeof_hdr" <- from$header$sizeofhdr
  nim@"data_type" <- from$header$datatype1
  nim@"db_name" <- from$header$dbname
  nim@"extents" <- from$header$extents
  nim@"session_error" <- from$header$sessionerror
  nim@"regular" <- from$header$regular
  nim@"dim_info" <- from$header$diminfo
  nim@"dim_" <- from$header$dimension
  nim@"intent_p1" <- from$header$intentp1
  nim@"intent_p2" <- from$header$intentp2
  nim@"intent_p3" <- from$header$intentp3
  nim@"intent_code" <- from$header$intentcode
  nim@"datatype" <- from$header$datatype
  nim@"bitpix" <- from$header$bitpix
  nim@"slice_start" <- from$header$slicestart
  nim@"pixdim" <- from$header$pixdim
  nim@"vox_offset" <- from$header$voxoffset
  nim@"scl_slope" <- from$header$sclslope
  nim@"scl_inter" <- from$header$sclinter
  nim@"slice_end" <- from$header$sliceend
  nim@"slice_code" <- as.numeric(from$header$slicecode)
  nim@"xyzt_units" <- as.numeric(from$header$xyztunits)
  nim@"cal_max" <- from$header$calmax
  nim@"cal_min" <- from$header$calmin
  nim@"slice_duration" <- from$header$sliceduration
  nim@"toffset" <- from$header$toffset
  nim@"glmax" <- from$header$glmax
  nim@"glmin" <- from$header$glmin
  nim@"descrip" <- from$header$describ ### bug in fmri? ###
  nim@"aux_file" <- from$header$auxfile
  nim@"qform_code" <- from$header$qform
  nim@"sform_code" <- from$header$sform
  nim@"quatern_b" <- from$header$quaternb
  nim@"quatern_c" <- from$header$quaternc
  nim@"quatern_d" <- from$header$quaternd
  nim@"qoffset_x" <- from$header$qoffsetx
  nim@"qoffset_y" <- from$header$qoffsety
  nim@"qoffset_z" <- from$header$qoffsetz
  nim@"srow_x" <- from$header$srowx
  nim@"srow_y" <- from$header$srowy
  nim@"srow_z" <- from$header$srowz
  nim@"intent_name" <- from$header$intentname
  nim@"magic" <- from$header$magic
  nim@"extender" <- from$header$extension
  ## convert voxel values from "raw"
  data <- readBin(from$ttt, "numeric", length(from$ttt)/4, 4)
  ## min/max values for visualization
  nim@"cal_max" <- as.numeric(max(data, na.rm=TRUE))
  nim@"cal_min" <- as.numeric(min(data, na.rm=TRUE))
  ## coerce voxel values into array
  dims <- 2:(1+nim@"dim_"[1])
  if (reorient) {
    nim@.Data <- reorient(nim, data, verbose)
    nim@"reoriented" <- TRUE
  } else {
    nim@.Data <- array(data, nim@"dim_"[dims])
  }
  ## Check validity
  validNIfTI <- getValidity(getClassDef("nifti"))
  validNIfTI(nim)
  if (getOption("niftiAuditTrail")) {
    if (is.null(call)) {
      call <- match.call()
    }
    nim <- niftiExtensionToAuditTrail(nim, workingDirectory=getwd(),
                                      filename="fmridata", call=call)
  }
  return(nim)
}
