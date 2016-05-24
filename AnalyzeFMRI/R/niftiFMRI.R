f.read.nifti.header <- function(file){
  #This function reads in the header information from a NIFTI format file
  
  is.nii <- substring(file, nchar(file) - 2,nchar(file)) # file extension

  file.name <- substring(file, 1, nchar(file) - 4) # file name without extension

  if (is.nii == "nii") {
    file.hdr <- paste(file.name, ".nii", sep = "")
    file.img <- paste(file.name, ".nii", sep = "")
  } 
  else {
    file.hdr <- paste(file.name, ".hdr", sep = "")
    file.img <- paste(file.name, ".img", sep = "")

    if(file.exists(file.hdr) == FALSE) stop(paste(file.hdr, "not found"))
  }

  
  

# Detect whether the data is big or little endian. The first part of a .hdr file is the size of the file which is always a C int (i.e. 4 bytes) and always has value 348.
# Therefore trying to read it in assuming little-endian will tell you if that is the correct mode.

    swap <- 0

    if(.C("swaptest_wrap_JM",
          ans = integer(1),
          file.hdr,
          PACKAGE="AnalyzeFMRI")$ans != 348) # $ans is sizeof_hdr
        swap <- 1


    
# A C function is used to read in all the components of the .hdr file
    a<-.C("read_nifti_header_wrap_JM",
          file.hdr,  # name of hdr file (with .hdr extension)         1
          as.integer(swap), # as defined above                        2
          integer(1), # sizeof_hdr                                    3
          paste(rep(" ", 10), sep = "", collapse = ""), # data_type   4
          paste(rep(" ", 18), sep = "", collapse = ""), # db_name     5
          integer(1), # extents                                       6
          integer(1), # session_error                                 7
          paste(rep(" ", 1), sep = "", collapse = ""), # regular      8
          paste(rep(" ", 1), sep = "", collapse = ""), # dim_info     9
          integer(8), # dim                                           10
          single(1), # intent_p1                                      11
          single(1), # intent_p2                                      12
          single(1), # intent_p3                                      13
          integer(1), # intent_code                                   14
          integer(1), # datatype                                      15
          integer(1), # bitpix                                        16
          integer(1), # slice_start                                   17
          single(8), # pixdim                                         18
          single(1), # vox_offset                                     19
          single(1), #         scl_slope                              20
          single(1), #        scl_inter                               21
          integer(1), #      slice_end                                22
          paste(rep(" ", 1), sep = "", collapse = ""), # slice_code   23       
          paste(rep(" ", 1), sep = "", collapse = ""), # xyzt_units   24
          single(1), # cal_max                                        25
          single(1), # cal_min                                        26
          single(1), # slice_duration                                 27
          single(1), #         toffset                                28
          integer(1), # glmax                                         29
          integer(1), # glmin                                         30
          paste(rep(" ", 80), sep = "", collapse = ""), # descrip     31
          paste(rep(" ", 24), sep = "", collapse = ""), # aux_file    32
          integer(1), #      qform_code                               33
          integer(1), #      sform_code                               34
          single(1), # quatern_b                                      35
          single(1), # quatern_c                                      36
          single(1), # quatern_d                                      37
          single(1), # qoffset_x                                      38
          single(1), # qoffset_y                                      39
          single(1), # qoffset_z                                      40
          single(4), # srow_x                                         41
          single(4), # srow_y                                         42
          single(4), # srow_z                                         43
          paste(rep(" ", 16), sep = "", collapse = ""), # intent_name 44
          paste(rep(" ", 4), sep = "", collapse = ""), # magic        45  
          rep(" ", 4),                                 # extension    46  
          PACKAGE="AnalyzeFMRI")

# A list (called L) is created containing all the components of the .hdr (or header part of .nii) file

    L <- list()
    L$file.name <- file.img
    L$swap <- a[[2]]
    L$sizeof.hdr <- a[[3]]
    L$data.type <- if (a[[4]] == "") rawToChar(raw(10)) else a[[4]]
    L$db.name <- if (a[[5]] == "") rawToChar(raw(18)) else a[[5]]
    L$extents <- a[[6]]
    L$session.error <- a[[7]]
    L$regular <- if (a[[8]] == "") rawToChar(raw(1)) else a[[8]]
    L$dim.info <- if (a[[9]] == "") rawToChar(raw(1)) else a[[9]]
    L$dim <- a[[10]]
    L$intent.p1 <- a[[11]]
    L$intent.p2 <- a[[12]]
    L$intent.p3 <- a[[13]]
    L$intent.code <- a[[14]]
    L$datatype <- a[[15]]
    L$bitpix <- a[[16]]
    L$slice.start <- a[[17]]
    L$pixdim <- a[[18]]
    L$vox.offset <- a[[19]]  
    L$scl.slope <- a[[20]]
    L$scl.inter <- a[[21]]
    L$slice.end <- a[[22]]
    L$slice.code <- if (a[[23]] == "") rawToChar(raw(1)) else a[[23]]
    L$xyzt.units <- if (a[[24]] == "") rawToChar(raw(1)) else a[[24]]
    L$cal.max <- a[[25]]
    L$cal.min <- a[[26]]
    L$slice.duration <- a[[27]]
    L$toffset <- a[[28]]
    L$glmax <- a[[29]]
    L$glmin <- a[[30]]
    L$descrip <- if (a[[31]] == "") rawToChar(raw(80)) else a[[31]] 
    L$aux.file <- if (a[[32]] == "") rawToChar(raw(24)) else a[[32]]
    L$qform.code <- a[[33]]
    L$sform.code <- a[[34]]
    L$quatern.b <- a[[35]]
    L$quatern.c <- a[[36]]
    L$quatern.d <- a[[37]]
    L$qoffset.x <- a[[38]]
    L$qoffset.y <- a[[39]]
    L$qoffset.z <- a[[40]]
    L$srow.x <- a[[41]]
    L$srow.y <- a[[42]]
    L$srow.z <- a[[43]]
    L$intent.name <- if (a[[44]] == "") rawToChar(raw(16)) else a[[44]]
    L$magic <- if (a[[45]] == "") rawToChar(raw(4)) else a[[45]]
    L$extension <- a[[46]]
    if (L$extension[1] == "")  tmp1 <- as.integer(0) else tmp1 <- as.integer(charToRaw(L$extension[1]))
    if (L$extension[2] == "")  tmp2 <- as.integer(0) else tmp2 <- as.integer(charToRaw(L$extension[2]))
    if (L$extension[3] == "")  tmp3 <- as.integer(0) else tmp3 <- as.integer(charToRaw(L$extension[3]))
    if (L$extension[4] == "")  tmp4 <- as.integer(0) else tmp4 <- as.integer(charToRaw(L$extension[4]))
    L$extension <- c(tmp1,tmp2,tmp3,tmp4)

    return(L)}



f.nifti.file.summary <- function(file){
#This function prints out a concise summary of the contents of a NIFTI .img/.hdr image pair (or .nii file)


  is.nii <- substring(file, nchar(file) - 2,nchar(file)) # file extension

  file.name <- substring(file, 1, nchar(file) - 4) # file name without extension
  
  if (is.nii == "nii") {
    file.img <- paste(file.name, ".nii", sep = "")
  } 
  else {
    file.img <- paste(file.name, ".img", sep = "")
  }


  hdr <- f.read.nifti.header(file)


  cat("\n")
  cat("       File name:", file.img, "\n")
  cat("  Data Dimension:", paste(hdr$dim[1], "-D", sep = ""), "\n")
  cat("     X dimension:", hdr$dim[2], "\n")
  cat("     Y dimension:", hdr$dim[3], "\n")
  cat("     Z dimension:", hdr$dim[4], "\n")
  cat("  Time dimension:", hdr$dim[5], "time points", "\n")
  cat("Voxel dimensions:", paste(hdr$pixdim[2], hdr$vox.units, "x",
                                 hdr$pixdim[3], hdr$vox.units, "x",
                                 hdr$pixdim[4], hdr$vox.units), "\n")
  cat("       Data type:", hdr$data.type, paste("(", hdr$bitpix, " bits per voxel)", sep = ""), "\n")
}


f.basic.hdr.nifti.list.create <- function(dim.mat, file){

#creates a basic list that can be used to write a .hdr file in NIFTI format (or the header part of a .nii file)

  is.nii <- substring(file, nchar(file) - 2,nchar(file)) # file extension

  file.name <- substring(file, 1, nchar(file) - 4) # file name without extension
  
  if (is.nii == "nii") {
    file.hdr <- paste(file.name, ".nii", sep = "")
  } 
  else {
    file.hdr <- paste(file.name, ".hdr", sep = "")
  }
  
  dim <- c(length(dim.mat), dim.mat, rep(0, 7 - length(dim.mat)))
  
  magic <- raw(4)

  if (is.nii == "nii") {magic[1:3] <- as.raw(charToRaw("n+1"))} else {magic[1:3] <- as.raw(charToRaw("ni1"))}

  magic <- rawToChar(magic)

  l <- list(file = file.hdr,
            sizeof.hdr = 348,
            data.type = rawToChar(raw(10)),
            db.name = rawToChar(raw(18)),
            extents = integer(1),
            session.error = integer(1),
            regular = character(1),
            dim.info = character(1),
            dim = as.integer(dim),
            intent.p1 = single(1),
            intent.p2 = single(1),
            intent.p3 = single(1),
            intent.code = integer(1),
            datatype = integer(1),
            bitpix = integer(1),
            slice.start = integer(1),
            pixdim = single(8),
            vox.offset = {if (is.nii == "nii") as.single(352) else as.single(0)},
            scl.slope = single(1),
            scl.inter = single(1),
            slice.end = integer(1),
            slice.code = character(1),
            xyzt.units = character(1),
            cal.max = single(1),
            cal.min = single(1),
            slice.duration = single(1),
            toffset = single(1),
            glmax = integer(1),
            glmin = integer(1),
            descrip = rawToChar(raw(80)),
            aux.file = rawToChar(raw(24)),
            qform.code = integer(1),
            sform.code = integer(1),
            quatern.b = single(1),
            quatern.c = single(1),
            quatern.d = single(1),
            qoffset.x = single(1),
            qoffset.y = single(1),
            qoffset.z = single(1),
            srow.x = single(4),
            srow.y = single(4),
            srow.z = single(4),
            intent.name = rawToChar(raw(16)),
            magic = magic) 
  return(l)
}


f.complete.hdr.nifti.list.create <- function(file,dim.info=character(1),dim,intent.p1=single(1),intent.p2=single(1),intent.p3=single(1),intent.code=integer(1),datatype=integer(1),bitpix=integer(1),slice.start=integer(1),pixdim=single(8),scl.slope=single(1),scl.inter=single(1),slice.end=integer(1),slice.code=character(1),xyzt.units=character(1),cal.max=single(1),cal.min=single(1),slice.duration=single(1),toffset=single(1),descrip=paste(rep(" ", 80), sep = "", collapse = ""),aux.file=paste(rep(" ", 24), sep = "", collapse = ""),qform.code=integer(1),sform.code=integer(1),quatern.b=single(1),quatern.c=single(1),quatern.d=single(1),qoffset.x=single(1),qoffset.y=single(1),qoffset.z=single(1),srow.x=single(4),srow.y=single(4),srow.z=single(4),intent.name=paste(rep(" ", 16), sep = "", collapse = "")){

#creates a complete list that can be used to write a .hdr file in NIFTI format (or the header part of a .nii file)

  if (nchar(descrip>80)) {
    descrip <- substr(descrip,1,80)
    warning("descrip is too long (>80) and has been truncated")
  }

  if (nchar(aux.file>24)) {
    aux.file <- substr(aux.file,1,24)
    warning("aux.file is too long (>24) and has been truncated")
  }
  
  is.nii <- substring(file, nchar(file) - 2,nchar(file)) # file extension

  file.name <- substring(file, 1, nchar(file) - 4) # file name without extension
  
  if (is.nii == "nii") {
    file.hdr <- paste(file.name, ".nii", sep = "")
  } 
  else {
    file.hdr <- paste(file.name, ".hdr", sep = "")
  }
  

  magic <- raw(4)

  if (is.nii == "nii") {magic[1:3] <- as.raw(charToRaw("n+1"))} else {magic[1:3] <- as.raw(charToRaw("ni1"))}

  magic <- rawToChar(magic)

  l <- list(file = file.hdr,
            sizeof.hdr = 348,
            data.type = rawToChar(raw(10)),
            db.name = rawToChar(raw(18)),
            extents = integer(1),
            session.error = integer(1),
            regular = "r",
            dim.info = dim.info,
            dim = as.integer(dim),
            intent.p1 = intent.p1,
            intent.p2 = intent.p2,
            intent.p3 = intent.p3,
            intent.code = intent.code,
            datatype = datatype,
            bitpix = bitpix,
            slice.start = slice.start,
            pixdim = pixdim,
            vox.offset = {if (is.nii == "nii") as.single(352) else as.single(0)},
            scl.slope = scl.slope,
            scl.inter = scl.inter,
            slice.end = slice.end,
            slice.code = slice.code,
            xyzt.units = xyzt.units,
            cal.max = cal.max,
            cal.min = cal.min,
            slice.duration = slice.duration,
            toffset = toffset,
            glmax = integer(1),
            glmin = integer(1),
            descrip = descrip,
            aux.file = aux.file,
            qform.code = qform.code,
            sform.code = sform.code,
            quatern.b = quatern.b,
            quatern.c = quatern.c,
            quatern.d = quatern.d,
            qoffset.x = qoffset.x,
            qoffset.y = qoffset.y,
            qoffset.z = qoffset.z,
            srow.x = srow.x,
            srow.y = srow.y,
            srow.z = srow.z,
            intent.name = intent.name,
            magic = magic) 
  return(l)
}



f.write.list.to.hdr.nifti <- function(L, file){

# To respect the length of some Nifti character fields
  strcomplete <- function(string,max.length) {
    string <- substr(string,1,max.length)
    as.character(paste(as.character(string),paste(rep(" ",max.length-nchar(string)),collapse=""),collapse="",sep=""))
  }
  

# Writes a list to a .hdr file in NIFTI format (and always in little-endian)
  a <- .C("write_nifti_header_wrap_JM",
          file,
          as.integer(L$sizeof.hdr),
          if (L$data.type != rawToChar(raw(10))) strcomplete(L$data.type,10) else rawToChar(raw(10)),
          if (L$db.name != rawToChar(raw(18))) strcomplete(L$db.name,18) else rawToChar(raw(18)),
          as.integer(L$extents),
          as.integer(L$session.error),
          if (L$regular != rawToChar(raw(1))) strcomplete(L$regular,1) else rawToChar(raw(1)),
          if (L$dim.info != rawToChar(raw(1))) strcomplete(L$dim.info,1) else rawToChar(raw(1)),
          as.integer(L$dim),
          as.single(L$intent.p1),
          as.single(L$intent.p2),
          as.single(L$intent.p3),
          as.integer(L$intent.code),
          as.integer(L$datatype),
          as.integer(L$bitpix),
          as.integer(L$slice.start),
          as.single(L$pixdim),
          as.single(L$vox.offset),
          as.single(L$scl.slope),
          as.single(L$scl.inter),
          as.integer(L$slice.end),
          if (L$slice.code != rawToChar(raw(1))) strcomplete(L$slice.code,1) else rawToChar(raw(1)),
          if (L$xyzt.units != rawToChar(raw(1))) strcomplete(L$xyzt.units,1) else rawToChar(raw(1)),
          as.single(L$cal.max),
          as.single(L$cal.min),
          as.single(L$slice.duration),
          as.single(L$toffset),
          as.integer(L$glmax),
          as.integer(L$glmin),
          if (L$descrip != rawToChar(raw(80))) strcomplete(L$descrip,80) else rawToChar(raw(80)),
          if (L$aux.file != rawToChar(raw(24))) strcomplete(L$aux.file,24) else rawToChar(raw(24)),
          as.integer(L$qform.code),
          as.integer(L$sform.code),
          as.single(L$quatern.b),
          as.single(L$quatern.c),
          as.single(L$quatern.d),
          as.single(L$qoffset.x),
          as.single(L$qoffset.y),
          as.single(L$qoffset.z),
          as.single(L$srow.x),
          as.single(L$srow.y),
          as.single(L$srow.z),
          if (L$intent.name != rawToChar(raw(16))) strcomplete(L$intent.name,16) else rawToChar(raw(16)),
          as.character(L$magic),
          PACKAGE="AnalyzeFMRI")
}






f.read.nifti.slice <- function(file, slice, tpt){
  #Reads in a .img file into an array


  is.nii <- substring(file, nchar(file) - 2,nchar(file)) # file extension

  file.name <- substring(file, 1, nchar(file) - 4) # file name without extension

  if (is.nii == "nii") {
    file.hdr <- paste(file.name, ".nii", sep = "")
    file.img <- paste(file.name, ".nii", sep = "")
    toadd <- 352 # because the image data begins at offset 0+348+4=352 in a .nii file
  } 
  else {
    file.hdr <- paste(file.name, ".hdr", sep = "")
    file.img <- paste(file.name, ".img", sep = "")
    toadd <- 0
    if(file.exists(file.img) == FALSE) stop(paste(file.img, "not found"))
  }


  hdr <- f.read.nifti.header(file.hdr)

  dim <- hdr$dim[2:3]

  num.data.pts <- dim[1] * dim[2]
  if(tpt < 1 || tpt > hdr$dim[5]) stop("tpt is not in range")
  if(slice < 1 || slice > hdr$dim[4]) stop("slice is not in range")
  
  offset <- (tpt - 1) * hdr$dim[2] * hdr$dim[3] * hdr$dim[4] + (slice - 1) * hdr$dim[2] * hdr$dim[3]

  if(hdr$datatype == 2){
    
    vol <- .C("readchar_v1_JM",
              mat = integer(num.data.pts),
              file.img,
              as.integer(hdr$swap),
              as.integer(num.data.pts),
              as.integer(offset * 1 + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
    vol <- array(vol$mat, dim = dim)
  #this works because array fills itself with the left most subscript moving fastest
  }
  if(hdr$datatype == 4){
    
    vol <- .C("read2byte_v1_JM",
              mat = integer(num.data.pts),
              file.img,
              as.integer(hdr$swap),
              as.integer(num.data.pts),
              as.integer(offset * 2 + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
    vol <- array(vol$mat, dim = dim)
#this works because array fills itself with the left most subscript moving fastest
  }
  if(hdr$datatype == 8){
    vol <- .C("read4byte_v1_JM",
              mat = integer(num.data.pts),
              file.img,
              as.integer(hdr$swap),
              as.integer(num.data.pts),
              as.integer(offset * 4 + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
    vol <- array(vol$mat, dim = dim)
    #this works because array fills itself with the left most subscript moving fastest
  }
  if(hdr$datatype == 16){
    vol <- .C("readfloat_v1_JM",
              mat = single(num.data.pts),
              file.img,
              as.integer(hdr$swap),
              as.integer(num.data.pts),
              as.integer(offset * 4 + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
    vol <- array(vol$mat, dim=dim)
#this works because array fills itself with the left most subscript moving fastest
  }
  if(hdr$datatype == 64){
    vol <- .C("readdouble_v1_JM",
              mat = numeric(num.data.pts),
              file.img,
              as.integer(hdr$swap),
              as.integer(num.data.pts),
              as.integer(offset * 8 + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
    vol <- array(vol$mat, dim = dim)
    #this works because array fills itself with the left most subscript moving fastest
  }
  
  if(hdr$datatype == 0 || hdr$datatype == 1 || hdr$datatype == 32 || hdr$datatype == 128 || hdr$datatype == 255) print(paste("The format", hdr$data.type, "is not supported yet. Please contact me if you want me to extend the functions to do this (marchini@stats.ox.ac.uk)"), quote=FALSE)
  
  return(vol)}



f.read.nifti.tpt <- function(file, tpt){
#Reads in one timepoint of a .img file into an array

  is.nii <- substring(file, nchar(file) - 2,nchar(file)) # file extension

  file.name <- substring(file, 1, nchar(file) - 4) # file name without extension

  if (is.nii == "nii") {
    file.hdr <- paste(file.name, ".nii", sep = "")
    file.img <- paste(file.name, ".nii", sep = "")
    toadd <- 352 # because the image data begins at offset 0+348+4=352 in a .nii file
  } 
  else {
    file.hdr <- paste(file.name, ".hdr", sep = "")
    file.img <- paste(file.name, ".img", sep = "")
    toadd <- 0
    if(file.exists(file.img) == FALSE) stop(paste(file.img, "not found"))
  }



  hdr <- f.read.nifti.header(file.hdr)

  dim <- hdr$dim[2:4]
  
  num.data.pts <- dim[1] * dim[2] * dim[3]
  if(tpt < 1 || tpt > hdr$dim[5]) stop("tpt is not in range")
  
  offset <- (tpt - 1) * hdr$dim[2] * hdr$dim[3] * hdr$dim[4]
  
  if(hdr$datatype == 2){
    
    vol <- .C("readchar_v1_JM",
                  mat = integer(num.data.pts),
              file.img,
              as.integer(hdr$swap),
              as.integer(num.data.pts),
              as.integer(offset * 1 + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
    vol <- array(vol$mat, dim = dim)
                                        #this works because array fills itself with the left most subscript moving fastest
  }
  if(hdr$datatype == 4){
    vol <- .C("read2byte_v1_JM",
              mat = integer(num.data.pts),
              file.img,
              as.integer(hdr$swap),
              as.integer(num.data.pts),
              as.integer(offset * 2 + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
    vol <- array(vol$mat, dim = dim)
    #this works because array fills itself with the left most subscript moving fastest
  }
  if(hdr$datatype == 8){
    vol <- .C("read4byte_v1_JM",
              mat = integer(num.data.pts),
              file.img,
              as.integer(hdr$swap),
              as.integer(num.data.pts),
              as.integer(offset * 4 + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
    vol <- array(vol$mat, dim = dim)
    #this works because array fills itself with the left most subscript moving fastest
  }
  if(hdr$datatype == 16){
    vol <- .C("readfloat_v1_JM",
              mat = single(num.data.pts),
              file.img,
              as.integer(hdr$swap),
              as.integer(num.data.pts),
              as.integer(offset * 4 + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
    vol <- array(vol$mat, dim = dim)
    #this works because array fills itself with the left most subscript moving fastest
  }
  if(hdr$datatype == 64){
    vol <- .C("readdouble_v1_JM",
              mat = numeric(num.data.pts),
              file.img,
              as.integer(hdr$swap),
              as.integer(num.data.pts),
              as.integer(offset * 8 + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
    vol <- array(vol$mat, dim = dim)
    #this works because array fills itself with the left most subscript moving fastest
  }
  
  if(hdr$datatype == 0 || hdr$datatype == 1 || hdr$datatype == 32 || hdr$datatype == 128 || hdr$datatype == 255) print(paste("The format", hdr$data.type, "is not supported yet. Please contact me if you want me to extend the functions to do this (marchini@stats.ox.ac.uk)"), quote = FALSE)
  
  return(vol)}



f.read.nifti.slice.at.all.timepoints <- function(file, slice){
  #Reads in a slice of a .img file at all time points into an array


  is.nii <- substring(file, nchar(file) - 2,nchar(file)) # file extension

  file.name <- substring(file, 1, nchar(file) - 4) # file name without extension

  if (is.nii == "nii") {
    file.hdr <- paste(file.name, ".nii", sep = "")
    file.img <- paste(file.name, ".nii", sep = "")
    toadd <- 352 # because the image data begins at offset 0+348+4=352 in a .nii file
  } 
  else {
    file.hdr <- paste(file.name, ".hdr", sep = "")
    file.img <- paste(file.name, ".img", sep = "")
    toadd <- 0
    if(file.exists(file.img) == FALSE) stop(paste(file.img, "not found"))
  }

  
  hdr <- f.read.nifti.header(file.hdr)

  dim <- hdr$dim[2:3]
  
  num.data.pts <- dim[1] * dim[2]
  if(slice < 1 || slice > hdr$dim[4]) stop("slice is not in range")
  
  vl <- array(0, dim = hdr$dim[c(2, 3, 5)])
  
  if(hdr$datatype == 2){
    for(i in 1:hdr$dim[5]){
      offset <- (i - 1) * hdr$dim[2] * hdr$dim[3] * hdr$dim[4] + (slice - 1) * hdr$dim[2] * hdr$dim[3]
      vol <- .C("readchar_v1_JM",
                mat = integer(num.data.pts),
                file.img,
                as.integer(hdr$swap),
                as.integer(num.data.pts),
                as.integer(offset * 1 + toadd),
                as.integer(1), PACKAGE="AnalyzeFMRI")
      vol <- array(vol$mat, dim = dim)
      vl[, , i] <- vol
    }
  }
  if(hdr$datatype == 4){
    for(i in 1:hdr$dim[5]){
      offset <- (i - 1) * hdr$dim[2] * hdr$dim[3] * hdr$dim[4] + (slice - 1) * hdr$dim[2] * hdr$dim[3]
      vol <- .C("read2byte_v1_JM",
                mat = integer(num.data.pts),
                file.img,
                as.integer(hdr$swap),
                as.integer(num.data.pts),
                as.integer(offset * 2 + toadd),
                as.integer(1), PACKAGE="AnalyzeFMRI")
      vol <- array(vol$mat, dim = dim)
      vl[, , i] <- vol
    }
  }
  
  if(hdr$datatype == 8){
    for(i in 1:hdr$dim[5]){
      offset <- (i - 1) * hdr$dim[2] * hdr$dim[3] * hdr$dim[4] + (slice - 1) * hdr$dim[2] * hdr$dim[3]
      vol <- .C("read4byte_v1_JM",
                mat = integer(num.data.pts),
                file.img,
                as.integer(hdr$swap),
                as.integer(num.data.pts),
                as.integer(offset * 4 + toadd),
                as.integer(1), PACKAGE="AnalyzeFMRI")
      vol <- array(vol$mat, dim = dim)
      vl[, , i] <- vol
    }
  }
  
  if(hdr$datatype == 16){
    for(i in 1:hdr$dim[5]){
      offset <- (i - 1) * hdr$dim[2] * hdr$dim[3] * hdr$dim[4] + (slice - 1) * hdr$dim[2] * hdr$dim[3]
      vol <- .C("readfloat_v1_JM",
                mat = single(num.data.pts),
                file.img,
                as.integer(hdr$swap),
                as.integer(num.data.pts),
                as.integer(offset * 4 + toadd),
                as.integer(1), PACKAGE="AnalyzeFMRI")
      vol <- array(vol$mat, dim = dim)
      vl[, , i] <- vol
    }
  }
  
  if(hdr$datatype == 64){
    for(i in 1:hdr$dim[5]){
      offset <- (i - 1) * hdr$dim[2] * hdr$dim[3] * hdr$dim[4] + (slice - 1) * hdr$dim[2] * hdr$dim[3]
      vol <- .C("readdouble_v1_JM",
                mat = numeric(num.data.pts),
                file.img,
                as.integer(hdr$swap),
                as.integer(num.data.pts),
                as.integer(offset * 8 + toadd),
                as.integer(1), PACKAGE="AnalyzeFMRI")
      vol <- array(vol$mat, dim = dim)
      vl[, , i] <- vol
    }}
  
  if(hdr$datatype == 0 || hdr$datatype == 1 || hdr$datatype == 32 || hdr$datatype == 128 || hdr$datatype == 255) print(paste("The format", hdr$data.type, "is not supported yet. Please contact me if you want me to extend the functions to do this (marchini@stats.ox.ac.uk)"), quote = FALSE)
  
  return(vl)}



f.read.nifti.ts <- function(file, x, y, z){
  #Reads in a .img file into an array

  is.nii <- substring(file, nchar(file) - 2,nchar(file)) # file extension

  file.name <- substring(file, 1, nchar(file) - 4) # file name without extension

  if (is.nii == "nii") {
    file.hdr <- paste(file.name, ".nii", sep = "")
    file.img <- paste(file.name, ".nii", sep = "")
    toadd <- 352 # because the image data begins at offset 0+348+4=352 in a .nii file
  } 
  else {
    file.hdr <- paste(file.name, ".hdr", sep = "")
    file.img <- paste(file.name, ".img", sep = "")
    toadd <- 0
    if(file.exists(file.img) == FALSE) stop(paste(file.img, "not found"))
  }


  
  hdr <- f.read.nifti.header(file.hdr)

  if(x < 1 || x > hdr$dim[2]) stop("x is not in range")
  if(y < 1 || y > hdr$dim[3]) stop("y is not in range")
  if(z < 1 || z > hdr$dim[4]) stop("z is not in range")
  
  offset.start <- (z - 1) * hdr$dim[2] * hdr$dim[3] + (y - 1) * hdr$dim[2] + (x - 1)
  offset.add <- hdr$dim[2] * hdr$dim[3] * hdr$dim[4]
  
  vol <- 1:hdr$dim[5]
  
  if(hdr$datatype == 2){
    for(i in 1:hdr$dim[5]){
      
      v <- .C("readchar_v1_JM",
              mat = integer(1),
              file.img,
              as.integer(hdr$swap),
              as.integer(1),
              as.integer(offset.start * 1 + 1 * (i - 1) * offset.add + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
      vol[i] <- v$mat}
  }
  
  if(hdr$datatype  == 4){
    for(i in 1:hdr$dim[5]){
      
      v <- .C("read2byte_v1_JM",
              mat = integer(1),
              file.img,
              as.integer(hdr$swap),
              as.integer(1),
              as.integer(offset.start * 2 + 2 * (i - 1) * offset.add + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
      vol[i] <- v$mat}
  }
  
  if(hdr$datatype  == 8){
    for(i in 1:hdr$dim[5]){
      
      v <- .C("read4byte_v1_JM",
              mat = integer(1),
              file.img,
              as.integer(hdr$swap),
              as.integer(1),
              as.integer(offset.start * 4 + 4 * (i - 1) * offset.add + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
      vol[i] <- v$mat}
  }
  
  if(hdr$datatype == 16){
    for(i in 1:hdr$dim[5]){
      
      v <- .C("readfloat_v1_JM",
              mat = single(1),
              file.img,
              as.integer(hdr$swap),
              as.integer(1),
              as.integer(offset.start * 4 + 4 * (i - 1) * offset.add + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
      vol[i] <- v$mat}
  }
  
  if(hdr$datatype == 64){
    for(i in 1:hdr$dim[5]){
      
      v <- .C("readdouble_v1_JM",
              mat = numeric(1),
              file.img,
              as.integer(hdr$swap),
              as.integer(1),
              as.integer(offset.start * 8 + 8 * (i - 1) * offset.add + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
      vol[i] <- v$mat}
  }
  
  if(hdr$datatype == 0 || hdr$datatype == 1 || hdr$datatype == 32 || hdr$datatype == 128 || hdr$datatype == 255) print(paste("The format", hdr$data.type, "is not supported yet. Please contact me if you want me to extend the functions to do this (marchini@stats.ox.ac.uk)"), quote = FALSE)
  
  return(vol)}


f.read.nifti.volume <- function(file){
  #Reads in a .img file into an array

  is.nii <- substring(file, nchar(file) - 2,nchar(file)) # file extension

  file.name <- substring(file, 1, nchar(file) - 4) # file name without extension

  if (is.nii == "nii") {
    file.hdr <- paste(file.name, ".nii", sep = "")
    file.img <- paste(file.name, ".nii", sep = "")
    toadd <- 352 # because the image data begins at offset 0+348+4=352 in a .nii file
  } 
  else {
    file.hdr <- paste(file.name, ".hdr", sep = "")
    file.img <- paste(file.name, ".img", sep = "")
    toadd <- 0
    if(file.exists(file.img) == FALSE) stop(paste(file.img, "not found"))
  }



  
  hdr <- f.read.nifti.header(file.hdr)
  num.dim <- hdr$dim[1]
  dim <- hdr$dim[1:num.dim + 1]
  
  if (num.dim < 4) dim <- c(dim,1)

  num.data.pts <- prod(dim)  
  
  if(hdr$datatype == 2){
    vol <- .C("readchar_v1_JM",
              mat = integer(num.data.pts),
              file.img,
              as.integer(hdr$swap),
              as.integer(num.data.pts),
              as.integer(0 + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
    vol <- array(vol$mat, dim = dim)
    #this works because array fills itself with the left most subscript moving fastest
  }
  if(hdr$datatype == 4){
    vol <- .C("read2byte_v1_JM",
              mat = integer(num.data.pts),
              file.img,
              as.integer(hdr$swap),
              as.integer(num.data.pts),
              as.integer(0 + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
    vol <- array(vol$mat, dim = dim)
#this works because array fills itself with the left most subscript moving fastest
  }
  if(hdr$datatype == 8){
    vol <- .C("read4byte_v1_JM",
              mat = integer(num.data.pts),
              file.img,
              as.integer(hdr$swap),
              as.integer(num.data.pts),
              as.integer(0 + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
    vol <- array(vol$mat, dim = dim)
#this works because array fills itself with the left most subscript moving fastest
  }
  if(hdr$datatype == 16){
    vol <- .C("readfloat_v1_JM",
              mat = single(num.data.pts),
              file.img,
              as.integer(hdr$swap),
              as.integer(num.data.pts),
              as.integer(0 + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
    vol <- array(vol$mat, dim = dim)
    #this works because array fills itself with the left most subscript moving fastest
  }
  if(hdr$datatype == 64){
    vol <- .C("readdouble_v1_JM",
              mat = numeric(num.data.pts),
              file.img,
              as.integer(hdr$swap),
              as.integer(num.data.pts),
              as.integer(0 + toadd),
              as.integer(1), PACKAGE="AnalyzeFMRI")
    vol <- array(vol$mat, dim = dim)
    #this works because array fills itself with the left most subscript moving fastest
  }

  if(hdr$datatype == 0 || hdr$datatype == 1 || hdr$datatype == 32 || hdr$datatype == 128 || hdr$datatype == 255) print(paste("The format", hdr$data.type, "is not supported yet. Please contact me if you want me to extend the functions to do this (marchini@stats.ox.ac.uk)"), quote = FALSE)
  
  return(vol)}




f.spectral.summary.nifti <- function(file, mask.file, ret.flag = FALSE)
{
  #for a NIFTI .img file the periodogram of the time series are divided by a flat spectral estimate using the median periodogram ordinate. The resulting values are then combined within each Fourier frequency and quantiles are plotted against freequency. This provides a fast look at a fMRI dataset to identify any artefacts that reside at single frequencies.

########################
#get info about dataset
########################
  is.nii <- substring(file, nchar(file) - 2,nchar(file)) # file extension

  file.name <- substring(file, 1, nchar(file) - 4) # file name without extension

  if (is.nii == "nii") {
    file.hdr <- paste(file.name, ".nii", sep = "")
    file.img <- paste(file.name, ".nii", sep = "")
    toadd <- 352 # because the image data begins at offset 0+348+4=352 in a .nii file
  } 
  else {
    file.hdr <- paste(file.name, ".hdr", sep = "")
    file.img <- paste(file.name, ".img", sep = "")
    toadd <- 0
    if(file.exists(file.img) == FALSE) stop(paste(file.img, "not found"))
  }


  hdr <- f.read.nifti.header(file.hdr)
  nsl <- hdr$dim[4]
  nim <- hdr$dim[5]
  pxdim <- hdr$pixdim[2:4]

#####################
#read in/create mask
#####################

  f.mask.create <- function(dat, pct = .1, slices = c(0)) {
  #function that creates a mask for an fMRI dataset by thresholding the mean of the pixel time series at a percentage point of the maximum intensity of the dataset
    is.nii <- substring(file, nchar(file) - 2,nchar(file)) # file extension
    
    file.name <- substring(file, 1, nchar(file) - 4) # file name without extension
    
    if (is.nii == "nii") {
      file.hdr <- paste(file.name, ".nii", sep = "")
      file.img <- paste(file.name, ".nii", sep = "")
      toadd <- 352 # because the image data begins at offset 0+348+4=352 in a .nii file
    } 
    else {
      file.hdr <- paste(file.name, ".hdr", sep = "")
    file.img <- paste(file.name, ".img", sep = "")
      toadd <- 0
      if(file.exists(file.img) == FALSE) stop(paste(file.img, "not found"))
    }
    
    hdr <- f.read.nifti.header(file.hdr)
    nsl <- hdr$dim[4]
    xc <- hdr$dim[2]
    yc <- hdr$dim[3]
    if(slices[1] == 0){slices <- seq(1, nsl)}
    mask <- array(0, dim = c(xc, yc, length(slices)))
    
    max.int <- 0
    for(k in 1:length(slices)){
      slice <- f.read.nifti.slice.at.all.timepoints(dat$file, slices[k])
      if(max(slice)>max.int){max.int <- max(slice)}
    }
    
    for(k in 1:length(slices)){
      slice <- f.read.nifti.slice.at.all.timepoints(dat$file, slices[k])
      
      for(i in 1:(xc * yc)){
        a <- (i - 1) %/% xc + 1
        b <- i - (a - 1) * xc
        mask[b, a, k] <- mean(slice[b, a, ])
        
        if(mask[b, a, k] >= (pct * max.int)){mask[b, a, k] <- 1}
        else{mask[b, a, k] <- 0}
      }
    }
    
    return(mask)
    
  }

  if(mask.file!=FALSE){mask <- f.read.nifti.volume(mask.file)}
  else{
    dat <- list(file = file, mask.file = mask.file)
    mask <- f.mask.create(dat = dat)}
  dim(mask) <- hdr$dim[2:4]

###############
#set constants
###############

  n <- floor(nim / 2) + 1


##########################
#initialise storage arrays
##########################

  res <- array(NA, dim = c(dim(mask), n))

#####################
#main evaluation loop
#####################
  cat("Processing slices...")
  for(l in 1:nsl){
    cat(" [", l, "]", sep = "")
    
    slice <- f.read.nifti.slice.at.all.timepoints(file, l)
    
    
    for(i in 1:(dim(slice)[1] * dim(slice)[2])){
      a <- (i - 1) %/% dim(slice)[1] + 1
      b <- i - (a - 1) * dim(slice)[1]
      if(mask[b, a, l] == 1){
        t <- Mod(fft(slice[b, a, ]) / sqrt(2 * pi * nim))[1:n]
        s <- median(t)
        res[b, a, l, ] <- t / s
      }
    }
  }
  cat("\n")
  
  b <- apply(res, 4, FUN = quantile, probs = seq(.5, 1, .05), na.rm = TRUE)
  plot(c(0, n - 1), c(0, 30), type = "n", xlab = "", ylab = "", axes = FALSE)
  axis(1, at = seq(0, n, 5))
  axis(2, at = seq(0, 30, 5))
  for(i in 0:(n - 1)){
    points(rep(i, 11), b[, i + 1])}
  if(ret.flag == TRUE)  return(b)
}




f.read.header <- function(file){
  #This function reads in the header information from a NIFTI (.nii or .hdr) or ANALYZE (.hdr) file depending on the magic field

	is.nii <- substring(file, nchar(file) - 2,nchar(file))

	file.name <- substring(file, 1, nchar(file) - 4)

	if (is.nii == "nii") {file.hdr <- paste(file.name, ".nii", sep = "")} else {file.hdr <- paste(file.name, ".hdr", sep = "")}

	if(file.exists(file.hdr) == FALSE) stop(paste(file.hdr, "not found"))

#Detect whether the data is big or little endian. The first part of a .hdr file is the size of the file which is always a C int (i.e. 4 bytes) and always has value 348. Therefore trying to read it in assuming little-endian will tell you if that is the correct mode

    swap <- 0

    if(.C("swaptest_wrap_JM",
          ans = integer(1),
          file.hdr,
          PACKAGE="AnalyzeFMRI")$ans != 348) # $ans is sizeof_hdr
        swap <- 1
    
    a<-substr(.C("read_nifti_magic_wrap",
          file.hdr,  
          as.integer(swap), 
	  magic =  paste(rep(" ", 4), sep = "", collapse = ""),
	PACKAGE="AnalyzeFMRI")$magic,start=1,stop=2)

if (a == "ni" | a == "n+") {res <- f.read.nifti.header(file.hdr)}
else if (a == "") {res <- f.read.analyze.header(file.hdr)}
else {stop("Problem in your magic field!")}

return(res)

}

f.read.volume <- function(file){
  #This function reads in the volume image file using header information from a NIFTI (.nii or .hdr) or ANALYZE (.hdr) file depending on the magic field

	is.nii <- substring(file, nchar(file) - 2,nchar(file))

	file.name <- substring(file, 1, nchar(file) - 4)

	if (is.nii == "nii") {file.hdr <- paste(file.name, ".nii", sep = "")} else {file.hdr <- paste(file.name, ".hdr", sep = "")}

	if(file.exists(file.hdr) == FALSE) stop(paste(file.hdr, "not found"))

#Detect whether the data is big or little endian. The first part of a .hdr file is the size of the file which is always a C int (i.e. 4 bytes) and always has value 348. Therefore trying to read it in assuming little-endian will tell you if that is the correct mode

    swap <- 0

    if(.C("swaptest_wrap_JM",
          ans = integer(1),
          file.hdr,
          PACKAGE="AnalyzeFMRI")$ans != 348) # $ans is sizeof_hdr
        swap <- 1
    
    a<-substr(.C("read_nifti_magic_wrap",
          file.hdr,  
          as.integer(swap), 
	  magic =  paste(rep(" ", 4), sep = "", collapse = ""),
	PACKAGE="AnalyzeFMRI")$magic,start=1,stop=2)

if (a == "ni" | a == "n+") {res <- f.read.nifti.volume(file.hdr)}
else if (a == "") {res <- f.read.analyze.volume(file.hdr)}
else {stop("Problem in your magic field!")}

return(res)

}





f.write.nifti <- function(mat, file, size = "float", L = NULL, nii = FALSE){
# Creates a NIFTI .img/.hdr pair of files or a .nii file from a given array


  if(max(mat) == "NA") stop("NA values in array not allowed. Files not written.")

  extension <- substring(file, nchar(file) - 2,nchar(file))
  is.nii <- extension
  
  if (extension == "nii" | extension == "img" | extension == "hdr") {

    if (is.nii == "nii") {
      if (!nii) stop("If you want to create a .nii file, you also shoud put argument nii to TRUE")
      file.hdr <- file
      file.img <- file
    }
    else {
      file.name <- substring(file, 1, nchar(file) - 4)
      file.hdr <- paste(file.name, ".hdr", sep = "")
      file.img <- paste(file.name, ".img", sep = "")
    }
  }

  else {

    if (nii) {
      file.hdr <- paste(file, ".nii", sep = "")
      file.img <- paste(file, ".nii", sep = "")
      
    }
    else {
      file.hdr <- paste(file, ".hdr", sep = "")
      file.img <- paste(file, ".img", sep = "")
    }
  } 
  
  if (is.null(L)) {

    L <- f.basic.hdr.nifti.list.create(dim(mat), file.hdr)
  }

  if (L$datatype == 16) size <- "float"
  if (L$datatype == 4) size <- "int"
  if (L$datatype == 2) size <- "char"

  
  if(size == "float"){
    L$datatype <- 16
    L$bitpix <- 32
    L$data.type <- "float"
    if (nii) {
      f.write.nii.array.to.img.float(mat, L, file.img)
    }
    else {
      f.write.array.to.img.float(mat, file.img)
      f.write.list.to.hdr.nifti(L, file.hdr)
    }
  }
  
  if(size == "int"){
    if(max(mat[!is.nan(mat)])>32767 || min(mat[!is.nan(mat)]) < ( -32768)) stop("Values are outside integer range. Files not written.")
    L$datatype <- 4
    L$bitpix <- 16
    L$data.type <- "signed sho" # signed short
    if (nii) {
      f.write.nii.array.to.img.2bytes(mat, L, file.img)
    }
    else {
      f.write.array.to.img.2bytes(mat, file.img)
      f.write.list.to.hdr.nifti(L, file.hdr)
    }
  }
  
  if(size == "char"){
    if(max(mat[!is.nan(mat)])>255 || min(mat[!is.nan(mat)]) < 0) stop("Values are outside integer range. Files not written.")
    L$datatype <- 2
    L$bitpix <- 8
    L$data.type <- "unsignchar" # unsigned char
    if (nii) {
      f.write.nii.array.to.img.8bit(mat, L, file.img)
    }
    else {
      f.write.array.to.img.8bit(mat, file.img)
      f.write.list.to.hdr.nifti(L, file.hdr)
    }
  }
  
}


f.write.nii.array.to.img.2bytes <- function(mat, L, file){
  #writes an array into a .img file of 2 byte integers
  # and add at the begining of the file the NIFTI header part

  f.write.list.to.hdr.nifti(L, file)
  
  dm <- dim(mat)
  dm.ln <- length(dm)
  num.data.pts <- prod(dm)
  
  null <- .C("write2byteappend_JM",
     as.integer(mat),
     file,
     as.integer(num.data.pts),NAOK=TRUE, PACKAGE="AnalyzeFMRI")
  # NAOK=TRUE is necessary here to be able to pass NaN values. This is important because SPM uses NaN values as markers for voxels for which it has not calculated any statistics.
  
}

f.write.nii.array.to.img.8bit <- function(mat, L, file){
#writes an array into a .img file of 8 bit (1 byte) integers
  # and add at the begining of the file the NIFTI header part

  f.write.list.to.hdr.nifti(L, file)

  dm <- dim(mat)
  dm.ln <- length(dm)
  num.data.pts <- prod(dm)
  
  null <- .C("write8bitappend_JM",
     as.integer(mat),
     file,
     as.integer(num.data.pts),NAOK=TRUE, PACKAGE="AnalyzeFMRI")
  # NAOK=TRUE is necessary here to be able to pass NaN values. This is important because SPM uses NaN values as markers for voxels for which it has not calculated any statistics.
  
}



f.write.nii.array.to.img.float <- function(mat, L, file){
  #writes an array into a .img file of 4 byte flotas
  # and add at the begining of the file the NIFTI header part

  f.write.list.to.hdr.nifti(L, file)

  dm <- dim(mat)
  dm.ln <- length(dm)
  num.data.pts <- prod(dm)
  
  null <- .C("writefloatappend_JM",
     as.single(mat),
     file,
     as.integer(num.data.pts),NAOK=TRUE,PACKAGE="AnalyzeFMRI")
  # NAOK=TRUE is necessary here to be able to pass NaN values. This is important because SPM uses NaN values as markers for voxels for which it has not calculated any statistics.
  
}

diminfo2fps <- function(dim.info) {

  # Extract freq_dim, phase_dim and slice_dim fields from the one byte dim.info field
  
  z1 <- z2 <- z3 <- raw(8)
  z3[3:4] <- as.raw(1)
  z2[5:6] <- as.raw(1)
  z1[7:8] <- as.raw(1)

  freq.dim <- as.integer(packBits(rev(z1 & rev(rawToBits(charToRaw(dim.info))))))
  phase.dim <- as.integer(rawShift(packBits(rev(z2 & rev(rawToBits(charToRaw(dim.info))))),-2))
  slice.dim <- as.integer(rawShift(packBits(rev(z3 & rev(rawToBits(charToRaw(dim.info))))),-4))
  
  return(list(freq.dim=freq.dim,phase.dim=phase.dim,slice.dim=slice.dim))
  
}


fps2diminfo <- function(freq.dim,phase.dim,slice.dim) {

  # Encode freq_dim, phase_dim and slice_dim fields into the one byte dim.info field
  
  res <- raw(8)

  if (length(freq.dim) != 0) {
    if (freq.dim == 1) res[8] <- as.raw(1)
    if (freq.dim == 2) res[7] <- as.raw(1)
    if (freq.dim == 3) res[7] <- res[8] <- as.raw(1)
  }
  
  if (length(phase.dim) != 0) {
    if (phase.dim == 1) res[6] <- as.raw(1)
    if (phase.dim == 2) res[5] <- as.raw(1)
    if (phase.dim == 3) res[5] <- res[6] <- as.raw(1)
  }
  
  if (length(slice.dim) != 0) {
    if (slice.dim == 1) res[4] <- as.raw(1)
    if (slice.dim == 2) res[3] <- as.raw(1)
    if (slice.dim == 3) res[3] <- res[4] <- as.raw(1)
  }
  
  return(list(dim.info=rawToChar(packBits(rev(res)))))
  
}

xyzt2st <- function(xyzt.units) {

  # Extract space and time dimensions fields from the one byte xyzt.units field

  z1 <- z2 <- raw(8)
  z2[3:5] <- as.raw(1)
  z1[6:8] <- as.raw(1)

  space <- as.integer(packBits(rev(z1 & rev(rawToBits(charToRaw((xyzt.units)))))))
  pr.space <- if (space == 0) "unknown" else if (space == 1) "meters" else if (space == 2) "millimeters" else if (space == 3) "micrometers"
  time <- as.integer(packBits(rev(z2 & rev(rawToBits(charToRaw((xyzt.units)))))))
  pr.time <- if (time == 8) "seconds" else if (time == 16) "milliseconds" else if (time == 24) "microseconds" else if (time == 32) "Hertz" else if (time == 40) "ppm" else if (time == 48) "radians per second"
  cat(paste("space: ",pr.space,"\t","time: ",pr.time,"\n"))
  
  return(list(space=space,time=time))

}


st2xyzt <- function(space,time) {

  # Encode space and time dimensions fields into the one byte xyzt.units field

  res <- raw(8)

  if (space == 1) res[8] <- as.raw(1)
  if (space == 2) res[7] <- as.raw(1)
  if (space == 3) res[7] <- res[8] <- as.raw(1)
 
  if (time == 8) res[5] <- as.raw(1)
  if (time == 16) res[4] <- as.raw(1)
  if (time == 24) res[4] <- res[5] <- as.raw(1)

  if (time == 32) res[3] <- as.raw(1)
  if (time == 40) res[3] <- res[5] <- as.raw(1)
  if (time == 48) res[3] <- res[4] <- as.raw(1)
  
  return(list(xyzt.units=rawToChar(packBits(rev(res)))))

  
}

magicfield <- function(file) {
#  Determine the type of a file : NIFIT .nii format, NIFTI .hdr/.img pair format, Analyze format

  hdr <- f.read.header(file)
  magic <- hdr$magic
  if (is.null(magic)) magic <- ""
  dim <- hdr$dim[5]

  
  if (substr(magic,start=1,stop=2) == "ni") {
    cat(paste("NIFTI ",magic," - .hdr/.img pair with ",dim," image(s)\n",sep=""))
  }
  
  else if (substr(magic,start=1,stop=2) == "n+") {
    cat(paste("NIFTI ",magic," - one .nii file with ",dim," image(s)\n",sep=""))
  }
  
  else if (magic == "") {
    cat(paste("ANALYZE .hdr/.img pair with ",dim," image(s)\n",sep=""))
  }
  else {
    stop("Problem in your magic field!")
  }

  return(list(magic=magic,dim=dim))
}


analyze2nifti <- function(file.in,path.in=".",path.out=".",file.out=NULL,is.nii=TRUE,qform.code=2,sform.code=2,data.type=rawToChar(raw(10)),db.name=rawToChar(raw(18)),dim.info=rawToChar(raw(1)),dim=NULL,TR=0,slice.code=rawToChar(raw(1)),xyzt.units=rawToChar(raw(1)),descrip=NULL,aux.file=rawToChar(raw(24)),intent.name=rawToChar(raw(16))) {
# Passage du format Analyze 4D (resp. 3D) vers le format Nifti 4D (resp. 3D)
# C'est un équivalent de la fonction SPM5: spm_write_vol

# [q/s]form.code values :
# [q/s]form.code=0 : Arbitrary coordinates (Method 1).
# [q/s]form.code=1 : Scanner-based anatomical coordinates
# [q/s]form.code=2 : Coordinates aligned to another file's, or to anatomical "truth".
# [q/s]form.code=3 : Coordinates aligned to Talairach-Tournoux Atlas; (0,0,0)=AC, etc.
# [q/s]form.code=4 : MNI 152 normalized coordinates. 

path.in <- if (substr(path.in,nchar(path.in),nchar(path.in)) != "/") paste(path.in,"/",sep="") else path.in
path.out <- if (substr(path.out,nchar(path.out),nchar(path.out)) != "/") paste(path.out,"/",sep="") else path.out


removeext <- function(filename) {res <- unlist(strsplit(x=filename,split="\\.")) ; n <- max(2,length(res)) ; if (nchar(res[n]) == 3) res <- paste(res[-n],collapse=".") else res <- filename ; return(res)}

if (is.null(file.out)) file.out <- removeext(file.in)

# Lecture des images ...
data <- f.read.volume(paste(path.in,file.in,".img",sep=""))
header.4D <- f.read.header(paste(path.in,file.in,".img",sep=""))

# ... puis traduction du header de l'Analyze vers le header du Nifti (en y intégrant les infos du .mat)
dim.mat <- header.4D$dim[2:5]
if (is.nii) {file <- paste(file.out,".nii",sep="")} else {file <- paste(file.out,".hdr",sep="")}
Lnifti <- f.basic.hdr.nifti.list.create(dim.mat,file)
# Lnifti$sizeof.hdr is already set to 348 in the f.basic.hdr.nifti.list.create function
Lnifti$data.type <- data.type # UNUSED
Lnifti$db.name <- db.name # UNUSED
# Lnifti$extents <- header.4D$extents # UNUSED
# Lnifti$session.error <- header.4D$session.error # UNUSED
Lnifti$regular <- header.4D$regular # UNUSED in NIFTI-1 but i put this to "r" to be conformant with SPM
Lnifti$dim.info <- dim.info # en fait, ici il devrait s'agir de l'information suivante codée sur 1 byte : (freq_dim,phase_dim,slice_dim)
                                        # chacune des 3 composantes étant codée sur 2 bits et pouvant prendre les valeurs 1, 2 ou 3 (pour x,y ou z respectivement)
                                        # These fields encode which spatial dimension (1,2, or 3) corresponds to which acquisition dimension for MRI data.

if (!is.null(dim)) Lnifti$dim <- dim
# Lnifti$dim a été créé par la fonction f.basic.hdr.nifti.list.create mais il devrait pouvoir être modifié en cas de besoin
# Uses of dimensions 6 (dim[7]) and 7 (dim[8]) are also specified in NIFTI-1 format.
# Pour l'instant    dim <- c(length(dim.mat), dim.mat, rep(0, 7 - length(dim.mat)))
# Mais on peut encore changer la 5eme dimension (dim[6]) : dimension 5 is for storing multiple values at each spatiotemporal voxel.
# Si dim[1]=5 et dim[6]>1 alors : The 5th dimension of the dataset contains multiple values (e.g., a vector) to be stored at each spatiotemporal location.
# Il faudra dans ce cas, probablement aussi modifier intent_code (et possiblement intent_p1, intent_p2, and intent_p3.)


Lnifti$intent.p1 <- as.single(0) # on devrait pouvoir changer cela
Lnifti$intent.p2 <- as.single(0) # on devrait pouvoir changer cela
Lnifti$intent.p3 <- as.single(0) # on devrait pouvoir changer cela
Lnifti$intent.code <- as.integer(0) # on devrait pouvoir changer cela

Lnifti$datatype <- as.integer(header.4D$datatype)
Lnifti$bitpix <- as.integer(header.4D$bitpix)
Lnifti$slice.start <- header.4D$dim.un0  # Indicates the start of the slice acquisition pattern, when slice_code is nonzero.  These values are present to allow
                                         # for the possible addition of "padded" slices at either end of the volume, which don't fit into the slice timing pattern.
                                         # If there are no padding slices, then slice_start=0 and slice_end=dim[slice_dim]-1 are the correct values.
                                         # For these values to be meaningful, slice_start must be non-negative and slice_end must be greater than slice_start.
Lnifti$pixdim <- header.4D$pixdim # Note : je le remodifie plus bas ...
Lnifti$pixdim[5] <- TR
Lnifti$vox.offset <- { if (is.nii) as.single(352) else as.single(0) }
Lnifti$scl.slope <- 1 # Voir le code de la fonction spm_write_vol qui calcule cette valeur qui peut etre differente de 1 !!!
Lnifti$scl.inter <- 0
Lnifti$slice.end <- as.integer(header.4D$funused3) # Indicates the end of the slice acquisition pattern, when slice_code is nonzero.  These values are present to allow
                                         # for the possible addition of "padded" slices at either end of the volume, which don't fit into the slice timing pattern.
                                         # If there are no padding slices, then slice_start=0 and slice_end=dim[slice_dim]-1 are the correct values.
                                         # For these values to be meaningful, slice_start must be non-negative and slice_end must be greater than slice_start.

Lnifti$slice.code <- slice.code   # If this is nonzero, AND if slice_dim is nonzero, AND if slice_duration is positive, indicates the timing pattern of the slice acquisition.
                          # The following codes are defined:
                          #                   NIFTI_SLICE_UNKNOWN  0
                          #                   NIFTI_SLICE_SEQ_INC  1
                          #                   NIFTI_SLICE_SEQ_DEC  2
                          #                   NIFTI_SLICE_ALT_INC  2
                          #                   NIFTI_SLICE_ALT_DEC  4

Lnifti$xyzt.units <- xyzt.units # SPM change ce code mais je ne sais pas ou il le prend, a voir avec Cecile !!!!
Lnifti$cal.max <- header.4D$cal.max
Lnifti$cal.min <- header.4D$cal.min
Lnifti$slice.duration <-  header.4D$compressed # If this is positive, AND if slice_dim is nonzero, indicates the amount of time used to acquire 1 slice.
                                               # slice_duration*dim[slice_dim] can be less than pixdim[4] with a clustered acquisition method, for example.
Lnifti$toffset <- header.4D$verified
# Lnifti$glmax <- header.4D$glmax # UNUSED
# Lnifti$glmin <- header.4D$glmin # UNUSED
Lnifti$descrip <- if (is.null(descrip)) header.4D$descrip else descrip
Lnifti$aux.file <- aux.file

# Il reste à changer les derniers champs grâce à la lecture du .mat (s'il existe!!!)

if (!is.na(file.info(paste(path.in,file.in,".mat",sep=""))[1])) { # le .mat existe vraiment

  if (qform.code == 0) qform.code <- 2 # pas sur ici!! 
  if (sform.code == 0) sform.code <- 2 # pas sur ici!!

##require("R.matlab")
M <- readMat(paste(path.in,file.in,".mat",sep=""))$M
M[1,] <- - M[1,] # Je ne sais pas s'il faut lire le $M ou le $mat . Dès fois le $mat n'est pas présent. La première ligne du $mat est (-1)*$M[1,]
############################ ATTENTION!!! Gros problème ici car les fichiers 1801_* ont tous des .mat différents!! Il
# faut s'assurer d'avoir le bon .mat pour le fichier d'entrée (surtout s'il a été obtenu en regroupant plusieurs fichiers 3D avec des .mat differents. Il faudra alors
# utiliser un reslice des fichiers Analyze d'origine pour qu'ils aient tous le même .mat avant de les merger avec threeDto4D.R
# prevoir de faire verifier tout cela par la fonction analyze2nifti


#     | R_xx R_xy R_xz T_x |
#     | R_yx R_yy R_yz T_y |
# M = | R_zx R_zy R_zz T_z |
#     | S_x  S_y  S_z   1  | 

#Translation and shearing can be described by each single vectors.
#If you want them in matrix form just take the translation/shear
#submatrix and fill the rest with the identity matrix (all 0
#except the diagonal).

#Rotation and Scaling are both described by the same 3x3
#submatrix. A rotation matrix has the property of being
#orthogonal. And if the scale is 1, then it is orthonormal. So to
#determine the scale you just have to determine the normalisation
#factor(s) of the rotation matrix.

#The rotational axis is the eigenvector of R (there's one and only
#one eigenvector if R is a rotation matrix - you might get two,
#but they're just antiparallel and have the same eigenvalue).


#An affine transformation is of the form Y = M*X + T, where X is the
#3x1 input, Y is the 3x1 output, M is a 3x3 matrix, and T is a 3x1 vector.
#T is the translation and is easy to identify in a 4x4 homogenous matrix.
#The essence of your question, though, is how to factor M into rotation
#and scales. The best you can do is the "singular value decomposition".
#You can factor M = L*S*R, where L and R are orthogonal matrices and
#S is a diagonal matrix whose diagonal entries are nonnegative. The
#factorization of M is generally not unique. Consider M = I, the identity
#matrix. In this case, S = I, R is any rotation matrix, and L is the inverse
#of R.

#It is not always possible to factor M = S*R or M = L*S.

#The "polar decomposition" factors M = A*P, where P is orthogonal and
#A is symmetric. You can think of A as "scales" but with respect to
#some coordinate system, not necessarily the one M applies to. Any
#symmetric matrix can be decomposed using eigenvalues/eigenvectors,
#A = L*S*Inverse(L), where L is orthogonal and S is a diagonal matrix.
#Notice that M = A*P = L*S*Inverse(L)*P = L*S*R, where
#R = Inverse(L)*P. Thus, the polar decomposition and singular value
#decomposition are related.

#The additional twist that 'fungus' mentions is when M is itself a rotation
#matrix and you want to factor it into products of coordinate-axis
#rotations. This is the topic of "Euler angles". Such a factorization
#requires you to specify the order of the rotations, and even here it is
#possible there is not a unique factorization. 



# Translation matrix = | 1 0 0 T_x |
#                      | 0 1 0 T_y |
#                      | 0 0 1 T_z |
#                      | 0 0 0  1  |

# Shear matrix = | 1   0   0   0 |  Pas sur de celle-là ????
#                | 0   1   0   0 |
#                | 0   0   1   0 |
#                | S_x S_y S_z 1 | 


translate111 <- diag(4) ; translate111[,4] <- 1 # translation de +1 cran en x, +1 cran en y et +1 cran en z
M <- M %*% translate111 # Convert from first voxel at (1,1,1) to first voxel at (0,0,0)  (SPM compte les voxels à partir de (1,1,1) mais NIFTI à partir de (0,0,0))

# Translations
Tmat <- M[1:3,4]

# Rotations and zooms
R <- M[1:3,1:3]
vx <- sqrt(apply(R^2,FUN=sum,MARGIN=2))
vx[vx==0] <- 1
R <- R %*% diag(1/vx)

# Ensure that R is O(3)
res.svd <- svd(R)
R <- res.svd$u %*% t(res.svd$v)



Lnifti$pixdim[2:4] <- vx

# see below for these two fields
#Lnifti$qform.code <- as.integer(qform.code)
#Lnifti$sform.code <- as.integer(sform.code)




if (det(R) > 0) qfac <- 1 else qfac <- -1 # A voir comment on fixe la valeur de qfac ... je pense que c'est avec le déterminant de R qui devrait être égal à 1, sinon ...???

R <- R %*% diag(c(1,1,qfac))
  
  
Lnifti$pixdim[1] <- qfac 




# Voir les fichiers htm que j'ai sur les quaternions ...



Q <- R2Q(R) # translate rotation matrix into quaternions


Lnifti$quatern.b <- as.single(Q[1])  
Lnifti$quatern.c <- as.single(Q[2]) 
Lnifti$quatern.d <- as.single(Q[3])  


Lnifti$qoffset.x <- as.single(Tmat[1])  
Lnifti$qoffset.y <- as.single(Tmat[2])  
Lnifti$qoffset.z <- as.single(Tmat[3])  
Lnifti$srow.x <- as.single(M[1,])
Lnifti$srow.y <- as.single(M[2,])
Lnifti$srow.z <- as.single(M[3,])

}

Lnifti$intent.name <- intent.name # name or meaning of data

magic <- raw(4)
if (is.nii) {magic[1:3] <- as.raw(charToRaw("n+1")) # cas du NIFTI.nii (one file)
           } else {magic[1:3] <- as.raw(charToRaw("ni1"))} # cas du NIFTI en paire}
Lnifti$magic <- rawToChar(magic)

Lnifti$qform.code <- as.integer(qform.code)
Lnifti$sform.code <- as.integer(sform.code)


# Ecriture du nifti .hdr/.img ou .nii
f.write.nifti(mat=data,file=paste(path.out,file.out,sep=""),size="int",L=Lnifti,nii=is.nii)

}


threeDto4D <- function(outputfile,path.in=NULL,prefix=NULL,regexp=NULL,times=NULL,list.of.in.files=NULL,path.out=NULL,is.nii.pair=FALSE,hdr.number=1) {


# Description:
# ------------
# To read tm functionnal images in ANALYZE or NIFTI format, and concatenate them to obtain one (hdr/img pair) 4D image file in Analyze or Nifti format which is written on disk

# Arguments:
# ----------
# outputfile: character. Name of the outputfile without extension
# path.in:   character with the path to the directory containing the image files
# prefix: character. common prefix to each file
# regexp: character. Regular expression to get all the files
# times:  vector. numbers of the image files to retrieve
# list.of.in.files: names of img files to concatenate (with full path)
# path.out: where to write the output hdr/img pair files. Will be taken as path.in if not provided.
  
# Values:
# -------
# None
  
# Example:
# --------
# path.fonc <- "/network/home/lafayep/Stage/Data/map284/functional/MondrianApril2007/preprocessing/1801/smoothed/"
# threeDto4D("essai",path.in=path.fonc,prefix="su1801_",regexp="????.img",times=1:120)


if (is.null(path.out)) path.out <- path.in

if (is.null(list.of.in.files)) {
if (is.null(path.in)) stop("Argument path.in must be provided")
if (is.null(prefix)) stop("Argument prefix must be provided")
if (is.null(regexp)) stop("Argument regexp must be provided")
if (is.null(times)) stop("Argument times must be provided")
list.of.in.files <- list.files(path=path.in,pattern=glob2rx(paste(prefix,regexp,sep="")))[times]
}
else {
times <- length(list.of.in.files)
}


# Reading of hdr from the hdr.number-th file
L <- f.read.header(paste(path.in,list.of.in.files[hdr.number],sep=""))


# L$dim[1] contient 3 ou 4 suivant qu'il y a un temps ou pas en L$dim[5]

if (L$dim[1] == 4 & L$dim[5] != 1) stop("This function should be used to read image files with time dimension not greater than 1")

# We create the final hdr list with correct time information
if (is.null(list.of.in.files)) {
 L$dim[5] <- length(times) 
}
else {
 L$dim[5] <- length(list.of.in.files)
}


L$dim[1] <- 4 # because we create 4D files

# We check if the headers of all files are the same (but a few unimportant fields: file.name, db.name, descrip)
if (is.null(L$magic)) rem <- c(1,5,28) else rem <- c(1,5,31)
if (length(times) > 1) {
  for (filename in list.of.in.files) {  
    file <- paste(path.in,filename,sep="")
    Ltemp <- f.read.header(file)
    Ltemp$dim[1] <- 4 # because we create 4D files
    somme <- 0
    for (i in (1:45)[-rem]) somme <- somme + sum(L[[i]] != Ltemp[[i]]) 
    if (somme != 1) stop("hdr part should be the same for all files")
  }
}

if (substr(path.out,nchar(path.out),nchar(path.out)) != "/") path.out <- paste(path.out,"/",sep="")

removeext <- function(filename) {res <- unlist(strsplit(x=filename,split="\\.")) ; n <- max(2,length(res)) ; if (nchar(res[n]) == 3) res <- paste(res[-n],collapse=".") else res <- filename ; return(res)}

# L$magic gives image format 
# NULL : Analyze hdr/img pair
# ni1  : NIFTI hdr/img pair
# n+1  : single NIFTI nii file
if (is.null(L$magic)) {
  f.write.list.to.hdr(L,paste(path.out,outputfile,".hdr",sep=""))
  if (file.exists(paste(path.in,removeext(list.of.in.files[1]),".mat",sep=""))) file.copy(from=paste(path.in,removeext(list.of.in.files[1]),".mat",sep=""), to=paste(path.out,outputfile,".mat",sep=""), overwrite = FALSE)
  
#  M <- readMat(paste(path.in,removeext(list.of.in.files[1]),".mat",sep=""))
#  writeMat(M,paste(path.out,outputfile,".mat",sep=""))
} else {
  if (is.nii.pair) {
    f.write.list.to.hdr.nifti(L,paste(path.out,outputfile,".hdr",sep=""))
  } else {
    f.write.list.to.hdr.nifti(L,paste(path.out,outputfile,".nii",sep=""))
  }
}

if (is.nii.pair | is.null(L$magic)) {
  # Opening of the img file to be written
  if (file.exists(paste(path.out,outputfile,".img",sep=""))) file.remove(paste(path.out,outputfile,".img",sep=""))
  con1 <- file(description = paste(path.out,outputfile,".img",sep=""), open = "ab")

} else {
  # Opening of the nii file to be written
#  if (file.exists(paste(path.out,outputfile,".nii",sep=""))) file.remove(paste(path.out,outputfile,".nii",sep=""))
  con1 <- file(description = paste(path.out,outputfile,".nii",sep=""), open = "ab")
  writeBin(raw(4),con=con1) # the 4 bytes (extension field) after the first 348 bytes of the header
}

# To force writing of the 4 bytes
close(con1)


swap <- if (is.nii.pair | is.null(L$magic)) f.read.header(paste(path.out,outputfile,".hdr",sep=""))$swap else f.read.nifti.header(paste(path.out,outputfile,".nii",sep=""))$swap

if (swap != L$swap) {
#  si le swap du header ecrit est different de l'ancien swap on reecrit le hdr et on change le swap

  if (is.nii.pair | is.null(L$magic)) {
  # Opening of the img file to be written
  tmp <- readBin(con=paste(path.out,outputfile,".hdr",sep=""),what="raw",n=348,size=1)
  con2 <- file(description = paste(path.out,outputfile,".hdr",sep=""), open = "wb")
  writeBin(tmp[(2:1)+rep((0:173)*2,each=2)],con=con2)

} else {
  # Opening of the nii file to be written
  tmp <- readBin(con=paste(path.out,outputfile,".nii",sep=""),what="raw",n=352,size=1)
  con2 <- file(description = paste(path.out,outputfile,".nii",sep=""), open = "wb")
  writeBin(tmp[(2:1)+rep((0:175)*2,each=2)],con=con2)
}

close(con2)

}


if (is.nii.pair | is.null(L$magic)) {
  # Reopening of the img file to be written
  con1 <- file(description = paste(path.out,outputfile,".img",sep=""), open = "ab")

} else {
  # Reopening of the nii file to be written
  con1 <- file(description = paste(path.out,outputfile,".nii",sep=""), open = "ab")
}

datatype <- L$datatype
if (datatype == 0) stop("datatype unknown!")
else if (datatype == 2) datatype <- 1
else if (datatype == 4) datatype <- 2
else if (datatype == 8) datatype <- 4
else if (datatype == 16) datatype <- 4
else if (datatype == 32) datatype <- 8
else if (datatype == 64) datatype <- 8
else if (datatype == 128) datatype <- 3
else if (datatype == 256) datatype <- 1
else if (datatype == 512) datatype <- 2
else if (datatype == 768) datatype <- 4
else if (datatype == 1024) datatype <- 8
else if (datatype == 1280) datatype <- 8
else if (datatype == 1536) datatype <- 16
else if (datatype == 1792) datatype <- 16
else if (datatype == 2048) datatype <- 32
else stop("datatype not in the list of codes given by NIFTI team!")


fourD <- c()
for (filename in list.of.in.files) {
  
  file <- paste(path.in,filename,sep="")
  
  if (is.null(L$magic) | (L[[45]] == "ni1")) {
    fourD <- readBin(con=file,what="raw",n=prod(L$dim[2:4])*datatype,size=1)
  }
  else if (L$magic == "n+1") { # we remove the header present at the begining of each nii file being read
    fourD <- readBin(con=file,what="raw",n=prod(L$dim[2:4])*datatype+352,size=1)[-(1:352)]
  }
  
  writeBin(fourD,con=con1)
  
  
  gc(FALSE)
}

close(con1)

}


twoDto4D <- function(x.2d, dim) {

# Description:
# ------------
# This function transform a 2D matrix of size tm x vm containing images in each row into a 4D array image
  
# Arguments:
# ----------
# x.2d: a 2D matrix to be transformed
# dim: vector of length 4 containing the dimensions of the array. dim[1:3] are the space dimensions. dim[4] is the time dimension
  
# Values:
# -------
# volume.4d: a 4D array image

# Example:
# --------


volume.4d <- array(t(x.2d),dim=dim)

return(volume.4d=volume.4d)

}

fourDto2D <- function(volume.4d, tm) {

# Description:
# ------------
# This function transform a 4D image in a 2D image matrix
  
# Arguments:
# ----------
# volume.4d: a 4D array to be transformed
# tm: number of time dimensions
  
# Values:
# -------
# x.2d:  matrix of size tm x vm which contains the tm images

# Example:
# --------


x.2d <- matrix(volume.4d,nrow=tm,byrow=TRUE)

return(x.2d=x.2d)

}

ijk2xyz <- function(ijk=c(1,1,1),method=2,L) {

#   There are 3 different methods by which continuous coordinates can be
#   attached to voxels.  The discussion below emphasizes 3D volumes, and
#   the continuous coordinates are referred to as (x,y,z).  The voxel
#   index coordinates (i.e., the array indexes) are referred to as (i,j,k),
#   with valid ranges:
#     i = 0 .. dim[1]-1
#     j = 0 .. dim[2]-1  (if dim[0] >= 2)
#     k = 0 .. dim[3]-1  (if dim[0] >= 3)
#   The (x,y,z) coordinates refer to the CENTER of a voxel.  In methods
#   2 and 3, the (x,y,z) axes refer to a subject-based coordinate system,
#   with
#     +x = Right  +y = Anterior  +z = Superior.
#   This is a right-handed coordinate system.  However, the exact direction
#   these axes point with respect to the subject depends on qform_code
#   (Method 2) and sform_code (Method 3).

#   N.B.: The i index varies most rapidly, j index next, k index slowest.
#    Thus, voxel (i,j,k) is stored starting at location
#      (i + j*dim[1] + k*dim[1]*dim[2]) * (bitpix/8)
#    into the dataset array.

#   N.B.: The ANALYZE 7.5 coordinate system is
#      +x = Left  +y = Anterior  +z = Superior
#    which is a left-handed coordinate system.  This backwardness is
#    too difficult to tolerate, so this NIFTI-1 standard specifies the
#    coordinate order which is most common in functional neuroimaging.

#   N.B.: The 3 methods below all give the locations of the voxel centers
#    in the (x,y,z) coordinate system.  In many cases, programs will wish
#    to display image data on some other grid.  In such a case, the program
#    will need to convert its desired (x,y,z) values into (i,j,k) values
#    in order to extract (or interpolate) the image data.  This operation
#    would be done with the inverse transformation to those described below.

#   N.B.: Method 2 uses a factor 'qfac' which is either -1 or 1; qfac is
#    stored in the otherwise unused pixdim[0].  If pixdim[0]=0.0 (which
#    should not occur), we take qfac=1.  Of course, pixdim[0] is only used
#    when reading a NIFTI-1 header, not when reading an ANALYZE 7.5 header.

#   N.B.: The units of (x,y,z) can be specified using the xyzt_units field.


    if (method != 1 & method != 2 & method != 3) {stop("method should be 1, 2 or 3")}

  # ijk is a matrix. Each column of ijk should contain a voxel index coordinates (i,j,k) to be mapped to its (x,y,z) real coordinates in some other space
  
    ijk <- as.matrix(ijk)
    
    ijk <- ijk - 1 # because the formulas given in nifti1.h are valid for voxel data indices beginning at 0 (as in C or C++) and not at 1 (as in R)
    
    nbpts <- ncol(ijk)
    
    if (method == 1) {
      
      xyz <- diag(L$pixdim[2:4]) %*% ijk
      
    }
    
    if (method == 2) {
      
      if (L$qform.code == 0) stop("Method 2 should be used for qform.code > 0")
      
      qfac <- L$pixdim[1]
      b <- L$quatern.b
      c <- L$quatern.c
      d <- L$quatern.d
      a <- sqrt(1.0-(b*b+c*c+d*d))
      R <- matrix(c(a*a+b*b-c*c-d*d,2*b*c-2*a*d,2*b*d+2*a*c,2*b*c+2*a*d,a*a+c*c-b*b-d*d,2*c*d-2*a*b,2*b*d-2*a*c,2*c*d+2*a*b,a*a+d*d-c*c-b*b),byrow=TRUE,nrow=3,ncol=3)
      xyz <- R %*% diag(c(L$pixdim[2:3],qfac*L$pixdim[4])) %*% ijk + c(L$qoffset.x,L$qoffset.y,L$qoffset.z)
      
    }
    
    if (method == 3) {
      
      if (L$sform.code == 0) stop("Method 3 should be used for sform.code > 0")
      
      M <- matrix(c(L$srow.x,L$srow.y,L$srow.z),byrow=TRUE,nrow=3,ncol=4)
      xyz <- M %*% rbind(ijk,rep(1,nbpts))
      
    }
    
    return(list(xyz=xyz))
    
  }


xyz2ijk <- function(xyz=c(1,1,1),method=2,L) {

#   There are 3 different methods by which continuous coordinates can be
#   attached to voxels.  The discussion below emphasizes 3D volumes, and
#   the continuous coordinates are referred to as (x,y,z).  The voxel
#   index coordinates (i.e., the array indexes) are referred to as (i,j,k),
#   with valid ranges:
#     i = 0 .. dim[1]-1
#     j = 0 .. dim[2]-1  (if dim[0] >= 2)
#     k = 0 .. dim[3]-1  (if dim[0] >= 3)
#   The (x,y,z) coordinates refer to the CENTER of a voxel.  In methods
#   2 and 3, the (x,y,z) axes refer to a subject-based coordinate system,
#   with
#     +x = Right  +y = Anterior  +z = Superior.
#   This is a right-handed coordinate system.  However, the exact direction
#   these axes point with respect to the subject depends on qform_code
#   (Method 2) and sform_code (Method 3).

#   N.B.: The i index varies most rapidly, j index next, k index slowest.
#    Thus, voxel (i,j,k) is stored starting at location
#      (i + j*dim[1] + k*dim[1]*dim[2]) * (bitpix/8)
#    into the dataset array.

#   N.B.: The ANALYZE 7.5 coordinate system is
#      +x = Left  +y = Anterior  +z = Superior
#    which is a left-handed coordinate system.  This backwardness is
#    too difficult to tolerate, so this NIFTI-1 standard specifies the
#    coordinate order which is most common in functional neuroimaging.

#   N.B.: The 3 methods below all give the locations of the voxel centers
#    in the (x,y,z) coordinate system.  In many cases, programs will wish
#    to display image data on some other grid.  In such a case, the program
#    will need to convert its desired (x,y,z) values into (i,j,k) values
#    in order to extract (or interpolate) the image data.  This operation
#    would be done with the inverse transformation to those described below.

#   N.B.: Method 2 uses a factor 'qfac' which is either -1 or 1; qfac is
#    stored in the otherwise unused pixdim[0].  If pixdim[0]=0.0 (which
#    should not occur), we take qfac=1.  Of course, pixdim[0] is only used
#    when reading a NIFTI-1 header, not when reading an ANALYZE 7.5 header.

#   N.B.: The units of (x,y,z) can be specified using the xyzt_units field.


    if (method != 1 & method != 2 & method != 3) {stop("method should be 1, 2 or 3")}

  # xyz is a matrix. Each column of xyz should contain a real set of coordinates (x,y,z) to be mapped to its (i,j,k) voxel index coordinates in the dataset
  
    xyz <- as.matrix(xyz)
    
    nbpts <- ncol(xyz)
    
    if (method == 1) {
      
      ijk <- diag(1/L$pixdim[2:4]) %*% xyz
      
    }
    
    if (method == 2) {
      
      if (L$qform.code == 0) stop("Method 2 should be used for qform.code > 0")
      
      qfac <- L$pixdim[1]
      b <- L$quatern.b
      c <- L$quatern.c
      d <- L$quatern.d
      a <- sqrt(1.0-(b*b+c*c+d*d))
      R <- matrix(c(a*a+b*b-c*c-d*d,2*b*c-2*a*d,2*b*d+2*a*c,2*b*c+2*a*d,a*a+c*c-b*b-d*d,2*c*d-2*a*b,2*b*d-2*a*c,2*c*d+2*a*b,a*a+d*d-c*c-b*b),byrow=TRUE,nrow=3,ncol=3)
      ijk <- diag(1/c(L$pixdim[2:3],qfac*L$pixdim[4])) %*% t(R) %*% (xyz - c(L$qoffset.x,L$qoffset.y,L$qoffset.z))  
      
    }
    
    if (method == 3) {
      
      if (L$sform.code == 0) stop("Method 3 should be used for sform.code > 0")
      
      M <- matrix(c(L$srow.x,L$srow.y,L$srow.z,0,0,0,1),byrow=TRUE,nrow=4,ncol=4)
      Minv <- matrix(NA,nrow=4,ncol=4)
        
      r11 <- M[1,1]
      r12 <- M[1,2]
      r13 <- M[1,3]
      r21 <- M[2,1]
      r22 <- M[2,2]
      r23 <- M[2,3]
      r31 <- M[3,1]
      r32 <- M[3,2]
      r33 <- M[3,3]
      v1 <- M[1,4]
      v2 <- M[2,4]
      v3 <- M[3,4]

      deti <- r11*r22*r33-r11*r32*r23-r21*r12*r33+r21*r32*r13+r31*r12*r23-r31*r22*r13 

      if ( deti != 0.0 ) deti <- 1.0 / deti 

      Minv[1,1] <- deti*( r22*r33-r32*r23)
      Minv[1,2] <- deti*(-r12*r33+r32*r13)
      Minv[1,3] <- deti*( r12*r23-r22*r13)
      Minv[1,4] <-  deti*(-r12*r23*v3+r12*v2*r33+r22*r13*v3-r22*v1*r33-r32*r13*v2+r32*v1*r23)

      Minv[2,1] <-  deti*(-r21*r33+r31*r23)
      Minv[2,2] <-  deti*( r11*r33-r31*r13)
      Minv[2,3] <-  deti*(-r11*r23+r21*r13)
      Minv[2,4] <-  deti*( r11*r23*v3-r11*v2*r33-r21*r13*v3+r21*v1*r33+r31*r13*v2-r31*v1*r23)

      Minv[3,1] <-  deti*( r21*r32-r31*r22)
      Minv[3,2] <-  deti*(-r11*r32+r31*r12)
      Minv[3,3] <-  deti*( r11*r22-r21*r12)
      Minv[3,4] <-  deti*(-r11*r22*v3+r11*r32*v2+r21*r12*v3-r21*r32*v1-r31*r12*v2+r31*r22*v1)

      Minv[4,1] <- Minv[4,2] <- Minv[4,3] <- 0.01
      Minv[4,4] <- if (deti == 0.0)  0.0 else 1.0

      
      ijk <- Minv %*% rbind(xyz,rep(1,nbpts))
      
    }
    
    ijk <- ijk + 1 # because the formulas given in nifti1.h are valid for voxel data indices beginning at 0 (as in C or C++) and not at 1 (as in R)
    
    return(list(ijk=ijk))
    
  }



R2Q <- function(R,qfac=NULL) {

# Convert from rotation matrix to quaternion form

  if (is.null(qfac) && det(R)<0) {
    qfac <- -1
    R[, 3] <- qfac * R[, 3]
  }
  R <- R[1:3,1:3]
  d <- diag(R)
  traceval <- sum(d) + 1
  if (traceval > 0.5) {
    s <- sqrt(traceval)*2;
    Q <- c((R[3,2]-R[2,3])/s,(R[1,3]-R[3,1])/s,(R[2,1]-R[1,2])/s,0.25*s)
  } else {
    traceval <- which.max(d)
    switch(traceval,
           {s <- 2*sqrt(1 + R[1,1] - R[2,2] - R[3,3]); Q <- c(0.25*s,(R[1,2]+R[2,1])/s,(R[3,1]+R[1,3])/s,(R[3,2]-R[2,3])/s)},
           {s <- 2*sqrt(1 + R[2,2] - R[1,1] - R[3,3]); Q <- c((R[1,2]+R[2,1])/s,0.25*s,(R[2,3]+R[3,2])/s,(R[1,3]-R[3,1])/s)},
           {s <- 2*sqrt(1 + R[3,3] - R[1,1] - R[2,2]); Q <- c((R[3,1]+R[1,3])/s,(R[2,3]+R[3,2])/s,0.25*s,(R[2,1]-R[1,2])/s)})
  }
  
  if (Q[4]<0) Q <- -Q
  return(Q)


#   The DICOM attribute (0020,0037) "Image Orientation (Patient)" gives the
#   orientation of the x- and y-axes of the image data in terms of 2 3-vectors.
#   The first vector is a unit vector along the x-axis, and the second is
#   along the y-axis.  If the (0020,0037) attribute is extracted into the
#   value (xa,xb,xc,ya,yb,yc), then the first two columns of the R matrix
#   would be
#              [ -xa  -ya ]
#              [ -xb  -yb ]
#              [  xc   yc ]
#   The negations are because DICOM's x- and y-axes are reversed relative
#   to NIFTI's.  The third column of the R matrix gives the direction of
#   displacement (relative to the subject) along the slice-wise direction.
#   This orientation is not encoded in the DICOM standard in a simple way;
#   DICOM is mostly concerned with 2D images.  The third column of R will be
#   either the cross-product of the first 2 columns or its negative.
#   a=(a1,a2,a3) ; b=(b1,b2,b3) ; a x b = (a2b3-a3b2,a3b1-a1b3,a1b2-a2b1)
#   Here is the third column of R: ±(-xb*yc+xc*yb,-xc*ya+xa*yc,xa*yb-xb*ya)
#   It is possible to infer the sign of the 3rd column by examining the coordinates
#   in DICOM attribute (0020,0032) "Image Position (Patient)" for successive
#   slices.  However, this method occasionally fails for reasons that I
#   (RW Cox) do not understand.
 
  
}


Q2R <- function(Q,qfac) {

# Generate a rotation matrix from a quaternion q=a+ib+jc+kd,
# where Q = [b c d], and a = 1-b^2-c^2-d^2.

  Q <- Q[1:3] 
  a <- sqrt(1 - sum(Q^2))
  b <- Q[1]
  c <- Q[2]
  d <- Q[3]
  if (a<10^(-7)) {
    a <- 1/sqrt(b*b+c*c+d*d)
    b <- b*a
    c <- c*a
    d <- d*a
    a <- 0
  }
  bb <- b*b; cc <- c*c; dd <- d*d; aa <- a*a; bc <- b*c; bd <- b*d; ab <- a*b; cd <- c*d; ac <- a*c; ad <- a*d;
  R <- matrix(c(1-2*cc-2*dd,2*(bc-ad),2*(bd+ac),2*(bc+ad),1-2*bb-2*dd,2*(cd-ab),2*(bd-ac),2*(cd+ab),1-2*cc-2*bb),nrow=3,ncol=3,byrow=TRUE)
  R[,3] <- qfac * R[,3]
  
  return(R)
  
}

nifti.quatern.to.mat44 <- function(L) {


  qb <- L$quatern.b
  qc <- L$quatern.c
  qd <- L$quatern.d

  qx <- L$qoffset.x
  qy <- L$qoffset.y
  qz <- L$qoffset.z

  dx <- L$pixdim[2]
  dy <- L$pixdim[3]
  dz <- L$pixdim[4]

  qfac <- L$pixdim[1]


  b <- qb; c <- qc; d <- qd

  R <- matrix(NA,nrow=4,ncol=4)
  
# last row is always [ 0 0 0 1 ] 

   R[4,1] <- R[4,2] <- R[4,3] <- 0.0 ; R[4,4] <- 1.0

# compute a parameter from b,c,d 

   a <- 1.0 - (b*b + c*c + d*d) 
   if ( a < 10^(-7) ) {                   # special case 
     a <- 1.0 / sqrt(b*b+c*c+d*d) 
     b <- b*a ; c <- c*a ; d <- d*a         # normalize (b,c,d) vector 
     a <- 0.0                         # a = 0 ==> 180 degree rotation 
   } else {
     a <- sqrt(a)                      # angle = 2*arccos(a) 
   }


# load rotation matrix, including scaling factors for voxel sizes 

   xd <- if (dx > 0.0) dx else 1.0        # make sure are positive 
   yd <- if (dy > 0.0) dy else 1.0 
   zd <- if (dz > 0.0) dz else 1.0 


  if ( qfac < 0.0 ) zd <- -zd          # left handedness? 

   R[1,1] <- (a*a+b*b-c*c-d*d) * xd 
   R[1,2] <- 2.0 * (b*c-a*d) * yd 
   R[1,3] <- 2.0 * (b*d+a*c) * zd 
   R[2,1] <- 2.0 * (b*c+a*d) * xd 
   R[2,2] <- (a*a+c*c-b*b-d*d) * yd 
   R[2,3] <- 2.0 * (c*d-a*b) * zd 
   R[3,1] <- 2.0 * (b*d-a*c) * xd 
   R[3,2] <- 2.0 * (c*d+a*b) * yd 
   R[3,3] <- (a*a+d*d-c*c-b*b) * zd 

# load offsets 

   R[1,4] <- qx ; R[2,4] <- qy ; R[3,4] <- qz 

   return(R )
 

}

mat34.to.TRSZ <- function(M) {
# Voir la fonction decompose_aff dans le fichier miscmaths.cc de FSL

if (nrow(M) == 3) M <- rbind(M,c(0,0,0,1))


# decomposes using the convention: mat = transl * rotmat * skew * scale
# order of parameters is 3 rotation + 3 translation + 3 scales + 3 skews
# angles are in radians

centre <- as.matrix(rep(0,3))
transl <- M[1:3,1:3]%*%centre+M[1:3,4]-centre
translation <- diag(rep(1,4))
translation[1:3,4] <- transl

  
x <- M[1:3,1]
y <- M[1:3,2]
z <- M[1:3,3]

sx <- sqrt(sum(x^2))
sy <- sqrt(sum(y^2)-(sum(x*y))^2/(sx)^2)
a <- sum(x*y)/(sx*sy)
x0 <- x/sx
y0 <- y/sy-a*x0
sz <- sqrt(sum(z^2)-(sum(x0*z))^2-(sum(y0*z))^2)

scales <- c(sx,sy,sz)

b <- sum(x0*z)/sz
c <- sum(y0*z)/sz

skews <- c(a,b,c)

skew <- matrix(c(1,a,b,0, 0,1,c,0, 0,0,1,0, 0,0,0,1),byrow=TRUE,nrow=4)


aff3 <- M[1:3,1:3]

rotmat <- diag(rep(1,4))

rotmat[1:3,1:3] <- aff3 %*% solve(diag(scales)) %*% solve(skew[1:3,1:3])

M <- rotmat

if (det(M) < 0) {

  warning("Rotation Rot is improper")
  Ref <- diag(c(-1,1,1,1))
  M <- Ref%*%M

} else {

  Ref <- diag(rep(1,4))

}

# Get Y Rotation and check that we aren't up a creek without a paddle
sy <- -M[3,1]
if (abs(sy) == 1) stop("cos Y = 0. Not solved yet")
# C'est compliqué. Ca fait appel a des equations du 3eme, 4eme et 5eme degré!

if (sy == 0) {

 cy <- -1

} else {

 cy <- -abs(cos(asin(sy)))

}

cx <- M[3,3]/cy
sx <- M[3,2]/cy
cz <- M[1,1]/cy
sz <- M[2,1]/cy


RotX <- diag(rep(1,4))
RotX[2,2] <- cx
RotX[3,2] <- sx
RotX[2,3] <- -sx
RotX[3,3] <- cx

RotY <- diag(rep(1,4))
RotY[1,1] <- cy
RotY[1,3] <- sy 
RotY[3,1] <- -sy 
RotY[3,3] <- cy

RotZ <- diag(rep(1,4))
RotZ[1,1] <- cz
RotZ[1,2] <- -sz
RotZ[2,1] <- sz
RotZ[2,2] <- cz


return(list(T=translation,Z=diag(c(scales,1)),S=skew,Rot=rotmat,RotX=RotX,RotY=RotY,RotZ=RotZ,Ref=Ref))

}


mat34.to.TZSR <- function(M) {
# Voir le fichier xfmdecomp.pl

# decomposes using the convention: mat = transl * scale * skew * rotation (rotation=Rz*Ry*Rx*Ref where Ref 
# is Reflexion if rotation is improper or Identity if rotation is proper)
# order of parameters is 3 rotation + 3 translation + 3 scales + 3 skews
# angles are in radians


if (nrow(M) == 3) M <- rbind(M,c(0,0,0,1))
  

# TRANSLATIONS : [M] = [T][Z][S][R]
# As of yet I am assuming the center of rotation is (0,0,0) as
# I am not exactly sure as to what SPM does here.
centre <- as.matrix(rep(0,3))
transl <- M[1:3,1:3]%*%centre+M[1:3,4]-centre
translation <- diag(rep(1,4))
translation[1:3,4] <- transl

  
# SCALES : [M] = inv[T][T][Z][S][R] = [Z][S][R]
# Here we use an identical method to Louis Collins's in mni_autoreg
# Namely multiply a unit vector in each direction and measure the length
# after the transformation.


M <- solve(translation) %*% M
Minv <- solve(M)

Sx <- sqrt(sum((Minv %*% as.matrix(c(1,0,0,0)))^2))
Sy <- sqrt(sum((Minv %*% as.matrix(c(0,1,0,0)))^2))
Sz <- sqrt(sum((Minv %*% as.matrix(c(0,0,1,0)))^2))

sx <- 1/Sx
sy <- 1/Sy
sz <- 1/Sz

scales <- c(sx,sy,sz)



# SHEARS : [M] = inv[T][T]inv[Z][Z][S][R] = [S][R]
# We assume the shear matrix:      S  [ 1 0 0 0 ]
# where x' = x                        [ a 1 0 0 ]
#       y' = ax + y                   [ b c 1 0 ]
#       z' = bx + cy + z              [ 0 0 0 1 ]
# 
# However as M at this point is in fact [S][R]
# we can extract a, b and c as such:
# 
# let [ M1 ]
#     [ M2 ]  =  [S][R]
#     [ M3 ] 
#
# thus:
# 
# a = (M2 . M1) / |M1|^2
# b = (M3 . M1) / |M1|^2
# c = (M3 . T)  / |T|^2     where T = M2 - (a . M1)
#
# We could also use the determinant to determine whether we have an 
# Orthogonal matrix and thus don;t have shears, but we don't do this yet....

M <- solve(diag(c(scales,1))) %*% M

a <- sum(M[2,]*M[1,]) / sum(M[1,]*M[1,])
b <- sum(M[3,]*M[1,]) / sum(M[1,]*M[1,])
tmp <- M[2,] - a*M[1,]
c <- sum(M[3,]*tmp) / sum(tmp*tmp)

skews <- c(a,b,c)

skew <- matrix(c(1,a,b,0, 0,1,c,0, 0,0,1,0, 0,0,0,1),byrow=FALSE,nrow=4)


# ROTATIONS : [M] = inv[T][T]inv[Z][Z]inv[S][S][R] = [R]
# We assume cy is positive to ensure we get one of the 2 possible solutions
# where rotations are between -pi and pi.
#
# Here we deduce Rx, Ry and Rz by virtue of the fact that the rotation
# matrix is as follows. 
#
# R = [ cos(Ry)*cos(Rz)  sin(Rx)*sin(Ry)*cos(Rz)-cos(Rx)*sin(Rz)  cos(Rx)*sin(Ry)*cos(Rz)+sin(Rx)*sin(Rz)  0 ]
#     [ cos(Ry)*sin(Rz)  sin(Rx)*sin(Ry)*sin(Rz)+cos(Rx)*cos(Rz)  cos(Rx)*sin(Ry)*sin(Rz)-sin(Rx)*cos(Rz)  0 ]
#     [ -sin(Ry)                    sin(Rx)*cos(Ry)                           cos(Rx)*cos(Ry)              0 ]
#     [ 0                                  0                                         0                     1 ]
#
# Then the quadrant of the angle must be deduced by the sign of
# cos and sin for the particular rotation.

M <- solve(skew) %*% M

Rot <- M

if (det(M) < 0) {

  warning("Rotation Rot is improper")
  Ref <- diag(c(-1,1,1,1))
  M <- Ref%*%M

} else {

  Ref <- diag(rep(1,4))

}

# Get Y Rotation and check that we aren't up a creek without a paddle
sy <- -M[3,1]
if (abs(sy) == 1) stop("cos Y = 0. Not solved yet") 
# C'est compliqué. Ca fait appel a des equations du 3eme, 4eme et 5eme degré!

if (sy == 0) {

 cy <- -1

} else {

 cy <- -abs(cos(asin(sy)))

}

cx <- M[3,3]/cy
sx <- M[3,2]/cy
cz <- M[1,1]/cy
sz <- M[2,1]/cy


RotX <- diag(rep(1,4))
RotX[2,2] <- cx
RotX[3,2] <- sx
RotX[2,3] <- -sx
RotX[3,3] <- cx

RotY <- diag(rep(1,4))
RotY[1,1] <- cy
RotY[1,3] <- sy 
RotY[3,1] <- -sy 
RotY[3,3] <- cy

RotZ <- diag(rep(1,4))
RotZ[1,1] <- cz
RotZ[1,2] <- -sz
RotZ[2,1] <- sz
RotZ[2,2] <- cz


return(list(T=translation,Z=diag(c(scales,1)),S=skew,Rot=Rot,RotX=RotX,RotY=RotY,RotZ=RotZ,Ref=Ref))

}


f.icast.fmri.gui <- function(){

  #starts GUI that allows user apply Spatial or Temporal ICA to an fMRI dataset

  path <- path.package(package = "AnalyzeFMRI")
  path.gui <- paste(path, "ICAst.gui.R", sep = .Platform$file.sep)
  source(path.gui)
}

f.icast.fmri <- function(foncfile,maskfile,is.spatial,n.comp.compute=TRUE,n.comp=0,hp.filter=TRUE) {

    #function for performing Spatial or Temporal ICA on an fMRI dataset
    #The function of Marchini avoids reading the dataset into R to minimise the memory used : see what can be done here !!


if ((!n.comp.compute) & (n.comp<1)) stop("You should provide n.comp")
  
#################################
# Lecture des données
#################################


volume.fonc <- f.read.volume(foncfile)

fonc.hdr <- f.read.header(foncfile)


nb.scans <- fonc.hdr$dim[5]
gc(FALSE)
volume.fonc <- fourDto2D(volume.fonc,nb.scans) # On passe en 2D (matrix of size tm x vm) car l'ACP et l'ACI fonctionnent avec des matrices 2D (le dépliage est arbitraire!)
gc(FALSE)


dimensions <- fonc.hdr$dim[2:4]
pixdim <- fonc.hdr$pixdim[2:4]
vox.units <- fonc.hdr$vox.units
cal.units <- fonc.hdr$cal.units

cat("Data have been read\n")
gc(FALSE)


#############################
# Pré-traitements
#############################

# masquage
# --------


mask.img <- f.read.volume(maskfile)
masque <- as.vector(fourDto2D(mask.img,tm=1))  # the vector of 0 and 1 values containing the mask (of length vm)
mask.img <- as.vector(mask.img)
aa <- (1:length(masque))[masque == 1]
X.masked <- volume.fonc[,aa] # matrix of size tm x length(aa)
dim.resmask <- length(masque)

gc(FALSE)


# centrage
# --------

if (is.spatial) {
Xcentred <- centering(X.masked,col.first=TRUE)$Xcentred  # Cas spatial
gc(FALSE)
} else {
Xcentred <- centering(X.masked,col.first=FALSE)$Xcentred  # Cas temporel
gc(FALSE)
}


# Phase des ondelettes (pour enlever la tendance ... à voir ...)
# require(waveslim)
#Xcentred <- apply(Xcentred[1:(nrow(Xcentred)-nrow(Xcentred)%%2),],MARGIN=1,myfunc<-function(x){dwt(x)[[3]]})
#gc(FALSE)


if (n.comp.compute) {

# réduction
# ---------

if (is.spatial) {
Xcr <- reduction(Xcentred,row.red=TRUE)$Xred  # Cas spatial
gc(FALSE)
} else {
Xcr <- reduction(Xcentred,row.red=FALSE)$Xred  # Cas temporel
gc(FALSE)
}



# Calcul des valeurs propres (ACP)
# --------------------------------

valpcr <- eigenvalues(Xcr,draw=TRUE)$eigenvalues
rm(Xcr)
gc(FALSE)

# selection du nombre de valeurs propres
# --------------------------------------

n.comp <- sum(valpcr>1)

cat("Eigenvalues have been computed\n")
cat(paste(n.comp," components will now be extracted using ICA\n",sep=""))

}


#############################
# Analyse statistique: ICA
###########################

if (is.spatial) { # spatial ICA
resICA <- ICAspat(Xcentred,n.comp,alg.typ="parallel",centering=FALSE,hp.filter)
gc(FALSE)
} else { # temporal ICA
resICA <- ICAtemp(t(Xcentred),n.comp,alg.typ="parallel",centering=FALSE,hp.filter)
gc(FALSE)
}

time.series <- resICA$time.series # tm x n.comp (spatial case) or n.comp x tm (temporal case)
spatial.components <- resICA$spatial.components # n.comp x vm (spatial case) or vm x n.comp (temporal case)



# Recréation des volumes 3D
littlefunc <- function(x,dim.resmask,dimensions,mask.img) {
  volume <- matrix(NA,nrow=dim.resmask,ncol=1)
  volume[mask.img!=0,] <- x
  volume <- array(volume,dim=dimensions)
  return(volume)
}

spatial.components.vol <- vector("list", n.comp)
for (i in 1:n.comp) {
  if (is.spatial) { # spatial ICA
    spatial.components.vol[[i]] <- littlefunc(spatial.components[i,],dim.resmask,dimensions,mask.img)
  } else { # temporal ICA
    spatial.components.vol[[i]] <- littlefunc(spatial.components[,i],dim.resmask,dimensions,mask.img)    
  }
}

cat("ICA has been performed\n")

#############################
# Ecriture de données
#################################

# on écrit les décours temporels

if (is.spatial) { # spatial ICA
  write.table(time.series,file=paste(substr(foncfile,1,nchar(foncfile)-4),"-ICAs-time-series.dat",sep=""))
} else { # temporal ICA
  write.table(time.series,file=paste(substr(foncfile,1,nchar(foncfile)-4),"-ICAt-time-series.dat",sep=""))
}


# on écrit toutes les images

Lbis <- fonc.hdr
Lbis$dim[5] <- n.comp
Lbis$intent.code <- 1001


name.prefix.fonc <- paste(substr(foncfile,1,nchar(foncfile)-4),"_ICA",sep="")

if (is.spatial) { # spatial ICA
  Lbis$descrip <- "Spatial ICA estimated (source) components"
  name.prefix.fonc <- paste(name.prefix.fonc,"s",sep="")
} else { # temporal ICA
  Lbis$descrip <- "Temporal ICA estimated (mixing) components"
  name.prefix.fonc <- paste(name.prefix.fonc,"t",sep="")
}

spatial.component.vol <- array(0,dim=Lbis$dim[2:5])
for (i in 1:n.comp) {
spatial.component.vol[,,,i] <- spatial.components.vol[[i]]
}

spatial.component.vol[is.na(spatial.component.vol)] <- 0
f.write.nifti(mat=spatial.component.vol,file=name.prefix.fonc,size="int",L=Lbis,nii=TRUE)




if (n.comp.compute) dev.off()
gc(FALSE)

}


reduction <- function(X,row.red=TRUE) {
# Description:
# ------------
# This function reduces the data in the row or col dimension
  
# Arguments:
# ----------
# X:       a matrix of size tm x vm which contains the functionnal images  
# row.red: Logical. Reduces the columns or the rows 
  
# Values:
# -------
# Xred: the reduced matrix

# Example:
# --------
# Xcr <- reduction(Xcentred,row.red=TRUE)$Xred


  vm <- ncol(X)
  tm <- nrow(X)

  if (row.red) {
    
    sd.X <- apply(X, 1, sd)*sqrt((vm-1)/vm)
    
    Xred <- sweep(X, 1,sd.X,FUN="/") # divide by standard deviation
  }
  
  else { 

    sd.X <- apply(X, 2, sd)*sqrt((tm-1)/tm)
  
    Xred <- sweep(X, 2,sd.X,FUN="/") # divide by standard deviation
  }
  
  
  return(list(Xred=Xred))

}

centering <- function(X,col.first=TRUE) {

# Description:
# ------------
# This function center the data in the two dimensions, the first dimension being indicated by col.first argument
  
# Arguments:
# ----------
# X:         a matrix of size tm x vm which contains the functionnal images  
# col.first: Logical. Center the columns or the rows first
  
# Values:
# -------
# Xcentred: the double centered matrix

# Example:
# --------
# Xcentred <- centering(X.masked,col.first=TRUE)$Xcentred

  
if (col.first) {

  mean.X <- apply(X, 2, mean)
  Xcentred <- sweep(X, 2, mean.X) # subtract the column mean for X
  mean.Xcentred <- apply(Xcentred, 1, mean)
  Xcentred <- sweep(Xcentred, 1, mean.Xcentred) # subtract the row mean
  
} else {
  
  mean.X <- apply(X, 1, mean)
  Xcentred <- sweep(X, 1, mean.X) # subtract the row mean for X
  mean.Xcentred <- apply(Xcentred, 2, mean)
  Xcentred <- sweep(Xcentred, 2, mean.Xcentred) # subtract the column mean
  
}

return(list(Xcentred=Xcentred))

}

eigenvalues <- function(X,draw=FALSE) {
  
# Description:
# ------------
# This function computes the eigenvalues of the centered and reduced data matrix
  
# Arguments:
# ----------
# X:    a matrix of size tm x vm which contains the functionnal images centered and reduced 
# draw: Logical. Should we plot the eigenvalues
  
# Values:
# -------
# eigenvalues: vector of the eigenvalues

# Example:
# --------
# valpcr <- eigenvalues(Xcr,draw=T)$eigenvalues


  tm <- nrow(X)
  vm <- ncol(X)

# On calcule les valeurs propres de la matrice des corrélations de X (donc X doit être centrée ET réduite!!)
  valp <- eigen(X%*%t(X)/vm,only.values=TRUE,symmetric=TRUE)$values

  if (draw) {
# Inspection visuelle du pouvoir explicatif des valeurs propres
    par(mfrow=c(2,1))
    plot(valp/sum(valp),type="l",col="blue",ylim=c(0,1),main=paste("Percentage of explained variance\nThreshold: eigenvalues > 1/tm=",round(1/tm,3)," are selected",sep=""),xlab="Index of sorted (decreasing) eigenvalues")
    abline(h=1/tm) # A VOIR!! pourquoi pas 1/vm ??
    plot(cumsum(valp/sum(valp)),type="l",col="red",main="Cumulated inertia",ylim=c(0,1),xlab="Index of sorted (decreasing) eigenvalues")
    x <- valp/sum(valp)>(1/tm) # A VOIR!! pourquoi pas 1/vm ??
    abline(v=which(x)[length(which(x))])
    
  }

  return(list(eigenvalues=valp))

}



ICAspat <- function(X,n.comp,alg.typ="parallel",centering=TRUE,hp.filter=TRUE) {

# On effectue l'ICA spatiale (avec blanchiment préalable)

# X: tm x vm
  
  if (centering)  X <- centering(X,col.first=TRUE)$Xcentred # Xcentred de taille tm x vm centrée d'abord en colonnes puis en lignes

  svd.Xc.spatial <- svd(X,nu=0,nv=n.comp) 
  gc(FALSE)

  Vr <- svd.Xc.spatial$v # vm x n.comp
  Drinv <- diag(1/svd.Xc.spatial$d[1:n.comp]) # n.comp x n.comp
  Ur <- X %*% Vr %*% Drinv # tm x n.comp
  
  vm <- ncol(X)
  Zr <- sqrt(vm)*t(Vr)  # de taille n.comp x vm


# high-pass filter
# ----------------
# A FAIRE!!! Réflechir comment!!!
  
  w.init <- matrix(rnorm(n.comp^2), n.comp, n.comp)
  maxit <- 200
  tol <- 1e-04
  verbose <- FALSE
  alpha <- 1
  fun <- c("logcosh")
  filter <- TRUE
# le prog estime la matrice de séparation monWZ donc monWz%*%monZ == monSest
  Wzr <- ica.R.def(Zr, n.comp, tol = tol, fun = fun,alpha = alpha, maxit = maxit, verbose = verbose,w.init = w.init) # n.comp x n.comp
  Sest <- Wzr%*%Zr # n.comp x vm
  gc(FALSE)

  Wxr <- sqrt(vm)*Wzr%*%Drinv%*%t(Ur) # n.comp  x tm
  Axr <- t(Wxr) %*% solve(Wxr %*% t(Wxr)) # de taille tm x n.comp
  

  gc(FALSE)
  return(list(time.series=Axr,spatial.components=Sest))
  
}

ICAtemp <- function(X,n.comp,alg.typ="parallel",centering=TRUE,hp.filter=TRUE) {

# On effectue l'ICA temporelle (avec blanchiment préalable)

# X: vm x tm



  
  X <- t(X) # now X: tm x vm

  if (centering)  X <- centering(X,col.first=FALSE)$Xcentred # Xcentred de taille tm x vm mais centrée d'abord en lignes puis en colonnes


  svd.Xc.spatial <- svd(X,nu=n.comp,nv=0) 
  
  gc(FALSE)
  Ur <- svd.Xc.spatial$u # Ur: tm x n.comp
  Drinv <- diag(1/svd.Xc.spatial$d[1:n.comp]) # n.comp x n.comp
  
# V <- svd(Xtcentree,nu=n.comp,nv=0)$u
# identique à
  Vr <- t(X)%*%Ur%*%Drinv # vm x n.comp
  
  tm <- nrow(X)
  Z <- sqrt(tm)*t(Ur) # n.comp x tm
  
# high-pass filter
# ----------------

  if (hp.filter) {
    M <- diag(rep(1,tm))
    M[col(M) == row(M)+1] <- -1
    Zstar <- Z%*%M
  } else {Zstar <- Z}

  
  w.init <- matrix(rnorm(n.comp^2), n.comp, n.comp)
  maxit <- 200
  tol <- 1e-04
  verbose <- FALSE
  alpha <- 1
  fun <- c("logcosh")
  filter <- TRUE

# le prog estime la matrice de séparation monWZ donc monWz%*%monZ == monSest
  Wz <- ica.R.def(Zstar, n.comp, tol = tol, fun = fun,alpha = alpha, maxit = maxit, verbose = verbose,w.init = w.init) # n.comp x n.comp
  Sestpseudo <- Wz%*%Z # n.comp x tm

  gc(FALSE)


  Wxt <- sqrt(tm)*Wz%*%Drinv%*%t(Vr) # n.comp x vm
# Contient, dans les colonnes, les cartes spatiales associées
  Axt <- t(Wxt) %*% solve(Wxt%*%t(Wxt)) # vm x n.comp

  gc(FALSE)


  return(list(time.series=Sestpseudo,spatial.components=Axt))


}





f.plot.volume.gui <- function(){

  #starts GUI that allows user apply Spatial or Temporal ICA to an fMRI dataset

  path <- path.package(package = "AnalyzeFMRI")
  path.gui <- paste(path, "plot.volume.gui.R", sep = .Platform$file.sep)
  source(path.gui)
}
