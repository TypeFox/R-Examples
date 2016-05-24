f.read.analyze.header <- function(file){
  #This function reads in the information from an ANALYZE format .hdr file
    file.name <- substring(file, 1, nchar(file) - 4)
    file.hdr <- paste(file.name, ".hdr", sep = "")
    file.img <- paste(file.name, ".img", sep = "")

    if(file.exists(file.img) == FALSE) return(paste(file.img, "not found"))
    if(file.exists(file.hdr) == FALSE) return(paste(file.hdr, "not found"))

#Detect whether the data is big or little endian. The first part of a .hdr file is the size of the file which is always a C int (i.e. 4 bytes) and always has value 348. Therefore trying to read it in assuming little-endian will tell you if that is the correct mode

    swap <- 0

    if(.C("swaptest_wrap_JM",
          ans = integer(1),
          file.hdr,
          PACKAGE="AnalyzeFMRI")$ans != 348) # $ans is sizeof_hdr
        swap <- 1


    
# A C function is used to read in all the components of the .hdr file
    a<-.C("read_analyze_header_wrap_JM",
          file.hdr,  # name of hdr file (with .hdr extension)       1
          as.integer(swap), # as defined above                      2
          integer(1), # sizeof_hdr                                  3
          paste(rep(" ", 10), sep = "", collapse = ""), # data_type 4
          paste(rep(" ", 18), sep = "", collapse = ""), # db_name   5
          integer(1), # extents                                     6
          integer(1), # session_error                               7
          paste(rep(" ", 1), sep = "", collapse = ""), # regular    8
          paste(rep(" ", 1), sep = "", collapse = ""), # hkey_un0   9
          integer(8), # dim                                         10
          paste(rep(" ", 4), sep = "", collapse = ""), # vox_units  11
          paste(rep(" ", 8), sep = "", collapse = ""), # cal_units  12
          integer(1), # unused1                                     13
          integer(1), # datatype                                    14
          integer(1), # bitpix                                      15
          integer(1), # dim_un0                                     16
          single(8), # pixdim                                       17
          single(1), # vox_offset                                   18
          single(1), # funused1                                     19
          single(1), # funused2                                     20
          single(1), # funused3                                     21
          single(1), # cal_max                                      22
          single(1), # cal_min                                      23
          single(1), # compressed                                   24
          single(1), # verified                                     25
          integer(1), # glmax                                       26
          integer(1), # glmin                                       27
          paste(rep(" ", 80), sep = "", collapse = ""), # descrip   28
          paste(rep(" ", 24), sep = "", collapse = ""), # aux_file  29
          paste(rep(" ", 1), sep = "", collapse = ""), # orient     30
         integer(5), # originator                                   31
#          paste(rep(" ", 10), sep = "", collapse = ""),
          paste(rep(" ", 10), sep = "", collapse = ""), # generated 32
          paste(rep(" ", 10), sep = "", collapse = ""), # scannum   33
          paste(rep(" ", 10), sep = "", collapse = ""),# patient_id 34
          paste(rep(" ", 10), sep = "", collapse = ""), # exp_date  35
          paste(rep(" ", 10), sep = "", collapse = ""), # exp_time  36
          paste(rep(" ",3 ), sep = "",collapse = ""), # hist_un0    37
          integer(1), # views                                       38
          integer(1), # vols_added                                  39
          integer(1), # start_field                                 40
          integer(1), # field_skip                                  41
          integer(1), # omax                                        42
          integer(1), # omin                                        43
          integer(1), # smax                                        44
          integer(1), # smin                                        45
          PACKAGE="AnalyzeFMRI")

# A list (called L) is created containing all the components of the .hdr file

    L <- list()
    L$file.name <- file.img
    L$swap <- a[[2]]
    L$sizeof.hdr <- a[[3]]
    L$data.type <- a[[4]]
    L$db.name <- a[[5]]
    L$extents <- a[[6]]
    L$session.error <- a[[7]]
    L$regular <- a[[8]]
    L$hkey.un0 <- a[[9]]
    L$dim <- a[[10]]
    L$vox.units <- a[[11]]
    L$cal.units <- a[[12]]
    L$unused1 <- a[[13]]
    L$datatype <- a[[14]]
    if(L$datatype == 0 ) L$data.type <- "unknown"
    if(L$datatype == 1) L$data.type <- "binary"
    if(L$datatype == 2) L$data.type <- "unsigned char"
    if(L$datatype == 4) L$data.type <- "signed short"
    if(L$datatype == 8) L$data.type <- "signed int"
    if(L$datatype == 16) L$data.type <- "float"
    if(L$datatype == 32) L$data.type <- "complex"
    if(L$datatype == 64) L$data.type <- "double precision"
    if(L$datatype == 128) L$data.type <- "rgb data"
    if(L$datatype == 255) L$data.type <- "all"
    L$bitpix <- a[[15]]
    L$dim.un0 <- a[[16]]
    L$pixdim <- a[[17]]
    L$vox.offset <- a[[18]]
    L$funused1 <- a[[19]]  # SPM extends the Analyze format by using a scaling factor for the image from the header.
    L$funused2 <- a[[20]]
    L$funused3 <- a[[21]]
    L$cal.max <- a[[22]]
    L$cal.min <- a[[23]]
    L$compressed <- a[[24]]
    L$verified <- a[[25]]
    L$glmax <- a[[26]]
    L$glmin <- a[[27]]
    L$descrip <- a[[28]]
    L$aux.file <- a[[29]]
    L$orient <- a[[30]]
    L$originator <- a[[31]] # SPM uses one of the Analyze header fields in an unorthodox way. This is the Originator field
    L$generated <- a[[32]]
    L$scannum <- a[[33]]
    L$patient.id <- a[[34]]
    L$exp.date <- a[[35]]
    L$exp.time <- a[[36]]
    L$hist.un0 <- a[[37]]
    L$views <- a[[38]]
    L$vols.added <- a[[39]]
    L$start.field <- a[[40]]
    L$field.skip <- a[[41]]
    L$omax <- a[[42]]
    L$omin <- a[[43]]
    L$smax <- a[[44]]
    L$smin <- a[[45]]

    return(L)}


f.analyze.file.summary <- function(file){
#This function prints out a concise summary of the contents of a .img/.hdr image pair

    file.name <- substring(file ,1 , nchar(file) - 4)
    file.hdr <- paste(file.name, ".hdr", sep = "")
    file.img <- paste(file.name, ".img", sep = "")
    hdr <- f.read.analyze.header(file.hdr)
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


f.basic.hdr.list.create <- function(X, file.hdr){

#creates a basic list that can be used to write a .hdr file

    dim <- dim(X)
    dim <- c(length(dim), dim, rep(0, 7 - length(dim)))

    l <- list(file = file.hdr,
              sizeof.hdr = 348,
              data.type = paste(rep(" ", 10), sep = "", collapse = ""),
              db.name = paste(rep(" ", 18), sep = "", collapse = ""),
              extents = 0,
              session.error = 0,
              regular = character(1),
              hkey.un0 = character(1),
              dim = as.integer(dim),
              vox.units = "mm",
              cal.units = "voxels",
              unused1 = 0,
              datatype = 0,
              bitpix = 0,
              dim.un0 = 0,
              pixdim = single(8),
              vox.offset = single(1),
              funused1 = single(1),
              funused2 = single(1),
              funused3 = single(1),
              cal.max = single(1),
              cal.min = single(1),
              compressed = single(1),
              verified = single(1),
              glmax = 0,
              glmin = 0,
              descrip = paste(rep(" ", 80), sep = "", collapse = ""),
              aux.file = paste(rep(" ", 24), sep = "", collapse = ""),
              orient = paste(rep(" ", 1), sep = "", collapse = ""),
        originator =  integer(5),
#              originator = paste(rep(" ", 10), sep = "", collapse = ""),
              generated = paste(rep(" ", 10), sep = "", collapse = ""),
              scannum = paste(rep(" ", 10), sep = "", collapse = ""),
              patient.id = paste(rep(" ", 10), sep = "", collapse = ""),
              exp.date = paste(rep(" ", 10), sep = "", collapse = ""),
              exp.time = paste(rep(" ", 10), sep = "", collapse = ""),
              hist.un0 = paste(rep(" ", 3), sep = "", collapse = ""),
              views = integer(1),
              vols.added = integer(1),
              start.field = integer(1),
              field.skip = integer(1),
              omax = integer(1),
              omin = integer(1),
              smax = integer(1),
              smin = integer(1) )
    return(l)
}

f.read.analyze.slice <- function(file, slice, tpt){
  #Reads in a .img file into an array
    file.name <- substring(file ,1 ,nchar(file) - 4)
    file.hdr <- paste(file.name ,".hdr", sep = "")
    file.img <- paste(file.name, ".img", sep = "")
    hdr <- f.read.analyze.header(file.hdr)
    #f.analyze.file.summary(file)

  #num.dim<-hdr$dim[1]
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
                  as.integer(offset * 1),
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
                  as.integer(offset * 2),
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
                  as.integer(offset * 4),
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
                  as.integer(offset * 4),
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
                  as.integer(offset * 8),
                  as.integer(1), PACKAGE="AnalyzeFMRI")
        vol <- array(vol$mat, dim = dim)
    #this works because array fills itself with the left most subscript moving fastest
    }

    if(hdr$datatype == 0 || hdr$datatype == 1 || hdr$datatype == 32 || hdr$datatype == 128 || hdr$datatype == 255) print(paste("The format", hdr$data.type, "is not supported yet. Please contact me if you want me to extend the functions to do this (marchini@stats.ox.ac.uk)"), quote=FALSE)

    return(vol)}


f.read.analyze.tpt <- function(file, tpt){
#Reads in one timepoint of a .img file into an array
    file.name <- substring(file, 1, nchar(file) - 4)
    file.hdr <- paste(file.name, ".hdr", sep = "")
    file.img <- paste(file.name, ".img", sep = "")
    hdr <- f.read.analyze.header(file.hdr)
  #f.analyze.file.summary(file)

  #num.dim <- hdr$dim[1]
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
                  as.integer(offset * 1),
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
                  as.integer(offset * 2),
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
                  as.integer(offset * 4),
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
                  as.integer(offset * 4),
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
                  as.integer(offset * 8),
                  as.integer(1), PACKAGE="AnalyzeFMRI")
        vol <- array(vol$mat, dim = dim)
    #this works because array fills itself with the left most subscript moving fastest
    }

    if(hdr$datatype == 0 || hdr$datatype == 1 || hdr$datatype == 32 || hdr$datatype == 128 || hdr$datatype == 255) print(paste("The format", hdr$data.type, "is not supported yet. Please contact me if you want me to extend the functions to do this (marchini@stats.ox.ac.uk)"), quote = FALSE)

    return(vol)}


f.read.analyze.slice.at.all.timepoints <- function(file, slice){
  #Reads in a slice of a .img file at all time points into an array
    file.name <- substring(file, 1, nchar(file) - 4)
    file.hdr <- paste(file.name, ".hdr", sep = "")
    file.img <- paste(file.name, ".img", sep = "")
    hdr <- f.read.analyze.header(file.hdr)
  #f.analyze.file.summary(file)

  #num.dim <- hdr$dim[1]
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
                      as.integer(offset * 1),
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
                      as.integer(offset * 2),
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
                      as.integer(offset * 4),
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
                      as.integer(offset * 4),
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
                      as.integer(offset * 8),
                      as.integer(1), PACKAGE="AnalyzeFMRI")
            vol <- array(vol$mat, dim = dim)
            vl[, , i] <- vol
        }}

    if(hdr$datatype == 0 || hdr$datatype == 1 || hdr$datatype == 32 || hdr$datatype == 128 || hdr$datatype == 255) print(paste("The format", hdr$data.type, "is not supported yet. Please contact me if you want me to extend the functions to do this (marchini@stats.ox.ac.uk)"), quote = FALSE)

    return(vl)}

f.read.analyze.ts <- function(file, x, y, z){
  #Reads in a .img file into an array
    file.name <- substring(file, 1, nchar(file) - 4)
    file.hdr <- paste(file.name, ".hdr", sep = "")
    file.img <- paste(file.name, ".img", sep = "")
    hdr <- f.read.analyze.header(file.hdr)
  #f.analyze.file.summary(file)

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
                    as.integer(offset.start * 1 + 1 * (i - 1) * offset.add),
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
                    as.integer(offset.start * 2 + 2 * (i - 1) * offset.add),
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
                    as.integer(offset.start * 4 + 4 * (i - 1) * offset.add),
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
                    as.integer(offset.start * 4 + 4 * (i - 1) * offset.add),
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
                    as.integer(offset.start * 8 + 8 * (i - 1) * offset.add),
                    as.integer(1), PACKAGE="AnalyzeFMRI")
            vol[i] <- v$mat}
    }

    if(hdr$datatype == 0 || hdr$datatype == 1 || hdr$datatype == 32 || hdr$datatype == 128 || hdr$datatype == 255) print(paste("The format", hdr$data.type, "is not supported yet. Please contact me if you want me to extend the functions to do this (marchini@stats.ox.ac.uk)"), quote = FALSE)

    return(vol)}



f.read.analyze.volume <- function(file){
  #Reads in a .img file into an array
    file.name <- substring(file, 1, nchar(file) - 4)
    file.hdr <- paste(file.name, ".hdr", sep = "")
    file.img <- paste(file.name, ".img", sep = "")
    hdr <- f.read.analyze.header(file.hdr)
  #f.analyze.file.summary(file)
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
                  as.integer(0),
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
                  as.integer(0),
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
                  as.integer(0),
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
              as.integer(0),
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
              as.integer(0),
              as.integer(1), PACKAGE="AnalyzeFMRI")
    vol <- array(vol$mat, dim = dim)
    #this works because array fills itself with the left most subscript moving fastest
}

if(hdr$datatype == 0 || hdr$datatype == 1 || hdr$datatype == 32 || hdr$datatype == 128 || hdr$datatype == 255) print(paste("The format", hdr$data.type, "is not supported yet. Please contact me if you want me to extend the functions to do this (marchini@stats.ox.ac.uk)"), quote = FALSE)

return(vol)}



f.spectral.summary <- function(file, mask.file, ret.flag = FALSE)
{
  #for an analyze .img file the periodogram of the time series are divided by a flat spectral estimate using the median periodogram ordinate. The resulting values are then combined within each Fourier frequency and quantiles are plotted against freequency. This provides a fast look at a fMRI dataset to identify any artefacts that reside at single frequencies.

########################
#get info about dataset
########################
    file.name <- substring(file, 1, nchar(file) - 4)
    file.hdr <- paste(file.name, ".hdr", sep = "")
    file.img <- paste(file.name, ".img", sep = "")
    hdr <- f.read.analyze.header(file.hdr)
    nsl <- hdr$dim[4]
    nim <- hdr$dim[5]
    pxdim <- hdr$pixdim[2:4]

#####################
#read in/create mask
#####################

    f.mask.create <- function(dat, pct = .1, slices = c(0)) {
  #function that creats a mask for an fMRI dataset by thresholding the mean of the pixel time series at a percentage point of the maximum intensity of the dataset
        file.name <- substring(dat$file, 1, nchar(dat$file) - 4)
        file.hdr <- paste(file.name, ".hdr", sep = "")
        hdr <- f.read.analyze.header(file.hdr)
        nsl <- hdr$dim[4]
        xc <- hdr$dim[2]
        yc <- hdr$dim[3]
        if(slices[1] == 0){slices <- seq(1, nsl)}
        mask <- array(0, dim = c(xc, yc, length(slices)))

        max.int <- 0
        for(k in 1:length(slices)){
            slice <- f.read.analyze.slice.at.all.timepoints(dat$file, slices[k])
            if(max(slice)>max.int){max.int <- max(slice)}
        }

        for(k in 1:length(slices)){
            slice <- f.read.analyze.slice.at.all.timepoints(dat$file, slices[k])

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

    if(mask.file!=FALSE){mask <- f.read.analyze.volume(mask.file)}
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

        slice <- f.read.analyze.slice.at.all.timepoints(file, l)


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


f.write.list.to.hdr <- function(L, file){

# To respect the length of some Analyze character fields
  strcomplete <- function(string,max.length) {
    string <- substr(string,1,max.length)
    as.character(paste(as.character(string),paste(rep(" ",max.length-nchar(string)),collapse=""),collapse="",sep="")) 
  }

  
# Writes a list to a .hdr file
    a <- .C("write_analyze_header_wrap_JM",
            file,
            as.integer(L$sizeof.hdr),
            strcomplete(L$data.type,10),
            strcomplete(L$db.name,18),
            as.integer(L$extents),
            as.integer(L$session.error),
            strcomplete(L$regular,1),
            strcomplete(L$hkey.un0,1),
            as.integer(L$dim),
            strcomplete(L$vox.units,4),
            strcomplete(L$cal.units,8),
            as.integer(L$unused1),
            as.integer(L$datatype),
            as.integer(L$bitpix),
            as.integer(L$dim.un0),
            as.single(L$pixdim),
            as.single(L$vox.offset),
            as.single(L$funused1),
            as.single(L$funused2),
            as.single(L$funused3),
            as.single(L$cal.max),
            as.single(L$cal.min),
            as.single(L$compressed),
            as.single(L$verified),
            as.integer(L$glmax),
            as.integer(L$glmin),
            strcomplete(L$descrip,80),
            strcomplete(L$aux.file,24),
            strcomplete(L$orient,1),
 #           as.character(L$originator),
           as.integer(L$originator),
            strcomplete(L$generated,10),
            strcomplete(L$scannum,10),
            strcomplete(L$patient.id,10),
            strcomplete(L$exp.date,10),
            strcomplete(L$exp.time,10),
            strcomplete(L$hist.un0,3),
            as.integer(L$views),
            as.integer(L$vols.added),
            as.integer(L$start.field),
            as.integer(L$field.skip),
            as.integer(L$omax),
            as.integer(L$omin),
            as.integer(L$smax),
            as.integer(L$smin),
            PACKAGE="AnalyzeFMRI")
}


f.write.analyze <- function(mat, file, size = "float", pixdim = c(4, 4, 6), vox.units = "mm", cal.units = "pixels", originator = rep(0,5)){

  #Creates a .img and .hdr pair of files from a given array

    if(max(mat) == "NA") return("NA values in array not allowed. Files not written.")

    file.img <- paste(file, ".img", sep = "")
    file.hdr <- paste(file, ".hdr", sep = "")

    L <- f.basic.hdr.list.create(mat, file.hdr)


        L$vox.units <- vox.units
        L$cal.units <- cal.units
        L$pixdim <- c(4, pixdim, 0, 0, 0, 0)
        L$glmax <- as.integer(floor(max(mat)))
        L$glmin <- as.integer(floor(min(mat)))
        L$originator <- as.integer(originator)

    
    if(size == "float"){
        L$datatype <- 16
        L$bitpix <- 32
       L$data.type <- "float"
        f.write.array.to.img.float(mat, file.img)
    }
    if(size == "int"){
        if(max(mat)>32767 || min(mat) < ( -32768)) return("Values are outside integer range. Files not written.")
        L$datatype <- 4
        L$bitpix <- 16
        L$data.type <- "signed sho" # signed short
        f.write.array.to.img.2bytes(mat, file.img)
    }
    if(size == "char"){
        if(max(mat)>255 || min(mat) < 0) return("Values are outside integer range. Files not written.")
        L$datatype <- 2
        L$bitpix <- 8
        L$data.type <- "unsignchar" # unsigned char
        f.write.array.to.img.8bit(mat, file.img)
    }

    f.write.list.to.hdr(L, file.hdr)
}


f.write.array.to.img.2bytes <- function(mat, file){
  #writes an array into a .img file of 2 byte integers

    dm <- dim(mat)
    dm.ln <- length(dm)
    num.data.pts <- length(mat)

    .C("write2byte_JM",
       as.integer(mat),
       file,
       as.integer(num.data.pts), PACKAGE="AnalyzeFMRI")

}

f.write.array.to.img.8bit <- function(mat, file){
#writes an array into a .img file of 8 bit (1 byte) integers

    dm <- dim(mat)
    dm.ln <- length(dm)
    num.data.pts <- length(mat)

    .C("write8bit_JM",
       as.integer(mat),
       file,
       as.integer(num.data.pts), PACKAGE="AnalyzeFMRI")

}



f.write.array.to.img.float <- function(mat, file){
  #writes an array into a .img file of 4 byte flotas

    dm <- dim(mat)
    dm.ln <- length(dm)
    num.data.pts <- length(mat)

    .C("writefloat_JM",
       as.single(mat),
       file,
       as.integer(num.data.pts), PACKAGE="AnalyzeFMRI")
}

f.analyzeFMRI.gui <- function(){

  #starts GUI that allows user to explore an fMRI dataset stored in an ANALYZEfile

    path <- path.package(package = "AnalyzeFMRI")
    path.gui <- paste(path, "AnalyzeFMRI.gui.R", sep = .Platform$file.sep)
    source(path.gui)}


f.ica.fmri <- function(file.name, n.comp, norm.col = TRUE, fun = "logcosh", maxit = 1000, alg.type = "parallel", alpha = 1, tol = 0.0001, mask.file.name = NULL, slices = NULL){

  #function for performing Spatial ICA on an fMRI dataset
  #The function avoids reading the dataset into R to minimise the memory used

    hdr <- f.read.analyze.header(file.name)

    if(length(slices) == 0) slices  <-  2:(hdr$dim[4] - 1)
    if(slices == "all") slices  <-  1:(hdr$dim[4])
    if(any(slices < 1 || slices>hdr$dim[4])) {
        return("some of selected slices out of allowable range")}

    ns <- hdr$dim[2] * hdr$dim[3] * hdr$dim[4] * n.comp
    na <- hdr$dim[5] * n.comp

    mask.flag <- 1
    if(length(mask.file.name) == 0){
        mask.flag <- 0
        mask.file.name <- ""}


    col.flag <- 1
    if(norm.col!=TRUE) col.flag <- 0

    fun.flag <- 1
    if(fun == "exp") fun.flag <- 2

    def.flag <- 0
    if(alg.type == "deflation") def.flag <- 1

    W <- matrix(rnorm(n.comp * n.comp), n.comp, n.comp)

    a <- .C("ica_fmri_JM",
            as.character(file.name),
            as.single(t(W)),
            as.integer(n.comp),
            as.integer(1),
            as.integer(col.flag),
            as.integer(fun.flag),
            as.integer(maxit),
            as.integer(def.flag),
            as.single(alpha),
            as.single(tol),
            as.integer(mask.flag),
            as.character(mask.file.name),
            as.integer(slices),
            as.integer(length(slices)),
            S = single(ns),
            A = single(na),
            PACKAGE="AnalyzeFMRI")

    S <- array(a$S, dim = c(hdr$dim[2], hdr$dim[3], hdr$dim[4], n.comp))

    A <- matrix(a$A, hdr$dim[5], n.comp, byrow = TRUE)

    return(list(A = A, S = S, file = file.name, mask = mask.file.name))
}

f.plot.ica.fmri.jpg <- function(ica.obj, file = "./ica", cols = heat.colors(100), width = 700,  height = 700){

    for(i in 1:dim(ica.obj$S)[4]){

        jpeg(filename = paste(file, ".comp.", i, ".jpeg", sep = ""), width = width, height = height)
        f.plot.ica.fmri(ica.obj,  i,  cols = cols)
        dev.off()
    }
    return()
}


f.plot.ica.fmri <- function (obj.ica,  comp,  cols = heat.colors(100))
{
    r  <-  range(obj.ica$S[,  ,  ,  comp],  na.rm = TRUE)
    tmp  <-  1000 * (obj.ica$S[,  ,  ,  comp]  == 0) + (obj.ica$S[,  ,  ,  comp]) * (obj.ica$S[,  ,  ,  comp] != 0)
    ncomp  <-  dim(obj.ica$S)[4]
    nsl  <-  dim(obj.ica$S)[3]
    t  <-  nrow(obj.ica$A)
    im  <-  floor(sqrt(nsl + 3)) + 1
    par(mfrow = c(im,  im),  mar = c(.5,  .5,  .5,  .5))
    plot(c(0, 1), c(0, 1), typ = "n", axes = FALSE, xlab = "", ylab = "")

    text(0, .9, "Spatial ICA", pos = 4, cex = 1.5)
    text(0, .75, paste("Component ", comp, sep = ""), pos = 4, cex = 1.5)
    l <- floor(nchar(obj.ica$file) / 20) + 1
    text(0, .6, paste("file: ", substring(obj.ica$file, 1, 20), sep = ""), pos = 4)
    for(i in 2:l){
        text(0, .6 - .07 * (i - 1), paste("       ", substring(obj.ica$file, 20 * (i - 1) + 1, 20 * i), sep = ""), pos = 4)
    }
    text(0, .6 - .07 * (l + 1), paste("Date: ", date(), sep = ""), pos = 4)

    for (i in 1:nsl) {
        image(tmp[,  ,  i],  zlim = r,  axes = FALSE,  col = cols)
        text(.5, .98, paste("slice", i), pos = 1)
        box()
    }
    r <- range(obj.ica$A[,  comp])
    plot(obj.ica$A[,  comp],  typ = "l",  axes = FALSE, ylim = r * 1.5)
    text(length(obj.ica$A[,  comp]) / 2, 1.5 * r[2], "Time course", pos = 1)
    box()
    s  <-  fft(obj.ica$A[,  comp]) / sqrt(2 * pi * t)
    s  <-  Mod(s[2:(floor(t / 2) + 1)])^2
    r <- c(0, max(s))
    plot(s,  axes = FALSE, typ = "l", ylim = r * 1.5)
    text(length(s) / 2, 1.5 * r[2], "Periodogram", pos = 1)
    box()
    par(mfrow = c(1,  1),  mar = c(5,  4,  4,  2))
}

f.ica.fmri.gui <- function(){

  #starts GUI that allows user apply Spatial ICA to an fMRI dataset

    path <- path.package(package = "AnalyzeFMRI")
    path.gui <- paste(path, "ICA.gui.R", sep = .Platform$file.sep)
    source(path.gui)}




