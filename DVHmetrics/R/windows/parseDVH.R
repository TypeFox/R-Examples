parseDVH <- function(x, type=c("Eclipse", "Cadplan", "Masterplan",
                               "Pinnacle", "Monaco", "HiArt", "RayStation")) {
    type <- match.arg(type)

    ## name them using patient IDs
    getPatID <- function(txt) {
        if(type == "Monaco") {
            header <- unlist(strsplit(txt[1], " [|] "))
            trimWS(sub("^Patient ID: .+[~](.+)$", "\\1", header[1]))
        } else if(type == "RayStation") {
            IDline <- txt[grep("^(#PatientId):.+", txt)]
            IDres  <- sub("^.+?:[[:blank:]]*([[:alnum:][:punct:][:blank:]]+$)", "\\1", IDline, perl=TRUE)
            collWS(trimWS(IDres, side="both"))
        } else if(type == "HiArt") {
            gsub("[^a-z0-9]", "\\1", tempfile(pattern="", tmpdir=""))
        } else if(type != "Pinnacle") {
            IDline <- txt[grep("^(Patient ID|Case)[[:blank:]]*:", txt)]
            IDres  <- sub("^.+?:[[:blank:]]+([[:alnum:][:punct:][:blank:]]+$)", "\\1", IDline, perl=TRUE)
            collWS(trimWS(IDres, side="both"))
        } else {
            pInf <- paste0(txt,"/PatInfo.csv")
            if(file.exists(pInf)) {
                PatInfo <- read.csv(pInf, header=TRUE)
                removeWS(PatInfo$MedicalRecordNumber)
            } else {
                randName <- basename(tempfile(pattern="", tmpdir=""))
                warning(c("No patient information found - set to:", randName))
                randName
            }
        }
    }

    DVHraw <- if(type != "Pinnacle") {
        files <- if(!missing(x)) {
            Sys.glob(x)
        } else if(interactive()) {  # && (.Platform$OS.type == "windows"))
            ## are we are in interactive mode AND under Windows?
        	## we are under Windows since this sits in a platform-specific directory
            ## choose files interactively
            myFilt <- rbind(Filters, txtCsvDat=c("Data files (*.txt, *.csv, *.dat)",
                                                 "*.txt;*.csv;*.dat"))
            choose.files(filters=myFilt[c("txtCsvDat", "All"), ], index=1)
        } else {
            character(0)
        }
    
        files <- Filter(function(y) file_test(op="-f", y), files)
        if(length(files) >= 1L) {
            ## read in files into a list of character vectors
            lapply(files, readLines)
        } else {
            character(0)
        }
    } else {
        ## Pinnacle data is spread across files in 1 dir/patID
        ## just return directory names
        if(!missing(x)) {
            if(length(x) > 1L) {
                warning(c("Will only use x=", x[1]))
                x <- x[1]
            }

            filesAndDirs <- if(file_test(op="-f", x)) {
                ## single file - assume zip file
                ## random sub-directory
                tmpf <- gsub("^[^0-9](.+)$", "\\1", tempfile(pattern="", tmpdir=""))
                tmpd <- normalizePath(paste0(tempdir(), "/pinnacle_", tmpf), mustWork=FALSE)
                
                ## check the hierarchy in the zip file
                zipFD <- unzip(x, list=TRUE)
                dirs  <- if(all(dirname(zipFD$Name) %in% c("Data", "."))) {
                    ## zip file was created within single patient directory
                    unzip(x, exdir=tmpd)
                    tmpd
                } else {
                    ## zip file is some superordinate directory
                    unzip(x, exdir=tmpd)
                    ## only use directories
                    allDirs <- zipFD[zipFD$Length == 0L, "Name"]
                    ## chop off trailing /
                    allDirs <- gsub("^(.+)/$", "\\1", allDirs)
                    ## chop off trailing /Data
                    dirs <- unique(gsub("^(.+)/Data$", "\\1", allDirs))
                    normalizePath(paste0(tmpd, "/", dirs))
                }
            } else if(file_test(op="-d", x)) {
                ## single directory - assume it's the patient directory
                x
            } else {
                ## assume globbing pattern
                Sys.glob(x)
            }

            ## only directories we can read from
            Filter(function(y) file_test(op="-d", y) & file_test(op="-x", y), filesAndDirs)
        } else if(interactive()) {
            x <- choose.dir(getwd(), "Choose a suitable folder")
            ## check if this is a patient directory or a superordinate directory
            filesAndDirs <- Sys.glob(paste0(x, "/*"))
            if(any(grepl("^.+PatInfo\\.csv$", filesAndDirs))) {
                ## patient directory
                filesAndDirs <- unique(dirname(filesAndDirs))
            }
            Filter(function(y) file_test(op="-d", y) & file_test(op="-x", y), filesAndDirs)
        } else {
            character(0)
        }
    }
    
    ## patient id's as names for list components
    names(DVHraw) <- vapply(DVHraw, getPatID, character(1))
    DVHraw
}
