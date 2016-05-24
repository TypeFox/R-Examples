# History Mar 28 2008  Allow infile to be a connection in scanFile
#                      In readTable, call createDummy using the factors
#                      variables.
#         Mar 31 2008  Change snpNames option    
#         Apr 01 2008  Add option in loadData for file.type = 3 to call
#                      scanFile instead of readTable.  
#         Apr 03 2008  Change the option creatDummy to make.dummy in readTable.
#         Apr 04 2008  Add in.miss to readTable
#                      Remove makeIdsNames from readTable.
#         Apr 11 2008  Add getNcols function.
#                      Use scanFile in table2Format
#         Apr 15 2008  Change in sas.list
#         May 12 2008  Add remove.miss to readTable
#                      Add baseline option for readTable
#         Jun 27 2008  Add function getColumns
#         Jun 27 2008  Use the tempfile function for generating temporary
#                      file names
#         Jun 28 2008  Do not extract snp names in loadData
#         Jun 30 2008  Add code for streamed input
#         Jul 02 2008  Add code for type 5 (.zip) files
#         Jul 03 2008  Add code for type 6, 8
#                      Add getNrows function
#         Jul 08 2008  Add changeLevels and subsetData options to readTable
#                      Unfactor all variables in the beginning of readTable
#         Jul 09 2008  Remove baseline option from readTable and
#                      keep all factors
#         Sep 17 2008  In loadData, check for outfile when it is needed.
#         Sep 24 2008  In getNcols, add option to return the columns
#         Mar 05 2009  Generalize getColumns function
#         Apr 20 2009  Change getNrows function
#         Jun 12 2009  Change default delimiter in loadData to "\t" 
#         Jun 15 2009  In loadData, method=2 file.type=3,6,8 change
#                      default returnMatrix to 1    
#         Jun 15 2009  For snpNames option in loadData, use the substring
#                      option in extractByStr   
#         Jul 20 2009  Fix bug in getNrows. Set delimiter to "\n".   
#         Aug 03 2009  Change in readTable: set as.is=FALSE in readTable
#                      Remove changeLevels option    
#         Aug 10 2009  Pass snp.list option snpNames.keep into extractByStr function  
#         Aug 19 2009  Add checkVars function   
#         Aug 27 2009  Set as.is=TRUE in readTable
#         Oct 16 2009  Changes for set up:
#                      readTable
#         Dec 24 2009  In readTable, call unfactor.all for file.type=1
#         Dec 28 2009  Add function getFileType
#                      Use getFileType in loadData 
#         Dec 29 2009  Add getFileDelim function
#         Jan 06 2010  Call getFileType in getFID
#         Jan 19 2010  Add function loadData.table
#         Mar 23 2010  Update getFileDelim
#         Apr 08 2010  Add function to read data from IMPUTE2 file
#         Apr 12 2010  Add option to load.impute2
#         Jun 09 2010  Return imputed flag for each snp in load.impute2
#         Jul 22 2010  Update getFileDelim for data passed in
#         Oct 11 2010  Update getFileDelim
#         Jan 04 2011  Add getFileHeader function
#         Jul 08 2011  Modify load.impute2 for missing genotypes (0 0 0)
#         Dec 22 2011  Allow for genetic.model = 4 (heterozygous) in load.impute2
#         Jan 13 2012  Add function to read a compressed file for a C function

# Function to load and return the data object (or subset)
loadData <- function(infile, p) {

  # infile         The R object file containing the data. Only 1 object should
  #                be in the file. The load() function is used to read in
  #                the data.
  #                No default
  #######################################################################
  # p         A list with the folowing names:
  #
  # file.type      1-8  Type 1 is for a file created with the save()
  #                command. Type 2 is for other flat files that have row 1
  #                as subject ids and then the snps as rows.
  #                Type 3 is for tables where column 1 is the subject id
  #                and the remining columns are snps. See delimiter,
  #                transpose and method below.
  #                Type 4 is a sas data set. See sas.list below.
  #                No default
  # start.row      The starting row. The default value is 1
  # stop.row       The end row or -1 to return all the data from row
  #                start.row onwards.
  #                The default is -1
  # include.row1   0 or 1 to include the first row of the data. The
  #                first row may be a header containing variable names.
  #                This option overrides start.row
  # snpNames       Character vector of the snp names to use. If NULL, then
  #                start.row and stop.row will be used.
  #                The default is NULL.
  # snpNames.keep  0 or 1 to remove or keep the snps in snpNames
  #                The default is 1 (keep).
  # read.n         Integer   This option is only used with file.type=3.
  #                It specifies the maximum number of lines to read in
  #                each iteration.
  #                The default value is -1, so that the entire table will
  #                be read in all at once.
  # delimiter      This option is only for file.type = 3. 
  #                The default is "|".
  # transpose      0 or 1 to transpose the data. Set to 1 for the snp data.
  #                This option is only for file.type = 3 and 4.
  #                The default is 0.
  # method         1 or 2  This option is only for file.type = 3 and 
  #                transpose = 0.
  #                1: call readTable
  #                2: call scanFile
  #                The default is 1.
  # stream         0 or 1 for streamed input
  #                The default is 0
  # outfile        Name of the output (temporary) file if one needs it to
  #                be created with stream = 1.
  #                The default is NULL. 
  # id.var         For type 3 and 4 only
  #                No default
  # For file.type = 4
  # temp.list
  # sas.list       Only for file.type = 4. A list with the following names:
  #         sas.exe
  #         sas.file
  #         shell
  # subject.list   Only for file.type = 9, 10.
  #                File containing subject ids
  #########################################################################

  # Set defaults
  p <- default.list(p, c("start.row", "stop.row", 
         "include.row1", "read.n", "delimiter", "transpose", 
         "stream", "snpNames.keep"), 
       list(1, -1, 1, -1, "\t", 0, 0, 1), 
       error=c(0, 0, 0, 0, 0, 0, 0, 0))

  type <- p[["file.type", exact=TRUE]]
  if (is.null(type)) type <- getFileType(infile, default=-1)
  if (type == -1) stop("ERROR in loadData: file.type must be specified")
  p$file.type <- type

  # Check for errors with start.row and stop.row
  if ((!p$include.row1) && (p$start.row == 1)) p$start.row <- 2

  # Check for errors
  if (p$start.row < 1) stop("ERROR: with start.row")
  if ((p$stop.row > 0) && (p$start.row > p$stop.row)) {
    stop("ERROR: with start.row and/or stop.row")
  }

  type      <- p$file.type
  stream    <- p$stream
  outfile   <- p[["outfile", exact=TRUE]]
  p$outfile <- NULL

  # Set default for some file type
  temp <- p[["method", exact=TRUE]]
  if (is.null(temp)) {
    if (type %in% c(6, 8)) {
      p$method <- 2
    } else {
      p$method <- 1
    }
  } 
  if ((type %in% c(3, 6, 8)) & (p$method == 2)) {
    p <- default.list(p, c("returnMatrix", "what"), list(1, "character"))
  }

  # Load the data
  if (type == 1) {
    # R object file (variables are rows)
    dat <- loadFile(infile, p)

  } else if (type == 2) {
    # Flat file (variables are rows)
    # If stream, then return file id
    if (stream) return(file(infile, "r"))

    p$delimiter <- "\n"
    dat <- scanFile(infile, p)

  } else if (type == 3) {
    # Flat file (variables are columns)
    if (p$transpose) {
      p$include.header <- p$include.row1
      p$vars           <- p[["snpNames", exact=TRUE]]
      p$start.col      <- p$start.row
      p$stop.col       <- p$stop.row
      dat              <- table2Format(infile, p)
    } else {
      if (p$method == 2) {
        dat <- scanFile(infile, p)
      } else {
        dat <- readTable(infile, p)
      }
    }
  } else if (type == 4) {
    # SAS data set
    p$vars            <- p$snpNames
    p$start.col       <- p$start.row
    p$stop.col        <- p$stop.row
    p$sas.list$id.var <- p$id.var

    dat <- getSASData(infile, p, transpose=p$transpose) 

  } else if (type == 5) {
    # Type 2 compressed zip file
    dat <- readFile.zip2(infile, p)
    if (stream) return(dat)

  } else if (type == 7) {
    # Type 2 compressed gz file
    dat <- readFile.gz2(infile, p)
    if (stream) return(dat)

  } else if (type %in% c(6, 8)) {
    # Type 3 compressed file
    dat <- readFile.gzip3(infile, p)
    if (stream) return(dat)
  } else if (type %in% c(9, 10)) {
    # Returned object is a list
    return(load.impute2(infile, p))

  } else {
    # GLU format
    dat <- glu.transform(infile, inFormat=type)
    return(dat)
  } 

  # Extract snps
  temp <- p[["snpNames", exact=TRUE]]
  if (!is.null(temp)) {
    op <- list(include.row1=p$include.row1, substr.vec=c(1, 15), 
               keep=p$snpNames.keep, exact=1)
    dat <- extractByStr(dat, temp, op=op)
  }

  # For streamed input
  if (stream) {
    if (is.null(outfile)) stop("ERROR: outfile must be specified")
    writeLines(dat, con=outfile)
    dat <- file(outfile, "r")
  }

  dat

} # END: loadData

# Function to load an object from a file that was created with the 
#  save() function
loadFile <- function(infile, p) {

  # infile       Path to the file to load
  ############################################################
  # p      A list of options with the names:
  # name         The name of the object in the file
  #              For now, the file should only contain 1 object
  #              The default value is "dat"
  # start.row    The default is 1
  # stop.row     The default is -1
  # include.row1 The default is 1

  p <- default.list(p, c("start.row", "stop.row", "include.row1", 
                        "name", "stream"), 
       list(1, -1, 1, "dat", 0))

  # Initialize 
  dat <- NULL

  # Load the object
  temp <- load(infile)
 
  if (length(temp) > 1) stop("ERROR: infile contains more than 1 object")

  # Rename the object
  if (temp != p$name) {
    dat <- eval(parse(text=temp))
    eval(parse(text=paste("rm(", temp, ")", sep="")))
    temp <- gc(verbose=FALSE)
  }

  # Get the number of rows
  n <- length(dat)

  rows <- NULL
  if ((p$stop.row < 0) || (p$stop.row > n)) p$stop.row <- n
  if ((p$include.row1) && (p$start.row <= 2)) {
    p$start.row    <- 1
    p$include.row1 <- 0
  } 

  if ((p$stop.row < n) || (p$start.row > 1)) {
    rows <- (p$start.row):(p$stop.row)
  }
  if (p$include.row1) rows <- c(1, rows)
  temp <- length(rows)
  if ((temp) && (temp != n)) dat <- dat[rows]

  dat

} # END: loadFile

# Function to read a data set using scan()
scanFile <- function(infile, p) {

  # infile        File to read. Infile can also be a connection.
  #               If it is a connection, then include.row1 gets set to
  #               0 if start.row > 2 and start.row gets set to 1.
  #               No default
  #################################################################
  # p      A list with the following names
  # start.row     Starting row to read.
  #               The default is 1.
  # stop.row      The stopping row to read.
  #               The default is -1, so that up to the end of the file
  #               will be read.
  # include.row1  0 or 1 to include row1 
  #               The default is 1.
  # what          Type of data to read
  #               The default is "character"
  # delimiter     The delimiter
  #               The default is "\n"
  # returnMatrix  0 or 1 to return a matrix. If set to 1 and include.row1=1,
  #               then it is assumed that row 1 is a header and the column
  #               names will be assigned.
  #               The default is 0 
  # ncol          If returnMatrix = 1 and infile is a connection, then
  #               ncol must be set to the number of columns in the returned
  #               matrix. ncol does not need to be specified otherwise.
  #               The default is NULL
  ##################################################################

  p <- default.list(p, c("start.row", "stop.row", "what", 
         "include.row1", "delimiter", "returnMatrix"), 
       list(1, -1, "character", 1, "\n", 0), 
       error=c(0, 0, 0, 0, 0, 0))
    
  start.row    <- floor(p$start.row)
  stop.row     <- floor(p$stop.row)
  include.row1 <- p$include.row1
  connFlag     <- 0
 
  if (stop.row > 0) {
    if ((start.row < 0) || (start.row > stop.row)) start.row <- 1
    if ((stop.row == 1) && (include.row1)) include.row1 <- 0
  }

  row1Flag <- 0
  temp <- (p$what == "character")
  if (!length(temp)) temp <- 0
  if (temp) {
    if (include.row1) {
      row1Flag <- 1
      if (start.row == 1) start.row <- 2
    }
    #if ((include.row1) && (start.row == 2)) start.row <- 1
    #if ((include.row1) && (start.row > 1)) row1Flag <- 1
  } else {
    if ((include.row1) && (p$returnMatrix)) row1Flag <- 1
    if ((include.row1) && (p$returnMatrix) && (start.row == 1)) start.row <- 2
  }

  # See if infile is a connection
  if ("connection" %in% class(infile)) {
    connFlag <- 1
    if ((p$returnMatrix) && (is.null(p$ncol))) stop("ERROR: ncol must be set")
    if (include.row1) {
      if (start.row <= 2) {
        row1Flag  <- 1
        start.row <- 1
      } else {
        row1Flag     <- 0
        include.row1 <- 0
      }
    }
    skip   <- start.row - 1
    nlines <- stop.row - start.row + 1
    if ((row1Flag) && (start.row == 1)) nlines <- nlines - 1
    if (nlines == 0) nlines <- 1
  } else {
    skip   <- start.row - 1
    nlines <- stop.row - start.row + 1 
  } 

  # Get row1
  if (row1Flag) {
    row1 <- scan(file=infile, what="character", nlines=1, 
                 sep=p$delimiter, quiet=TRUE)
  } else {
    row1 <- NULL
  }

  # Get the rest
  dat <- scan(file=infile, what=p$what, skip=skip, quiet=TRUE,
                 nlines=nlines, sep=p$delimiter)

  if (p$returnMatrix) {
    # Get the number of columns
    if (!is.null(p$nc)) {
      nc <- p$nc
    } else {
      if (!include.row1) {
        nc <- length(scan(file=infile, what=p$what, nlines=1, 
                     sep=p$delimiter, quiet=TRUE))
      } else {
        nc <- length(row1)
      }
    }

    # Convert to a matrix
    dat <- matrix(data=dat, ncol=nc, byrow=TRUE)

    # Assign column names
    if (include.row1) colnames(dat) <- row1
  } else {
    if (!is.null(row1)) dat <- c(row1, dat)
  }

  dat

} # END: scanFile

# Function to read in a table and convert data to the proper format
table2Format <- function(infile, p) {

  # infile         Path to the file. The file must have a header.
  #                infile can also be a data frame or matrix
  #                No default
  ##################################################################
  # p        A list with the following names:
  # delimiter      Delimiter in infile (infile is a file)
  #                The default is " "
  # out.delimiter  Delimiter to use in transposed data
  #                The default is "|"
  # id.var         Name or column number of the id variable
  #                No default
  # include.header 0 or 1 to include the header as the first column
  #                in the transposed data.
  #                The default is 1.
  # read.n         Number of rows to read each iteration.
  #                The default is that all rows will be read.
  # vars           Vector of variable names or column numbers to read.
  #                If column numbers, then the vector must be numeric.
  #                The default is NULL
  # start.col      Starting column to read
  #                The default is 1
  # stop.col       Ending column to read.
  #                The default is -1.
  # outfile        Output file
  #                The default is NULL
  # out.type       1 = compressed R object file (using save() function)
  #                2 = flat file  
  #                The default is 1.  
  ###################################################################

  # Local function
  f2 <- function(col) {
    paste(col, sep="", collapse=delimiter)
  }

  p <- default.list(p, c("id.var", "start.col", "stop.col", "read.n", 
         "include.header", "out.delimiter", "out.type", "delimiter"), 
       list("ERROR", 1, -1, -1, 1, "\t", 1, " "),
        error=c(1, 0, 0, 0, 0, 0, 0, 0))

  # Initialize 
  include.header <- p$include.header
  delimiter      <- p$out.delimiter
  start.col      <- p$start.col
  stop.col       <- p$stop.col
  vars           <- p$vars
  stop           <- 0
  index          <- 1
  ids            <- NULL
  rflag          <- 0
  readFlag       <- 0
  idFlag         <- 0
  fid            <- NULL
  tlist          <- list()

  if ((!is.data.frame(infile)) && (!is.matrix(infile))) {
    # Get the number of columns
    temp <- getNcols(infile, p)

    # Open a connection
    fid <- file(infile, "r")

    # Define a list for scanFile
    tlist <- list(what="character", returnMatrix=1, delimiter=p$delimiter,
                  include.row1=1, start.row=1, stop.row=p$read.n,
                  ncol=temp)
  }

  while (!stop) {
    if (!is.null(fid)) {
      # Read in the data
      infile <- scanFile(fid, tlist)
    
      # Check for error
      if (class(infile) == "try-error") {
        if (readFlag) {
          break
        } else {
          # Set read.n and continue
          if (read.n < 1) {
            read.n <- 1000
          } else {
            read.n <- ceiling(read.n/2)
          }
          next
        }
      } else {
        if (!length(infile)) break
      }
    } # END: if (!is.null(fid))
    readFlag <- 1
    tlist$include.row1 <- 0

    if (index == 1) {
      nc <- ncol(infile)
      cnames <- colnames(infile)
      if (nc < 2) stop("ERROR: the table has too few columns")

      # Get the name of the id variable
      if (is.numeric(p$id.var)) {
        if ((p$id.var < 1) || (p$id.var > nc)) stop("ERROR with id.var")
        idname <- cnames[p$id.var]
      } else {
        idname   <- p$id.var
        p$id.var <- match(idname, cnames)
        if (is.na(p$id.var)) stop("ERROR: with id.var") 
      }

      # Get the column numbers for the vars
      if (!is.null(vars)) {
        rflag <- 1
        if (!is.numeric(vars)) {
          vars <- match(vars, cnames)

          # Remove NAs
          vars <- vars[!is.na(vars)]
          
          if (!length(vars)) {
            if (!is.null(fid)) close(fid)
            return("")
          }
 
        } else {
          if ((any(vars < 1)) || (any(vars > nc))) 
            stop("ERROR: column numbers are not correct")
        }

        # Set start row
        start.col <- min(vars)
      } else {
        if (start.col < 1) start.col <- 1
        if ((start.col == 2) && (include.header)) start.col <- 1
        if (stop.col < 1) stop.col   <- nc
        if (stop.col > nc) stop.col  <- nc
        #if ((start.col==1) && (stop.col==1)) include.header <- 1

        if ((stop.col < nc) || (start.col > 1)) {
          rflag <- 1
          vars  <- start.col:stop.col 
        }
      }
      # See if the id variable is the only variable
      if ((length(vars) == 1) && (vars == p$id.var)) include.header <- 1

    } # END if (index == 1)

    # Get the ids  
    ids <- c(ids, as.character(infile[, p$id.var]))
  
    # Keep columns
    if (rflag) {
      infile <- removeOrKeepCols(infile, vars, which=1)
      if (index == 1) {
        cnames <- colnames(infile)
        nc     <- ncol(infile)
        if (idname %in% cnames) {
          idFlag <- 1
          idPos  <- match(idname, cnames)
        }
      }
    } else {
      idFlag <- 1
      idPos  <- p$id.var
    }
 
    # Remove the id column if it was not removed before
    if (idFlag) {
      infile <- removeOrKeepCols(infile, idPos, which=-1)
      if (index == 1) {
        cnames <- colnames(infile)
        nc     <- length(cnames)
      }
    }

    # Get rows 
    temp <- try(apply(infile, 2, f2), silent=TRUE)

    # Check for error
    if (class(temp) == "try-error") {
      # It failed, so loop over each column
      temp <- character(nc)

      for (i in 1:nc) temp[i] <- f2(infile[, i])
    } 

    if (index == 1) {
      data <- temp
    } else {
      data <- paste(data, temp, sep=delimiter)
    }

    if (is.null(fid)) stop <- 1
    index <- index + 1

  } # END: while (!stop)

  # Close the file
  if (!is.null(fid)) close(fid)

  rm(infile, vars, tlist, fid)
  temp <- gc(verbose=FALSE)

  # First row are the subject ids
  if (include.header) {
    data   <- c(f2(ids), data)
    cnames <- c(idname, cnames) 
  }

  # Add the snp ids
  data <- paste(cnames, data, sep=delimiter)

  # Save data
  temp <- getListName(p, "outfile")
  if (!is.null(temp)) {
    if (p$out.type == 1) {
      save(data, file=temp)
    } else {
      writeLines(data, con=temp)
    }
  } 

  data

} # END: table2Format

# Function to return a text file from a sas data set. 
# Missing values will be represented as -9. The delimiter will be "|".
getSASData <- function(sas.data, p, transpose=0) {

  # sas.data     The complete path to the SAS data set.
  #              No default
  # p            A list with the names "sas.list", "temp.list" and
  #              other names as in scanFile (for transpose = 1) or
  #              readTable (for transpose = 0)
  #              No default
  # transpose    0 or 1  1 is for retrieving the snp data.
  #              0 is for reading a sas table with the read.table()
  #              function.
  #              The default is 0.

  # Check the lists
  if (is.null(p$sas.list)) stop("ERROR: SAS options list is NULL")
  p$sas.list <- default.list(p$sas.list, 
   c("sas.exe", "sas.file", "id.var", "shell"),
   list("ERROR", "ERROR", "ERROR", "bash"),
    error=c(1, 1, 1, 0))
  p$sas.list$delimiter <- "|"

  p$sas.list$temp.list <- check.temp.list(p$sas.list$temp.list)

  # Initialize
  dir    <- p$sas.list$temp.list$dir
  id     <- p$sas.list$temp.list$id
  delete <- p$sas.list$temp.list$delete
  str    <- paste("_sas", id, "_", sep="")

  # Define the file names
  p$sas.list$outfile     <- getTempfile(dir, prefix=c(str, "out"), ext=".txt")
  p$sas.list$temp.file   <- getTempfile(dir, prefix=c(str, "tmp"), ext=".sas") 
  p$sas.list$temp.script <- getTempfile(dir, prefix=c(str, "bat"), ext=".bat") 
  p$sas.list$idFile      <- getTempfile(dir, prefix=c(str, "ids"), ext=".txt") 
  p$sas.list$delete      <- delete 

  if (transpose) {
    # Create the text files
    p$sas.list$out.miss <- -9
    ret <- SAS2Format(sas.data, p)

    # Read in the text file
    p$what         <- "character"
    p$delimiter    <- "\n"
    p$returnMatrix <- 0
    p$start.row    <- 1
    p$stop.row     <- -1
    p$snpNames     <- NULL
    dat            <- scanFile(p$sas.list$outfile, p)

    # Read in the file of subject ids
    ids <- scanFile(p$sas.list$idFile, p)
    ids <- paste(ids, sep="", collapse="|")

    # Combine
    dat <- c(ids, dat)

    # Delete file
    if (delete) file.remove(p$sas.list$idFile)

  } else {
    # Create the text file
    ret <- exportSAS(sas.data, p$sas.list)

    # Set the list
    p$file.type <- 3
    p$header    <- 1
    p$delimiter <- "|"
    p$in.miss   <- c(p$in.miss, "", " ", "  ")

    # Read the table
    dat <- readTable(p$sas.list$outfile, p) 
  }

  # Delete file
  if (delete) file.remove(p$sas.list$outfile)

  dat

} # END: getSASData

# Function to export a SAS data set
exportSAS <- function(sas.data, p) {

  # See SAS2Format for the options  

  p <- default.list(p, 
   c("sas.exe", "sas.file", "outfile", "temp.file", "delimiter", 
     "temp.script", "delete", "shell"),
   list("ERROR", "ERROR", "outfile.txt", "tempfile.sas", " ", 
        "temp.bat", 1, "bash"),
    error=c(1, 1, 0, 0, 0, 0, 0, 0))

  # Check the existence of the sas data set
  if (check.files(sas.data)) stop()

  # Check the existence of the sas source file
  if (check.files(p$sas.file)) stop()

  # Get the path to the library
  libname <- dirname(sas.data)

  # Get the name of the sas data set
  temp <- strsplit(basename(sas.data), ".", fixed=TRUE)[[1]]
  data <- temp[length(temp)-1] 

  # Open a connection to the temporary file
  fid <- try(file(p$temp.file, "w"), silent=TRUE)
  if (class(fid)[1] == "try-error") stop("ERROR: temp.file is invalid")
  
  temp <- paste('libname _tmp4125 "', libname, '"; \n', sep="")
  cat(temp, file=fid)

  # % include file
  temp <- paste('%include "', p$sas.file, '"; \n', sep="")
  cat(temp, file=fid)

  # Define macro variables
  temp <- paste("%let data = _tmp4125.", data, "; \n", sep="")
  cat(temp, file=fid)

  temp <- paste('%let outfile = "', p$outfile, '"; \n', sep="")
  cat(temp, file=fid)

  temp <- paste('%let delimiter = "', p$delimiter, '"; \n', sep="")
  cat(temp, file=fid)

  # Define the macro call
  temp <- "%export(data=&data, outfile=&outfile, \n"
  cat(temp, file=fid)
  temp <- "delimiter=&delimiter); \n"
  cat(temp, file=fid)

  # Close the file
  close(fid)

  # Create the batch file and call sas
  ret <- createScript(p$temp.script, p$sas.exe, p$temp.file, shell=p$shell,
                 delete=p$delete) 

  # Delete file
  if (p$delete) file.remove(p$temp.file)

  # Check the return code
  if (ret) stop("ERROR: calling exportSAS")

  ret

} # END: exportSAS

# Function to transform a sas data set to a text file, where the columns
# of the sas data set will be the rows.
SAS2Format <- function(sas.data, p) {

  # sas.data         The complete path to the SAS data set
  #                  No default
  ###################################################################
  # p    A list with the following names
  # start.col
  # stop.col
  # vars
  ################################################################
  # sas.list  A list with the following names
  # sas.exe          The complete path to the executable file to 
  #                  start SAS. If there is a command (eg sas) to start
  #                  SAS from any directory, then use the command instead.
  #                  No default
  # sas.file         The file containing the SAS macros needed to
  #                  convert the SAS data set.
  #                  No default
  # outfile          The output file that will contain the data in
  #                  the new format.
  #                  The default is "outfile.txt"
  # temp.file        The temporary file that will contain SAS commands
  #                  to run the SAS macro %table2Format.
  #                  The default is "tempfile.sas"
  # id.var           The id variable on the SAS data set
  #                  The default is NULL
  # delimiter        The delimiter to be used in the outfile.
  #                  The default is "|".
  # temp.script      The temporary script file that will call SAS.
  #                  The default is "temp.bat"
  # delete           0 or 1 to delete the temporary files temp.file and
  #                  temp.script.
  #                  The default is 1.
  # shell            The default is "bash"
  # out.miss         The numeric value to denote missing values.
  #                  Use "." for a regular SAS missing value
  #                  The default is -9999
  # idFile           Output file containing the SNP ids.
  #                  The default is "id.txt"
  ##################################################################

  p$sas.list <- default.list(p$sas.list, 
   c("sas.exe", "sas.file", "outfile", "temp.file", "delimiter", 
     "temp.script", "delete", "shell", "out.miss", "idFile"),
   list("ERROR", "ERROR", "outfile.txt", "tempfile.sas", "|", 
        "temp.bat", 1, "bash", -9999, "id.txt"),
    error=c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0))
 
  p <- default.list(p, c("start.col", "stop.col"), list(1, -1))
  s <- p$sas.list

  # Check the existence of the sas data set
  if (check.files(sas.data)) stop()

  # Check the existence of the sas source file
  if (check.files(s$sas.file)) stop()

  # Convert vars to a string
  if (!is.null(p$vars)) p$vars <- paste(p$vars, collapse=" ")

  # Get the path to the library
  libname <- dirname(sas.data)

  # Get the name of the sas data set
  temp <- strsplit(basename(sas.data), ".", fixed=TRUE)[[1]]
  data <- temp[length(temp)-1] 

  # Open a connection to the temporary file
  fid <- try(file(s$temp.file, "w"), silent=TRUE)
  if (class(fid)[1] == "try-error") stop("ERROR: temp.file is invalid")
  
  temp <- paste('libname _tmp4125 "', libname, '"; \n', sep="")
  cat(temp, file=fid)

  # % include file
  temp <- paste('%include "', s$sas.file, '"; \n', sep="")
  cat(temp, file=fid)

  # Define macro variables
  temp <- paste("%let data = _tmp4125.", data, "; \n", sep="")
  cat(temp, file=fid)
  temp <- paste('%let outfile = "', s$outfile, '"; \n', sep="")
  cat(temp, file=fid)
  temp <- paste("%let outMiss = ", as.numeric(s$out.miss), "; \n", sep="")
  cat(temp, file=fid)
  temp <- paste("%let idVar = ", s$id.var, "; \n", sep="")
  cat(temp, file=fid)
  temp <- paste('%let delimiter = "', s$delimiter, '"; \n', sep="")
  cat(temp, file=fid)
  temp <- paste("%let startCol = ", p$start.col, "; \n", sep="")
  cat(temp, file=fid)
  temp <- paste("%let stopCol = ", p$stop.col, "; \n", sep="")
  cat(temp, file=fid)
  temp <- paste("%let vars = ", p$vars, "; \n", sep="")
  cat(temp, file=fid)
  temp <- paste('%let idFile = "', s$idFile, '"; \n', sep="")
  cat(temp, file=fid)

  # Define the macro call
  temp <- "%table2Format(data=&data, idVar=&idVar, outfile=&outfile, \n"
  cat(temp, file=fid)
  temp <- "outMiss=&outMiss, delimiter=&delimiter, idFile=&idFile,\n"
  cat(temp, file=fid)
  temp <- "startCol=&startCol, stopCol=&stopCol, vars=&vars); \n"
  cat(temp, file=fid)

  # Close the file
  close(fid)

  # Create the batch file and call sas
  ret <- createScript(s$temp.script, s$sas.exe, s$temp.file, shell=s$shell,
                 delete=s$delete) 

  # Delete file
  if (s$delete) file.remove(s$temp.file)

  # Check the return code
  if (ret) stop("ERROR: calling SAS2Format")

  ret
  
} # END: SAS2Format 

# Function to create and execute a script file
createScript <- function(script.file, exe, run.file, shell="bash",
                 delete=1) {

  # script.file   The batch file to create.
  #               No default
  # exe           The executable file
  #               No default
  # run.file      The file that the executable program will run
  #               No default
  # shell         Shell to use (UNIX)
  #               The default is "bash"
  # delete        0 or 1 to delete script.file
  #               The default is 1

  # Open the script file
  fid <- try(file(script.file, "w"), silent=TRUE)
  if (class(fid)[1] == "try-error") stop("ERROR: script file is invalid")

  # Get the path seperator
  sep <- .Platform$file.sep

  # Get the directory names of the executable program
  dirs <- strsplit(dirname(exe), sep, fixed=TRUE)[[1]]
  file <- basename(exe)

  # Determine the platform
  os      <- .Platform$OS.type
  winFlag <- (os == "windows")

  if (winFlag) {
    temp <- paste(dirs[1], " \n", sep="")
    cat(temp, file=fid)
    cat("cd \\ \n", file=fid)
    dirs <- dirs[-1]
  }
  for (dir in dirs) {
    temp <- paste("cd ", dir, " \n", sep="") 
    cat(temp, file=fid)
  }
  
  # Define the command to call sas
  temp <- paste(file, ' "', run.file, '"', sep='')
  cat(temp, file=fid)

  # Close the file
  close(fid)

  # Call the script file
  if (winFlag) {
    ret <- shell(script.file)
  } else {
    ret <- system(paste(shell, script.file, sep=" "))
  }

  # Delete the file
  if (delete) file.remove(script.file)
 
  # Check the return code
  if (ret > 1) stop("ERROR: with system call")

  0

} # END: createScript

# Function to read in a table (data frame)
readTable <- function(file, p) {

  # file           Data file. 
  #                No default.
  #################################################################
  # p is a list with the folowing fields:
  # file.type      1 or 3  1 is for an R object file created with the
  #                save() function. 3 is for a table that will be read in
  #                with read.table()
  #                The default is 3
  # header         0 or 1 . Set to 0 if the file does not contain a header.
  #                The default is 1.
  # delimiter      The delimiter used in the file.
  #                The default is "".
  # id.var         Name or column number of the id variable. 
  #                Use NULL if there is no id variable.
  #                No default.
  # keep.vars      Vector of variable names or column numbers to keep. 
  #                The default is NULL, so that all variables will be kept.
  # remove.vars    Vector of variable names or column numbers to remove. 
  #                The default is NULL, so that all variables will be kept.
  #                Both use.vars and remove.vars cannot be specified. 
  # factor.vars    Vector of variable names or column numbers to convert
  #                into factors.
  #                The default is NULL.
  # make.dummy     0 or 1 to create dummy variables for the factor variables.
  #                The default is 0.
  # keep.ids       Vector of ids or row numbers to keep.
  #                The default is NULL.
  # remove.ids     Vector of ids or row numbers to remove.
  #                Both keep.ids and remove.ids cannot be specified
  #                The default is NULL.
  # in.miss        Character vector of strings to define missing values
  #                in read.table (na.strings)
  #                The default is "NA"
  # remove.miss    0 or 1 to remove rows with missing values.
  #                Rows are removed after keep.vars or remove.vars has
  #                been applied.
  #                The default is 0.
  # subsetData     List of sublists. Each sublist must contain the names
  #                "var", "operator", and "value". 
  #                See the function subsetData.list in the file wga_util.R
  #                Ex: list(list(var="GENDER", operator"==", value="MALE"))
  #                The default is NULL
  ##################################################################

  # Set defaults
  p <- default.list(p, 
      c("header", "file.type", "delimiter", "make.dummy", "in.miss",
        "remove.miss", "is.the.data"),
      list(1, 3, "", 0, "NA", 0, 0))

  # Check for errors
  if ((!is.null(p$keep.vars)) && (!is.null(p$remove.vars))) {
    stop("ERROR: keep.vars and remove.vars cannot both be specified.")
  }
  if ((!is.null(p$keep.ids)) && (!is.null(p$remove.ids))) {
    stop("ERROR: keep.ids and remove.ids cannot both be specified.")
  }

  # Read in the data
  if (p$is.the.data) {
    x <- p$data
    p$data <- NULL
  } else if (p$file.type == 1) {
    x <- loadFile(file)
  } else {
    x <- read.table(file, header=p$header, sep=p$delimiter, 
          na.strings=p$in.miss, as.is=TRUE)
  }

  # See if there is an id variable
  idFlag <- 1 - is.null(p$id.var)

  # Check the variables
  nc     <- ncol(x)
  cnames <- colnames(x)
  temp   <- list(maxValue=nc, minValue=1, checkList=cnames)
  check.vec(p$keep.vars, "keep.vars", temp)
  check.vec(p$remove.vars, "remove.vars", temp)
  check.vec(p$factor.vars, "factor.vars", temp)
  temp$len <- 1

  if (idFlag) check.vec(p$id.var, "id.var", temp)

  # Get the rows to keep
  if ((!is.null(p$remove.ids)) || (!is.null(p$keep.ids))) {
  
    if (!is.null(p$remove.ids)) {
      temp.ids <- p$remove.ids
      flag     <- -1
    } else {
      temp.ids <- p$keep.ids
      flag     <- 1
    }

    temp.ids <- as.character(unique(temp.ids))
    nrows    <- nrow(x)
    if (idFlag) {
      id <- as.character(x[, p$id.var])
    } else {
      id <- as.character(1:nrows)
    }
  
    temp <- id %in% temp.ids
    temp[is.na(temp)] <- FALSE
    x <- removeOrKeepRows(x, temp, which=flag)
  }

  # Subset data by variables
  sub <- getListName(p, "subsetData")
  if (!is.null(sub)) {

    # Get var names and unfactor
    for (i in 1:length(sub)) {
      temp <- getListName(sub[[i]], "var")
      x[, temp] <- unfactor(x[, temp], fun=NULL)
    } 

    x <- subsetData.list(x, sub)
  }

  # Remove variables
  temp <- getListName(p, "remove.vars")
  if (!is.null(temp)) {
    temp <- unique(temp)
    x    <- removeOrKeepCols(x, temp, which=-1) 
  }

  # Keep variables
  temp <- getListName(p, "keep.vars")
  if (!is.null(temp)) {
    temp <- unique(temp)
    x    <- removeOrKeepCols(x, temp, which=1) 
  }

  # Un-factor all factors
  if (p$file.type == 1) x <- unfactor.all(x)

  # Factor variables
  if (!is.null(p$factor.vars)) {
    factorFlag <- 1
    p$factor.vars <- unique(p$factor.vars)
    for (var in p$factor.vars) x[, var] <- factor(x[, var])

    if ((p$make.dummy) && (is.numeric(p$factor.vars))) {
      # If factor.vars is numeric, then get the variable names
      factor.char <- cnames[p$factor.vars]
    } else {
      factor.char <- p$factor.vars
    }
  } else {
    factorFlag <- 0
  }

  # Remove missing values
  if (p$remove.miss) x <- removeMiss(x)

  # Create dummy variables
  if ((p$make.dummy) && (factorFlag)) {
    n <- length(factor.char)
    b <- rep("x_!!_@#^&#*", times=n)
    d <- rep(1, times=n)
    x <- createDummy(x, vars=factor.char, baseline=b, keep.factor=d)$data
  }

  x

} # END: readTable

# Function to get the number of columns in a file or data object
getNcols <- function(infile, p) {

  p <- default.list(p, c("delimiter", "return.cols"), list("\t", 0))

  # See if infile is a data frame or matrix
  temp <- dim(infile)
  if (!is.null(temp)) return(temp[2])

  p$open <- "r"
  fid    <- getFID(infile, p)
  temp   <- scan(file=fid, what="character", sep=p$delimiter, 
               nlines=1, quiet=TRUE)
  close(fid)

  if (p$return.cols) {
    return(temp)
  } else {
    return(length(temp))
  }

} # END: getNcols

# Function to read in a file and return certain columns
getColumns <- function(file.list, vars, temp.list=NULL, op=NULL) {

  # file.list   (See file.list.wordpad)
  # vars        Vector of variable names or column numbers
  #############################################################
  # op         List with names
  #   return.type   0, 1, 2 0=list, 1=matrix, 2=data frame

  # Returns a list with each element is a name form vars, 
  # or if vars is numeric the position in the list.

  file.list <- default.list(file.list, 
                            c("file", "file.type", "delimiter", "header"), 
                            list("ERROR", 3, "\t", 1), error=c(1, 0, 0, 0))
  op <- default.list(op, c("return.type"), list(0))

  type <- getListName(file.list, "file.type")
  file <- getListName(file.list, "file")
  vars <- unique(vars)

  if (type == 1) {
    data <- loadFile(file)
  } else if (type == 3) {
    file.list$include.row1 <- file.list$header
    file.list$what         <- "character"
    file.list$returnMatrix <- 1

    data <- scanFile(file, file.list)

  } else if (type == 4) {
    temp.list <- check.temp.list(temp.list)
    file.list$sas.list$delimiter <- "|"
    file.list$sas.list$temp.list <- temp.list
   
    # Define an id variable to prevent an error
    file.list$sas.list$id.var <- "id"
    data <- getSASData(file, file.list, transpose=0)

  } else if (type %in% c(6, 8)) {
    file.list$method       <- 2
    file.list$include.row1 <- file.list$header
    file.list$what         <- "character"
    file.list$returnMatrix <- 1

    data <- readFile.gzip3(file, file.list)
  }

  # Check for errors
  oplist <- list(checkList=colnames(data), minValue=1, maxValue=ncol(data))
  check.vec(vars, "vars", oplist)

  return.type <- op$return.type
  if (return.type == 0) {
    ret <- list()
    for (var in vars) ret[[var]] <- as.vector(data[, var])
  } else {
    ret <- removeOrKeepCols(data, vars, which=1)
    if (return.type == 1) {
      ret <- as.matrix(ret)
    } else {
      ret <- data.frame(ret)
    }
  }

  ret
  
} # END: getColumns

# Function for type 2 zip files
readFile.zip2 <- function(infile, p) {

  p <- default.list(p, c("zipFile", "start.row", "stop.row", "stream"),
       list("ERROR", 1, -1, 0), error=c(1, 0, 0, 0))

  fid <- unz(infile, p$zipFile, open="r")
  if (p$stream) {
    return(fid)
  } else {
    return(readFile.type2(fid, p))
  }

} # END: readFile.zip2

# Function for type 2 zip files
readFile.gz2 <- function(infile, p) {

  p <- default.list(p, c("start.row", "stop.row", "stream"),
       list(1, -1, 0), error=c(0, 0, 0))

  fid <- gzfile(infile, open="r")
  if (p$stream) {
    return(fid)
  } else {
    return(readFile.type2(fid, p))
  }

} # END: readFile.gz2

# Function for readFile.zip2 and readFile.gz2
readFile.type2 <- function(fid, p) {

  # Read row 1
  row1 <- scan(file=fid, what="character", nlines=1, sep="\n", quiet=TRUE)

  # Update start.row and stop.row 
  p$start.row    <- max(p$start.row - 1, 1)
  if (p$stop.row > 0) p$stop.row <- max(p$stop.row - 1, 1)
  p$what         <- "character"
  p$delimiter    <- "\n"
  p$include.row1 <- 0
  p$returnMatrix <- 0
  dat            <- scanFile(fid, p)
  close(fid)
  dat            <- c(row1, dat)
  dat

} # END: readFile.type2

# Function for opening a file
getFID <- function(infile, p) {

  p <- default.list(p, c("open"), list("r"))
  if (is.null(p[["file.type", exact=TRUE]])) p$file.type <- getFileType(infile)

  type <- p$file.type
  if (type %in% c(2, 3, 9)) {
    fid <- file(infile, open=p$open)
  } else if (type %in% c(5, 6)) {
    fid <- unz(infile, p$zipFile, open=p$open)
  } else if (type %in% c(7, 8, 10)) {
    fid <- gzfile(infile, open=p$open)
  }

  fid

} # END: getFID

# Function for type 3 zip or gz files
readFile.gzip3 <- function(infile, p) {

  p <- default.list(p, c("file.type", "method", "delimiter"),
       list("ERROR", 1, "\t"), error=c(1, 0, 0),
       checkList=list(c(6, 8), 1:2, NA))

  # Call scanFile or readTable
  # NOTE: read.table automatically closes the connection
  if (p$method == 2) {
    # Set up list for scanFile
    p$open         <- "r"
    p$returnMatrix <- 1
    p$ncol         <- getNcols(infile, p)
    fid            <- getFID(infile, p) 
    dat            <- scanFile(fid, p)
    close(fid)
  } else {
    p$open         <- ""
    fid            <- getFID(infile, p) 
    dat            <- readTable(fid, p)
  }
  
  dat

} # END: readFile.gzip3

# Function to get the number of rows in a data set
getNrows <- function(infile, file.type=3, delimiter="\n") {

  fid <- getFID(infile, list(file.type=file.type, open="r"))
  i   <- 0
  while (1) {
    temp <- scan(fid, what="character", sep="\n", nlines=1, 
                 n=1, quiet=TRUE)
    if (!length(temp)) break
    i <- i + 1
  }
  close(fid)
  return(i)

} # END: getNrows

# Function to check variables in a data set
checkVars <- function(flist, vars) {

  flist <- default.list(flist, c("file", "file.type", "delimiter", "header"), 
                        list("ERROR", 3, "\t", 1), error=c(1, 0, 0, 0))
  flist$return.cols <- 1
  f <- flist$file
  cols <- getNcols(f, flist) 
  if (is.numeric(vars)) {
    op <- list(minValue=1, maxValue=length(cols))
  } else {
    op <- list(checkList=cols)
  }
  
  ret <- check.vec(vars, paste("Variables for ", f, sep=""), op) 
  if (ret) stop("ERROR in checkVars")

  0

} # END: checkVars

# Function to return the file type
getFileType <- function(f, default=3) {

  ret  <- default
  vec  <- getVecFromStr(f, delimiter=".")
  n    <- length(vec)
  ext1 <- ""
  if (n) {
    ext0 <- vec[n]
    if (n > 1) ext1 <- vec[n-1]
    if (ext0 == "gz") {
      if (ext1 %in% c("txt", "xls", "csv", "dat", "data", "sdat", "def")) {
        ret <- 8
      } else if (ext1 %in% c("ldat")) {
        ret <- 7
      }
    } else if (ext0 %in% c("txt", "xls", "csv", "dat", "data", "sdat", "def")) {
      ret <- 3
    } else if (ext0 %in% c("rda")) {
      ret <- 1
    } else if (ext0 %in% c("ldat")) {
      ret <- 2
    } else if (ext0 %in% c("zip")) {
      if (ext1 %in% c("txt")) {
        ret <- 6
      } else if (ext1 %in% c("ldat")) {
        ret <- 5
      }
    } else if (ext0 %in% c("lbat", "sbat")) {
      ret <- ext0
    }
  }

  ret

} # END: getFileType

# Function for guessing what the delimiter in a file is
getFileDelim <- function(f, type=NULL, default="\t") {

  delim <- default
  if (length(f) == 1) {
    if (is.null(type)) type <- getFileType(f)
    if (type %in% c("lbat", "sbat")) return(delim)
    fid <- getFID(f, list(file.type=type))
    x <- scan(fid, what="character", nlines=2, sep="\n")
    close(fid)
    if (length(x) < 2) return(delim)
  } else {
    # Assume f is the data
    x <- f
  }

  x1 <- getVecFromStr(x[1], delimiter="")
  x1 <- sort(table(x1, exclude=NULL), decreasing=TRUE)
  x2 <- getVecFromStr(x[2], delimiter="")
  x2 <- sort(table(x2, exclude=NULL), decreasing=TRUE)
  n1 <- names(x1)
  n2 <- names(x2)

  temp <- (n1 %in% n2) & (x1 %in% x2)
  if (!any(temp)) return("\n")
  x1   <- x1[temp]
  n1   <- n1[temp]
  temp <- (n2 %in% n1) & (x2 %in% x1)
  if (!any(temp)) return("\n")
  x2   <- x2[temp]
  n2   <- n2[temp]

  check <- c("\t", " ", "\n", "|", ",")
  len1  <- length(n1)
  for (j in 1:len1) {
    temp <- n1[j]
    if ((temp %in% check)  && (x1[j] == x2[temp])) {
      delim <- temp
      break
    }
  }

  if (!(delim %in% check)) delim <- "\n"

  delim 

} # END: getFileDelim

# Function to load a table
loadData.table <- function(flist) {

  if (!is.list(flist)) {
    if (!is.character(flist)) stop("ERROR in loadData.table: incorrect type of input argument")
    temp <- flist
    flist <- list(file=temp)
  }
  flist <- check.file.list(flist)
  flag  <- 0
  type  <- flist$file.type
  if (type == 1) {
    # Load the object
    temp <- load(flist$file)
    name <- flist[["name", exact=TRUE]]
    if (is.null(name)) name <- "data"
    if ((length(temp) > 1) || (name != "data")) {
    
      # Rename the object
      data <- eval(parse(text=name))
      eval(parse(text=paste("rm(", temp, ")", sep="")))
      temp <- gc(verbose=FALSE)
    }
  } else if (type == 3) {
    method <- flist[["method", exact=TRUE]]
    if (is.null(method)) method <- 1
    if (method == 2) {
      flist$returnMatrix <- 1
      data <- scanFile(flist$file, flist)
    } else {
      data <- readTable(flist$file, flist)
      flag <- 1
    }
  } else if (type == 4) {
    # SAS data set
    flist$sas.list$id.var <- flist$id.var
    data <- getSASData(flist$file, flist, transpose=0) 
  } else if (type %in% c(6, 8)) {
    # Type 3 compressed file
    method <- flist[["method", exact=TRUE]]
    if (is.null(method)) flist$method <- 2
    data <- readFile.gzip3(flist$file, flist)
  } else {
    stop("ERROR in loadData.table: file.type must be 1, 3, 4, 6, or 8")
  } 

  # Apply other arguments
  if (!flag) {
    flist$is.the.data <- 1
    flist$data <- data
    data <- readTable(flist$file, flist)
  }

  data

} # END: loadData.table

# Function to load data from impute2 output
load.impute2 <- function(infile, op) {

  # op                List of options
  #   genetic.model   0 = Return expected value
  #                   1 = Dominant
  #                   2 = Recessive
  #                   3 = Return all probs delimited by "|"
  #                   The default is 0
  #   subject.list    List of subject ids

  # Imputed snps have "---" as first column

  gmodel <- function(v1, v2, v3) {

    if (retType == 0) {
      ret <- v2 + 2*v3
    } else if (retType == 1) {
      ret <- v2 + v3
    } else if (retType == 2) {
      ret <- v3
    } else if (retType == 3) {
      ret <- paste(v1, v2, v3, collapse="\t", sep="|") 
    } else if (retType == 4) {
      ret <- v2
    } 

    ret

  } # END: gmodel

  op <- default.list(op, c("genetic.model", "stop.row", "file.type"), 
                    list(0, -1, 10))

  # Adjust options
  op$delimiter <- "\n"
  if (op$stop.row > 1) op$stop.row <- op$stop.row - 1
  if (op$file.type == 9) {
    x <- scanFile(infile, op)
  } else {
    x <- readFile.gz2(infile, op)
  }

  # Get the subject ids
  subs     <- getIdsFromFile(op$subject.list, id.vec=NULL) 
  nsubs    <- length(subs)
  subs     <- paste(subs, collapse="\t", sep="")
  subs     <- paste("ldat\t", subs, sep="")
  nperLine <- 3*nsubs + 5
  s1       <- seq(from=1, to=nperLine-5, by=3)
  s2       <- seq(from=2, to=nperLine-5, by=3)
  s3       <- seq(from=3, to=nperLine-5, by=3)
  rem      <- 1:5
  x        <- c(subs, x)
  nx       <- length(x)
  snpNames <- character(nx-1)
  alleles  <- character(nx-1) # First should be major
  MAF      <- numeric(nx-1)
  imputed  <- rep.int(1, times=nx-1) 
  retType  <- op$genetic.model

  for (i in 1:(nx-1)) {
    vec <- getVecFromStr(x[i+1], delimiter=" ")
    if (length(vec) != nperLine) {
      temp <- paste("ERROR in load.impute2: with line ", i+1, sep="")
      stop(temp)
    }
    if (vec[1] != "---") imputed[i] <- 0
    snpNames[i] <- vec[2]
    a1    <- vec[4]
    a2    <- vec[5]
    vec   <- as.numeric(vec[-rem])
    svec1 <- vec[s1]
    svec2 <- vec[s2]
    svec3 <- vec[s3]

    # Check for missing values
    temp <- svec1 + svec2 + svec3 == 0
    temp[is.na(temp)] <- TRUE
    if (any(temp)) {
      svec1[temp] <- NA
      svec2[temp] <- NA
      svec3[temp] <- NA
    }
    
    sum1  <- sum(svec1, na.rm=TRUE)
    sum2  <- sum(svec2, na.rm=TRUE)
    sum3  <- sum(svec3, na.rm=TRUE)

    if (sum1 >= sum3) {
      alleles[i] <- paste(a1, a2, sep="")
      MAF[i]     <- (sum2 + 2*sum3)/(2*(sum1 + sum2 + sum3))
      value      <- gmodel(svec1, svec2, svec3)
    } else {
      alleles[i] <- paste(a2, a1, sep="")
      MAF[i]     <- (sum2 + 2*sum1)/(2*(sum1 + sum2 + sum3))
      value      <- gmodel(svec3, svec2, svec1)
    }
    
    temp <- paste(value, collapse="\t", sep="")
    x[i+1] <- paste(snpNames[i], "\t", temp, sep="") 
  }

  rm(s1, s2, s3, svec2, svec3)
  temp <- gc()

  vec <- op[["snpNames", exact=TRUE]]
  if (!is.null(vec)) {
    temp <- snpNames %in% vec
    if (any(temp)) {
      x        <- x[c(TRUE, temp)]
      snpNames <- snpNames[temp]
      alleles  <- alleles[temp]
      MAF      <- MAF[temp]
    } else {
      temp <- paste("NOTE: No SNPs found in file ", infile, sep="")
      print(temp)
    }
  }

  list(snpData=x, alleles=alleles, MAF=MAF, snpNames=snpNames, imputed=imputed)

} # END: load.impute2

# Function to determine if a file has a header
getFileHeader <- function(flist) {

  ret <- 1

  flist$return.cols <- 1
  cols <- getNcols(flist$file, flist)
  numcols <- as.numeric(cols)
  temp <- !is.na(numcols)
  if (any(temp)) ret <- 0

  ret  

} # END: getFileHeader

# Function to read a compressed file for a C function
# For this function, there must be a global variable GLOBAL_GZ_FPT declared 
#read_C_GZ_file <- function(filename) {
#
#  if (exists("GLOBAL_GZ_FPT")) {
#    if (is.null(GLOBAL_GZ_FPT)) {
#      GLOBAL_GZ_FPT <<- gzfile(filename, open="r")
#      ret <- scan(GLOBAL_GZ_FPT, nlines=1, what="character", sep="\n")
#    }
#  } else {
#    stop("GLOBAL_GZ_FPT must be declared as a global variable")
#  }
#  ret <- scan(GLOBAL_GZ_FPT, nlines=1, what="character", sep="\n")
#  ret
#
#} # END: read_C_GZ_file



