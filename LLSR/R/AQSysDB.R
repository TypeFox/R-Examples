####################################################################################################################
options(digits = 14)
require(XLConnect)
require(digest)
####################################################################################################################
#' @import XLConnect
#' @import digest
####################################################################################################################
#' @rdname AQSysDB
#' @name AQSysDB
#' @title AQSysDB
#' @description Import DB data from an Excel Worksheet.
#' @export
#' @param path String containing the full path to the XLS or XLSX file.
#' @param order Defines how the data is organized in the Worksheet. Use "xy" whether the first column corresponds to the lower phase fraction and "yx" whether the opposite.
#' @param CAS The user has the option to identify the component's cells in the worksheet with the CAS (CAS = TRUE) or with the row number that matches a CAS entry in the CASDB worksheet (CAS = FALSE)
#' @examples
#' \dontrun{
#' AQSysDB("C:/data.xls", order = "xy", CAS = FALSE)
#'}
####################################################################################################################
# AQSysDB() is a simple approach that is ready to use any three-parameter equation
# and thus
#
AQSysDB <- function(path, order = "xy", CAS = FALSE) {
  # path must point to a xlsx or xls file
  if (grepl(".xlsx",path) | grepl(".xls", path)) {
    # load file
    wrbk <- loadWorkbook(path)
    # load sheets
    sheets <- getSheets(wrbk)
    # initiate refdb - dataframe that will get reference data
    refdb <- data.frame()
    # find a sheet with the "REFDB" fragment in its name and load it
    refdb <-
      readWorksheet(wrbk, grep("REFDB", sheets), header = FALSE)
    # initiate its second column
    refdb[, 2] <- NA
    # define refdb headers
    names(refdb) <- c("REF.NAME", "REF.MD5")
    # encrypt entries found in the file using md5 and store it in refdb
    refdb[, 2] <- sapply(refdb[, 1], digest, algo = "md5")
    # initiate casdb - dataframe that will get CAS data
    casdb <- data.frame()
    # find a sheet with the "CASDB" fragment in its name and load it
    casdb <-
      readWorksheet(wrbk, grep("CASDB", sheets), header = FALSE)
    #define casdb headers
    names(casdb) <- c("CAS.CODE", "CHEM.NAME", "CHEM.COMMON")
    #
    #CONSIDER FIND OTHERS SHEETS FOLLOWING A PATTERN AND EVALUATE ALL OF THEM
    #
    #Loop through worksheet index that have the datasource string in its name
    #and aggregate data under one data prior going got analysis.
    #
    #wsdt <- readWorksheet(wrbk, 3, header = FALSE)
    #
    #if (is.odd(ncol(wsdt))) AQSys.err("2")
    # all references and cas numbers used in the analysis must be stored in its
    # respective sheer but given the nature of experimental system's data it is
    # useful to have the possibility to split data among several sheets.
    # AQSys.merge creates this possibility. Check AQSysUtils.R for details.
    wsdt <- AQSys.merge(wrbk, sheets)
    # Each system have two columns, thus the total number of columns divided by two
    # gives the number of systems
    nSys <- ncol(wsdt) / 2
    # set llsrb as a datafram which data are not converted automatically to factors
    llsrdb <- data.frame(stringsAsFactors = FALSE)
    # System's info location starts in the row below
    db.info <- 1
    # System's data location starts in the row below
    db.data <- 6
    #
  } else {
    # if an invalid path is loaded, it triggers an error
    # (check AQSys.err.R for details)
    AQSys.err("1")
  }
  # Just giving user an output on R prompt, showing what system is under analysis
  cat('\014')
  cat(paste("Analysing ", nSys," systems. \n\n", collapse = NULL))
  # the experimental phase diagram data fetched in the lines above will be used
  # to calculate the nonlinear parameters for all equations in AQSysList()
  for (i in AQSysList()) {
    # print in the prompt which equations is being used in the systems
    cat(i,": ", sep = "")
    # verify if llsrdb already have data and how many lines already were added
    db.n.row <- nrow(llsrdb)
    # If this is the first run for the dataset under analysis, initialize llsrdb
    if (db.n.row == 0) {
      db.first.row <- 1
      db.last.row <- nSys
      # if not, shift first and last rows indexes
    }else{
      db.first.row <- db.n.row + 1
      db.last.row <- db.n.row + nSys
    }
    # populate llsrdb based in the number of systems loaded from file
    for (j in db.first.row:db.last.row) {
      # c1 e c2 are the index for the systems unders analysis at the momment
      c1 <- 2 * (j - db.n.row) - 1
      c2 <- c1 + 1
      # get the data length of system under analysis
      lSys <- length(wsdt[, c1])
      # add md5 encoded ref to llsrdb
      llsrdb[j, 1] <- refdb[wsdt[db.info + 3, c1], 2]
      # if cas field in sysdb is filled with the cas
      if (CAS == TRUE) {
        # add Component's CAS directly to llsrdb
        llsrdb[j, 2] <- wsdt[db.info + 2, c1]
        llsrdb[j, 3] <- wsdt[db.info + 2, c2]
        # if cas field in sysdb in filled with an index refering to casdb
      } else{
        # Cross reference indexes and add only the Component's CAS
        llsrdb[j,2] <- casdb[wsdt[db.info + 2, c1], 1]
        llsrdb[j,3] <- casdb[wsdt[db.info + 2, c2], 1]
      }
      # populate db with system's pH, additive, additive conc and temperature
      llsrdb[j,4] <- wsdt[db.info, c1]
      llsrdb[j,5] <- wsdt[db.info, c2]
      llsrdb[j,6] <- wsdt[db.info + 1, c2]
      llsrdb[j,7] <- wsdt[db.info + 1, c1]
      # select phase diagram's data only
      rawSys <- wsdt[db.data:lSys, c1:c2]
      # remove NA entries and convert to dataframe
      db.Sys <- as.data.frame(na.exclude(rawSys))
      # checkinf if all entries are numbers and uses dot as decimal
      if (is.numeric(db.Sys[1,1])) {
        
      } else{
        db.first.col <- as.numeric(sub(",", ".", db.Sys[,1], fixed = TRUE))
        db.second.col <-
          as.numeric(sub(",", ".", db.Sys[,2], fixed = TRUE))
      }
      # When calling AQSysDB the user can use as input the order the data is in the sheet
      # if XY, analyse it normally
      if (tolower(order) == "xy") {
        resSys <-
          summary(AQSys(LLSRxy(db.first.col,db.second.col),mathDesc = i))
      }
      else{
        # if YX, invert columns
        resSys <-
          summary(AQSys(LLSRxy(db.second.col,db.first.col),mathDesc = i))
      }
      # populate llsrdb with the appropriated parameters from the nonlinear regression
      #
      llsrdb[j,8] <- resSys$parameters[1,1]
      llsrdb[j,9] <- resSys$parameters[2,1]
      llsrdb[j,10] <- resSys$parameters[3,1]
      llsrdb[j,11] <- resSys$sigma
      llsrdb[j,12] <- sum(resSys$residuals ^ 2)
      llsrdb[j,13] <- resSys$parameters[1,2]
      llsrdb[j,14] <- resSys$parameters[2,2]
      llsrdb[j,15] <- resSys$parameters[3,2]
      llsrdb[j,16] <- resSys$parameters[1,3]
      llsrdb[j,17] <- resSys$parameters[2,3]
      llsrdb[j,18] <- resSys$parameters[3,3]
      llsrdb[j,19] <- resSys$convInfo$finTol
      llsrdb[j,20] <- length(db.first.col)
      llsrdb[j,21] <- i
    }
    # At the end of the analysis for each equation,  return an OK
    cat("[OK] \n")
    #
  }
  # set llsrdb header
  names(llsrdb) <-
    c(
      "REF.MD5", "UP.Rich", "LP.Rich", "Sys.pH", "Sys.Additive",
      "Sys.Additive.Conc", "Sys.Temp","P1", "P2", "P3",  "Res.Std.Err",
      "SSR",	"P1.Std.Err",	"P2.Std.Err",	"P3.Std.Err",
      "P1.tValue",	"P2.tValue",	"P3.tValue", "Ach.Conv.Tol",
      "n.Points", "math.Desc"
    )
  # return silently all data obtained from the worksheet in a list of three dataframes
  invisible(list(
    "db.ref" = refdb, "db.sys" = llsrdb, "db.cas" = casdb
  ))
}
