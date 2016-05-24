

createBinaryID <- function(ssn, o.write) {

  if (file.exists("binaryID.db") == TRUE) {
    if (o.write == TRUE) {
      unlink("binaryID.db")
    } else {
      cat("binaryID.db already exists - no changes were made to binaryID.db table\n")
      mm<- T}
  }

  options(show.error.messages = FALSE)
  m <- try(if (file.exists("binaryID.db") == FALSE) {
    mm<- F
    # DEFINE DATABASE DRIVER
    driver <- RSQLite::SQLite() 

    # DEFINE CONNECTION and create a new table
    db.name <- "binaryID.db"

    connect <- dbConnect(SQLite(), db.name)

    # get number of networks from observed sites attribute table...
    #net.no <- as.numeric(levels(ssn@obspoints@SSNPoints[[1]]@point.data$netID))
    net.no <- as.numeric(levels(ssn@network.line.coords$NetworkID))

############################################################################
# change this so that it doesn't always start with 1
# use file.list and stuff in sub.ssn function to do this....
############################################################################
    # read data into SQLite directly from file

    for (i in 1:length(net.no)) {

      network <- paste("net",net.no[i], sep = "")
      file.name <- paste("netID", net.no[i], ".dat", sep = "")

      if(dbExistsTable(connect, network)){
        dbRemoveTable(connect, network)
      }

      dbWriteTable(connect, network, read.table(file = file.name, header = T, sep = ",",
        colClasses = c("numeric","character")),overwrite = T, row.names = F)

      ## dbWriteTable(connect, network, read.table(file = file.name, header = T, sep = ",", colClasses = c("numeric","character")),
      ##   overwrite = T, row.names = F, col.names = T)

      # Check to ensure binary files were imported to SQLite database
      if (i == length(net.no)) {
        if (length(dbListTables(connect)) != length(net.no)) {
          ##sqliteCloseConnection(connect)
          dbDisconnect(connect)
          ##sqliteCloseDriver(driver)
          stop("ERROR: binary tables did not import to SQLite database properly")
        }
      }}

    # close the connection and driver
    dbDisconnect(connect)
    ##sqliteCloseConnection(connect)
    ##sqliteCloseDriver(driver)
  }, silent = TRUE)
  options(show.error.messages = TRUE)

  if (mm != T)  {
    if (m !=T) {
      dbDisconnect(connect)
      ##sqliteCloseConnection(connect)
      ##sqliteCloseDriver(driver)
      stop("ERROR: binary tables did not import to SQLite database properly")}}

}
