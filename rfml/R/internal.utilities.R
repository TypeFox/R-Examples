###############################################################################
# Internal functions used in the package
###############################################################################

# verify that the database has everything neccessary in order to use the package
# if it exists it returns TRUE if not FALSE and if any other status_code than 404
# it will add a warning
.check.database <- function(mlHost, username, password) {

  mlURL <- paste(mlHost, "/v1/resources/rfml.check", sep="")
  response <- GET(mlURL,authenticate(username, password, type="digest"), accept_json())
  rContent <- content(response)
  if (response$status_code != 200){
    return(FALSE)
  }
  curRfmlVer <- as.character(packageVersion("rfml"))
  return(curRfmlVer == rContent$rfmlVersion)
}

# internal function to add rest interface used
.insert.lib <- function(mlHost, username, password, mlLibName)  {

  mlLibFile <- paste(mlLibName, ".sjs", sep='')
  file <- system.file("lib",mlLibFile ,package = "rfml")
  lib <- upload_file(file, "application/vnd.marklogic-javascript")
  mlURL <- paste(mlHost, "/v1/ext/rfml/", mlLibFile, sep="")
  # add or replace search options to the database
  response <- PUT(mlURL, authenticate(username, password, type="digest"), body=lib, accept_json())
  status_code <- response$status_code
  if (status_code != 201 && status_code != 204) {
    rContent <- content(response)
    errorMsg <- paste("statusCode: ",
                      rContent$errorResponse$statusCode,
                      ", status: ", rContent$errorResponse$status,
                      ", message: ", rContent$errorResponse$message, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))

  }
  # return the name of the search options
  return(TRUE)
}

# internal function to remove lib used
.remove.lib <- function(mlHost, username, password, mlLibName)  {

  mlLibFile <- paste(mlLibName, ".sjs", sep='')
  mlURL <- paste(mlHost, "/v1/ext/rfml/", mlLibFile, sep="")
  # add or replace search options to the database
  response <- DELETE(mlURL, authenticate(username, password, type="digest"), accept_json())

  if (response$status_code != 204) {

    rContent <- content(response)
    errorMsg <- paste("statusCode: ",
                      rContent$errorResponse$statusCode,
                      ", status: ", rContent$errorResponse$status,
                      ", message: ", rContent$errorResponse$message, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))
  }
  # return the name of the search options
  return(TRUE)
}

# internal function to add exstention to rest interface used
.insert.ext <- function(mlHost, username, password, mlExtName)  {

  mlExtFile <- paste(mlExtName, ".sjs", sep='')
  file <- system.file("ext",mlExtFile ,package = "rfml")
  ext <- upload_file(file, "application/vnd.marklogic-javascript")
  #  'http://localhost:8004/v1/config/resources/example'
  mlURL <- paste(mlHost, "/v1/config/resources/", mlExtName, sep="")
  # add or replace search options to the database
  response <- PUT(mlURL, authenticate(username, password, type="digest"), body=ext, accept_json())
  status_code <- response$status_code
  if (status_code != 201 && status_code != 204) {
    rContent <- content(response)
    errorMsg <- paste("statusCode: ",
                      rContent$errorResponse$statusCode,
                      ", status: ", rContent$errorResponse$status,
                      ", message: ", rContent$errorResponse$message, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))

  }
  # return the name of the search options
  return(TRUE)
}

# internal function to remove ext used
.remove.ext <- function(mlHost, username, password, mlExtName)  {

  ##mlLibFile <- paste(mlLibName, ".sjs", sep='')
  mlURL <- paste(mlHost, "/v1/config/resources/", mlExtName, sep="")
  # add or replace search options to the database
  response <- DELETE(mlURL, authenticate(username, password, type="digest"), accept_json())

  if (response$status_code != 204) {

    rContent <- content(response)
    errorMsg <- paste("statusCode: ",
                      rContent$errorResponse$statusCode,
                      ", status: ", rContent$errorResponse$status,
                      ", message: ", rContent$errorResponse$message, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))
  }
  # return the name of the search options
  return(TRUE)
}
.get.ml.metadata <- function(mlDf, nrows=0, searchOption=NULL) {

}
.get.ml.data <- function(mlDf, nrows=0, searchOption=NULL) {

  conn <- mlDf@.conn

  if (length(.rfmlEnv$key) < conn@.id) {
    stop("Need a valid connection. Use ml.connection to create one!")
  }
  key <- .rfmlEnv$key[[conn@.id]]
  password <- tryCatch(rawToChar(PKI::PKI.decrypt(conn@.password, key))
                       , error = function(err) stop("Need a valid connection. Use ml.connection to create one!"))
  username <- conn@.username
  queryComArgs <- mlDf@.queryArgs

  mlHost <- paste("http://", conn@.host, ":", conn@.port, sep="")
  mlSearchURL <- paste(mlHost, "/v1/resources/rfml.dframe", sep="")

  # need to pick start and end from mlDf...
  nStart=mlDf@.start
  if (nrows>0 && nrows<mlDf@.nrows) {
    nPageLength <- nrows
  } else {
    nPageLength <- mlDf@.nrows
  }
  queryArgs <- c(queryComArgs, 'rs:start'=nStart,'rs:pageLength'=nPageLength, 'rs:return'="data")
  # Need to check if extracted then we could have changed the rows...
  if (mlDf@.extracted) {
    # create a extfields parameter...
    extFields <- "{"
    for (i in 1:length(mlDf@.col.name)) {
      if (nchar(extFields) > 1) {
        extFields <- paste(extFields, ',', sep='')
      }
      extFields <- paste(extFields, '"', mlDf@.col.name[i],
                         '":{"fieldDef":"',mlDf@.col.name[i],
                         '","orgField":"',mlDf@.col.org_name[i],
                         '","orgPath":"',mlDf@.col.org_xpath[i],
                         '","orgFormat":"',mlDf@.col.format[i],'"}',sep='')
    }
    extFields <- paste(extFields, '}', sep='')
    queryArgs <- c(queryArgs,'rs:extfields'=extFields)
  }

  # create
  if (length(mlDf@.col.defs) > 0) {
    fields <- "{"
    for (i in 1:length(mlDf@.col.defs)) {
      if (nchar(fields) > 1) {
        fields <- paste(fields, ',', sep='')
      }
      fields <- paste(fields, '"', names(mlDf@.col.defs[i]), '":{"fieldDef":"',mlDf@.col.defs[[i]] ,'"}',sep='')
    }
    fields <- paste(fields, '}', sep='')
    queryArgs <- c(queryArgs, 'rs:fields'=fields)
  }
   # do a search
  response <- GET(mlSearchURL, query = queryArgs, authenticate(username, password, type="digest"), accept_json())
  # check that we get an 200
  rContent <- content(response, as = "text")
  if(response$status_code != 200) {
    errorMsg <- paste("statusCode: ",
                      rContent, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))
  }

  if (validate(rContent)) {
    return(fromJSON(rContent, simplifyDataFrame = TRUE)$results)
  } else {
    stop("The call to MarkLogic did not return valid data. The ml.data.frame data could be missing in the database.")
  }

}
# Internal used function that creats new documentsin MarkLogic based on a
# ml.data.frame object.
# Each line is added as a document, put in a collection named after
# myCollection value
.save.ml.data <- function(mlDf, myCollection, directory) {

  conn <- mlDf@.conn
  # get connection imformation
  key <- .rfmlEnv$key[[conn@.id]]
  password <- rawToChar(PKI::PKI.decrypt(conn@.password, key))
  username <- conn@.username
  mlHost <- paste("http://", conn@.host, ":", conn@.port, sep="")
  mlSearchURL <- paste(mlHost, "/v1/resources/rfml.dframe", sep="")

  rfmlCollection <- myCollection
  # generate the directory URI
  if (directory == "") {
    rfmlDirectory <- paste("/rfml/", username, "/", myCollection, "/", sep="")
  } else {
    rfmlDirectory <- directory
  }

  # need to pick start and end from mlDf...
  nStart=mlDf@.start
  nPageLength <- mlDf@.nrows
  queryComArgs <- mlDf@.queryArgs
  queryArgs <- c(queryComArgs, 'rs:start'=nStart,'rs:pageLength'=nPageLength, 'rs:saveDirectory'=rfmlDirectory, 'rs:saveCollection'=rfmlCollection)
  # Need to check if extracted then we could have changed the rows...
  if (mlDf@.extracted) {
    # create a extfields parameter...
    extFields <- "{"
    for (i in 1:length(mlDf@.col.name)) {
      if (nchar(extFields) > 1) {
        extFields <- paste(extFields, ',', sep='')
      }
      extFields <- paste(extFields, '"', mlDf@.col.name[i],
                         '":{"fieldDef":"',mlDf@.col.name[i],
                         '","orgField":"',mlDf@.col.org_name[i],
                         '","orgPath":"',mlDf@.col.org_xpath[i],
                         '","orgFormat":"',mlDf@.col.format[i],'"}',sep='')
    }
    extFields <- paste(extFields, '}', sep='')
    queryArgs <- c(queryArgs,'rs:extfields'=extFields)
  }

  # create
  if (length(mlDf@.col.defs) > 0) {
    fields <- "{"
    for (i in 1:length(mlDf@.col.defs)) {
      if (nchar(fields) > 1) {
        fields <- paste(fields, ',', sep='')
      }
      fields <- paste(fields, '"', names(mlDf@.col.defs[i]), '":{"fieldDef":"',mlDf@.col.defs[[i]] ,'"}',sep='')
    }
    fields <- paste(fields, '}', sep='')
    queryArgs <- c(queryArgs, 'rs:fields'=fields)
  }


  # do a search
  response <- PUT(mlSearchURL, query = queryArgs, authenticate(username, password, type="digest"), accept_json())
  # check that we get an 200
  rContent <- content(response, as = "text")
  if(response$status_code != 204) {
    errorMsg <- paste("statusCode: ",
                      rContent, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))
  }

  return(rfmlCollection)

}
# Internal used function that inserts data.frame data into MarkLogic.
# Each line is added as a document, put in a collection named after
# myCollection value
.insert.ml.data <- function(conn, myData, myCollection, format, directory) {

  # get connection imformation
  key <- .rfmlEnv$key[[conn@.id]]
  password <- rawToChar(PKI::PKI.decrypt(conn@.password, key))
  username <- conn@.username
  mlHost <- paste("http://", conn@.host, ":", conn@.port, sep="")

  mlPostURL <- paste(mlHost, "/v1/documents", sep="")

  rfmlCollection <- myCollection
  # generate the directory URI
  if (directory == "") {
    rfmlDirectory <- paste("/rfml/", username, "/", myCollection, "/", sep="")
  } else {
    rfmlDirectory <- directory
  }

  if (format == "XML") {
    bodyFile <- .generate.xml.body(myData, myCollection, rfmlDirectory)
  } else if (format == "json") {
    bodyFile <- .generate.json.body(myData, myCollection, rfmlDirectory)
  } else {
    stop("Unkown format")
  }

  response <- POST(mlPostURL,  body = upload_file(bodyFile, type = "multipart/mixed; boundary=BOUNDARY"), authenticate(username, password, type="digest"), encode = "multipart", accept_json())
  #suppressWarnings(closeAllConnections())
  suppressWarnings(unlink(bodyFile))

  # message(bodyFile)
  if(response$status_code != 200) {
    rContent <- content(response, as = "text")
    errorMsg <- paste("statusCode: ",rContent, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))
  }

  return(rfmlCollection)

}
# Internal used function that deletes data created by as.ml.data.frame.
.delete.ml.data <- function(myData, directory) {
  conn <- myData@.conn
  # get connection imformation
  key <- .rfmlEnv$key[[conn@.id]]
  password <- rawToChar(PKI::PKI.decrypt(conn@.password, key))
  username <- conn@.username
  mlHost <- paste("http://", .rfmlEnv$conn$host, ":", conn@.port, sep="")

  mlDelURL <- paste(mlHost, "/v1/resources/rfml.dframe", sep="")

  rfmlCollection <- myData@.queryArgs$`rs:collection`
  if (nchar(rfmlCollection) == 0) {
    stop("Can only delete data for a ml.data.frame that has been created using as.mld.data.frame!")
  }
  # generate the directory URI
  if (directory == "") {
    rfmlDirectory <- paste("/rfml/", username, "/", rfmlCollection, "/", sep="")
  } else {
    rfmlDirectory <- directory
  }
  queryArgs <- list('rs:collection'=rfmlCollection, 'rs:directory'=rfmlDirectory)

  response <- DELETE(mlDelURL, query = queryArgs, authenticate(username, password, type="digest"), accept_json())
  if(response$status_code != 204) {
    rContent <- content(response, as = "text")
    errorMsg <- paste("statusCode: ",rContent, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))
  }

  return(TRUE)

}
# executes a statistic function
.ml.stat.func <- function(mlDf, fields, func) {
  conn <- mlDf@.conn
  key <- .rfmlEnv$key[[conn@.id]]
  password <- rawToChar(PKI::PKI.decrypt(conn@.password, key))
  username <- conn@.username
  queryComArgs <- mlDf@.queryArgs

  mlHost <- paste("http://", conn@.host, ":", conn@.port, sep="")
  mlSearchURL <- paste(mlHost, "/v1/resources/rfml.stat", sep="")
  nPageLength <- mlDf@.nrows
  queryArgs <- c(queryComArgs, 'rs:pageLength'=nPageLength, 'rs:statfunc'=func,'rs:fields'=fields)

  response <- GET(mlSearchURL, query = queryArgs, authenticate(username, password, type="digest"), accept_json())
  rContent <- content(response) #, as = "text""
  if(response$status_code != 200) {
    errorMsg <- paste("statusCode: ",
                      rContent, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))
  }
  return(rContent)

}

# Get data for the summary function
.ml.matrix <- function(mlDf, matrixfunc) {
  conn <- mlDf@.conn
  key <- .rfmlEnv$key[[conn@.id]]
  password <- rawToChar(PKI::PKI.decrypt(conn@.password, key))
  username <- conn@.username
  queryComArgs <- mlDf@.queryArgs

  mlHost <- paste("http://", conn@.host, ":", conn@.port, sep="")
  mlSearchURL <- paste(mlHost, "/v1/resources/rfml.matrix", sep="")

  nStart=1
  nPageLength <- mlDf@.nrows

  queryArgs <- c(queryComArgs, 'rs:pageLength'=nPageLength, 'rs:matrixfunc'=matrixfunc)

  # create
  if (length(mlDf@.col.defs) > 0) {
    fields <- "{"
    for (i in 1:length(mlDf@.col.defs)) {
      if (nchar(fields) > 1) {
        fields <- paste(fields, ',', sep='')
      }
      fields <- paste(fields, '"', names(mlDf@.col.defs[i]), '":{"fieldDef":"',mlDf@.col.defs[[i]] ,'"}',sep='')
    }
    fields <- paste(fields, '}', sep='')
    queryArgs <- c(queryArgs, 'rs:fields'=fields)
  }


  # do a search
  response <- GET(mlSearchURL, query = queryArgs, authenticate(username, password, type="digest"), accept_json())
  # check that we get an 200
  rContent <- content(response)
  if(response$status_code != 200) {
    errorMsg <- paste("statusCode: ",
                      rContent, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))
  }

  return(rContent)
}

# executes a Moving Average function
.ml.movavg.func <- function(mlTs, fields, func, n) {
  conn <- mlTs@.conn
  key <- .rfmlEnv$key[[conn@.id]]
  password <- rawToChar(PKI::PKI.decrypt(conn@.password, key))
  username <- conn@.username
  queryComArgs <- mlTs@.queryArgs

  mlHost <- paste("http://", conn@.host, ":", conn@.port, sep="")
  mlSearchURL <- paste(mlHost, "/v1/resources/rfml.movavg", sep="")
  nPageLength <- mlTs@.nrows
  queryArgs <- c(queryComArgs, 'rs:pageLength'=nPageLength, 'rs:avgfunc'=func,'rs:fields'=fields)

  response <- GET(mlSearchURL, query = queryArgs, authenticate(username, password, type="digest"), accept_json())
  rContent <- content(response) #, as = "text""
  if(response$status_code != 200) {
    errorMsg <- paste("statusCode: ",
                      rContent, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))
  }
  return(rContent)

}

.generate.xml.body <- function(data, name, directory) {

  boundary <- "--BOUNDARY"
  contentType <- "Content-Type: application/xml"
  bodyText <- c(boundary,contentType)

  # add metadata
  bodyText <- c(bodyText, "Content-Disposition: inline; category=metadata")
  bodyText <- c(bodyText, "")
  bodyText <- c(bodyText, '<?xml version="1.0" encoding="UTF-8"?>')
  bodyText <- c(bodyText, '<rapi:metadata xmlns:rapi="http://marklogic.com/rest-api">')
  bodyText <- c(bodyText, paste('<rapi:collections><rapi:collection>', name, '</rapi:collection></rapi:collections>', sep=""))
  bodyText <- c(bodyText, '</rapi:metadata>')

  # start loop
  for (i in 1:nrow(data)) {
    bodyText <- c(bodyText,boundary,contentType, paste("Content-Disposition: attachment;filename=",directory,row.names(data[i,]), ".xml" ,sep=""), "")
    myXml <- xmlTree()
    myXml$addTag(name, close=FALSE)
    for (j in names(data)) {
      myXml$addTag(j, data[i, j])
    }
    myXml$closeTag()
    bodyText <- c(bodyText, saveXML(myXml,indent = FALSE,prefix = '<?xml version="1.0"?>'))
    #bodyText <- c(bodyText, saveXML(myXml,indent = FALSE, prefix = ''))
  }
  bodyText <- c(bodyText, "--BOUNDARY--", "")
  # add it
  tf <- tempfile()
  multipartBody <- file(tf, open = "wb")
  #writeLines(text = bodyText, con = multipartBody
  # need to have CRLF no matter of which platform it is running on...
  writeLines(text = bodyText, con = multipartBody, sep="\r\n")
  close(multipartBody)
  #writeBin(bodyText, multipartBody)
  #message(tf)
  return(tf)

}

.generate.json.body <- function(data, name, directory) {

  boundary <- "--BOUNDARY"
  contentType <- "Content-Type: application/json"
  bodyText <- c(boundary,contentType)

  # add metadata
  bodyText <- c(bodyText, "Content-Disposition: inline; category=metadata")
  bodyText <- c(bodyText, "")
  bodyText <- c(bodyText, paste('{"collections" : ["', name, '"] }', sep=""))

  if (is.data.frame(data)) {
    # start loop
    for (i in 1:nrow(data)) {
      bodyText <- c(bodyText,boundary,contentType, paste('Content-Disposition: attachment;filename="',directory,row.names(data[i,]), '.json"', sep=""), "")
      jsonData <- toJSON(data[i,])
      # need to remove the [ ] in the doc, before sending it
      jsonData <- gsub("\\]", "", gsub("\\[", "", jsonData))

      bodyText <- c(bodyText,jsonData)
    }
  } else if (is.character(data)) {
    uploadFiles <- list.files(system.file(data, package = "rfml"))
    for (i in 1:length(uploadFiles)) {
      bodyText <- c(bodyText,boundary,contentType, paste("Content-Disposition: inline;extension=json;directory=", directory,sep=""), "")
      fileName <- system.file(data, uploadFiles[i], package="rfml")
      jsonData <-readChar(fileName, file.info(fileName)$size) #toJSON(data[i,])

      # need to remove the [ ] in the doc, before sending it
      #jsonData <- gsub("\\]", "", gsub("\\[", "", jsonData))

      bodyText <- c(bodyText, jsonData)
    }
  }
  bodyText <- c(bodyText, "--BOUNDARY--", "")
  # add it
  tf <- tempfile()
  multipartBody <- file(tf, open = "wb")
  #writeLines(text = bodyText, con = multipartBody
  # need to have CRLF no matter of which platform it is running on...
  writeLines(text = bodyText, con = multipartBody, sep="\r\n")
  close(multipartBody)
  #writeBin(bodyText, multipartBody)
  #message(tf)
  return(tf)
}
