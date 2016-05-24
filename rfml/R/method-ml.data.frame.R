#' Creates a \link{ml.data.frame} object
#'
#' A ml.data.frame object is an abstraction layer of data stored in a MarkLogic Server database. It is created based
#' on the provided query, collection, directory and/or fiedlFilter parameters. For query and fieldFilter
#' parameters see details section. It present data in MarkLogic Server in a tabular format.
#' The ml.data.frame object enables many of the operations that can be used with a data.frame object.
#'
#' The query parameter are using the \href{http://docs.marklogic.com/guide/search-dev/string-query}{string search grammar}
#' for searching for data, all of the syntax is supported except contstraints. This enables searches such as
#' "dog AND cat"  or "dog NEAR cat". The search is always done on all fields in the data, for a more precise search
#' use the fieldFilter.
#'
#' fieldFilter enables filtering on a specific field using comparison operators can be used. For
#' the ">"  "<"  "!=" "<=" ">=" operators there muset exist a
#' \href{http://docs.marklogic.com/guide/admin/range_index#id_93351}{element range index}
#' on the source field or a error will be raised, element range index can be created using the \link{ml.add.index}
#' function. "==" operator will always work since it does not depend of range indexes.
#'
#' @seealso \code{\link[rfml]{as.data.frame}} for pulling data, \code{\link{as.ml.data.frame}} for uploading data, \code{\link{rm.ml.data.frame}} for delete uploaded data
#'
#' @param conn A \link{ml.conn-class} object created by \link{ml.connect}
#' @param query The query string used to define the result, see details for more information about syntax.
#' @param fieldFilter Field level filtering. Multiple field filters are separated by , See details for limitations.
#' @param ns A character with the namespace URI to be used with fieldFilter, default is none
#' @param collection A list of collection URI:s to filter on.
#' @param directory A list of directory URI:s to filter on.
#' @param relevanceScores TRUE/FALSE. If the result attributes score, confidence and fitness should be included. Default is FALSE
#' @param docUri TRUE/FALSE. If the uri of the documents in the results should be included. Default is FALSE.
#' @return A ml.data.frame object.
#' @examples
#' \dontrun{
#'  library(rfml)
#'  localConn <- ml.connect()
#'  # create a ml.data.frame based on a search
#'  mlIris <- ml.data.frame(localConn, "setosa")
#'  # using search and collection filtering
#'  mlIris <- ml.data.frame(localConn, "setosa", collection = "iris")
#'  # using field filter
#'  mlIris <- ml.data.frame(localConn, fieldFilter = "Species == setosa")
#' }
#' @export
ml.data.frame <- function (conn, query="", fieldFilter="", ns = "NA", collection = c(), directory = c(), relevanceScores = FALSE, docUri = FALSE)
{
  # if (length(.rfmlEnv$conn) != 4) {
  #   stop("Need create a connection object. Use ml.connect first.")
  # }

  if (class(conn) != "ml.conn" || missing(conn)) {
    stop("Need a valid ml.conn object. Use ml.connect to create one.")
  }
  # get data from ML
  # we need to create a "unique" name for the frame that we use to save the resultset
  dframe <- format(Sys.time(),"%Y%m%d%H%M%S")
  # need to check that the key exits...
  key <- .rfmlEnv$key[[conn@.id]]
  password <- rawToChar(PKI::PKI.decrypt(conn@.password, key))
  username <- conn@.username

  mlHost <- paste("http://", conn@.host, ":", conn@.port, sep="")
  mlSearchURL <- paste(mlHost, "/v1/resources/rfml.dframe", sep="")
  nPageLength=30

    # These are the arguments that are common to all calls to MarkLogic
  queryComArgs <- list('rs:q'=query, 'rs:relevanceScores'=relevanceScores, 'rs:docUri' = docUri)

  # fieldQuery
  # operators:  "==" ">"  "<"  "!=" "<=" ">=", only == is currently supported, all the rest require Range Indexs
  if (nchar(fieldFilter) > 0) {
    # seperated by ,
    fieldExprs <- unlist(strsplit(fieldFilter, ",", fixed = TRUE))
    fieldQuery <- "{"
    for (i in 1:length(fieldExprs)) {
      # get the values on both side of the operator
      fieldExpr <- unlist(strsplit(fieldExprs[i], "==|>|<|!=", perl = TRUE))
      if (length(fieldExpr) != 2) {
        stop("Need to provide a valid expression")
      }
      # get the operator
      opMatch <- regexpr("==|>|<|!=", fieldExprs[i],perl = TRUE)
      opIndex <- opMatch[[1]]
      opLength <- attr(opMatch, "match.length")
      op <- trimws(substr(fieldExprs[i], opIndex, opIndex+opLength))
      if (i > 1) {
        fieldQuery <- paste(fieldQuery, ',', sep='')
      }
      fieldQuery <- paste(fieldQuery, '"', trimws(fieldExpr[1]),
                          '":{"value":"',trimws(fieldExpr[2]),
                          '","operator":"', op ,'","orgPath":"","orgFormat":"","xmlns":"',
                          ns, '"}',sep='')
    }
    fieldQuery <- paste(fieldQuery, '}', sep='')
    queryComArgs <- c(queryComArgs,  'rs:fieldQuery'=fieldQuery)
  }
  # Collection and/or directory filtering
  if (length(collection) > 0) {
    strColl <- ''
    for (i in 1:length(collection)) {
      if (i>1) {
        strColl <- paste(strColl, ',', sep='')
      }
      strColl <- paste(strColl, collection[i], sep='')
    }
    queryComArgs <- c(queryComArgs, 'rs:collection'=strColl)
  }
  if (length(directory) > 0) {
    strDir <- ''
    for (i in 1:length(directory)) {
      if (i>1) {
        strDir <- paste(strDir, ',', sep='')
      }
      strDir <- paste(strDir, directory[i], sep='')
    }
    queryComArgs <- c(queryComArgs, 'rs:directory'=strDir)
  }

  queryArgs <- c(queryComArgs,'rs:start'=1, 'rs:pageLength'=nPageLength, 'rs:return'="meta")
  # do a search
  response <- GET(mlSearchURL, query = queryArgs, authenticate(username, password, type="digest"), accept_json())

  # get the content
  rContent <- content(response)

  if(response$status_code != 200) {
    errorMsg <- paste("statusCode: ",
                      rContent$errorResponse$statusCode,
                      ", status: ", rContent$errorResponse$status,
                      ", message: ", rContent$errorResponse$message, sep="")
    stop(paste("Ops, something went wrong.", errorMsg))
  }
  if (rContent$nrows == 0) {
    stop("Search did not produce any result");
  }

  res <- new("ml.data.frame")
  res@.name <- dframe
  res@.conn <- conn
  res@.queryArgs <- queryComArgs
  res@.nrows <- as.integer(rContent$nrows)
  res@.start <- 1L
  res@.extracted=FALSE
  fieldList <- rContent$dataFrameFields
  fieldNames <- c()
  fieldTypes <- c()
  fieldOrgNames <- c()
  fieldOrgXPaths <- c()
  fieldFormat <- c()
  fieldXmlns <- c()
  for (i in 1:length(fieldList)) {
    fieldNames[i] <-  as.character(attributes(fieldList[i]))
    fieldTypes[i] <- fieldList[[i]]$fieldType
    fieldOrgNames[i] <- fieldList[[i]]$orgField
    fieldOrgXPaths[i] <- fieldList[[i]]$orgPath
    fieldFormat[i] <- fieldList[[i]]$orgFormat
    if (!is.null(fieldList[[i]]$xmlns)) {
      fieldXmlns[i] <- fieldList[[i]]$xmlns
    }
  }
  res@.col.name <- fieldNames
  res@.col.data_type <- fieldTypes
  res@.col.org_name <- fieldOrgNames
  res@.col.org_xpath <- fieldOrgXPaths
  res@.col.format <- fieldFormat
  if (!is.null(fieldXmlns)) {
    res@.col.xmlns <- fieldXmlns
  }
  res@.col.defs <- list()
  return(res);

}

################ Generic methods for upload and download of data ############################

################ as.data.frame ############################
#' Pull data from MarkLogic server based on a \link{ml.data.frame} object and return it as a data.frame.
#'
#' @param x a  \link{ml.data.frame} object
#' @param max.rows maximum rows to return. Default all rows.
#' @param ... Not used.
#' @examples
#' \dontrun{
#'  library(rfml)
#'  localConn <- ml.connect()
#'  # create a ml.data.frame based on a search
#'  mlIris <- ml.data.frame(localConn, "setosa")
#'  lIris <- as.data.frame(mlIris)
#'  }
#' @aliases as.data.frame
#' @seealso \code{\link{ml.data.frame}}, \code{\link{as.ml.data.frame}} for uploading data, \code{\link{rm.ml.data.frame}} for delete uploaded data
#' @concept array
#' @export
setMethod("as.data.frame", signature(x="ml.data.frame"),
          function (x, max.rows=NULL, ...) {
            if (is.null(max.rows)) {
              max.rows <- 0
            }
            result <- return(.get.ml.data(x,max.rows))
          }
)
#' Upload data in a data.frame object or create data based on a \link{ml.data.frame} object
#'
#' The function will upload the data within a data.frame object or create data in MarkLogic Server
#' based on a \link{ml.data.frame} object. Data created based on \link{ml.data.frame} will be flat and
#' fields will have the same names as in the .col.name slot. See details for more information about how
#' data is created.
#'
#' When data is uploaded or created it will be stored as json documents default, the format parameter controls, and
#' Document URIs, the identifier of a document, is generated based on the string "rfml", the rowname if a data.frame
#' or a counter if it is a ml.data.frame, the loged in username and the name parameter, for example /rfml/admin/iris/.
#' The documents will also belong to a collection named after tne name parameter.
#'
#' @param conn A ml.conn object that has a valid connection to a MarkLogic Server
#' @param x a Data Frame or ml.data.frame object.
#' @param name The name of the object.
#' @param format The format od the documents that is created, json or XML. Default is json
#' @param directory The directory to save the documents, needs to start and end with a /. Default saved to /rfml/[username]/[name]/
#' @return A ml.data.frame object.
#' @examples
#' \dontrun{
#'  library(rfml)
#'  ml.connect()
#'  # create a ml.data.frame based on the iris data set
#'  mlIris <- as.ml.data.frame(iris, "iris")
#' }
#' @seealso \code{\link{ml.data.frame}}, \code{\link[rfml]{as.data.frame}} for pulling data, \code{\link{rm.ml.data.frame}} for delete uploaded data
#' @concept array
#' @export
as.ml.data.frame <- function (conn, x, name, format = "json", directory = "") {

  if (is.data.frame(x)) {
    if (class(conn) != "ml.conn" || missing(conn)) {
      stop("Need a valid ml.conn object. Use ml.connect to create one.")
    }
    suppressWarnings(rfmlCollection <- .insert.ml.data(conn, x, name, format, directory))
  } else if (is.ml.data.frame(x)) {
    #stop("Only objects of ml.data.frame type are supported!")
    rfmlCollection <- .save.ml.data(x, name, directory)
    conn <- x@.conn
  } else {
    stop("Only objects of data.frame or ml.data.frame type are supported!")
  }

  # create a ml.data.frame object based on a collection search
  return(ml.data.frame(conn, collection=c(rfmlCollection)));
}

#' Remove the data of a ml.data.frame object in MarkLogic server database.
#'
#' Removes the data that whas saved to MarkLogic server database using the \link{as.ml.data.frame} function.
#' If using a directory parameter it that call the same value needs to be provided for this function.
#' The function will also delete the x object form tne R environment.
#'
#' @param x a ml.data.frame object.
#' @param directory Optional. The directory where the data is stored, needs to start and end with a /.
#' @return A ml.data.frame object.
#' @examples
#' \dontrun{
#'  rm.ml.data.frame(mlIris)
#' }
#' @seealso \code{\link{ml.data.frame}}, \code{\link{as.ml.data.frame}} for uploading data, \code{\link[rfml]{as.data.frame}} for pulling data
#' @concept array
#' @export
rm.ml.data.frame <- function(x, directory = "" ){
  if (!is.ml.data.frame(x)) {
    stop("Only objects of ml.data.frame type are supported!")
  }
  call <- match.call()
  if(.delete.ml.data(x, directory)) {
    retMsg <- paste("Data for ", call$x, " has been deleted", sep="")
  } else {
    retMsg <- paste("Could not delete data for ", call$x, sep="")
  }
  message(retMsg)
  TRUE
}
################ [ ############################
#' Extract subsets of a ml.data.frame
#'
#' Extract subset of columns and/or rows of a ml.data.frame. When extracting rows a ml.col.def
#' referense can be used or a search text, see \link{ml.data.frame} for query string grammar.
#' See details for limitations when using a reference.
#' The row filtering will be used togheter with the existing query of the ml.data.frame
#'
#' When extracting rows using ml.col.def comparison operators can be used. For
#' the ">"  "<"  "!=" "<=" ">=" operators there muset exist a
#' \href{http://docs.marklogic.com/guide/admin/range_index#id_93351}{element range index}
#' on the source field or a error will be raised, element range index can be created using the \link{ml.add.index}
#' function. "==" operator will always work since it does not depend of range indexes.
#'
#' @param x a ml.data.frame from which to extract element(s).
#' @param i,j Indices specifying elements to extract. Indices are 'numeric' or 'character' vectors or empty (missing) or 'NULL'.
#' @param ... Not used.
#' @param drop Not used.
#' @return A ml.data.frame object is returned
#' @examples
#' \dontrun{
#'  library(rfml)
#'  localConn <- ml.connect()
#'  # create a ml.data.frame based on the iris data set
#'  mlIris <- as.ml.data.frame(localConn, iris, "iris")
#'  # select first three columns
#'  mlIris2 <- mlIris[1:3]
#'  # same
#'  mlIris2 <- mlIris[,1:3]
#'  # same
#'  mlIris2 <- mlIris[,c("Sepal.Length","Sepal.Width","Petal.Length")]
#'  # select first three columns for all rows with Spieces = setosa
#'  mlIris2 <- mlIris[mlIris$Species=="setosa", 1:3]
#'  # select all columns for all rows with Spieces = setosa
#'  mlIris2 <- mlIris[mlIris$Species=="setosa",]
#'  # select all columns for all rows with "setosa" in any column
#'  mlIris2 <- mlIris["setosa",]
#' }
#' @concept array
#' @aliases [,ml.data.frame-method
#' @export
setMethod("[", c(x = "ml.data.frame", i="ANY", j="ANY"),
          function (x, i, j,..., drop=NA)
          {
            colArg <- NULL
            rowArg <- NULL
            cols <- c()
            n <- nargs()
            if (n == 1) {
              stop("Argument is missing!")
            }
            if (n == 2) {
              # select columns -> mlDf[1]/mlDf[1:4]/mlDf["column"]/mlDf[c("col1","col2)]
              colArg <- i
            } else if (n == 3) {
              # several cases
              if (missing(i)) {
                # if i is missing -> mlDf[, 1:2]
                # select columns based on j
                colArg <- j
              } else if (missing(j)) {
                # if j is missing -> mlDf[1:2,] -> row selection...
                rowArg <- i
              } else {
                # else mlDf[1:2,1:4]
                # i is row selection ...
                colArg <- j
                rowArg <- i
              }
            }
            # rows selection
            if (!is.null(rowArg)) {
              newQueryArgs <- x@.queryArgs
              if (is.numeric(rowArg)) {
                stop("row numbering is not allowed")
              } else if (class(rowArg)=="ml.col.def") {
                if(rowArg@.type!='logical' && validate(toJSON(rowArg@.expr))) {
                  stop("Column expression must resolve into a boolean value for row selection.")
                }
                newQueryArgs <- c(newQueryArgs, 'rs:fieldQuery'=rowArg@.expr)
                # we should update the row count since we are changing the definition...
                # need a estimate function....
                #x@.nrows <- NA
              } else if (is.character(rowArg)) {
                qText <- newQueryArgs$`rs:q`
                if (nchar(qText) > 0) {
                  qText <- paste(qText, " AND ", sep="")
                }
                qText <- paste(qText, rowArg, sep="")
                newQueryArgs$`rs:q` <- qText
              } else {
                stop("row object does not specify a subset")
              }
              x@.queryArgs <- newQueryArgs
            }
            # column selection.
            # must verify that we handle added columns as well.
            if (!is.null(colArg)) {
              if (is.numeric(colArg)) { # column selection
                # using  x[1:2]
                cols <- c(cols,as.integer(colArg))
              } else if (!is.integer(colArg)) {
                # using column names
                if (is.character(colArg)){
                  for (n in colArg){
                    # get the index of columnname
                    cols <- c(cols, which(names(x)==n))
                    # need to check if the name exists...
                  }
                }
              }
              if (!is.null(x@.col.name) ) {
                x@.col.name <- x@.col.name[cols]
                x@.col.data_type <- x@.col.data_type[cols]
                x@.col.org_name <- x@.col.org_name[cols]
                x@.col.org_xpath <- x@.col.org_xpath[cols]
                x@.col.format <- x@.col.format[cols]
                x@.col.xmlns <- x@.col.xmlns[cols]
                # need to check if selected column is part of .col.defs
                colDefs <- list()
                if (length(x@.col.defs) > 0) {
                  for (i in 1:length(x@.col.name)) {
                    if (!is.null(x@.col.defs[[x@.col.name[i]]])) {
                      colDefs <- c(colDefs, x@.col.defs[x@.col.name[i]])
                    }
                  }
                }
                x@.col.defs <- colDefs
                x@.extracted <- TRUE
              }
            }
            return(x)
          }
)

################ $ ############################
#' Returns a \link{ml.data.frame} field as a \link{ml.col.def-class}
#'
#' @param x an \link{ml.data.frame} object
#' @param name field name
#' @return \link{ml.col.def-class} object
#' @export
setMethod("$", signature(x = "ml.data.frame"),
          function(x, name) {
            if(!(name %in% x@.col.name)) {
              stop("Column not found in ml.data.frame.")
            }
            # pickup the data type
            i <- which(x@.col.name %in% name)
            dataType <- x@.col.data_type[i]
            orgName <- x@.col.org_name[i]
            colFormat <- x@.col.format[i]
            # check if the column are a added or already existing
            if(is.null(x@.col.defs[[name]])) {
              return(new(Class="ml.col.def",.expr=paste("rfmlResult[\'",name, "\']", sep=''),
                         .name=name,.data_type=dataType, .org_name=orgName, .format=colFormat,
                         .parent=x,.type="field",.aggType="none"));
            } else {

              return(new(Class="ml.col.def",.expr=x@.col.defs[[name]],.name=name,.data_type=dataType,
                         .org_name=orgName, .format=colFormat,
                         .parent=x,.type="expr",.aggType="none"));
            }

          }
)
################ $<- ############################
#' Adds a new ml.data.frame field as a \link{ml.col.def-class}
#'
#' The fields  only exists within the object and are not created at the database side.
#'
#' @param x A ml.data.frame object
#' @param name Name of the new field
#' @param value The value for the new field. Typical a expression
#' @return \link{ml.col.def-class} object
#' @examples
#' \dontrun{
#'  library(rfml)
#'  locConn <- ml.connect()
#'  # create a ml.data.frame based on the iris data set
#'  mlIris <- as.ml.data.frame(locConn, iris, "iris")
#'  # create a field based on an existing
#'  mlIris$newField <- mlIris$Petal.Width
#'  # create a field based calculation on existing
#'  mlIris$newField2 <- mlIris$Petal.Width + mlIris$Petal.Length
#'  # create a field based on an previous created
#'  mlIris$newField3 <- mlIris$Petal.Width + 10
#'  mlIris$abs_width <- abs(mlIris$Petal.Width)
#' }
#' @export
setMethod("$<-", signature(x = "ml.data.frame"),
          function(x, name,value) {
             if(is.null(value)) {
               #remove col def
               if(!is.null(x@.col.defs[[name]])) {
                 x@.col.defs[[name]] <- NULL;
               }
               x@.col.name <- setdiff(x@.col.name,name)
             } else {
               if(!inherits(value,"ml.col.def"))
                 stop("Column definition is not valid for a ml.data.frame, please refer to the documentation of ml.data.frame for details on usage.")
               if(((value@.parent)@.name != x@.name))
                stop("Column defintions are only allowed on the same ml.data.frame.");

               if(value@.aggType!='none') {
                 stop("Cannot add column that contains aggregation term to ml.data.frame.")
               }
               # if we are refering an existing column ie x$field
               # then value is the defenition of that field.
               x@.col.defs[[name]]<-value@.expr;
               if(!(name %in% x@.col.name)) {
                 x@.col.name<-c(x@.col.name,name);
                 x@.col.data_type<-c(x@.col.data_type, value@.data_type)
               }
              }

            return(x);
          }
)

#' Check if an object is of type ml.data.frame
#'
#' This function checks if the input is of type \link{ml.data.frame}.
#'
#' @param x The input can be of any type.
#' @return True if it is a ml.data.frame object. False otherwise.
#' @export
is.ml.data.frame <-
  function(x) {
    return(inherits(x, "ml.data.frame"))
  }
################ dim ############################
#' Dimensions of an \link{ml.data.frame} object
#'
#' @param x an ml.data.frame object
#' @export
setMethod("dim", signature(x="ml.data.frame"),
          function(x) {
              return(c(x@.nrows, length(x@.col.name)))
          }
)
################ colnames ############################
#' Column Names of an \link{ml.data.frame} object
#'
#' @param x an ml.data.frame object
#' @export
setMethod("colnames", signature(x="ml.data.frame"),
          function(x) { x@.col.name }
)
################ head ############################
#' Return the First Part of an \link{ml.data.frame}
#'
#' @param x an ml.data.frame object
#' @param n a single integer. The number of rows to return, default is 6
#' @param ... not used
#' @export
setMethod("head", signature(x="ml.data.frame"),
          function(x, n = 6, ...) {
            #if (length(x@.conn) != 4) {
            #  stop("Need create a connection object. Use ml.connect first.")
            #}
            if (n >= 0) {
              return(.get.ml.data(x,n))
            } else {
              nr <- nrow(x)
              n <- abs(n)
              #ans <- idaQuery(idadf.query(x), " FETCH FIRST ", format(nr - n, scientific = FALSE), " ROWS ONLY")

              #if ((nr-n) != 0) rownames(ans) <- 1:(nr-n);
              #return(ans)
            }
          }
)

################ Basic ml.data.frame functions and methods ############################
################ print ############################
#' Prints information of a \link{ml.data.frame} object
#'
#' @param x an ml.data.frame object
#' @export
setMethod("print", signature(x="ml.data.frame"),
          function (x) {
            print(format(x@.queryArgs))
          }
)
################ show ############################
#' Prints information of a \link{ml.data.frame} object
#'
#' @param object an ml.data.frame object
#' @export
setMethod("show", signature(object="ml.data.frame"),
          function (object) {
            print(format(object@.queryArgs))
          }
)

################ names ############################
#' Shows field names of a \link{ml.data.frame} object
#'
#' @param x an ml.data.frame object
#' @export
setMethod("names", signature(x="ml.data.frame"),
          function(x) { return(x@.col.name) }
)

