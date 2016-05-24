## -----------------------------------------------------------------------
## Private variables of the package
## -----------------------------------------------------------------------

# A package specific enviroment, used to store the RSA key and connection info
.rfmlEnv <- new.env()
.rfmlEnv$key <- list()
# name of libs used
.rfmlEnv$mlLibs <- c("rfmlUtilities", "xml2json")
# name of exstentions used
.rfmlEnv$mlExts <- c("rfml.dframe", "rfml.lm", "rfml.stat", "rfml.matrix", "rfml.collection", "rfml.arules", "rfml.check")

#' An S4 class to represent a connection to a MarkLogic Server Database
#'
#' @slot .id A integer with the connection number.
#' @slot .host A string with the MarkLogic Server hostname or ip-adress
#' @slot .port A string with the port number to the HTTP server for the MarkLogic Databse used
#' @slot .mlversion A string with the version of the MarkLogic Server
#' @slot .username A string with username
#' @slot .password Encrypted password
setClass("ml.conn",
         slots=c( .id="integer",
                  .host="character",
                  .port="character",
                  .mlversion="character",
                  .username="character",
                  .password="raw")
         )

#' An S4 class to represent a ml.data.frame.
#'
#' @slot .name A string with the internal name for the ml.data.frame
#' @slot .conn The \link{ml.conn-class} object that was created with ml.connect
#' @slot .queryArgs A list with parameters used to query MarkLogic Server
#' @slot .start A integer with the index of the first result to get
#' @slot .nrows A integer with the number of rows in the result
#' @slot .extracted A logical value indicating if we have selected a subset of fields
#' @slot .col.name A character vector with the field names
#' @slot .col.data_type A character vector with the data types of the fields
#' @slot .col.org_name A character vector with the original names of fields in the source documents
#' @slot .col.org_xpath A character vector with the xpath to the original names in the source documents
#' @slot .col.format A character vector withthe  source document format XML/JSON
#' @slot .col.xmlns A character vector with the namespace for the source document
#' @slot .col.defs  A list of \link{ml.col.def-class} added fields
setClass("ml.data.frame",
         slots=c(
           .name="character",
           .conn="ml.conn",
           .queryArgs="list", # parameters used to query ML
           .start="integer", # the index of the first result to get
           .nrows="integer",  # the number of rows in the result or maximum number or result to get
           .extracted="logical", # if we have selected a subset of columns
           .col.name="character", # column names
           .col.data_type = "character", # column types
           .col.org_name = "character", # name of field in source document
           .col.org_xpath = "character",# xpath in the source document
           .col.format = "character", # source document format XML/JSON
           .col.xmlns = "character", # the namespace for the source document
           .col.defs = "list" # added columns
          )
    )

#' An S4 class to represent a ml.col.def.
#'
#' @slot .expr A string with expresion that define the ml.col.def
#' @slot .parent Pointer to the \link{ml.data.frame-class} object that the field belongs to
#' @slot .type A string with the type of field
#' @slot .name A string with name of the field
#' @slot .data_type A string with the data type of the field
#' @slot .org_name A string with the original names of field
#' @slot .format A string with the format of the source field
#' @slot .xmlns A string with the namespace of the source field
#' @slot .aggType A string
setClass("ml.col.def",
         slots=c(.expr="character",
                 .parent="ml.data.frame",
                 .type = "character", # column types
                 .name = "character",
                 .data_type = "character",
                 .org_name = "character",
                 .format = "character",
                 .xmlns = "character",
                 .aggType="character"
                )
         )

# setClass("ml.ts",
#          slots=c(
#            .start="integer", # the index of the first result to get
#            .end="integer",  # the number of rows in the result or maximum number or result to get
#            .frequency="integer",
#            .queryArgs="list", # parameters used to query ML
#            .field="character" # field used in Time series
#          )
# )
