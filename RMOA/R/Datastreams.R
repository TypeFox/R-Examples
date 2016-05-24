#' @title Datastream objects and methods
#' 
#' @description Reference object of class datastream. This is a generic class which holds general 
#' information about the data stream.\cr
#' Currently streams are implemented for data in table format (streams of read.table, read.csv, read.csv2,
#' read.delim, read.delim2), data in RAM (data.frame, matrix), data in ff (on disk).\cr
#' See the documentation of \code{\link{datastream_file}}, \code{\link{datastream_dataframe}}, \code{\link{datastream_matrix}},
#' and \code{\link{datastream_ffdf}}
#' 
#' @name datastream
#' @export datastream 
#' @docType class
#' @param description The name how the stream is labelled
#' @param args a list with arguments used to set up the stream and used in the datastream methods
#' @return A class of type datastream which contains
#' \describe{
#'   \item{description: }{character with the name how the stream is labelled.}
#'   \item{state: }{integer with the current state at which the stream will read new instances of data}
#'   \item{processed: }{integer with the number of instances already processed}
#'   \item{finished: }{logical indicating if the stream has finished processing all the instances}
#'   \item{args: }{list with arguments passed on to the stream when it is created (e.g. arguments of read.table)}
#' }
#' @seealso \code{\link{datastream_file}}
#' @examples
#' ## Basic example, showing the general methods available for a datastream object
#' x <- datastream(description = "My own datastream", args = list(a = "TEST"))
#' x
#' str(x)
#' try(x$get_points(x))
datastream <- setRefClass(Class="datastream", 
                          fields = list(
                            description = "character",
                            state = "integer",
                            processed = "integer",
                            finished = "logical",
                            args = "list"))
datastream$methods(
  initialize = function(description = "Default datastream object", args = list()){
    .self$description <- description 
    .self$state <- 1L
    .self$processed <- 0L
    .self$finished <- FALSE
    .self$args <- args
    .self
  },
  hasread = function(n){
    n <- as.integer(n)
    .self$state <- .self$state + n
    .self$processed <- .self$processed + n
  },  
  hasfinished = function(){
    .self$finished <- TRUE
  },
  isfinished = function(){
    .self$finished
  },
  reset = function(){
    .self$state <- 0L
    .self$processed <- 0L
    .self$finished <- FALSE
  },
  show = function() {
    cat("Reference object of class", classLabel(class(.self)), "\n")
    cat(" Current state:", state, "\n")
    cat(" Processed:", processed, "\n")
    cat(" Processing finished:", finished, "\n")
  },
  get_points = function(...){
    stop(sprintf("get_points not implemented for class '%s'", classLabel(class(.self))))
  }) 


#' @title File data stream
#' 
#' @description Reference object of class \code{datastream_file}.
#' This is a class which inherits from class \code{datastream} and which can be used to read in a stream
#' from a file. A number of file readers have been implemented, namely 
#' \strong{\code{datastream_table}}, \strong{\code{datastream_csv}}, \strong{\code{datastream_csv2}},
#' \strong{\code{datastream_delim}}, \strong{\code{datastream_delim2}}.\cr
#' See the examples.
#' 
#' @name datastream_file
#' @export datastream_file 
#' @aliases datastream_file datastream_table datastream_csv datastream_csv2 datastream_delim datastream_delim2
#' @docType class
#' @param description The name how the stream is labelled
#' @param FUN The function to use to read in the file. 
#' Defaults to \code{read.table} for \code{datastream_table}, \code{read.csv} for \code{datastream_csv}, 
#' \code{read.csv2} for \code{datastream_csv2}, \code{read.delim} for \code{datastream_delim}, 
#' \code{read.delim2} for \code{datastream_delim2}
#' @param columnnames optional character vector of column to overwrite the column names of the data read in with in \code{get_points}
#' @param file The file to read in. See e.g. \code{read.table}
#' @param ... parameters passed on to \code{FUN}. See e.g. \code{read.table}
#' @return A class of type \code{datastream_file} which contains
#' \describe{
#'   \item{FUN: }{The function to use to read in the file}
#'   \item{connection: }{A connection to the file}
#'   \item{columnnames: }{A character vector of column names to overwrite the column names with in \code{get_points}}
#'   \item{all fields of the datastream superclass: }{See \code{\link{datastream}}}
#' }
#' @section Methods:
#' \itemize{
#'   \item \code{get_points(n)} Get data from a datastream object.
#'      \describe{
#'        \item{n}{integer, indicating the number of instances to retrieve from the datastream}
#'      }
#' }
#' @seealso \code{\link{read.table}}, \code{\link{read.csv}}, \code{\link{read.csv2}}, \code{\link{read.delim}}, \code{\link{read.delim2}}
#' @examples
#' mydata <- iris
#' mydata$Species[2:3] <- NA
#' ## Example of a CSV file stream
#' myfile <- tempfile()
#' write.csv(iris, file = myfile, row.names=FALSE, na = "")
#' x <- datastream_csv(file = myfile, na.strings = "")
#' x
#' x$get_points(n=10)
#' x
#' x$get_points(n=10)
#' x
#' x$stop()
#' 
#' ## Create your own specific file stream
#' write.table(iris, file = myfile, row.names=FALSE, na = "")
#' x <- datastream_file(description="My file defintion stream", FUN=read.table, 
#'  file = myfile, header=TRUE, na.strings="")
#' x$get_points(n=10)
#' x
datastream_file <- setRefClass(Class="datastream_file", 
                          fields = list(connection = "ANY", FUN="function", columnnames = "character", file = "ANY"),
                          contains = "datastream")
datastream_file$methods(
  initialize = function(description="File stream", FUN=read.table, columnnames=character(0), file=character(0), ...){
    callSuper(description=description, args=list(...))
    .self$file <- file
    .self$FUN <- FUN    
    .self$columnnames <- columnnames
    if(length(.self$file) == 1){
      .self$connection <- file(.self$file)      
      open(.self$connection) 
    }
    .self
  },
  show = function() {
    callSuper()
    cat(.self$description, "\n")
    cat(" File:", .self$file, "\n")
    cat(" File size (bytes):", file.info(.self$file)$size, "\n")
  },
  stop = function(){
    close(con=.self$connection)
    .self$hasfinished()
  },
  get_points = function(n=1, ...) {
    fargs <- .self$args
    fargs$file <- .self$connection
    if(.self$processed > 0){
      fargs$header <- FALSE  
    }    
    fargs$nrows <- n
    x <- try(do.call(.self$FUN, args=fargs), silent=TRUE)
    if(inherits(x, 'try-error')){
      .self$hasfinished()
      close(con=.self$connection)
      return(invisible())
    }else{
      if(.self$processed == 0 & !length(.self$columnnames) > 0){
        .self$columnnames <- colnames(x)  
      }    
      .self$hasread(nrow(x))  
      stopifnot(length(.self$columnnames) == ncol(x))
      colnames(x) <- .self$columnnames
      return(x)
    }
  })

#' @name datastream_table
#' @rdname datastream_file
#' @export datastream_table
datastream_table <- setRefClass("datastream_table", contains = "datastream_file", methods = list(
  initialize = function(file=character(0), columnnames=character(0), ...){
    callSuper(description="table file stream", FUN=read.table, columnnames=columnnames, file=file, ...)
  }))
#' @name datastream_csv
#' @rdname datastream_file
#' @export datastream_csv
datastream_csv <- setRefClass("datastream_csv", contains = "datastream_file", methods = list(
  initialize = function(description="csv file stream", FUN=read.csv, columnnames=character(0), file=character(0), ...){
    callSuper(description=description, FUN=FUN, columnnames=columnnames, file=file, ...)
    }))
#' @name datastream_csv2
#' @rdname datastream_file
#' @export datastream_csv2
datastream_csv2 <- setRefClass("datastream_csv2", contains = "datastream_file", methods = list(
  initialize = function(description="csv2 file stream", FUN=read.csv2, columnnames=character(0), file=character(0), ...){
    callSuper(description=description, FUN=FUN, columnnames=columnnames, file=file, ...)
  }))
#' @name datastream_delim
#' @rdname datastream_file
#' @export datastream_delim
datastream_delim <- setRefClass("datastream_delim", contains = "datastream_file", methods = list(
  initialize = function(description="delim file stream", FUN=read.delim, columnnames=character(0), file=character(0), ...){
    callSuper(description=description, FUN=FUN, columnnames=columnnames, file=file, ...)
  }))
#' @name datastream_delim2
#' @rdname datastream_file
#' @export datastream_delim2
datastream_delim2 <- setRefClass("datastream_delim2", contains = "datastream_file", methods = list(
  initialize = function(description="delim2 file stream", FUN=read.delim2, columnnames=character(0), file=character(0), ...){
    callSuper(description=description, FUN=FUN, columnnames=columnnames, file=file, ...)
  }))





#' @title data streams on a data.frame
#' 
#' @description Reference object of class \code{datastream_dataframe}.
#' This is a class which inherits from class \code{datastream} and which can be used to read in a stream
#' from a data.frame.
#' 
#' @name datastream_dataframe
#' @export datastream_dataframe 
#' @aliases datastream_dataframe
#' @docType class
#' @param data a data.frame to extract data from in a streaming way
#' @return A class of type \code{datastream_dataframe} which contains
#' \describe{
#'   \item{data: }{The data.frame to extract instances from}
#'   \item{all fields of the datastream superclass: }{See \code{\link{datastream}}}
#' }
#' @section Methods:
#' \itemize{
#'   \item \code{get_points(n)} Get data from a datastream object.
#'      \describe{
#'        \item{n}{integer, indicating the number of instances to retrieve from the datastream}
#'      }
#' }
#' @seealso \code{\link{datastream}}
#' @examples
#' x <- datastream_dataframe(data=iris)
#' x$get_points(10)
#' x
#' x$get_points(10)
#' x
datastream_dataframe <- setRefClass(Class="datastream_dataframe", 
                          fields = list(data = "data.frame"),
                          contains = "datastream")
datastream_dataframe$methods(
  initialize = function(data){
    callSuper(description="data.frame stream")
    .self$data <- data     
    .self
  },
  show = function() {
    callSuper()
    cat(.self$description, "\n")
    cat(" NROWS:", nrow(.self$data), "\n")
  },
  get_points = function(n=1, ...) {
    if(.self$state > nrow(.self$data)){
      .self$hasfinished()
      return(invisible())
    }else{
      x <- .self$data[.self$state:min(.self$state+n-1, nrow(.self$data)), , drop=FALSE]
      .self$hasread(nrow(x))  
      return(x)
    }
  }) 


#' @title data streams on a matrix
#' 
#' @description Reference object of class \code{datastream_matrix}.
#' This is a class which inherits from class \code{datastream} and which can be used to read in a stream
#' from a matrix.
#' 
#' @name datastream_matrix
#' @export datastream_matrix 
#' @aliases datastream_matrix
#' @docType class
#' @param data a matrix to extract data from in a streaming way
#' @return A class of type \code{datastream_matrix} which contains
#' \describe{
#'   \item{data: }{The matrix to extract instances from}
#'   \item{all fields of the datastream superclass: }{See \code{\link{datastream}}}
#' }
#' @section Methods:
#' \itemize{
#'   \item \code{get_points(n)} Get data from a datastream object.
#'      \describe{
#'        \item{n}{integer, indicating the number of instances to retrieve from the datastream}
#'      }
#' }
#' @seealso \code{\link{datastream}}
#' @examples
#' data <- matrix(rnorm(1000*10), nrow = 1000, ncol = 10)
#' x <- datastream_matrix(data=data)
#' x$get_points(10)
#' x
#' x$get_points(10)
#' x
datastream_matrix <- setRefClass(Class="datastream_matrix", 
                                 fields = list(data = "matrix"),
                                 contains = "datastream")
datastream_matrix$methods(
  initialize = function(data){
    callSuper(description="matrix stream")
    .self$data <- data     
    .self
  },
  show = function() {
    callSuper()
    cat(.self$description, "\n")
    cat(" NROWS:", nrow(.self$data), "\n")
  },
  get_points = function(n=1, ...) {
    if(.self$state > nrow(.self$data)){
      .self$hasfinished()
      return(invisible())
    }else{
      x <- .self$data[.self$state:min(.self$state+n-1, nrow(.self$data)), , drop=FALSE]
      .self$hasread(nrow(x))  
      return(x)
    }
  }) 


#' @title data streams on an ffdf
#' 
#' @description Reference object of class \code{datastream_ffdf}.
#' This is a class which inherits from class \code{datastream} and which can be used to read in a stream
#' from a ffdf from the ff package.
#' 
#' @name datastream_ffdf
#' @export datastream_ffdf 
#' @aliases datastream_ffdf
#' @docType class
#' @param data a data.frame to extract data from in a streaming way
#' @return A class of type \code{datastream_ffdf} which contains
#' \describe{
#'   \item{data: }{The ffdf to extract instances from}
#'   \item{all fields of the datastream superclass: }{See \code{\link{datastream}}}
#' }
#' @section Methods:
#' \itemize{
#'   \item \code{get_points(n)} Get data from a datastream object.
#'      \describe{
#'        \item{n}{integer, indicating the number of instances to retrieve from the datastream}
#'      }
#' }
#' @seealso \code{\link{datastream}}
#' @examples
#' ## You need to load package ff before you can use datastream_ffdf
#' require(ff)
#' irisff <- as.ffdf(factorise(iris))
#' x <- datastream_ffdf(data=irisff)
#' x$get_points(10)
#' x
#' x$get_points(10)
#' x
datastream_ffdf <- setRefClass(Class="datastream_ffdf", 
                               fields = list(data = "ANY"),
                               contains = "datastream")
datastream_ffdf$methods(
  initialize = function(data){
    callSuper(description="ffdf stream")
    .self$data <- data     
    .self
  },
  show = function() {
    callSuper()
    cat(.self$description, "\n")
    cat(" NROWS:", nrow(.self$data), "\n")
  },
  get_points = function(n=1, ...) {
    if(.self$state > nrow(.self$data)){
      .self$hasfinished()
      return(invisible())
    }else{
      x <- .self$data[.self$state:min(.self$state+n-1, nrow(.self$data)), , drop=FALSE]
      .self$hasread(nrow(x))  
      return(x)
    }
  }) 



