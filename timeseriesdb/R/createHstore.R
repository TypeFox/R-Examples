#' Function to Create Hstore Key Value Pair Mapping 
#' 
#' This function creates a key value pair mapping from a time series object. 
#' It returns an hstore object that can be inserted to a PostgreSQL database relation field of type hstore. 
#' @author Matthias Bannert
#' @title Create Hstore
#' @param x a time series object, a two column data frame or object of S3 class
#' miro (meta information for R objects).
#' @param ... optional arguments, fct = TRUE create text expressions of hstore function calls.
#' also for data.frames key_pos and value_pos could be given if they are different from 1 and 2. 
#'  e.g. position of the key col and
#' pasition of the value col in a data.frame.
#' @examples
#' ts1 <- ts(rnorm(100),start = c(1990,1),frequency = 4)
#' createHstore(ts1)
#' 
#' @export
createHstore <- function(x,...) UseMethod("createHstore")

#' @rdname createHstore
#' @export
createHstore.ts <- function(x,...){
  tm <- time(x)
  paste(sprintf('"%s"=>"%s"',
                zooLikeDateConvert(tm),
                as.character(x)),
        collapse=",")
}

#' @rdname createHstore
#' @export
createHstore.data.frame <- function(x,...){
  # only allow to cols because its KEY => VALUE
  if(!exists('key_pos')) key_pos <- 1
  if(!exists('value_pos')) value_pos <- 2
  
  stopifnot(ncol(x) == 2)
  
  paste(sprintf('"%s"=>"%s"',
                as.character(x[,key_pos]),
                as.character(x[,value_pos])),
        collapse=",")
}



#' @rdname createHstore
#' @export
createHstore.list <- function(x,...){
  dot_args <- list(...)
  # check if list is more than 2 dim
  if(getListDepth(x) != 1) stop('Only key-value pairs are accepted,
                         this list has too many dimensions!') 
  
  if(is.null(names(x))) stop('Only named lists are accepted.')
  
  # the => operator is deprecated in 
  # Postgres so if you want to use the new version function
  # based version use fct = T
  # the operator will be kept alive as long as postgres does 
  # the same 
  deprecated_hstore_operator <- paste(sprintf('"%s"=>"%s"',
                                              names(x),
                                              as.character(unlist(x))),
                                      collapse=",")
  
  if(exists("fct",dot_args)){
    if(dot_args$fct){
      paste(sprintf("hstore('%s','%s')",
                    names(x),
                    as.character(unlist(x))),
            collapse="||")  
    } else{
      deprecated_hstore_operator
    }
    
  } else {
    deprecated_hstore_operator
  }
  
}

