#' Convert JSON to BSON Object
#'
#' Converts a JSON string to a mongo BSON object.
#'
#' @param JSON (string) A valid JSON string.
#' @param simplifyVector (FALSE) coerse JSON arrays containing only scalars into a vector.
#' @param ... additional parameters parsed to fromJSON
#'
#' @return A BSON object.
#'
#' @seealso \code{\link{mongo.find}},\cr \code{\link{mongo.bson.from.list}}, \cr \code{\link{mongo.bson}}
#' , \cr \code{\link{fromJSON}}.
#'
#' @examples
#' mongo.bson.from.JSON('{"name" : "Peter"}')
#' mongo.bson.from.JSON('{"_id" : 1}')
#'
#'
#' @export mongo.bson.from.JSON
mongo.bson.from.JSON <- function(JSON, simplifyVector=FALSE, ...){

  if( !jsonlite::validate(I(JSON)) ){
    stop("Not a valid JSON content: ", JSON)
  }

  json_list <- jsonlite::fromJSON(JSON, simplifyVector=simplifyVector, ...)

  if( length(json_list) == 0 ){
    bson <- mongo.bson.empty()
  } else
    bson <- mongo.bson.from.list( as.list( json_list  ) )

  return(bson)
}





#' Convert a mongo.bson object to an R object.
#'
#' Convert a \link{mongo.bson} object to an R object.
#'
#' Note that this function and \code{\link{mongo.bson.from.list}()} do not
#' always perform inverse conversions since \code{mongo.bson.to.Robject}() will
#' convert objects and subobjects to atomic vectors if possible.
#'
#' This function is somewhat schizophrenic depending on the types of the fields
#' in the mongo.bson object. If all fields in an object (or subobject/array)
#' can be converted to the same atomic R type (for example they are all strings
#' or all integer, you'll actually get out a vector of the atomic type with the
#' names attribute set.
#'
#' For example, if you construct a mongo.bson object like such:
#'
#' \preformatted{
#' b <- mongo.bson.from.JSON('{"First":"Joe", "Last":"Smith"}')
#' l <- mongo.bson.to.Robject(b)
#' }
#'
#' You'll get a vector of strings out of it which may be indexed by number,
#' like so:
#'
#' \code{print(l[1]) # display "Joe"}
#'
#' or by name, like so:
#'
#' \code{print(l[["Last"]]) # display "Smith"}
#'
#' If, however, the mongo.bson object is made up of disparate types like such:
#'
#' \preformatted{
#' b <- mongo.bson.from.JSON('{"First":"Joe Smith", "Last":21.5}')
#' l <- mongo.bson.to.Robject(b)
#' }
#'
#' You'll get a true list (with the names attribute set) which may be indexed
#' by number also:
#'
#' \code{print(l[1]) # display "Joe Smith"}
#'
#' or by name, in the same fashion as above, like so
#'
#' \code{print(l[["Name"]]) # display "Joe Smith"}
#'
#' \strong{but} also with the $ operator, like so:
#'
#' \code{print(l$age) # display 21.5}
#'
#' Note that \code{mongo.bson.to.Robject()} operates recursively on subobjects and
#' arrays and you'll get lists whose members are lists or vectors themselves.
#' See \code{\link{mongo.bson.value}()} for more information on the conversion
#' of component types.
#'
#' This function also detects the special wrapper as output by
#' \code{\link{mongo.bson.buffer.append.object}()} and will return an
#' appropriately attributed object.
#'
#' Perhaps the best way to see what you are going to get for your particular
#' application is to test it.
#'
#'
#' @param b (\link{mongo.bson}) The mongo.bson object to convert.
#' @return Best guess at an appropriate R object representing the mongo.bson
#' object.
#' @seealso \code{\link{mongo.bson.from.list}},\cr \code{\link{mongo.bson.to.list}},\cr \link{mongo.bson}.
#' @examples
#'
#' b <- mongo.bson.from.JSON('{"name":"Fred", "city":"Dayton"}')
#'
#' l <- mongo.bson.to.Robject(b)
#' print(l)
#'
#' @export mongo.bson.to.Robject
mongo.bson.to.Robject <- function(b)
  .Call(".mongo.bson.to.list", b)



#' Convert a mongo.bson object to an R list object.
#'
#' Convert a \link{mongo.bson} object to an R list object.
#'
#'
#' @param b (\link{mongo.bson}) The mongo.bson object to convert.
#' @param simplify \link{logical} (default: TRUE); should the bson arrays be simplified to a vectors if possible?
#' If types of values in bson array are heterogeneous or non-primitive, array will be converted into list.
#' @return an R object of the type list
#' @seealso \code{\link{mongo.bson.from.list}}, \code{\link{mongo.bson.to.Robject}},\cr \link{mongo.bson}.
#' @note Now arrays in bson document are 1) converted into unnamed lists 2) If simplify == TRUE,  function tries
#' to turn arrays of primitive types into R vectors.
#' Please see examples below;
#'
#' @examples
#' # arrays will be converted into unnamed lists without any symplifying:
#' l <- list(storageArray = list('value_1', 'value_2'))
#' # Here we construct bson of form {'storageArray':['value_1''value_2']}
#' b <- mongo.bson.from.list(l)
#' # simplify
#' print(mongo.bson.to.list(b, simplify = TRUE))
#' # not simplify
#' print(mongo.bson.to.list(b, simplify = FALSE))
#' # heterogeneous types of array values
#' print(mongo.bson.to.list(mongo.bson.from.list(list(x = list('a', 1))), simplify = TRUE))
#' # identical to call with simplify = F
#' print(mongo.bson.to.list(mongo.bson.from.list(list(x = list('a', 1))), simplify = FALSE))
#' @export mongo.bson.to.list
mongo.bson.to.list <- function(b, simplify = TRUE) {
  stopifnot(is.logical(simplify), (length(simplify) == 1))
  .Call("R_ConvertObject", b, simplify)
}


#' Convert a list to a mongo.bson object
#'
#' Convert a list to a \link{mongo.bson} object.
#'
#' This function permits the simple and convenient creation of a mongo.bson
#' object. This bypasses the creation of a \link{mongo.bson.buffer}, appending
#' fields one by one, and then turning the buffer into a mongo.bson object with
#' \code{\link{mongo.bson.from.buffer}()}.
#'
#' Note that this function and \code{\link{mongo.bson.to.list}()} perform inverse conversions.
#'
#' @param lst (list) The list to convert.
#'
#' This \emph{must} be a list, \emph{not} a vector of atomic types; otherwise,
#' an error is thrown; use \code{as.list()} as necessary.
#' @return (\link{mongo.bson}) A mongo.bson object serialized from \code{lst}.
#' @note Function converts unnamed R lists into bson arrays.
#' It is very easy to construct bson object of any form using this function and list.
#' @seealso \code{\link{mongo.bson.to.list}},\cr \link{mongo.bson},\cr
#' \code{\link{mongo.bson.destroy}}.
#' @examples
#'
#' lst <- list(name="John", age=32)
#' b <- mongo.bson.from.list(lst)
#' # the above produces a BSON object of the form:
#' # { "name" : "John", "age" : 32.0 }
#'
#' # Convert a vector of an atomic type to a list and
#' # then to a mongo.bson object
#' v <- c(president="Jefferson", vice="Burr")
#' b <- mongo.bson.from.list(as.list(v))
#' # the above produces a BSON object of the form:
#' # { "president" : "Jefferson", "vice" : "Burr" }
#' # Let's try to construct bson with array.
#' # This one
#' mongo.bson.from.list(list(fruits = list('apple', 'banana', 'orange')))
#' # will produce a BSON object of the form:
#' # {"fruits" : ["apple", "banana", "orange"]}
#' @export mongo.bson.from.list
mongo.bson.from.list <- function(lst){
  if(length(lst)){
    .Call(".mongo.bson.from.list", lst)
  } else {
    mongo.bson.empty()
  }
}

#' Convert a data.frame to a mongo.bson object
#'
#' Convert a data.frame to a \link{mongo.bson} object.
#'
#' This function permits the simple and convenient creation of a mongo.bson
#' object.  This bypasses the creation of a \link{mongo.bson.buffer}, appending
#' fields one by one, and then turning the buffer into a mongo.bson object with
#' \code{\link{mongo.bson.from.buffer}()}.
#'
#' @param df (data.frame) The data.frame to convert.
#'
#' This \emph{must} be a data.frame, \emph{not} a vector of atomic types; otherwise,
#' an error is thrown; use \code{as.data.frame()} as necessary.
#' @return (\link{mongo.bson}) A mongo.bson object serialized from \code{data.frame}.
#' @seealso \code{\link{mongo.bson.to.list}},\cr \link{mongo.bson},\cr
#' \code{\link{mongo.bson.destroy}}.
#'
#' @examples
#' df <- data.frame(name=c("John", "Peter"), age=c(32,18))
#' b <- mongo.bson.from.df(df)
#'
#' @export mongo.bson.from.df
mongo.bson.from.df <- function(df){
  # Put each row to a seperate list item
  #data_list = apply(dataframe,1,as.list) # Cannot be done like this because then the data type is coersed

  data_list <- vector("list", nrow(df))
  for( i in 1:nrow(df) ) data_list[[i]] <- as.list(df[i, ])
  # Iterate over the table and create the BSON object
  bson_data <- lapply(data_list,function(x){
    idx <- 1
    names <- names(x)
    buf <- mongo.bson.buffer.create()

    lapply(x,function(y) {
      mongo.bson.buffer.append(buf, names[idx], y)
      idx <<- idx+1
    })

    mongo.bson.from.buffer(buf)
  })


  return(bson_data)

}
