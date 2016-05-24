
#' Create an empty mongo.bson object
#' 
#' Returns an empty mongo.bson object. mongo.bson objects have "mongo.bson" as
#' their class and contain an externally managed pointer to the actual data.
#' This pointer is stored in the "mongo.bson" attribute of the object.
#' 
#' 
#' @return An empty mongo.bson object
#' @seealso \link{mongo.bson}
#' @examples
#' 
#' # Use an empty mongo.bson for the query object which matches everything.
#' # This happens to be the default value for the query
#' # parameter to mongo.count,  but we explicity use mongo.bson.empty()
#' # here for an example.
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo))
#'     print(mongo.count(mongo, "test.people", query=mongo.bson.empty()))
#' 
#' @export mongo.bson.empty
mongo.bson.empty <- function()
    .Call(".mongo.bson.empty")



#' Get the size of a mongo.bson object
#' 
#' Get the number of bytes taken up by the BSON data attached to the mongo.bson
#' object
#' 
#' 
#' @param b (\link{mongo.bson}) the mongo.bson object to examine.
#' @return (integer) the number of bytes taken up by the BSON data attached to
#' the mongo.bson object.
#' @seealso \link{mongo.bson}
#' @examples
#' 
#' # should report 5
#' print(mongo.bson.size(mongo.bson.empty()))
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append(buf, "name", "Fred")
#' mongo.bson.buffer.append(buf, "city", "Dayton")
#' y <- mongo.bson.from.buffer(buf)
#' # should report 37
#' print(mongo.bson.size(y))
#' 
#' @export mongo.bson.size
mongo.bson.size <- function(b)
    .Call(".mongo.bson.size", b)



#' Destroy a mongo.bson object
#' 
#' Releases the resources associated with a \link{mongo.bson} object. It is not
#' absolutely necessary to call this function since R's garbage collection will
#' eventually get around to doing it for you.
#' 
#' 
#' @param b A (\link{mongo.bson}) object.
#' @return NULL
#' @seealso \link{mongo.bson},\cr \code{\link{mongo.bson.from.list}},\cr
#' \code{\link{mongo.bson.from.buffer}}.
#' @examples
#' 
#' b <- mongo.bson.from.list(list(name="Cheryl", age=29))
#' print(b)
#' mongo.bson.destroy(b)
#' 
#' @export mongo.bson.destroy
mongo.bson.destroy <- function(b)
    .Call(".mongo.bson.destroy", b)



#' Display a mongo.bson object
#' 
#' Display formatted output of a mongo.bson object.
#' 
#' Output is tabbed (indented to show the nesting level of subobjects and
#' arrays).
#' 
#' 
#' @param x (\link{mongo.bson}) The mongo.bson object to display.
#' @param ... Parameters passed from generic.
#' @return The parameter is returned unchanged.
#' @seealso \link{mongo.bson}
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append(buf, "name", "Fred")
#' mongo.bson.buffer.append(buf, "city", "Dayton")
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # all display the same thing
#' mongo.bson.print(b)
#' print.mongo.bson(b)
#' print(b)
#' 
#' @export mongo.bson.print
mongo.bson.print <- function(x, ...) {
    if (typeof(x) == "list")
        x <- mongo.bson.from.list(x)
    invisible(.Call(".mongo.bson.print", x))
}



#' Display a mongo.bson object
#' 
#' Display formatted output of a mongo.bson object.
#' 
#' Output is tabbed (indented to show the nesting level of subobjects and
#' arrays).
#' 
#' This version is an alias of mongo.bson.print() so that print() will properly
#' handle the mongo.bson class.
#' 
#' 
#' @param x (\link{mongo.bson} The object to display.
#' @param ... Parameters passed from generic.
#' @return The parameter is returned unchanged.
#' @seealso \code{\link{mongo.bson.print}},\cr \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append(buf, "name", "Fred")
#' mongo.bson.buffer.append(buf, "city", "Dayton")
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # all display the same thing
#' print.mongo.bson(b)
#' mongo.bson.print(b)
#' print(b)
#' 
#' @export print.mongo.bson
#' @method print mongo.bson
print.mongo.bson <- function(x, ...)
    invisible(.Call(".mongo.bson.print", x))








#' Create a mongo.bson.iterator object
#' 
#' Create a \link{mongo.bson.iterator} object used to step through a given
#' \link{mongo.bson} object one field at a time.
#' 
#' 
#' @param b (\link{mongo.bson}) The mongo.bson object through which to iterate.
#' 
#' \code{b} may also be a mongo.bson.iterator and is expected to point to a
#' subobject or array. The iterator returned may be used to step through the
#' subobject or array.
#' @return (\link{mongo.bson.iterator}) An iterator initialized to 'before' the
#' start of the given mongo.bson object.
#' \code{\link{mongo.bson.iterator.next}()} should be used on the iterator
#' first to step to the first field.
#' @seealso \link{mongo.bson.iterator},\cr \code{\link{mongo.bson.find}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr
#' \code{\link{mongo.bson.iterator.key}},\cr
#' \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.value}},\cr \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' # Append a string
#' mongo.bson.buffer.append(buf, "name", "Joe")
#' # Append a date/time
#' mongo.bson.buffer.append(buf, "created", Sys.time())
#' # Append a NULL
#' mongo.bson.buffer.append(buf, "cars", NULL)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' iter <- mongo.bson.iterator.create(b)
#' while (mongo.bson.iterator.next(iter))
#'     if (mongo.bson.iterator.key(iter) == "created") {
#'         print(mongo.bson.iterator.value(iter))
#'         break
#'     }
#' 
#' # The above is given for illustrative purposes, but may be performed 
#' # much easier (and faster) by the following:
#' iter <- mongo.bson.find(b, "created")
#' print(mongo.bson.iterator.value(iter))
#' 
#' @export mongo.bson.iterator.create
mongo.bson.iterator.create <- function(b)
    .Call(".mongo.bson.iterator.create", b)



#' Find a field within a mongo.bson object by name
#' 
#' Find a field within a \link{mongo.bson} object by the name (key) of the
#' field\cr and return a \link{mongo.bson.iterator} pointing to that field.
#' 
#' The search parameter may also be a 'dotted' reference to a field in a
#' subobject or array. See examples.
#' 
#' 
#' @param b (\link{mongo.bson}) The object in which to find the field.
#' @param name (string) The name of the field to find.
#' @return (\link{mongo.bson.iterator}) An iterator pointing to the field found
#' if name was found among the names of the fields; otherwise, NULL.
#' @seealso \link{mongo.bson.iterator},\cr
#' \code{\link{mongo.bson.iterator.value}},\cr \link{mongo.bson}.
#' @examples
#' 
#' b <- mongo.bson.from.list(list(name="John", age=32L, 
#'      address=list(street="Vine", city="Denver", state="CO")))
#' iter <- mongo.bson.find(b, "age")
#' print(mongo.bson.iterator.value(iter)) # print 32
#' 
#' iter <- mongo.bson.find(b, "address.city")
#' print(mongo.bson.iterator.value(iter)) # print Denver
#' 
#' x <- c(1,1,2,3,5)
#' b <- mongo.bson.from.list(list(fib=x))
#' iter <- mongo.bson.find(b, "fib.3")  # BSON arrays are 0-based
#' print(mongo.bson.iterator.value(iter)) # print 3
#' 
#' @export mongo.bson.find
mongo.bson.find <- function(b, name)
    .Call(".mongo.bson.find", b, name)



#' Return the value of a mongo.bson field
#' 
#' Search a \link{mongo.bson} object for a field by name and retrieve its
#' value.
#' 
#' The search parameter may also be a 'dotted' reference to a field in a
#' subobject or array. See examples.
#' 
#' 
#' @param b A \link{mongo.bson} object.
#' @param name (string) The name of a field within \code{b}.
#' @return NULL if name is not found;\cr otherwise, the value of the found
#' field.
#' 
#' This function returns an appropriate R object depending on the type of the
#' field found. This mapping to values is as follows: \tabular{ll}{
#' \code{\link{mongo.bson.double}} \tab A double \cr
#' \code{\link{mongo.bson.string}} \tab A string \cr
#' \code{\link{mongo.bson.object}} \tab (See below).\cr
#' \code{\link{mongo.bson.array}} \tab (See below).\cr
#' \code{\link{mongo.bson.binary}} \tab A raw object. (See below).\cr
#' \code{\link{mongo.bson.undefined}} \tab A \link{mongo.undefined} object \cr
#' \code{\link{mongo.bson.oid}} \tab A \link{mongo.oid} object \cr
#' \code{\link{mongo.bson.bool}} \tab A logical \cr
#' \code{\link{mongo.bson.date}} \tab A "POSIXct" class object \cr
#' \code{\link{mongo.bson.null}} \tab NULL \cr \code{\link{mongo.bson.regex}}
#' \tab A \link{mongo.regex} object \cr \code{\link{mongo.bson.dbref}} \tab
#' Error! (deprecated -- see link) \cr \code{\link{mongo.bson.code}} \tab A
#' \link{mongo.code} object \cr \code{\link{mongo.bson.symbol}} \tab A
#' \link{mongo.symbol} object \cr \code{\link{mongo.bson.code.w.scope}} \tab A
#' \link{mongo.code.w.scope} object \cr \code{\link{mongo.bson.int}} \tab An
#' integer \cr \code{\link{mongo.bson.timestamp}} \tab A \link{mongo.timestamp}
#' object \cr \code{\link{mongo.bson.long}} \tab A double \cr }
#' 
#' Special handling:
#' 
#' \code{\link{mongo.bson.object}}: If the object is recognized as a complex
#' value (of the form { "r" : double, "i" : double }), a complex value is
#' returned. If the special wrapper as output by\cr
#' \code{\link{mongo.bson.buffer.append.object}()} is detected, an
#' appropriately attributed R object is returned; otherwise, a list is returned
#' containing the subfields.
#' 
#' \code{\link{mongo.bson.array}}: If all fields of the array are of the same
#' atomic type, a vector of that type is returned.  (Multidimensinal arrays are
#' detected and the \code{dims} attribute will be set accordingly.  Arrays of
#' complex values are also detected as above). Otherwise, a list is returned
#' containing the subfields.
#' 
#' \code{\link{mongo.bson.binary}}: If non-zero, the subtype of the binary data
#' is stored in the attribute "subtype". See
#' \code{\link{mongo.bson.buffer.append.raw}()}.
#' @seealso \link{mongo.bson.iterator.value},\cr \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' # Append a string
#' mongo.bson.buffer.append(buf, "name", "Joe")
#' # Append a date/time
#' mongo.bson.buffer.append(buf, "created", Sys.time())
#' # Append a NULL
#' mongo.bson.buffer.append(buf, "cars", NULL)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # Display the date appended above
#' print(mongo.bson.value(b, "created"))
#' 
#' b <- mongo.bson.from.list(list(name="John", age=32L, 
#'      address=list(street="Vine", city="Denver", state="CO")))
#' print(mongo.bson.value(b, "age")) # print 32
#' print(mongo.bson.value(b, "address.state")) # print CO
#' 
#' x <- c(1,1,2,3,5)
#' b <- mongo.bson.from.list(list(fib=x)) # BSON arrays are 0-based
#' print(mongo.bson.value(b, "fib.4")) # print 5
#' 
#' @export mongo.bson.value
mongo.bson.value <- function(b, name)
    .Call(".mongo.bson.value", b, name)




#' Advance an iterator to the first or next field
#' 
#' Advance a \link{mongo.bson.iterator} to the first or next field.
#' 
#' 
#' @param iter A \link{mongo.bson.iterator}.
#' @return (integer) The type of the next of the field pointed to by the
#' iterator as indicated by the folllowing constants: \itemize{
#' \item\link{mongo.bson.eoo} -- End of Object (0L)
#' \item\link{mongo.bson.double} \item\link{mongo.bson.string}
#' \item\link{mongo.bson.object} \item\link{mongo.bson.array}
#' \item\link{mongo.bson.binary} \item\link{mongo.bson.undefined}
#' \item\link{mongo.bson.oid} \item\link{mongo.bson.bool}
#' \item\link{mongo.bson.date} \item\link{mongo.bson.null}
#' \item\link{mongo.bson.regex} \item\link{mongo.bson.dbref} -- deprecated
#' (follow link for more info) \item\link{mongo.bson.code}
#' \item\link{mongo.bson.symbol} \item\link{mongo.bson.code.w.scope}
#' \item\link{mongo.bson.int} \item\link{mongo.bson.timestamp}
#' \item\link{mongo.bson.long} }
#' @seealso \link{mongo.bson.iterator},\cr
#' \code{\link{mongo.bson.iterator.create}},\cr
#' \code{\link{mongo.bson.find}},\cr \code{\link{mongo.bson.iterator.key}},\cr
#' \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.value}},\cr \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' # Append a string
#' mongo.bson.buffer.append(buf, "name", "Joe")
#' # Append a date/time
#' mongo.bson.buffer.append(buf, "created", Sys.time())
#' # Append a NULL
#' mongo.bson.buffer.append(buf, "cars", NULL)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' iter <- mongo.bson.iterator.create(b)
#' # Advance to the "cars" field
#' while (mongo.bson.iterator.next(iter) != mongo.bson.null)
#' { 
#'     # NOP
#' }
#' print(mongo.bson.iterator.value(iter))
#' 
#' # The above is given for illustrative purposes, but may be performed 
#' # much easier by the following:
#' iter <- mongo.bson.find(b, "cars")
#' print(mongo.bson.iterator.value(iter))
#' 
#' # iterate through all values and print them with their keys (names)
#' iter <- mongo.bson.iterator.create(b)
#' while (mongo.bson.iterator.next(iter)) { # eoo at end stops loop
#'     print(mongo.bson.iterator.key(iter))
#'     print(mongo.bson.iterator.value(iter))
#' }
#' 
#' @export mongo.bson.iterator.next
mongo.bson.iterator.next <- function(iter)
    .Call(".mongo.bson.iterator.next", iter)



#' Return the key (name) of the field pointed to by an iterator
#' 
#' Return the key (name) of the field pointed to by a
#' \link{mongo.bson.iterator}.
#' 
#' 
#' @param iter A \link{mongo.bson.iterator}.
#' @return (string) The key (name) of the field pointed to by \code{iter}
#' @seealso \link{mongo.bson.iterator},\cr
#' \code{\link{mongo.bson.iterator.create}},\cr
#' \code{\link{mongo.bson.find}},\cr \code{\link{mongo.bson.iterator.next}},\cr
#' \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.value}},\cr \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' # Append a string
#' mongo.bson.buffer.append(buf, "name", "Joe")
#' # Append a date/time
#' mongo.bson.buffer.append(buf, "created", Sys.time())
#' # Append a NULL
#' mongo.bson.buffer.append(buf, "cars", NULL)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # iterate through all values and print them with their keys (names)
#' iter <- mongo.bson.iterator.create(b)
#' while (mongo.bson.iterator.next(iter)) { # eoo at end stops loop
#'     print(mongo.bson.iterator.key(iter))
#'     print(mongo.bson.iterator.value(iter))
#' }
#' 
#' @export mongo.bson.iterator.key
mongo.bson.iterator.key <- function(iter)
    .Call(".mongo.bson.iterator.key", iter)



#' Get the type of data pointed to by an iterator
#' 
#' Return the type of the field currently pointed to by a
#' \link{mongo.bson.iterator}.
#' 
#' 
#' @param iter A \link{mongo.bson.iterator}.
#' @return (integer) The type of the field pointed to by the iterator as
#' indicated by the following constants: \itemize{ \item\link{mongo.bson.eoo}
#' -- End of Object (0L) \item\link{mongo.bson.double}
#' \item\link{mongo.bson.string} \item\link{mongo.bson.object}
#' \item\link{mongo.bson.array} \item\link{mongo.bson.binary}
#' \item\link{mongo.bson.undefined} \item\link{mongo.bson.oid}
#' \item\link{mongo.bson.bool} \item\link{mongo.bson.date}
#' \item\link{mongo.bson.null} \item\link{mongo.bson.regex}
#' \item\link{mongo.bson.dbref} -- deprecated (follow link for more info)
#' \item\link{mongo.bson.code} \item\link{mongo.bson.symbol}
#' \item\link{mongo.bson.code.w.scope} \item\link{mongo.bson.int}
#' \item\link{mongo.bson.timestamp} \item\link{mongo.bson.long} }
#' @seealso \link{mongo.bson.iterator},\cr
#' \code{\link{mongo.bson.iterator.create}},\cr
#' \code{\link{mongo.bson.find}},\cr \code{\link{mongo.bson.iterator.next}},\cr
#' \code{\link{mongo.bson.iterator.key}},\cr
#' \code{\link{mongo.bson.iterator.value}},\cr \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' # Append a string
#' mongo.bson.buffer.append(buf, "name", "Joe")
#' # Append a date/time
#' mongo.bson.buffer.append(buf, "created", Sys.time())
#' # Append a NULL
#' mongo.bson.buffer.append(buf, "cars", NULL)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' iter <- mongo.bson.iterator.create(b)
#' while (mongo.bson.iterator.next(iter))
#'     if (mongo.bson.iterator.type(iter) == mongo.bson.date) {
#'         print(mongo.bson.iterator.value(iter))
#'         break
#'     }
#' 
#' # The above is given for illustrative purposes, but may be performed 
#' # much easier by the following:
#' iter <- mongo.bson.find(b, "created")
#' print(mongo.bson.iterator.value(iter))
#' 
#' @export mongo.bson.iterator.type
mongo.bson.iterator.type <- function(iter)
    .Call(".mongo.bson.iterator.type", iter)



#' BSON data type constant for 'End Of Object'
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (0L) at
#' the end of the object when there are no more fields through which to
#' iterate.
#' 
#' 
#' @return 0L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson}.
#' @export mongo.bson.eoo
mongo.bson.eoo       <- 0L # end of object


#' BSON data type constant for a double value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (1L) to
#' indicate that the value pointer to by an iterator is a double.
#' 
#' 
#' @return 1L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson}.
#' @export mongo.bson.double
mongo.bson.double    <- 1L


#' BSON data type constant for a string value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (2L) to
#' indicate that the value pointer to by an iterator is a string.
#' 
#' 
#' @return 2L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson}.
#' @export mongo.bson.string
mongo.bson.string    <- 2L


#' BSON data type constant for a subobject value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (3L) to
#' indicate that the value pointer to by an iterator is a subobject.
#' 
#' 
#' @return 3L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson}.
#' @export mongo.bson.object
mongo.bson.object    <- 3L


#' BSON data type constant for an array
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (4L) to
#' indicate that the value pointer to by an iterator is an array (containing
#' child values).
#' 
#' 
#' @return 4L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson.iterator},\cr
#' \link{mongo.bson}.
#' @export mongo.bson.array
mongo.bson.array     <- 4L


#' BSON data type constant for a binary data value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (5L) to
#' indicate that the value pointer to by an iterator is binary data.
#' 
#' 
#' @return 5L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson.iterator},\cr
#' \link{mongo.bson}.
#' @export mongo.bson.binary
mongo.bson.binary    <- 5L


#' BSON data type constant for a undefined value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (6L) to
#' indicate that the value pointer to by an iterator is a undefined.
#' 
#' 
#' @return 6L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson}.
#' @export mongo.bson.undefined
mongo.bson.undefined <- 6L


#' BSON data type constant for a oid value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (7L) to
#' indicate that the value pointer to by an iterator is a oid (Object ID).
#' 
#' 
#' @return 7L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson}.
#' @export mongo.bson.oid
mongo.bson.oid       <- 7L


#' BSON data type constant for a bool value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (8L) to
#' indicate that the value pointer to by an iterator is a bool.
#' 
#' 
#' @return 8L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson.iterator},\cr
#' \link{mongo.bson}.
#' @export mongo.bson.bool
mongo.bson.bool      <- 8L


#' BSON data type constant for a date value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (9L) to
#' indicate that the value pointer to by an iterator is a date/time.
#' 
#' 
#' @return 9L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson}.
#' @export mongo.bson.date
mongo.bson.date      <- 9L


#' BSON data type constant for a null value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (10L) to
#' indicate that the value pointer to by an iterator is a null.
#' 
#' 
#' @return 10L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson}.
#' @export mongo.bson.null
mongo.bson.null      <- 10L


#' BSON data type constant for a regex value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (11L) to
#' indicate that the value pointer to by an iterator is a regular expression.
#' 
#' 
#' @return 11L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson}.
#' @export mongo.bson.regex
mongo.bson.regex     <- 11L


#' BSON data type constant for a dbref value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (12L) to
#' indicate that the value pointed to by an iterator is a dbref (database
#' reference).
#' 
#' Note that this BSON data type is deprecated and rmongodb provides no support
#' for it. Attempting to fetch the value of a dbref with
#' \code{\link{mongo.bson.to.list}()} or\cr
#' \code{\link{mongo.bson.iterator.value}()} will throw an error. The field
#' must be skipped by calling \code{\link{mongo.bson.iterator.next}()}.
#' 
#' 
#' @return 12L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson}.
#' @export mongo.bson.dbref
mongo.bson.dbref     <- 12L # deprecated


#' BSON data type constant for a code value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (13L) to
#' indicate that the value pointer to by an iterator is javascript code.
#' 
#' 
#' @return 13L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson}.
#' @export mongo.bson.code
mongo.bson.code      <- 13L


#' BSON data type constant for a symbol value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (14L) to
#' indicate that the value pointer to by an iterator is a symbol.
#' 
#' 
#' @return 14L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson}.
#' @export mongo.bson.symbol
mongo.bson.symbol    <- 14L


#' BSON data type constant for a code with scope value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (15L) to
#' indicate that the value pointer to by an iterator is a javascript with a
#' scope.
#' 
#' 
#' @return 15L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson}
#' @export mongo.bson.code.w.scope
mongo.bson.code.w.scope <- 15L


#' BSON data type constant for a integer value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (16L) to
#' indicate that the value pointer to by an iterator is a integer (32-bit).
#' 
#' 
#' @return 16L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson}.
#' @export mongo.bson.int
mongo.bson.int       <- 16L


#' BSON data type constant for a timestamp value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (17L) to
#' indicate that the value pointer to by an iterator is a timestamp.
#' 
#' 
#' @return 17L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson}.
#' @export mongo.bson.timestamp
mongo.bson.timestamp <- 17L


#' BSON data type constant for a long value
#' 
#' \code{\link{mongo.bson.iterator.type}()} and
#' \code{\link{mongo.bson.iterator.next}()} will return this constant (18L) to
#' indicate that the value pointer to by an iterator is a long integer (64
#' bits).
#' 
#' 
#' @return 18L
#' @seealso \code{\link{mongo.bson.iterator.type}},\cr
#' \code{\link{mongo.bson.iterator.next}},\cr \link{mongo.bson}.
#' @export mongo.bson.long
mongo.bson.long      <- 18L



#' BSON binary data subtype constant for standard binary data
#' 
#' BSON binary data subtype constant for standard binary data.
#' 
#' 
#' @return 0L
#' @seealso \code{\link{mongo.bson.buffer.append.raw}},\cr \link{mongo.bson}.
#' @export mongo.binary.binary
mongo.binary.binary  <- 0L 


#' BSON binary data subtype constant for function data
#' 
#' BSON binary data subtype constant for function data.
#' 
#' 
#' @return 1L
#' @seealso \code{\link{mongo.bson.buffer.append.raw}},\cr \link{mongo.bson}.
#' @export mongo.binary.function
mongo.binary.function <- 1L 


#' BSON binary data subtype constant for old format data
#' 
#' BSON binary data subtype constant for old format data (deprecated).
#' 
#' 
#' @return 2L
#' @seealso \code{\link{mongo.bson.buffer.append.raw}},\cr \link{mongo.bson}.
#' @export mongo.binary.old
mongo.binary.old     <- 2L 


#' BSON binary data subtype constant for uuid data
#' 
#' BSON binary data subtype constant for uuid data.
#' 
#' 
#' @return 4L
#' @seealso \code{\link{mongo.bson.buffer.append.raw}},\cr \link{mongo.bson}.
#' @export mongo.binary.uuid
mongo.binary.uuid    <- 3L 


#' BSON binary data subtype constant for md5 data
#' 
#' BSON binary data subtype constant for md5 data.
#' 
#' 
#' @return 5L
#' @seealso \code{\link{mongo.bson.buffer.append.raw}},\cr \link{mongo.bson}.
#' @export mongo.binary.md5
mongo.binary.md5     <- 5L


#' BSON binary data subtype constant for user data
#' 
#' BSON binary data subtype constant for user data.
#' 
#' 
#' @return 128L
#' @seealso \code{\link{mongo.bson.buffer.append.raw}},\cr \link{mongo.bson}.
#' @export mongo.binary.user
mongo.binary.user    <- 128L



#' Return the value of the field pointed to by an iterator
#' 
#' Return the value of the field pointed to by a \link{mongo.bson.iterator}.
#' 
#' 
#' @param iter A \link{mongo.bson.iterator}.
#' @return The value of the field pointed to by \code{iter}.
#' 
#' This function returns an appropriate R object depending on the type of the
#' field pointed to by the iterator. This mapping to values is as follows:
#' \tabular{ll}{ \code{\link{mongo.bson.eoo}} \tab 0L \cr
#' \code{\link{mongo.bson.double}} \tab A double \cr
#' \code{\link{mongo.bson.string}} \tab A string \cr
#' \code{\link{mongo.bson.object}} \tab (See below).\cr
#' \code{\link{mongo.bson.array}} \tab (See below).\cr
#' \code{\link{mongo.bson.binary}} \tab A raw vector.  (See below).\cr
#' \code{\link{mongo.bson.undefined}} \tab A \link{mongo.undefined} object \cr
#' \code{\link{mongo.bson.oid}} \tab A \link{mongo.oid} object \cr
#' \code{\link{mongo.bson.bool}} \tab A logical \cr
#' \code{\link{mongo.bson.date}} \tab A "POSIXct" class object \cr
#' \code{\link{mongo.bson.null}} \tab NULL \cr \code{\link{mongo.bson.regex}}
#' \tab A \link{mongo.regex} object \cr \code{\link{mongo.bson.dbref}} \tab
#' Error! (deprecated -- see link) \cr \code{\link{mongo.bson.code}} \tab A
#' \link{mongo.code} object \cr \code{\link{mongo.bson.symbol}} \tab A
#' \link{mongo.symbol} object \cr \code{\link{mongo.bson.code.w.scope}} \tab A
#' \link{mongo.code.w.scope} object \cr \code{\link{mongo.bson.int}} \tab An
#' integer \cr \code{\link{mongo.bson.timestamp}} \tab A \link{mongo.timestamp}
#' object \cr \code{\link{mongo.bson.long}} \tab A double \cr }
#' 
#' Special handling:
#' 
#' \code{\link{mongo.bson.object}}: If the object is recognized as a complex
#' value (of the form { "r" : double, "i" : double }), a complex value is
#' returned. If the special wrapper as output by\cr
#' \code{\link{mongo.bson.buffer.append.object}()} is detected, an
#' appropriately attributed R object is returned; otherwise, a list is returned
#' containing the subfields.
#' 
#' \code{\link{mongo.bson.array}}: If all fields of the array are of the same
#' atomic type, a vector of that type is returned.  (Multidimensinal arrays are
#' detected and the \code{dims} attribute will be set accordingly.  Arrays of
#' complex values are also detected as above). Otherwise, a list is returned
#' containing the subfields.
#' 
#' \code{\link{mongo.bson.binary}}: If non-zero, the subtype of the binary data
#' is stored in the attribute "subtype". See
#' \code{\link{mongo.bson.buffer.append.raw}()}.
#' @seealso \link{mongo.bson.iterator},\cr
#' \code{\link{mongo.bson.iterator.create}},\cr
#' \code{\link{mongo.bson.find}},\cr \code{\link{mongo.bson.iterator.next}},\cr
#' \code{\link{mongo.bson.iterator.key}},\cr
#' \code{\link{mongo.bson.iterator.type}},\cr \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' # Append a string
#' mongo.bson.buffer.append(buf, "name", "Joe")
#' # Append a date/time
#' mongo.bson.buffer.append(buf, "created", Sys.time())
#' # Append a NULL
#' mongo.bson.buffer.append(buf, "cars", NULL)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # iterate through all values and print them with their keys (names)
#' iter <- mongo.bson.iterator.create(b)
#' while (mongo.bson.iterator.next(iter)) { # eoo at end stops loop
#'     print(mongo.bson.iterator.key(iter))
#'     print(mongo.bson.iterator.value(iter))
#' }
#' 
#' @export mongo.bson.iterator.value
mongo.bson.iterator.value <- function(iter)
    .Call(".mongo.bson.iterator.value", iter)



#' Create a mongo.oid object ftom a string
#' 
#' Create from a 24-character hex string a mongo.oid object representing a
#' MongoDB Object ID.
#' 
#' See \url{http://www.mongodb.org/display/DOCS/Object+IDs}
#' 
#' 
#' @param hexstr (string) 24 hex characters representing the OID.
#' 
#' Note that although an error is thrown if the length is not 24, no error is
#' thrown if the characters are not hex digits; you'll get zero bits for the
#' invalid digits.
#' @return A \link{mongo.oid} object constructed from hexstr.
#' @seealso \link{mongo.oid},\cr \code{\link{mongo.oid.create}},\cr
#' \code{\link{as.character.mongo.oid}},\cr
#' \code{\link{mongo.oid.to.string}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.oid}},\cr
#' \code{\link{mongo.bson.buffer.append.list}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' oid <- mongo.oid.from.string("ABCD1234EFAB5678CDEF9012")
#' mongo.bson.buffer.append(buf, "_id", oid)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' @export mongo.oid.from.string
mongo.oid.from.string <- function(hexstr)
    .Call(".mongo.oid.from.string", hexstr)



#' Convert a mongo.oid object to a string
#' 
#' Convert a \link{mongo.oid} object to a string of 24 hex digits. This
#' performs the inverse operation of \code{\link{mongo.oid.from.string}()}.
#' 
#' This function is an alias of \code{\link{as.character.mongo.oid}()} which
#' you may perfer to use since the class mechanism of R allows that to be
#' called simply by \code{as.character(oid)}.
#' 
#' See \url{http://www.mongodb.org/display/DOCS/Object+IDs}
#' 
#' 
#' @param oid (\link{mongo.oid}) The OID to be converted.
#' @return (string) A string of 24 hex digits representing the bits of
#' \code{oid}.
#' @seealso \link{mongo.oid},\cr \code{\link{mongo.oid.create}},\cr
#' \code{\link{as.character.mongo.oid}},\cr
#' \code{\link{mongo.oid.from.string}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.oid}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' oid <- mongo.oid.create()
#' print(mongo.oid.to.string(oid))
#' print(as.character(oid))  # print same thing as above line
#' 
#' @export mongo.oid.to.string
mongo.oid.to.string <- function(oid)
    .Call(".mongo.oid.to.string", oid)



#' Convert a mongo.oid object to a string
#' 
#' Convert a \link{mongo.oid} object to a string of 24 hex digits. This
#' performs the inverse operation of \code{\link{mongo.oid.from.string}()}.
#' 
#' This function is an alias of \code{\link{mongo.oid.to.string}()} so that the
#' class mechanism of R allows it to be called simply by
#' \code{as.character(oid)}.
#' 
#' See \url{http://www.mongodb.org/display/DOCS/Object+IDs}
#' 
#' 
#' @param x (\link{mongo.oid}) The OID to be converted.
#' @param ... Parameters passed from generic.
#' @return (string) A string of 24 hex digits representing the bits of oid
#' \code{x}.
#' @seealso \link{mongo.oid},\cr \code{\link{mongo.oid.create}},\cr
#' \code{\link{as.character.mongo.oid}},\cr
#' \code{\link{mongo.oid.to.string}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.oid}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' oid <- mongo.oid.create()
#' print(as.character.mongo.oid(oid))
#' print(as.character(oid))  # print same thing as above line
#' 
#' @export as.character.mongo.oid
#' @method as character.mongo.oid
as.character.mongo.oid <- function(x, ...)
    .Call(".mongo.oid.to.string", x)



#' Create a mongo.oid object
#' 
#' Create a \link{mongo.oid} object for appending to a buffer with\cr
#' \code{\link{mongo.bson.buffer.append.oid}()} or
#' \code{\link{mongo.bson.buffer.append}()}, or for embedding in a list such
#' that \code{\link{mongo.bson.buffer.append.list}()} will properly insert an
#' Object ID value into a mongo.bson.buffer object.
#' 
#' See \url{http://www.mongodb.org/display/DOCS/Object+IDs}
#' 
#' 
#' @return A \link{mongo.oid} object that is reasonably assured of being
#' unique.
#' @seealso \link{mongo.oid},\cr \code{\link{mongo.oid.from.string}},\cr
#' \code{\link{as.character.mongo.oid}},\cr
#' \code{\link{mongo.oid.to.string}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.oid}},\cr
#' \code{\link{mongo.bson.buffer.append.list}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' oid <- mongo.oid.create()
#' mongo.bson.buffer.append(buf, "_id", oid)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' @export mongo.oid.create
mongo.oid.create <- function()
    .Call(".mongo.oid.create")



#' Get an Object ID's time
#' 
#' Get the 32-bit UTC time portion of an OID (Object ID).
#' 
#' See \url{http://www.mongodb.org/display/DOCS/Object+IDs}
#' 
#' 
#' @param oid (\link{mongo.oid}) The OID to be examined.
#' @return (integer) ("POSIXct") The time portion of the given \code{oid}.
#' @seealso \link{mongo.oid},\cr \code{\link{mongo.oid.create}},\cr
#' \code{\link{as.character.mongo.oid}},\cr
#' \code{\link{mongo.oid.to.string}},\cr
#' \code{\link{mongo.oid.from.string}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.oid}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' oid <- mongo.oid.create()
#' print(mongo.oid.time(oid))
#' 
#' @export mongo.oid.time
mongo.oid.time <- function(oid)
    .Call(".mongo.oid.time", oid)



#' Display a mongo.oid object
#' 
#' Display formatted output of a \link{mongo.oid} object.
#' 
#' This version is an alias of \code{\link{print.mongo.oid}()} which allows
#' \code{print()} to properly handle the mongo.oid class.
#' 
#' 
#' @param x \link{mongo.oid} The object to display.
#' @return The parameter is returned unchanged.
#' @seealso \code{\link{mongo.oid.print}},\cr
#' \code{\link{mongo.oid.to.string}},\cr \link{mongo.bson.oid},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' oid <- mongo.oid.create()
#' 
#' # all display the same thing
#' print.mongo.oid(oid)
#' mongo.oid.print(oid)
#' print(oid)
#' 
#' @export mongo.oid.print
mongo.oid.print <- function(x)
    invisible(.Call(".mongo.oid.print", x))



#' Display a mongo.oid object
#' 
#' Display formatted output of a \link{mongo.oid} object.
#' 
#' Output is tabbed (indented to show the nesting level of subobjects and
#' arrays).
#' 
#' This version is an alias of \code{\link{mongo.oid.print}()} so that
#' \code{print()} will properly handle the mongo.oid class.
#' 
#' 
#' @param x \link{mongo.oid} The object to display.
#' @param ... Parameters passed from generic.
#' @return The parameter is returned unchanged.
#' @seealso \code{\link{mongo.oid.print}},\cr
#' \code{\link{mongo.oid.to.string}},\cr \link{mongo.bson.oid},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' oid <- mongo.oid.create()
#' 
#' # all display the same thing
#' print.mongo.oid(oid)
#' mongo.oid.print(oid)
#' print(oid)
#' 
#' @export print.mongo.oid
#' @method print mongo.oid
print.mongo.oid <- function(x, ...)
    invisible(.Call(".mongo.oid.print", x))




#' Create a mongo.timestamp object
#' 
#' Create a \link{mongo.timestamp} object for appending to a buffer with\cr
#' \code{\link{mongo.bson.buffer.append.timestamp}()} or
#' \code{\link{mongo.bson.buffer.append}()}, or for embedding in a list such
#' that \code{\link{mongo.bson.buffer.append.list}()} will properly insert a
#' timestamp value into the mongo.bson.buffer object.
#' 
#' See \url{http://www.mongodb.org/display/DOCS/Timestamp+Data+Type}
#' 
#' 
#' @param time (integer) date/time value (milliseconds since UTC epoch).
#' 
#' This may also be a "POSIXct" or "POSIXlt" class object.
#' @param increment increment ordinal
#' @return A \link{mongo.timestamp} object
#' @seealso \link{mongo.timestamp},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.time}},\cr
#' \code{\link{mongo.bson.buffer.append.list}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' mongo <- mongo.create()
#' if (mongo.is.connected(mongo)) {
#'     buf <- mongo.bson.buffer.create()
#'     # special Null timestamp -- automatically filled in
#'     # if one of first two fields in a record
#'     ts <- mongo.timestamp.create(0,0)
#'     mongo.bson.buffer.append(buf, "InsertTime", ts)
#'     mongo.bson.buffer.append(buf, "name", "Joe")
#'     b <- mongo.bson.from.buffer(buf)
#'     mongo.insert(mongo, "test.people", b)
#' 
#'     # create using a POSIXlt
#'     ts <- mongo.timestamp.create(strptime("05-12-2012",
#'         "%m-%d-%Y"), increment=1)
#' }
#' 
#' @export mongo.timestamp.create
mongo.timestamp.create <- function(time, increment)
    .Call(".mongo.timestamp.create", time, increment)



#' Create a mongo.code object
#' 
#' Create a mongo.code object for appending to a buffer with
#' \code{\link{mongo.bson.buffer.append}()} or for embedding in a list such
#' that\cr \code{\link{mongo.bson.buffer.append.list}()} will properly insert a
#' code value into the mongo.bson.buffer object.
#' 
#' 
#' @param code (string) javascript code
#' @return A \link{mongo.code} object
#' @seealso \link{mongo.code},\cr \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.list}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' code <- mongo.code.create("y = x")
#' mongo.bson.buffer.append(buf, "Code", code)
#' lst <- list(c1 = code, One = 1)
#' mongo.bson.buffer.append.list(buf, "listWcode", lst)
#' mongo.bson.buffer.append.code(buf, "Code2", "a = 1")
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # the above will create a mongo.bson object of the following form:
#' # { "Code": (CODE) "y = x",
#' #   "listWcode" : { "c1"  : (CODE) "y = x", 
#' #                   "One" : 1 },
#' #   "Code2" : (CODE) "a = 1" }
#' 
#' @export mongo.code.create
mongo.code.create <- function(code)
    .Call(".mongo.code.create", code)



#' Create a mongo.code.w.scope object
#' 
#' Create a mongo.code.w.scope object for appending to a buffer with\cr
#' \code{\link{mongo.bson.buffer.append}()} or for embedding in a list such
#' that \code{\link{mongo.bson.buffer.append.list}()} will properly insert a
#' code value into the mongo.bson.buffer object.
#' 
#' 
#' @param code (string) javascript code
#' @param scope (\link{mongo.bson}) the scope object
#' @return A \link{mongo.code.w.scope} object
#' @seealso \link{mongo.code.w.scope},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.list}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append(buf, "scopevar", "scopevalue")
#' scope <- mongo.bson.from.buffer(buf)
#' codeWscope <- mongo.code.w.scope.create("y = x", scope)
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append(buf, "CodeWscope", codeWscope)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form:
#' # { "CodeWscope" : (CODEWSCOPE) "y = x"
#' #                  (SCOPE) { "scopevar" : "scopevalue" } }
#' 
#' @export mongo.code.w.scope.create
mongo.code.w.scope.create <- function(code, scope)
    .Call(".mongo.code.w.scope.create", code, scope)



#' Create a mongo.symbol object
#' 
#' Create a mongo.symbol object for appending to a buffer with\cr
#' \code{\link{mongo.bson.buffer.append}()} or for embedding in a list such
#' that \code{\link{mongo.bson.buffer.append.list}()} will properly insert a
#' symbol value into the mongo.bson.buffer object.
#' 
#' 
#' @param value (string) The value of the symbol
#' @return a \link{mongo.symbol} object
#' @seealso \link{mongo.symbol},\cr \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.list}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' sym <- mongo.symbol.create("Alpha")
#' mongo.bson.buffer.append(buf, "A", sym)
#' lst <- list(s1 = sym, One = 1)
#' mongo.bson.buffer.append.list(buf, "listWsym", lst)
#' mongo.bson.buffer.append.symbol(buf, "D", "Delta")
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # the above will create a mongo.bson object of the following form:
#' # { "A": (SYMBOL) "Alpha",
#' #   "listWsym" : { "a1"  : (SYMBOL) "Aplha",
#' #                  "One" : 1 }, 
#' #   "D" : (SYMBOL) "Delta" }
#' 
#' @export mongo.symbol.create
mongo.symbol.create <- function(value)
    .Call(".mongo.symbol.create", value)



#' Create a mongo.undefined object
#' 
#' Create a mongo.undefined object for appending to a buffer with\cr
#' \code{\link{mongo.bson.buffer.append}()} or for embedding in a list such
#' that \code{\link{mongo.bson.buffer.append.list}()} will properly insert an
#' undefined value into the mongo.bson.buffer object.
#' 
#' 
#' @return a \link{mongo.undefined} object
#' @seealso \link{mongo.undefined},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.list}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' undef <- mongo.undefined.create()
#' mongo.bson.buffer.append(buf, "Undef", undef)
#' l <- list(u1 = undef, One = 1)
#' mongo.bson.buffer.append.list(buf, "listWundef", l)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # the above will create a mongo.bson object of the following form:
#' # { "Undef": UNDEFINED, "listWundef" : { "u1" : UNDEFINED, "One" : 1 } }
#' 
#' @export mongo.undefined.create
mongo.undefined.create <- function()
    .Call(".mongo.undefined.create")



#' Create a mongo.regex object
#' 
#' Create a \link{mongo.regex} object for appending to a buffer with\cr
#' \code{\link{mongo.bson.buffer.append.regex}()} or
#' \code{\link{mongo.bson.buffer.append}()}, or for embedding in a list such
#' that \code{\link{mongo.bson.buffer.append.list}()} will properly insert a
#' regular expression value into a mongo.bson.buffer object.
#' 
#' See
#' \url{http://www.mongodb.org/display/DOCS/Advanced+Queries#AdvancedQueries-RegularExpressions}
#' 
#' 
#' @param pattern (string) The regular expression.
#' @param options (string) Options governing the parsing done with the pattern.
#' @return A \link{mongo.regex} object
#' @seealso \link{mongo.regex},\cr \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.buffer.append.regex}},\cr
#' \code{\link{mongo.bson.buffer.append.list}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' regex <- mongo.regex.create("acme.*corp", options="i")
#' mongo.bson.buffer.append.regex(buf, "MatchAcme", regex)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' @export mongo.regex.create
mongo.regex.create <- function(pattern, options="")
    .Call(".mongo.regex.create", pattern, options)




#' Create an new mongo.bson.buffer object
#' 
#' Returns a fresh mongo.bson.buffer object ready to have data appended onto
#' it.
#' 
#' mongo.bson.buffer objects are used to build mongo.bson objects.
#' 
#' 
#' @return A fresh \link{mongo.bson.buffer} object
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append(buf, "name", "Donna")
#' b <- mongo.bson.from.buffer(buf)
#' 
#' @export mongo.bson.buffer.create
mongo.bson.buffer.create <- function()
    .Call(".mongo.bson.buffer.create")



#' Convert a mongo.bson.buffer object to a mongo.bson object
#' 
#' Convert a \link{mongo.bson.buffer} object to a \link{mongo.bson} object.
#' 
#' Use this after appending data to a buffer to turn it into a mongo.bson
#' object for network transport.
#' 
#' No futher data may be appended to the buffer after calling this function.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer to convert.
#' @return A \link{mongo.bson} object as converted from the buffer parameter.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr \code{\link{mongo.bson.destroy}}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append(buf, "name", "Fred")
#' mongo.bson.buffer.append(buf, "city", "Dayton")
#' b <- mongo.bson.from.buffer(buf)
#' print(b)
#' mongo.bson.destroy(b)
#' 
#' @export mongo.bson.from.buffer
mongo.bson.from.buffer <- function(buf)
    .Call(".mongo.bson.from.buffer", buf)



#' Append an integer field onto a mongo.bson.buffer
#' 
#' Append an integer or vector of integers onto a \link{mongo.bson.buffer}.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param value (integer vector) The integer(s) to append to the buffer.
#' 
#' If value has a \code{dims} attribute of length > 1, any \code{names} or
#' \code{dimnames} attribute is ignored and a nested array is appended.\cr (Use
#' \code{\link{mongo.bson.buffer.append.object}()} if you want to preserve
#' \code{dimnames}).
#' 
#' If value has a names attribute, a subobject is appended and the subfields
#' are given the indicated names.
#' 
#' Otherwise, if more than one element is present in value it must be a vector
#' of integers and the integers are appended as a subarray.
#' 
#' In the last case, the single value must be coerible to an integer.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \code{\link{mongo.bson.buffer.append}}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.int(buf, "age", 23L)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # the above produces a BSON object of the form { "age" : 21 }
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.int(buf, "ages", c(21L, 19L, 13L))
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # the above produces a BSON object of the form { "ages" : [21, 19, 13] }
#' 
#' buf <- mongo.bson.buffer.create()
#' dim <- c(2L, 4L, 8L)
#' names(dim) <- c("width", "height", "length")
#' mongo.bson.buffer.append.int(buf, "board", dim)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # theabove produces a BSON object of the form:
#' # { "board" : { "width" : 2, "height" : 4, "length" : 8 } }
#' 
#' @export mongo.bson.buffer.append.int
mongo.bson.buffer.append.int <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append.int", buf, name, value)



#' Append a boolean field onto a mongo.bson.buffer
#' 
#' Append an logical (boolean) or vector of logical values onto a
#' \link{mongo.bson.buffer}.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param value (logical vector) the booleans(s) to append to the buffer.
#' 
#' If value has a \code{dims} attribute of length > 1, any \code{names} or
#' \code{dimnames} attribute is ignored and a nested array is appended.\cr (Use
#' \code{\link{mongo.bson.buffer.append.object}()} if you want to preserve
#' \code{dimnames}).
#' 
#' If value has a names attribute, a subobject is appended and the subfields
#' are given the indicated names.
#' 
#' Otherwise, if more than one element is present in value, the booleans are
#' appended as a subarray.
#' 
#' In the last case, a single as.boolean is appended as the value of the field.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \code{\link{mongo.bson.buffer.append}}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.bool(buf, "wise", TRUE)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form { "wise" : true }
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.bool(buf, "bools", c(TRUE, FALSE, FALSE))
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form:
#' # { "bools" : [true, false, false] }
#' 
#' buf <- mongo.bson.buffer.create()
#' flags <- c(FALSE, FALSE, TRUE)
#' names(flags) <- c("Tall", "Fat", "Pretty")
#' mongo.bson.buffer.append.bool(buf, "Looks", flags)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form:
#' # { "Looks" : { "Tall" : false, "Fat" : false, "Pretty" : true } }
#' 
#' @export mongo.bson.buffer.append.bool
mongo.bson.buffer.append.bool <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append.bool", buf, name, value)



#' Append a double field onto a mongo.bson.buffer
#' 
#' Append a double or vector of doubles onto a \link{mongo.bson.buffer}.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param value (double vector) The values(s) to append to the buffer.
#' 
#' If value has a \code{dims} attribute of length > 1, any \code{names} or
#' \code{dimnames} attribute is ignored and a nested array is appended.\cr (Use
#' \code{\link{mongo.bson.buffer.append.object}()} if you want to preserve
#' \code{dimnames}).
#' 
#' If value has a names attribute, a subobject is appended and the subfields
#' are given the indicated names.
#' 
#' Otherwise, if more than one element is present in value, the values are
#' appended as a subarray.
#' 
#' In the last case, a single as.double is appended as the value of the field.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \code{\link{mongo.bson.buffer.append}}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.double(buf, "YearSeconds",
#'      365.24219 * 24 * 60 * 60)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form:
#' # { "YearSeconds" : 31556925.2 }
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.double(buf, "dbls", 
#'     c(1.7, 87654321.123, 12345678.321))
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form:
#' # { "dbls" : [1.7, 87654321.123, 12345678.321] }
#' 
#' buf <- mongo.bson.buffer.create()
#' fractions <- c(0.5, 0.25, 0.333333)
#' names(fractions) <- c("Half", "Quarter", "Third")
#' mongo.bson.buffer.append.double(buf, "Fractions", fractions)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form:
#' # { "Fractions" : { "Half"    : 0.5,
#' #                   "Quarter" : 0.25,
#' #                   "Third"   : 0.333333 } }
#' 
#' @export mongo.bson.buffer.append.double
mongo.bson.buffer.append.double <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append.double", buf, name, value)



#' Append a double field onto a mongo.bson.buffer
#' 
#' Append a double or vector of doubles onto a \link{mongo.bson.buffer}.
#' 
#' Note that since BSON has no built-in complex type, R's complex values are
#' appended as subobjects with two fields: "r" : the real part and "i" : the
#' imaginary part.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param value (complex vector) The values(s) to append to the buffer.
#' 
#' If value has a \code{dims} attribute of length > 1, any \code{names} or
#' \code{dimnames} attribute is ignored and a nested array is appended.\cr (Use
#' \code{\link{mongo.bson.buffer.append.object}()} if you want to preserve
#' \code{dimnames}).
#' 
#' If value has a names attribute, a subobject is appended and the subfields
#' are given the indicated names.
#' 
#' Otherwise, if more than one element is present in value, the values are
#' appended as a subarray.
#' 
#' In the last case, a single complex is appended as the value of the field.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \code{\link{mongo.bson.buffer.append}}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.complex(buf, "Alpha", 3.14159 + 2i)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form:
#' # { "Alpha" : { "r" : 3.14159, "i" : 2 } }
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.complex(buf, "complexi", c(1.7 + 2.1i, 97.2))
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form:
#' # { "complexi" : [ { "r": 1.7, i : 2.1}, { "r": 97.2, "i" : 0} ] }
#' 
#' buf <- mongo.bson.buffer.create()
#' values <- c(0.5 + 0.1i, 0.25)
#' names(values) <- c("Theta", "Epsilon")
#' mongo.bson.buffer.append.complex(buf, "Values", values)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form:
#' # { "Values" : { "Theta"   : { "r" : 0.5, "i" : 0.1 },
#' #                "Epsilon" : { " r" : 0.25, "i" : 0 } } }
#' 
#' @export mongo.bson.buffer.append.complex
mongo.bson.buffer.append.complex <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append.complex", buf, name, value)



#' Append a long valued field onto a mongo.bson.buffer
#' 
#' Append a long value or vector of longs onto a \link{mongo.bson.buffer}.
#' 
#' Note that since R has no long (64-bit integer) type, doubles are used in R,
#' but are converted to 64-bit values when stored in the buffer; some loss of
#' precision may occur.
#' 
#' This is the only case in which \code{\link{mongo.bson.buffer.append}()}
#' cannot make the proper guess about what type to encode into the buffer.\cr
#' You must call \code{mongo.bson.buffer.append.long()} explicitly; otherwise,
#' doubles are appended.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param value (double vector) The values(s) to append to the buffer.
#' 
#' If value has a \code{dims} attribute of length > 1, any \code{names} or
#' \code{dimnames} attribute is ignored and a nested array is appended.\cr (Use
#' \code{\link{mongo.bson.buffer.append.object}()} if you want to preserve
#' \code{dimnames}; however, this can't append value as longs).
#' 
#' If value has a names attribute, a subobject is appended and the subfields
#' are given the indicated names.
#' 
#' Otherwise, if more than one element is present in value, the values are
#' appended as a subarray.
#' 
#' In the last case, a single long is appended as the value of the field.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \code{\link{mongo.bson.buffer.append}}.\cr
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.long(buf, "YearSeconds",
#'     365.24219 * 24 * 60 * 60)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form:
#' # { "YearSeconds" : 31556925 }
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.long(buf, "longs", 
#'     c(1, 9087654321, 1234567809))
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form:
#' # { "longs" : [1, 9087654321, 1234567809] }
#' 
#' buf <- mongo.bson.buffer.create()
#' distances <- c(473, 133871000, 188178313)
#' names(distances) <- c("Sol", "Proxima Centari", "Bernard's Star")
#' mongo.bson.buffer.append.long(buf, "Stars", distances)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form:
#' # { "Stars" : { "Sol" : 474,
#' #               "Proxima Centari" : 133871000, 
#' #               "Bernard's Star"  : 188178313 } }
#' 
#' @export mongo.bson.buffer.append.long
mongo.bson.buffer.append.long <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append.long", buf, name, value)



#' Append a double field onto a mongo.bson.buffer
#' 
#' Append a NULL value onto a \link{mongo.bson.buffer}.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \code{\link{mongo.bson.buffer.append}}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.null(buf, "Nil")
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form { "Nil" : NULL }
#' 
#' @export mongo.bson.buffer.append.null
mongo.bson.buffer.append.null <- function(buf, name)
    .Call(".mongo.bson.buffer.append.null", buf, name)



#' Append a undefined field onto a mongo.bson.buffer
#' 
#' Append a undefined value onto a \link{mongo.bson.buffer}.
#' 
#' BSON has a special field type to indicate an undefined value. This function
#' appends such an indicator as the value of a field.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.undefined},\cr \code{\link{mongo.undefined.create}},\cr
#' \code{\link{mongo.bson.buffer.append}}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.undefined(buf, "Undef")
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form { "Undef" : UNDEFINED }
#' 
#' # The same result can be produced by the following code:
#' buf <- mongo.bson.buffer.create()
#' undef <- mongo.undefined.create()
#' mongo.bson.buffer.append(buf, "Undef", undef)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' @export mongo.bson.buffer.append.undefined
mongo.bson.buffer.append.undefined <- function(buf, name)
    .Call(".mongo.bson.buffer.append.undefined", buf, name)



#' Append a string field onto a mongo.bson.buffer
#' 
#' Append an string or vector of strings onto a \link{mongo.bson.buffer}.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param value (string vector) The strings(s) to append to the buffer.
#' 
#' If value has a \code{dims} attribute of length > 1, any \code{names} or
#' \code{dimnames} attribute is ignored and a nested array is appended.\cr (Use
#' \code{\link{mongo.bson.buffer.append.object}()} if you want to preserve
#' \code{dimnames}).
#' 
#' If value has a names attribute, a subobject is appended and the subfields
#' are given the indicated names.
#' 
#' Otherwise, if more than one element is present in value, the strings are
#' appended as a subarray.
#' 
#' In the last case, a single string is appended as the value of the field.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \code{\link{mongo.bson.buffer.append}}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.string(buf, "name", "Joe")
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form { "name" : "Joe" }
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.string(buf, "names", c("Fred", "Jeff", "John"))
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form:
#' # { "names" : ["Fred", "Jeff", "John"] }
#' 
#' buf <- mongo.bson.buffer.create()
#' staff <- c("Mark", "Jennifer", "Robert")
#' names(staff) <- c("Chairman", "President", "Secretary")
#' mongo.bson.buffer.append.string(buf, "board", staff)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form:
#' # { "board" : { "Chairman"  : "Mark",
#' #               "President" : "Jennifer",
#' #               "Secretary" : "Robert" } }
#' 
#' @export mongo.bson.buffer.append.string
mongo.bson.buffer.append.string <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append.string", buf, name, value)



#' Append a time value into a mongo.bson.buffer
#' 
#' Append a date/time value into a \link{mongo.bson.buffer}.
#' 
#' BSON has a special field type to indicate a date/time; these are 64-bit
#' values.
#' 
#' However, R has a 'standard' object of class "POSIXct" used to represent
#' date/time values, such as that returned by Sys.time(). Internally these are
#' a 32-bit integer number of milliseconds since midnight January 1, 1970. On
#' January 19, 2038, 32-bit versions of the the Unix time stamp will cease to
#' work, as it will overflow the largest value that can be held in a signed
#' 32-bit number. At such time, many applications, including R and this driver,
#' will need to address that issue.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param time (integer) A time value.  This may also be an object of\cr class
#' "POSIXct", "POSIXlt" or "mongo.timestamp".
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.timestamp},\cr \code{\link{mongo.timestamp.create}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.time(buf, "Now", Sys.time())
#' b <- mongo.bson.from.buffer(buf)
#' 
#' @export mongo.bson.buffer.append.time
mongo.bson.buffer.append.time <- function(buf, name, time)
    .Call(".mongo.bson.buffer.append.time", buf, name, time)



#' Append a timestamp value into a mongo.bson.buffer
#' 
#' Append a timestamp value into a \link{mongo.bson.buffer}.
#' 
#' \link{mongo.timestamp} objects extend the "POSIXct" class to include an
#' attrubute "increment".
#' 
#' See \code{\link{mongo.bson.buffer.append.time}()}.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param value A (\link{mongo.timestamp}) value as created by
#' \code{\link{mongo.timestamp.create}()}.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.timestamp.create},\cr
#' \code{\link{mongo.bson.buffer.append.time}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr \link{mongo.bson},\cr
#' \link{mongo.bson.buffer}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.timestamp(buf, "Now-27",
#'     mongo.timestamp.create(Sys.time(), 27))
#' b <- mongo.bson.from.buffer(buf)
#' 
#' @export mongo.bson.buffer.append.timestamp
mongo.bson.buffer.append.timestamp <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append.timestamp", buf, name, value)



#' Append a code field onto a mongo.bson.buffer
#' 
#' Append a javascript code value onto a \link{mongo.bson.buffer}.
#' 
#' BSON has a special field type to indicate javascript code. This function
#' appends such an indicator as the type of a field with its value.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param value string
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.code},\cr \code{\link{mongo.code.create}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr \link{mongo.bson},\cr
#' \link{mongo.bson.buffer}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.code(buf, "SetXtoY", "x = y")
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form:
#' # { "SetXtoY" : (CODE) "x = y" }
#' 
#' # The same result can be produced by the following code:
#' buf <- mongo.bson.buffer.create()
#' code <- mongo.code.create("x = y")
#' mongo.bson.buffer.append(buf, "SetXtoY", code)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' @export mongo.bson.buffer.append.code
mongo.bson.buffer.append.code <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append.code", buf, name, value)



#' Append a code field with a scope onto a mongo.bson.buffer
#' 
#' Append a javascript code value with a scope object onto a
#' \link{mongo.bson.buffer}.
#' 
#' BSON has a special field type to indicate javascript code with a scope. This
#' function appends such an indicator as the type of a field with its value.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param value \link{mongo.code.w.scope} The scoped javascript code.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.code.w.scope},\cr
#' \code{\link{mongo.code.w.scope.create}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.from.list}},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.bson}.
#' @examples
#' 
#' scope <- mongo.bson.from.list(list(scopevar="scopevalue"))
#' buf <- mongo.bson.buffer.create()
#' codeWscope <- mongo.code.w.scope.create("y = x", scope)
#' mongo.bson.buffer.append.code.w.scope(buf, "CodeWscope1",
#'      codeWscope)
#' 
#' # mongo.bson.buffer.append() will give the same result
#' # as it can detect the mongo.code.w.scope object
#' mongo.bson.buffer.append(buf, "CodeWscope2", codeWscope)
#' 
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form:
#' # { "CodeWscope1" : (CODEWSCOPE) "y = x"
#' #        (SCOPE) { "scopevar" : "scopevalue" },
#' #   "CodeWscope2" : (CODEWSCOPE) "y = x"
#' #        (SCOPE) { "scopevar" : "scopevalue" } }
#' 
#' 
#' @export mongo.bson.buffer.append.code.w.scope
mongo.bson.buffer.append.code.w.scope <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append.code.w.scope", buf, name, value)



#' Append a symbol field onto a mongo.bson.buffer
#' 
#' Append a symbol value onto a \link{mongo.bson.buffer}.
#' 
#' BSON has a special field type to indicate a symbol. This function appends
#' such an indicator as the type of a field with its value.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param value (string) The value of the symbol.
#' 
#' Note that the value may simply be a string representing the symbol's value
#' and not necessarily a \link{mongo.symbol} object.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.symbol},\cr \code{\link{mongo.symbol.create}},\cr
#' \code{\link{mongo.bson.buffer.append}}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.symbol(buf, "A", "Alpha")
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # The above produces a BSON object of the form { "A" : (SYMBOL) "Alpha" }
#' 
#' # The same result can be produced by the following code:
#' buf <- mongo.bson.buffer.create()
#' sym <- mongo.symbol.create("Alpha")
#' mongo.bson.buffer.append(buf, "A", sym)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' @export mongo.bson.buffer.append.symbol
mongo.bson.buffer.append.symbol <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append.symbol", buf, name, value)



#' Append a timestamp value into a mongo.bson.buffer
#' 
#' Append a regular expression value into a \link{mongo.bson.buffer}.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param value (\link{mongo.regex}) A regular expression as created\cr by
#' \code{\link{mongo.regex.create}()}.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.regex.create},\cr
#' \code{\link{mongo.bson.buffer.append.regex}},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr \link{mongo.bson},\cr
#' \link{mongo.bson.buffer}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' regex <- mongo.regex.create("acme.*corp", options="i")
#' mongo.bson.buffer.append.regex(buf, "MatchAcme", regex)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' @export mongo.bson.buffer.append.regex
mongo.bson.buffer.append.regex <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append.regex", buf, name, value)



#' Append a raw (binary) field onto a mongo.bson.buffer
#' 
#' Append raw (binary) data onto a \link{mongo.bson.buffer}.
#' 
#' BSON has a special field type to indicate binary data. This function appends
#' such an indicator as the type of a field with its value.
#' 
#' If value has a \code{dims} attribute of length > 1, any \code{names} or
#' \code{dimnames} attribute is ignored and a nested array is appended.\cr (Use
#' \code{\link{mongo.bson.buffer.append.object}()} if you want to preserve
#' \code{dimnames}).
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param value (raw) the binary data.
#' @param subtype (as.integer) The binary data subtype.  If subtype == NULL,
#' the "subtype" attribute of the raw is used.  If this is not present,
#' mongo.binary.binary is used.  The following constants are defined: \itemize{
#' \item\code{\link{mongo.binary.binary}} (0L)
#' \item\code{\link{mongo.binary.function}} (1L)
#' \item\code{\link{mongo.binary.old}} (2L)
#' \item\code{\link{mongo.binary.uuid}} (3L)
#' \item\code{\link{mongo.binary.md5}} (5L)
#' \item\code{\link{mongo.binary.user}} (128L) }
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \code{\link{mongo.bson.buffer.append}},\cr \link{mongo.bson},\cr
#' \link{mongo.bson.buffer}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' bin <- raw(3)
#' for (i in 0:2)
#'     bin[i] <- as.raw(i * 3)
#' mongo.bson.buffer.append.raw(buf, "bin1", bin)
#' 
#' # Note that mongo.bson.buffer.append()
#' # will detect whether the value parameter 
#' # is a raw object and append the appropriate value.
#' 
#' mongo.bson.buffer.append(buf, "bin2", bin)  # gives same result
#' 
#' @export mongo.bson.buffer.append.raw
mongo.bson.buffer.append.raw <- function(buf, name, value, subtype=NULL)
    .Call(".mongo.bson.buffer.append.raw", buf, name, value, subtype)



#' Append a OID into a mongo.bson.buffer
#' 
#' Append a OID (Object ID) value into a \link{mongo.bson.buffer}.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param value (\link{mongo.oid}) An OID value.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \link{mongo.oid.create},\cr \code{\link{mongo.bson.buffer.append}}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.oid(buf, "Now", mongo.oid.create())
#' b <- mongo.bson.from.buffer(buf)
#' 
#' @export mongo.bson.buffer.append.oid
mongo.bson.buffer.append.oid <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append.oid", buf, name, value)



#' Append a mongo.bson object into a mongo.bson.buffer
#' 
#' Append a \link{mongo.bson} object into a \link{mongo.bson.buffer} as a
#' subobject.
#' 
#' Note that \code{\link{mongo.bson.buffer.append}()} will detect if its value
#' parameter is a mongo.bson object and perform the same action as this
#' function.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the subobject field appended to the
#' buffer.
#' @param value (\link{mongo.bson}) a mongo.bson object.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \code{\link{mongo.bson.from.list}},\cr
#' \code{\link{mongo.bson.buffer.append}}.
#' @examples
#' 
#' name <- mongo.bson.from.list(list(first="Joe", last="Smith"))
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.bson(buf, "name", name)
#' mongo.bson.buffer.append.string(buf, "city", "New York")
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # the above will create a mongo.bson object of the following form:
#' # { "name" : { "first" : "Joe", "last" : "Smith" }, "city" : "New York" }
#' 
#' @export mongo.bson.buffer.append.bson
mongo.bson.buffer.append.bson <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append.bson", buf, name, value)



#' Append a mongo.bson.iterator's element into a mongo.bson.buffer
#' 
#' Append a \link{mongo.bson.iterator}'s element into a
#' \link{mongo.bson.buffer}.
#' 
#' \code{\link{mongo.bson.buffer.append}()} will detect if its value parameter
#' is a mongo.bson.iterator object and perform the same action as this
#' function.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the subobject field appended to the
#' buffer.
#' 
#' If NULL, the name appended will come from the element pointed to by the
#' iterator.
#' @param value A (\link{mongo.bson.iterator}) object.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \code{\link{mongo.bson.find}},\cr \code{\link{mongo.bson.from.list}},\cr
#' \code{\link{mongo.bson.buffer.append}}.
#' @examples
#' 
#' name <- mongo.bson.from.list(list(first="Joe", last="Smith"))
#' iter <- mongo.bson.find(name, "last")
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.element(buf, "last", iter)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # the above will create a mongo.bson object (b) of the following form:
#' # { "last" : "Smith" }
#' 
#' @export mongo.bson.buffer.append.element
mongo.bson.buffer.append.element <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append.element", buf, name, value)



#' Append a name/value pair into a mongo.bson.buffer
#' 
#' Append a name/value pair into a \link{mongo.bson.buffer}.
#' 
#' This function is a generic version of many 'append' functions.  It will
#' detect the type of the \code{value} parameter and perform the same action as
#' the specific functions.  These functions are: \itemize{ \item
#' \code{\link{mongo.bson.buffer.append.int}()} \item
#' \code{\link{mongo.bson.buffer.append.string}()} \item
#' \code{\link{mongo.bson.buffer.append.bool}()} \item
#' \code{\link{mongo.bson.buffer.append.double}()} \item
#' \code{\link{mongo.bson.buffer.append.complex}()} \item
#' \code{\link{mongo.bson.buffer.append.null}()} \item
#' \code{\link{mongo.bson.buffer.append.undefined}()} \item
#' \code{\link{mongo.bson.buffer.append.symbol}()} \item
#' \code{\link{mongo.bson.buffer.append.code}()} \item
#' \code{\link{mongo.bson.buffer.append.code.w.scope}()} \item
#' \code{\link{mongo.bson.buffer.append.raw}()} \item
#' \code{\link{mongo.bson.buffer.append.time}()} \item
#' \code{\link{mongo.bson.buffer.append.timestamp}()} \item
#' \code{\link{mongo.bson.buffer.append.regex}()} \item
#' \code{\link{mongo.bson.buffer.append.oid}()} \item
#' \code{\link{mongo.bson.buffer.append.bson}()} \item
#' \code{\link{mongo.bson.buffer.append.element}()} \item
#' \code{\link{mongo.bson.buffer.append.list}()} }
#' 
#' \code{\link{mongo.bson.buffer.append.long}()} is missing from the above list
#' since R has no 64-bit long integer type.  If you wish a value to be stored
#' in the BSON data as a long you must explicity call that function.
#' 
#' All of the above functions will lose the attributes of the object other than
#' "names". When vectors of length > 1 are appended, "names" are preserved.\cr
#' \code{\link{mongo.bson.buffer.append.object}()} gets around this shortcoming
#' and allows most R objects to be stored in a database without loss of
#' attributes.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param value The value of the field.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' # Append a string
#' mongo.bson.buffer.append(buf, "name", "Joe")
#' # Append a date/time
#' mongo.bson.buffer.append(buf, "created", Sys.time())
#' # Append a NULL
#' mongo.bson.buffer.append(buf, "cars", NULL)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' @export mongo.bson.buffer.append
mongo.bson.buffer.append <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append", buf, name, value)



#' Append a list onto a mongo.bson.buffer
#' 
#' Append a list onto a \link{mongo.bson.buffer}.
#' 
#' Note that the value parameter must be a true list, not an vector of a single
#' atomic type.
#' 
#' Also note that this function is recursive and will append items that are
#' lists themselves as subobjects.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param value (list) The list to append to the buffer as a subobject.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \code{\link{mongo.bson.buffer.append}}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' l <- list(fruit = "apple", hasSeeds = TRUE)
#' mongo.bson.buffer.append.list(buf, "item", l)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # this produces a BSON object of the form:
#' # { "item" : { "fruit" : "apple", "hasSeeds" : true } }
#' 
#' @export mongo.bson.buffer.append.list
mongo.bson.buffer.append.list <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append.list", buf, name, value)



#' Append an R object onto a mongo.bson.buffer
#' 
#' Append an R object onto a \link{mongo.bson.buffer}.
#' 
#' This function allows you to store higher level R objects in the database
#' without losing their attribute information.  It will correctly handle data
#' frames, matrices and arrays for instance; although, empty objects, such as a
#' data frame with no rows, are not permitted.
#' 
#' Note that the names attribute will not be preserved if the object is
#' multidimensional (although dimnames will be).
#' 
#' The object's value will look like this in the buffer: \preformatted{ { ...
#' name : { R_OBJ : true, value : xxx, attr : { attr1 : yyy, attr2 : zzz } }
#' ...  } }
#' 
#' \code{name} will be substituted with the value of the \code{name}
#' parameter.\cr \code{xxx} will be substituted with the low level value of the
#' object (as would be appended by
#' \code{\link{mongo.bson.buffer.append}()}).\cr \code{attr1} and \code{attr2}
#' will be substituted with the names of attributes.\cr \code{yyy} and
#' \code{zzz} will be substituted with the values of those attributes.\cr
#' 
#' Note that it is inadvised to construct this wrapper manually as
#' \code{\link{mongo.bson.value}()} and
#' \code{\link{mongo.bson.iterator.value}()} bypass the special checking and
#' handling that is done by R code that set attributes.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the field appended to the buffer.
#' @param value (object) The object to append to the buffer as a subobject.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \code{\link{mongo.bson.buffer.append}},\cr
#' \code{\link{mongo.bson.value}},\cr \code{\link{mongo.bson.iterator.value}}
#' @examples
#' 
#' age <- c(5, 8)
#' height <- c(35, 47)
#' d <- data.frame(age=age, height=height)
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append.object(buf, "table", d)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # this produces a BSON object of the form:
#' # { "table" : { "R_OBJ" : true,
#' #               "value" : {
#' #                    "age"    : [ 5, 8 ],
#' #                    "height" : [35, 47 ]
#' #               },
#' #               "attr" : {
#' #                  "row.names" : [ -2147483648, -2 ],
#' #                  "class" : "data.frame"
#' #               }
#' #             }
#' # }
#' # row.names is stored in the compact form used for integer row names.
#' 
#' @export mongo.bson.buffer.append.object
mongo.bson.buffer.append.object <- function(buf, name, value)
    .Call(".mongo.bson.buffer.append.object", buf, name, value)



#' Start a subobject within a mongo.bson.buffer
#' 
#' BSON documents may themselves contain nested documents.  Call this function
#' to start a subobject within a \link{mongo.bson.buffer}.
#' 
#' \code{\link{mongo.bson.buffer.finish.object}()} must be called when finsihed
#' appending subfields.\cr (\code{mongo.bson.buffer.start.object()},
#' \code{mongo.bson.buffer.start.array()})\cr and
#' \code{mongo.bson.buffer.finish.object()} may be called in a stackwise (LIFO)
#' order to further nest documents and arrays.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the subobject to be appended to the
#' buffer.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \code{\link{mongo.bson.buffer.finish.object}},\cr
#' \code{\link{mongo.bson.buffer.start.array}},\cr
#' \code{\link{mongo.bson.buffer.append}}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.start.object(buf, "name")
#' mongo.bson.buffer.append(buf, "first", "Jeff")
#' mongo.bson.buffer.append(buf, "last", "Davis")
#' mongo.bson.buffer.finish.object(buf)
#' mongo.bson.buffer.append(buf, "city", "Toronto")
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # the above produces a BSON object of the form:
#' # { "name" : { "first" : "Jeff", "last" : "Davis" }, "city" : "Toronto" }
#' 
#' @export mongo.bson.buffer.start.object
mongo.bson.buffer.start.object <- function(buf, name)
    .Call(".mongo.bson.buffer.start.object", buf, name)



#' Start an array within a mongo.bson.buffer
#' 
#' Call this function to start an array within a \link{mongo.bson.buffer}.\cr
#' \code{\link{mongo.bson.buffer.finish.object}()} must be called when finished
#' appending the elements of the array.
#' 
#' (\code{mongo.bson.buffer.start.object()},
#' \code{mongo.bson.buffer.start.array()}) and\cr
#' \code{mongo.bson.buffer.finsih.object()} may be called in a stackwise (LIFO)
#' order to further nest arrays and documents.
#' 
#' The names of the elements appended should properly be given sequentially
#' numbered strings.
#' 
#' Note that arrays will be automatically appended by the 'append' functions
#' when appending vectors (containing more than one element) of atomic types.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object to which to append.
#' @param name (string) The name (key) of the array to be appended to the
#' buffer.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \code{\link{mongo.bson.buffer.finish.object}},\cr
#' \code{\link{mongo.bson.buffer.start.array}},\cr
#' \code{\link{mongo.bson.buffer.append}}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.start.array(buf, "Fibonacci")
#' x <- 0
#' mongo.bson.buffer.append.int(buf, "0", x)
#' y <- 1
#' mongo.bson.buffer.append.int(buf, "1", y)
#' for (i in 2:8) {
#'     z <- x + y
#'     mongo.bson.buffer.append.int(buf, as.character(i), z)
#'     x <- y
#'     y <- z
#' }
#' mongo.bson.buffer.finish.object(buf)
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # the above produces a BSON object of the form:
#' # { "Fibonacci" : [ 0, 1, 1, 2, 3, 5, 8, 13, 21 ] }
#' 
#' @export mongo.bson.buffer.start.array
mongo.bson.buffer.start.array <- function(buf, name)
    .Call(".mongo.bson.buffer.start.array", buf, name)



#' Finish a subobject or array within a mongo.bson.buffer
#' 
#' BSON documents may themselves contain nested documents.  Call this function
#' to finish a subobject within a \link{mongo.bson.buffer}.
#' 
#' \code{\link{mongo.bson.buffer.start.object}()} and
#' \code{mongo.bson.buffer.finish.object()} may be called in a stackwise (LIFO)
#' order to further nest documents.
#' 
#' This function must also be called to finish arrays.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) The buffer object on which to finish a
#' subobject.
#' @return TRUE if successful; otherwise, FALSE if an error occured appending
#' the data.
#' @seealso \link{mongo.bson},\cr \link{mongo.bson.buffer},\cr
#' \code{\link{mongo.bson.buffer.start.object}},\cr
#' \code{\link{mongo.bson.buffer.start.array}},\cr
#' \code{\link{mongo.bson.buffer.append}}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.start.object(buf, "name")
#' mongo.bson.buffer.append(buf, "first", "Jeff")
#' mongo.bson.buffer.append(buf, "last", "Davis")
#' mongo.bson.buffer.finish.object(buf)
#' mongo.bson.buffer.append(buf, "city", "Toronto")
#' b <- mongo.bson.from.buffer(buf)
#' 
#' # the above produces a BSON object of the form:
#' # { "name" : { "first" : "Jeff", "last" : "Davis" }, "city" : "Toronto" }
#' 
#' @export mongo.bson.buffer.finish.object
mongo.bson.buffer.finish.object <- function(buf)
    .Call(".mongo.bson.buffer.finish.object", buf)



#' Get the size of a mongo.bson.buffer object
#' 
#' Get the number of bytes which would be taken up by the BSON data when the
#' buffer is converted to a mongo.bson object with
#' \code{\link{mongo.bson.from.buffer}()}.
#' 
#' 
#' @param buf (\link{mongo.bson.buffer}) the mongo.bson.buffer object to
#' examine.
#' @return (integer) the number of bytes which would be taken up by the BSON
#' data with the buffer is converted to a mongo.bson object with
#' \code{\link{mongo.bson.from.buffer}()}.
#' @seealso \link{mongo.bson.buffer},\cr \link{mongo.bson}.
#' @examples
#' 
#' buf <- mongo.bson.buffer.create()
#' mongo.bson.buffer.append(buf, "name", "Fred")
#' mongo.bson.buffer.append(buf, "city", "Dayton")
#' # both should report 37
#' print(mongo.bson.buffer.size(buf))
#' y <- mongo.bson.from.buffer(buf)
#' print(mongo.bson.size(y))
#' 
#' @export mongo.bson.buffer.size
mongo.bson.buffer.size <- function(buf)
    .Call(".mongo.bson.buffer.size", buf)

