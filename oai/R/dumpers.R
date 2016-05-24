#' Result dumpers
#'
#' Result dumpers are functions allowing to handle the chunks of results from
#' OAI-PMH service "on the fly". Handling can include processing, writing to
#' files, databases etc.
#'
#' Often the result of a request to a OAI-PMH service are so large that it is
#' split into chunks that need to be requested separately using
#' \code{resumptionToken}. By default functions like
#' \code{\link{list_identifiers}} or \code{\link{list_records}} request these
#' chunks under the hood and return all concatenated in a single R object. It is
#' convenient but insufficient when dealing with large result sets that might
#' not fit into RAM. A result dumper is a function that is called on each result
#' chunk. Dumper functions can write chunks to files or databases, include initial
#' pre-processing or extraction, and so on.
#'
#' A result dumper needs to be function that accepts at least the arguments:
#' \code{res}, \code{args}, \code{as}. They will get values by the enclosing
#' function internally. There may be additional arguments, including \code{...}.
#' Dumpers should return \code{NULL} or a value that will
#' be collected and returned by the function calling the dumper (e.g.
#' \code{\link{list_records}}).
#'
#' Currently result dumpers can be used with functions:
#' \code{\link{list_identifiers}},
#' \code{\link{list_records}}, and
#' \code{\link{list_sets}}.
#' To use a dumper with one of these functions you need to:
#' \enumerate{
#' \item{Pass it as an additional argument \code{dumper}}
#' \item{Pass optional addtional arguments to the dumper function in a list as the \code{dumper_args} argument}
#' }
#' See Examples. Below we provide more details on the dumpers currently implemented.
#'
#'
#'
#' @param res results, depends on \code{as}, not to be specified by the user
#' @param args list, query arguments, not to be specified by the user
#' @param as character, type of result to return, not to be specified by the
#'   user
#' @param ... arguments passed to/from other functions
#'
#' @return Dumpers should return \code{NULL} or a value that will be collected
#'   and returned by the function using the dumper.
#'
#' @references OAI-PMH specification
#'   \url{https://www.openarchives.org/OAI/openarchivesprotocol.html}
#'
#' @seealso Functions supporting the dumpers:
#' \code{\link{list_identifiers}},
#' \code{\link{list_sets}}, and
#' \code{\link{list_records}}
#'
#' @aliases dumpers
#' @name dumpers
#' @example man-roxygen/dumpers.R
NULL



#' @rdname dumpers
#'
#' @details
#'
#' \code{dump_raw_to_txt} writes raw XML to text files. It requires
#' \code{as=="raw"}. File names are created using \code{\link{tempfile}}. By
#' default they are written in the current working directory and have a format
#' \code{oaidump*.xml} where \code{*} is a random string in hex.
#'
#' @param file_pattern,file_dir,file_ext character respectively: initial part of
#'   the file name, directory name, and file extension used to create file
#'   names. These arguments are passed to \code{\link{tempfile}} arguments
#'   \code{pattern}, \code{tmpdir}, and \code{fileext} respectively.
#'
#' @return \code{dump_raw_to_txt} returns the name of the created file.
#'
#' @export
dump_raw_to_txt <- function(res, args, as, file_pattern="oaidump", file_dir=".",
                            file_ext=".xml") {
  stopifnot(as == "raw")
  filename <- tempfile(pattern=file_pattern, tmpdir=file_dir, fileext = file_ext)
  cat( as.character(res), file=filename )
  filename
}


#' @rdname dumpers
#'
#' @details
#'
#' \code{dump_to_rds} saves results in an \code{.rds} file via \code{\link{saveRDS}}.
#' Type of object being saved is determined by the \code{as} argument. File names
#' are generated in the same way as by \code{dump_raw_to_txt}, but with default
#' extension \code{.rds}.
#'
#' @return \code{dump_to_rds} returns the name of the created file.
#'
#' @export
dump_to_rds <- function(res, args, as, file_pattern="oaidump", file_dir=".", file_ext=".rds") {
  filename <- tempfile(pattern=file_pattern, tmpdir=file_dir, fileext = file_ext)
  saveRDS(res, file=filename)
  filename
}





#' @rdname dumpers
#'
#' @details
#'
#' \code{dump_xml_to_db} writes raw XML to a single text column of a table in a
#' database. Requires \code{as == "raw"}. Database connection \code{dbcon}
#' should be a connection object as created by \code{\link[DBI]{dbConnect}} from
#' package \pkg{DBI}. As such, it can connect to any database supported by
#' \pkg{DBI}. The records are written to a field \code{field_name} in a table
#' \code{table_name} using \code{\link[DBI]{dbWriteTable}}. If the table does not
#' exist, it is created. If it does, the records are appended. Any additional
#' arguments are passed to \code{\link[DBI]{dbWriteTable}}.
#'
#' @param dbcon \pkg{DBI}-compliant database connection
#' @param table_name character, name of the database table to write into
#' @param field_name character, name of the field in database table to write
#'   into
#'
#' @return \code{dump_xml_to_db} returns \code{NULL}
#'
#' @export
dump_raw_to_db <- function(res, args, as, dbcon, table_name, field_name, ...) {
  stopifnot( as == "raw" )
  dframe <- data.frame(x=as.character(res), stringsAsFactors = FALSE)
  names(dframe) <- field_name
  DBI::dbWriteTable(con=dbcon, name=table_name, value=dframe, row.names=FALSE, append=TRUE, ...)
  invisible(NULL)
}









# Checking dumper arguments
valid_dumper <- function(dumper, dumper_args) {
  rval <- NULL
  stopifnot(is.function(dumper))

  # Arguments of dumper function
  dargs <- formals(dumper)
  darg_has_default <- !sapply(dargs, is.symbol)
  d_has_ddd <- "..." %in% names(dargs)

  # Dumper has minimal obligatory args
  a <- c("res", "args", "as")
  has_a <- a %in% names(dargs)
  if(!all(has_a))
    rval <- c(rval,
              paste("dumper misses obligatory arguments:",
                    paste(a[has_a], collapse=", ") )
    )

  # Additional user arguments for dumper: dumper_args

  # Dumper requires user arguments which are not supplied
  aa <- setdiff( names(dargs[!darg_has_default]), c(a, "..."))
  z <- aa %in% names(dumper_args)
  if(any(!z))
    rval <- c(rval,
              paste("required user arguments (dumper_args) missing:",
                    paste(aa[!z], collapse=", ") )
              )



  if(!is.null(dumper_args)) {
    uargs <- names(dumper_args)

    # User supplied res/args/as argument(s)
    z <- a %in% names(dumper_args)
    if( any(z) )
      rval <- c(rval, paste("dumper_args not allowed:",
                            paste(a[z], collapse=", ") ))

    # User arguments that will not be accepted by dumper
    z <- uargs %in% names(dargs)
    if( any(!z) && !d_has_ddd )
      rval <- c(rval, paste("user arguments (dumper_args) not accepted by dumper:",
                            paste(names(dumper_args)[!z], collapse=", ") )
      )
  }

  # Returning
  if(is.null(rval)) {
    return(TRUE)
  } else {
    stop( paste(rval, collapse="\n"))
  }
}
