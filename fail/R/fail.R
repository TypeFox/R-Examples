#' Create a file abstraction interface layer (FAIL) object.
#'
#' The general idea is to not bother about file path joining or file extensions.
#' Instead, FAIL offers a key-value like interface to RData files in a specified directory.
#' The filename (without extension) acts as the key while the stored R objects are the values.
#' Fail provides an interface to the basic file system actions: listing, reading / loading,
#' writing / saving, removing and applying functions on files.
#' An implemented cache mechanism can be used to avoid repeated disk reads.
#'
#' @param path [\code{character(1)}]\cr
#'   Path to work in, will be created if it does not exists.
#' @param extension [\code{character(1)}]\cr
#'   File extension to work with.
#'   Default is \dQuote{RData}.
#' @param all.files [\code{logical(1)}]\cr
#'   Also include hidden files, i.e. files whose name start with a dot (\dQuote{.}).
#'   Default is \code{FALSE}.
#' @param use.cache [\code{logical(1)}]\cr
#'   Use a memory cache per global default.
#'   Global option which can locally be overwritten in most functions.
#'   Default is \code{FALSE}
#' @param simplify [\code{character(1)}]\cr
#'   If only one object is stored in a R data file,
#'   should the return value be simplified?
#'   If set to \code{TRUE},
#'   instead of a list containing one element the object itself will be returned.
#' @return Object of class \code{fail}. See details.
#' @details
#'   For a quick introduction on the usage, see \url{https://github.com/mllg/fail}.
#'
#'   An object with the following functions is returned by \code{fail}:
#'   \describe{
#'     \item{\code{ls(pattern=NULL)}}{
#'       Function to list keys in directory \code{path} matching a regular expression pattern \code{pattern}.
#'       Returns a character vector of keys.
#'     }
#'     \item{\code{get(key, use.cache)}}{
#'       Function to load a file identified by \code{key} from directory \code{path}.
#'       To load many objects at once, use \code{as.list}, \code{assign} or \code{get} together with \code{\link[base]{lapply}}.
#'       Argument \code{use.cache} can be set to temporarily overwrite the global \code{use.cache} flag.
#'     }
#'     \item{\code{put(..., li, keys, use.cache)}}{
#'       Function to save objects to directory \code{path}.
#'       Names for objects provided via \code{...} will be looked up or can be provided using a \code{key = value} syntax.
#'       More objects can be passed as a named list using the argument \code{li}: Each list item will be saved to a separate file.
#'       If you provide \code{keys} as a character vector, these names will be taken for the arguments passed via \code{...}.
#'       Argument \code{use.cache} temporarily overwrites the global \code{use.cache} flag.
#'       Returns a character vector of stored keys.
#'     }
#'     \item{\code{remove(keys)}}{
#'       Function to remove files identified by \code{keys} from directory \code{path}.
#'       Returns a character vector of deleted keys.
#'     }
#'     \item{\code{apply(FUN, ..., keys, use.cache, simplify=FALSE, use.names=TRUE)}}{
#'       Apply function \code{FUN} on files identified by \code{keys}.
#'       \code{keys} defaults to all keys available and will be used to name the returned list.
#'       The loaded R objects will be past unnamed as first argument.
#'       Use \code{...} for additional function arguments.
#'       Argument \code{use.cache} can be set to temporarily overwrite the global \code{use.cache} flag.
#'       For arguments \code{simplify} and \code{use.names}, see \code{\link[base]{lapply}}.
#'    }
#'     \item{\code{mapply(FUN, ..., keys, use.cache, moreArgs = NULL, simplify=FALSE, use.names=TRUE)}}{
#'       Apply function \code{FUN} on files identified by \code{keys}.
#'       \code{keys} defaults to all keys available and will be used to name the returned list.
#'       The function \code{FUN} must have the formal arguments \dQuote{key} and \dQuote{value}.
#'       Both key and value will be passed named.
#'       Use \code{...} and/or \code{moreArgs} for additional function arguments.
#'       Argument \code{use.cache} can be set to temporarily overwrite the global \code{use.cache} flag.
#'       For arguments \code{moreArgs}, \code{simplify} and \code{use.names}, see \code{\link[base]{mapply}}.
#'    }
#'    \item{\code{as.list(keys, use.cache)}}{
#'       Return a named list of objects identified by \code{keys}. \code{keys} defaults to all keys available.
#'       Argument \code{use.cache} can be set to temporarily overwrite the global \code{use.cache} flag.
#'    }
#'    \item{\code{assign(keys, envir=parent.frame(), use.cache)}}{
#'       Assigns all objects identified by the character vector \code{keys} in the environment \code{envir}.
#'       Argument \code{use.cache} can be set to temporarily overwrite the global \code{use.cache} flag.
#'       Returns a character vector of assigned keys.
#'    }
#'    \item{\code{clear(keys)}}{
#'       Clear the cache to free memory. \code{keys} defaults to all keys available.
#'       Returns a character vector of cleared keys.
#'    }
#'    \item{\code{cached()}}{
#'       Returns a character vector of keys of cached objects.
#'    }
#'    \item{\code{size(keys, unit="b")}}{
#'       Get the file size in Bytes of the files identified by \code{keys}. \code{keys} defaults to all keys available.
#'       Argument \code{unit} accepts \dQuote{b}, \dQuote{Kb}, \dQuote{Mb} and \dQuote{Gb} and can be used to convert Bytes to KiloBytes, MegaBytes or GigaBytes, respectively.
#'    }
#'    \item{\code{info()}}{
#'       Returns a named list with \code{path}, \code{extension} and \code{use.cache}.
#'       Internally used for the \code{\link[base]{print}} method with a much nicer summary of the FAIL object.
#'    }
#'   }
#'   Furthermore, the package provides S3 methods for \code{\link[base]{print}} and \code{\link[base]{as.list}}.
#'
#'   Be aware of the following restriction regarding file names and keys:
#'   The package performs some basic checks for illegal characters on the key names.
#'   In principle all characters matching the pattern \dQuote{[a-zA-Z0-9._-]} are allowed and should work on most or all file systems.
#'   But be careful with key names which are not compatible with R's variable naming restrictions, e.g. using the minus character or
#'   key names starting with a number: these provoke unwanted side effects and will result in errors if used with \code{assign}.
#'
#'   If two files would collide on case-insensitive file systems like Windows' NTFS, the package will throw warnings. Best practice
#'   is to not rely on case sensitivity.
#'
#' @export
#' @examples
#' # initialize a FAIL in a temporary directory
#' path <- tempfile("")
#' files <- fail(path)
#'
#' # save x and y, vectors of random numbers
#' x <- runif(100)
#' files$put(x, y = runif(100))
#'
#' # save columns of the iris data set as separate files
#' files$put(li = as.list(iris))
#'
#' # load all RData files in a named list as a one-liner
#' as.list(fail(path))
#'
#' # load a single object from the file system
#' files$get("Species")
#' files$as.list(c("x", "y"))
#'
#' # remove an object (and related file)
#' files$remove("Species")
#'
#' # apply a function over files
#' files$apply(mean)
#' files$mapply(function(key, value) sprintf("%s -> %f", key, mean(value)), simplify = TRUE)
#'
#' # show file size informations
#' files$size(unit = "Mb")
#'
#' # get an object and cache it
#' files$get("x", use.cache = TRUE)
#' files$cached()
#' files$clear()
#' files$cached()
#'
#' # assign variables in the current environment
#' files$assign("y")
#' mean(y)
fail = function(path = getwd(), extension = "RData", all.files = FALSE, use.cache = FALSE, simplify = TRUE) {
  ### argument checks
  .self = list(
    path = checkPath(path),
    extension = checkExtension(extension),
    all.files = asFlag(all.files),
    use.cache = asFlag(use.cache),
    simplify = asFlag(simplify, na.ok = TRUE),
    cache = Cache(),
    loadFun = loadRData,
    saveFun = saveRData
  )
  checkCollision(Ls(.self))
  setClasses(makeObject(.self), "fail")
}
