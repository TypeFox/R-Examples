#############################################################################
# io.R
#############################################################################

#' @include style.R
NULL

#' Initialises a pubprint object
#'
#' \code{pubprint} returns an empty pubprint S3 object.
#'
#' This function initialises an empty pubprint S3 object and returns it.
#' This is mandatory for using the pull and push functions of the pubprint
#' package.
#'
#' @seealso See \code{\link{pubprint-package}} for package documentation.
#'
#' @examples
#' ppo <- pubprint()
#' ppo
#'
#' @export
pubprint <- function()
{
    x <- list(pipe = list(),
              memory = list())

    class(x) <- "pubprint"

    return(x)
}

#' Adds another item to an object
#'
#' \code{push<-} is a generic function and is used to put another item to a stack
#' or pipe.
#'
#' There is no default function, so you have to see a specific \code{push<-}
#' function for further information.
#'
#' @param x an object used to select a method.
#' 
#' @param ... further arguments passed to or from other methods.
#'
#' @param value an item pushed to \code{x}.
#'
#' @return The updated object.
#' 
#' @export
`push<-` <- function(x, ..., value) UseMethod("push<-")

#' Adds another item to a pubprint object
#'
#' \code{push<-.pubprint} adds the given item to named memory or pipe of a
#' pubprint object.
#'
#' No further details.
#'
#' @param x a pubprint object to which \code{value} is added.
#'
#' @param item numeric or character. If \code{item} is a numeric, \code{value}
#' is added to pipe. If \code{item} is a character, \code{value} is added to
#' named memory.  A warning is thrown, if an existing item is overwritten.
#'
#' @param add logical, indicating if \code{value} is added to an existing item.
#' If \code{item} is specified, \code{value} is added to this item, else
#' argument \code{n} is used.
#'
#' @param n numeric. If \code{item} is missing and \code{add} is true,
#' \code{n} indicates to which pipe position (backwards) \code{value} is added.
#' Therefore, \code{n = 1} adds the item to the last pipe item, \code{n = 2}
#' to the second last item and so on.
#' 
#' @param ... further arguments passed to or from other methods.
#'
#' @param value an item pushed to \code{x}.
#'
#' @return The updated pubprint object.
#'
#' @seealso See \code{\link{push<-}} for the generic function,
#' \code{\link{pull}} to extract items again.
#'
#' @examples
#' ppo <- pubprint()
#' push(ppo) <- t.test(1:100, 2:101)
#' push(ppo, add = TRUE) <- .8123 # add d value to last pipe item
#' push(ppo, item = "i1") <- t.test(1:30, 2:31)
#'
#' pull(ppo)
#' pull(ppo, item = "i1")
#' 
#' @export
`push<-.pubprint` <- function(x, 
                                 item, 
                                 add = FALSE,
                                 n = 1,
                                 ..., 
                                 value)
{
    if (add)
    {
        if (missing(item))
        {
            mypos <- length(x$pipe) + 1 - n
            x$pipe[[mypos]] <- c(x$pipe[[mypos]], list(value))
        }
        else if ("numeric" == class(item))
        {
            x$pipe[[item]] <- c(x$pipe[[item]], list(value))
        }
        else
        {
            x$memory[[item]] <- c(x$memory[[item]], list(value))
        }
    }
    else
    {
        if (missing(item))
            x$pipe <- c(x$pipe, list(list(value)))
        else 
        {
            if ("numeric" == class(item))
                x$pipe[[item]] <- list(value)
            else
            {
                if (!is.null(x$memory[[item]]))
                    warning("overwriting item")
                x$memory[[item]] <- list(value)
            }
        }
    }

    return(x) 
}

#' Pulls an item from an object
#'
#' \code{pull} is a generic function and is used to pull an item from a stack
#' or pipe.
#'
#' There is no default function, so you have to see a specific \code{pull}
#' function for further information.
#'
#' @param x an object used to select a method.
#' 
#' @param ... further arguments passed to or from other methods.
#'
#' @return The updated object.
#' 
#' @export
pull <- function(x, ...) UseMethod("pull")

#' Pulls an item from a pubprint object 
#'
#' \code{pull.pubprint} is used to pull an item from the pipe or the named
#' memory of a pubprint object.
#'
#' No further details.
#'
#' @param x a pubprint object
#'
#' @param item the item to pull. If item is numeric, pipe and if it is a
#' character, named memory is chosen.
#'
#' @param remove either a logical, \code{"pipe"} or \code{"memory"}. If
#' \code{remove} is \code{TRUE}, every returned item is removed from pipe or
#' memory. If it is \code{"pipe"} (or \code{"memory"}), only accessed pipe (or
#' memory) items will be removed.
#' 
#' @param ... further arguments passed to \code{\link{pprint}} or the internal
#' style functions.
#'
#' @return The updated object.
#'
#' @seealso See \code{\link{pull}} for the generic function,
#' \code{\link{push<-}} to put items to pipe or named memory.
#'
#' @examples
#' ppo <- pubprint()
#' push(ppo) <- t.test(1:100, 2:101)
#' push(ppo, add = TRUE) <- .8123 # add d value to last pipe item
#' push(ppo, item = "i1") <- t.test(1:30, 2:31)
#'
#' pull(ppo)
#' pull(ppo, item = "i1")
#' pull(ppo, item = "i1", remove = TRUE) # removes item as well
#' 
#' @export
pull.pubprint <- function(x, 
                             item = 1, 
                             remove = pp_opts$get("removeItems"),
                             ...)
{
    objName <- deparse(substitute(x))

    if ("numeric" == class(item))
    {
        if (!length(x$pipe) || length(x$pipe) < item) 
            stop("subscript out of bounds")

        ret <- x$pipe[[item]]
        if ((is.logical(remove) && remove) || "pipe" == remove) 
            x$pipe <- x$pipe[-item]
    }    
    else
    {
        if (!length(x$memory) || ! item %in% names(x$memory)) 
            stop("item \"", item, "\" not available")

        ret <- x$memory[[item]]
        if ((is.logical(remove) && remove) || "memory" == remove) 
            x$memory <- x$memory[item != names(x$memory)]
    }

    ret <- pprint(ret, ...)

    assign(objName, x, envir = parent.frame())
    return(ret)
}


#' Prints a pubprint object
#' 
#' Prints the contents of a pubprint object
#'
#' Prints contents of named memory and pipe of a pubprint object. 
#'
#' @param x object of class \code{pubprint}.
#' 
#' @param ... further arguments. Ignored.
#' 
#' @examples
#' ppo <- pubprint()
#' push(ppo) <- t.test(1:10)
#' print(ppo)
#' 
#' @export
print.pubprint <- function(x, ...)
{
    cat("Values in unnamed register (pipe):\n")
    if (length(x$pipe)) 
        print(lapply(x$pipe, pprint, format = "object"))
    else 
        cat("empty\n")
    cat("\n")

    cat("Values in named register (memory):\n")
    if (length(x$memory)) 
        print(lapply(x$memory, pprint, format = "object"))
    else 
        cat("empty\n")
}


#' Prints an object in a publishable manner
#' 
#' \code{pprint} formats the output of the given object in a specified way
#'
#' This function calls internal style functions (depending on specified output
#' format) to convert the output of the object into the specified publication
#' style. It offers options to put a math mode and surrounding characters
#' around the (concatenated) output.
#'
#' If argument \code{format} is missing, a function tries to determine a
#' default format specifier. Can be specified to simple return the input
#' object (\code{"object"}). It is possible to set it to any internal style
#' function, the selected style supports.
#' 
#' @param x object which output should be printed. Can be a list to deliver
#' additional information to internal functions.
#' 
#' @param format optional format specifier. Character vector, see details.
#'
#' @param ... optional arguments passed to internal style functions. See their
#' help files for more information.
#'
#' @param concat logical, whether returned result is in a single character or
#' a character vector with parts of the statistical output.
#'
#' @param mmode logical indicating if the returned result should be set in
#' math mode (depends on output format).
#'
#' @param separator character string specifying the surrounding characters.
#'
#' @return Simply the unmodified object \code{x} in a list if \code{format} is
#' \code{"object"}, else a character vector.
#'
#' @seealso See \code{\link{pp_opts_style}} for setting publication style and
#' \code{\link{pp_opts_out}} for setting output format.
#' 
#' @examples
#' pprint(t.test(1:30))
#' pprint(t.test(1:30, 2:31))
#' pprint(t.test(1:30, 2:31), format = "object")
#' pprint(t.test(1:30, 2:31), mmode = FALSE, separator = NULL)
#' pprint(list(t.test(1:30), .843))
#'
#' @export
pprint <- function(x,
                   format,
                   ...,
                   concat = TRUE,
                   mmode = pp_opts$get("mmode"),
                   separator = pp_opts$get("separator"))
{
    if ("list" != class(x))
        x <- list(x)

    if (missing(format))
        format <- utils.get.format(x[[1]])

    # format == "object" the whole list will be returned
    if ("object" != format)
    {
        x <- pp_opts_style$get(format)(x, ...)

        if (concat)
            x <- out.concat(x)

        x <- out.math(x, mmode = mmode)

        if (!is.null(separator))
        {
            if (separator %in% c("brackets", "delimiter"))
                x <- out.bracket(x, brackets = pp_opts$get(separator))
            else
                x <- out.bracket(x, brackets = separator)
        }
    }

    return(x)
}
