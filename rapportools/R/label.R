#' Variable Name
#'
#' This function returns character value previously stored in variable's \code{name} attribute. If none found, the function defaults to object's name.
#' @param x an R (atomic or data.frame/list) object to extract names from
#' @return a character value with variable's label
#' @examples \dontrun{
#' name(mtcars$am)
#' x <- 1:10
#' name(x)
#' }
#' @export
#' @aliases name rp.name
name <- function(x){

    if (base::missing(x))
        stop('variable not provided')

    if (is.atomic(x)){
        n <- attr(x, which = 'name', exact = TRUE)
        if (is.null(n)) {
            return (tail(as.character(substitute(x)), 1)) # return variable name if no label
        } else {
            if (length(n) > 1)
                warning('variable name is not a length-one vector, only the first element is displayed')
            return(attr(x, 'name'))                       # return variable label
        }
    }

    if (is.recursive(x)){
        n <- sapply(x, attr, which = 'name', exact = TRUE)
        n.nil <- sapply(n, is.null)

        ## no labels found
        if (all(n.nil)){
            n <- names(n)
        } else
            n[n.nil] <- names(n)[n.nil]

        return(n)
    }

    stop('Wrong R object type provided!')
}
#' @export
rp.name <- name


#' Get Variable Label
#'
#' This function returns character value previously stored in variable's \code{label} attribute. If none found, and \code{fallback} argument is set to \code{TRUE} (default), the function returns object's name (retrieved by \code{deparse(substitute(x))}), otherwise \code{NA} is returned with a warning notice.
#' @param x an R object to extract labels from
#' @param fallback a logical value indicating if labels should fallback to object name(s)
#' @param simplify coerce results to a vector (\code{TRUE} by default), otherwise, a \code{list} is returned
#' @return a character vector with variable's label(s)
#' @examples \dontrun{
#' x <- rnorm(100)
#' label(x)         # returns "x"
#' label(x, FALSE)  # returns NA and issues a warning
#'
#' label(mtcars$hp) <- "Horsepower"
#' label(mtcars)         # returns "Horsepower" instead of "hp"
#' label(mtcars, FALSE)  # returns NA where no labels are found
#' label(sleep, FALSE)   # returns NA for each variable and issues a warning
#' }
#' @export
#' @aliases label rp.label
label <- function(x, fallback = TRUE, simplify = TRUE){

    if (base::missing(x))
        stop('Variable not provided.')

    if (is.null(x))
        return (NULL)

    if (is.atomic(x)){
        lbl <- attr(x, which = 'label', exact = TRUE)
        if (is.null(lbl)){
            if (fallback){
                lbl <- tail(as.character(substitute(x)), 1)
            } else {
                warning('Atomic object has no labels.')
                lbl <- NA
            }
        } else {
            if (length(lbl) > 1){
                warning('Variable label is not a length-one vector, only first element is returned.')
                lbl <- lbl[1]
            }
        }
    } else {
        lbl <- sapply(x, attr, which = 'label', exact = TRUE)
        lbl.nil <- sapply(lbl, is.null)

        ## no labels found
        if (all(lbl.nil)){
            if (fallback){
                lbl <- structure(names(lbl), .Names = names(lbl))
            } else {
                warning('No labels found in recursive object,')
                lbl[lbl.nil] <- NA
            }
        } else {
            if (fallback)
                lbl[lbl.nil] <- names(lbl)[lbl.nil]
            else
                lbl[lbl.nil] <- NA
        }
    }

    if (simplify)
        lbl <- unlist(lbl)

    return(lbl)
}
#' @export
rp.label <- label


#' Set Variable Label
#'
#' This function sets a label to a variable, by storing a character string to its \code{label} attribute.
#' @param var a variable (see \code{\link{is.variable}} for details)
#' @param value a character value that is to be set as variable label
#' @usage label(var) <- value
#' @seealso \code{\link{label}}
#' @examples \dontrun{
#' label(mtcars$mpg) <- "fuel consumption"
#' x <- rnorm(100)
#' (label(x) <- "pseudo-random normal variable")
#' }
#' @export
#' @aliases label<- rp.label<-
`label<-` <- function(var, value){

    if (base::missing(var) | base::missing(value))
        stop('Both variable name and label should be provided.')

    if (!is.variable(var))
        stop('Label can only be assigned to a variable.')

    if (!is.string(value))
        stop('Only a character string can be assigned to a variable label.')

    attr(var, 'label') <- value

    return (var)
}
#' @export
`rp.label<-` <- `label<-`
