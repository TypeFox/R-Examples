#' Internal functions for EEM package
#' 
#' Internal functions for EEM package
#' @name EEM-misc
#' @aliases generatePoint
#' @aliases generateColor
#' 
#' @param n number
#' @param pch Either an integer specifying a symbol or a single character to be used as the default in plotting points.
#' @inheritParams drawEEM
#' @param string string or vector of strings
#' @inheritParams base:::round
#' 
#' @details `generatePoint` and `generateColor` are used to create point and color 
#' vector from specified number (n) and palette. 
#' 

#' @export
generatePoint <- function(n, pch = NULL){
    if (is.null(pch)){
        pch <- c(16, 17, 15, 1, 2, 4, 6, 8)
    }
    times <- (n %/% length(pch)) +1
    pointType <- rep(pch, times)[1:n]
    return(pointType)
}

#' @describeIn EEM-misc generate colors
#' @export
generateColor <- function(n, color.palette = NULL){
    if (is.null(color.palette)) {
        color.palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "black", 
                      "#A65628", "#F781BF", "#999999", "#1B9E77", "#D95F02", "#7570B3", 
                      "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
    }
    times <- (n %/% length(color.palette)) +1
    colorType <- rep(color.palette, times)[1:n]
    return(colorType)
}

#' @describeIn EEM-misc get EX value
#' @export
getEX <- function(string, digits = NULL){
    pat <- "(EX|.*EX)(.*)(EM)(.*)"
    ex <- as.numeric(sub(pat, "\\2", string))
    if (!is.null(digits)){
        ex <- round(ex, digits)
    }
    return(ex)
}

#' @describeIn EEM-misc get EM value
#' @export
getEM <- function(string, digits = NULL){
    pat <- "(EX|.*EX)(.*)(EM)(.*)"
    ex <- as.numeric(sub(pat, "\\4", string))
    if (!is.null(digits)){
        ex <- round(ex, digits)
    }
    return(ex)
}

#' @export
c.EEM <- function(..., recursive = FALSE)  {
    dots <- list(...)
    res <- structure(unlist(dots, recursive = FALSE), class = "EEM")
    return (res)
}