
## ->  sigma  generic function
#' Extract scale parameter sigma from a model fit
#'
#' This function is generic; method functions can be written to handle specific classes of objects.
#'
#' @param object an object for which scale parameter can be extracted.
#' @param \dots some methods for this generic function may require additional arguments.
#' @return numeric scalar value.
#' @author Andrzej Galecki and Tomasz Burzykowski
#' @examples
#'  ## sigma (fm1)
#' @export

sigma <-  function(object, ...) UseMethod("sigma")


##
#' @S3method sigma default
sigma.default <- function(object, ...) object$sigma

## -> missPat function
#'
#' Extract pattern of missing data
#'
#' This function allows to compactly present pattern of missing data in a given vector/matrix/data frame or combination of thereof.
#'
#' @export
#' @param \dots one or more vectors/matrices/data frames. They need to be compatible for columnwise binding.
#' @param symbols vector containing two single characters used to indicate NA and remaining values. By defualt it has values: \code{X} and \code{-}.
#' @param collapse an optional character string. It is used in the internal call to \code{paste()} function to separate the results. Rarely used. By default set to NULL.
#' @param missData logical. If \code{TRUE} data frame with pattern of missing values is saved in \code{missData} attribute of the vector returned by this function.
#' @return character vector with as many elements as length of vectors(s)/number of rows in matrices and/or data frames in \code{\dots{}} argument(s).
#'   Attribute \code{cnames} contains names of vectors/columns/variables. 
#'   Optional attribute \code{missData} contains data frame with missing pattern.
#' @author Andrzej Galecki and Tomasz Burzykowski
#' @examples
#'
#' dtf <- subset(armd.wide, 
#'              select = c(visual12, visual24, visual52))
#' missPat(dtf, symbols = c("?","+"))
#'
missPat <- function(..., symbols = c("X","-"), collapse = "", missData = FALSE){
     .functionName <- "missPat"                                # Function name (recommended)
     .traceR <- if (is.null(options()$traceR)) function(...){} else options()$.traceR       

     .traceR(1, "-> missPat STARTS")
     args  <- as.list(substitute(list(...)))[-1]
     argsL <- lapply(args, eval)
     dt <- data.frame(argsL)
     nms <- lapply(args, FUN= function(el){
       elx <- eval(el)
       if (is.null(colnames(elx))) {
        nc <- ncol(elx)
        if (is.null(nc)) as.character(el) else paste(el, 1:nc, sep = ":")
     }   
         else colnames(elx)
     })
     cx1 <- symbols[1]
     cx2 <- symbols[2]
     miss.frame <- as.data.frame(ifelse(is.na(dt), cx1, cx2))
     names(miss.frame) <- unlist(nms)
     res <- apply(miss.frame, 1, paste, collapse = collapse)
     attr(res, "cnames") <- unlist(nms)
     if (missData) attr(res, "missData") <- miss.frame
     .traceR(1, "missPat ENDS <-")
     res
}

## -> runScript function
#' Executes scripts from GB book
#'
#' Default call of the function without arguments, prints a list of available scripts.
#'
#' @param script character string containing name of the script to be executed. By default is set to NA.
#' @param package character string containing package name. By default nlmeU.
#' @param subdir subdirectory containing scripts. By default: scriptsR2.15.0.
#' @param echo logical. Used by source function. By default set to TRUE.
#' @return Script is executed and results are printed.
#' @author Andrzej Galecki and Tomasz Burzykowski
#' @export
#' @examples runScript()
#'

runScript <- function(script= NA,  package = "nlmeU", subdir = "scriptsR2.15.0", 
    echo = TRUE){
    scriptsDir <- system.file(subdir, package = package)
    scriptsList <- list.files(scriptsDir, pattern = "[[:alnum:]][.][R]$")
    scriptFile <- file.path(scriptsDir, script)
    if (!(script %in% scriptsList)) {
    if (is.na(script)) {
            errFun <- message
            errMsg <- ""
        }
        else {
            errFun <- stop
            errMsg <- paste("Example", example, "does not exist. ")
        }
        errFun(errMsg, "Scripts in ", scriptsDir, " are: \n", paste("\"",scriptsList, 
            collapse = "\", \n", sep=""), "\"")
        if (subdir == "scriptsR2.15.0")
            cat ('\n Scripts employing lme4.0 package are stored in: \n',
            file.path(scriptsDir,"lme4.0"), ' directory',
            ' and can be found by issuing:\n   runScript(subdir = "scriptsR2.15.0/lme4.0") command \n')         
    }
    else {
        sourceText <- source(scriptFile, echo=echo)
        sourceText
    }
}



