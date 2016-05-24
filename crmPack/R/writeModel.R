#####################################################################################
## Author: Daniel Sabanes Bove [sabanesd *a*t* roche *.* com]
## Project: crmPack
##
## Time-stamp: <[writeModel.R] by DSB Mon 05/01/2015 17:38>
##
## Description:
## This is the write.model functionality from R2WinBUGS. We only need this,
## therefore we do not want to require the whole package, which in turn would
## require WinBUGS installation.
##
## History:
## 08/12/2014   file creation
#####################################################################################

##' Creating a WinBUGS model file
##'
##' Convert R function to a \pkg{WinBUGS} model file. BUGS models follow
##' closely S syntax. It is therefore possible to write most BUGS models as R
##' functions.
##' As a difference, BUGS syntax allows truncation specification like this:
##' \code{dnorm(...) I(...)}  but this is illegal in R. To overcome this
##' incompatibility, use dummy operator \code{\%_\%} before \code{I(...)}:
##' \code{dnorm(...) \%_\% I(...)}. The dummy operator \code{\%_\%} will be
##' removed before the BUGS code is saved.
##' In S-PLUS, a warning is generated when the model function is defined if the
##' last statement in the model is an assignment. To avoid this warning, add the
##' line \code{invisible()} to the end of the model definition. This line will be
##' removed before the BUGS code is saved.
##'
##' @param model R function containing the BUGS model in the BUGS
##' model language, for minor differences see Section Details.
##' @param con passed to \code{\link{writeLines}} which actually writes the
##' model file
##' @param digits number of significant digits used for \pkg{WinBUGS}
##' input, see \code{\link{formatC}}
##' @return Nothing, but as a side effect, the model file is written
##'
##' @export
##' @author original idea by Jouni Kerman, modified by Uwe Ligges, DSB removed S
##' Plus part
writeModel <- function(model, con = "model.bug", digits = 5)
{
    if (is.R()){
        model.text <- c("model", replaceScientificNotationR(body(model),
                                                            digits = digits))
                                        # "[\+\-]?\d*\.?[Ee]?[\+\-]?\d*"
    } else {
        ## In S-PLUS the source code of a function can be obtained with
        ## as.character(function_name).
        ## This omits the "function_name <- function()" piece
        model.text <- paste("model", as.character(model))
    }
    model.text <- gsub("%_%", "", model.text)
    writeLines(model.text, con = con)
}

##' @keywords internal
replaceScientificNotationR <- function(bmodel, digits = 5){
    env <- new.env()
    assign("rSNRidCounter", 0, envir=env)
    replaceID <- function(bmodel, env, digits = 5){
        for(i in seq_along(bmodel)){
            if(length(bmodel[[i]]) == 1){
                if(as.character(bmodel[[i]]) %in% c(":", "[", "[[")) return(bmodel)
                if((typeof(bmodel[[i]]) %in% c("double", "integer")) && ((abs(bmodel[[i]]) < 1e-3) || (abs(bmodel[[i]]) > 1e+4))){
                    counter <- get("rSNRidCounter", envir=env) + 1
                    assign("rSNRidCounter", counter, envir=env)
                    id <- paste("rSNRid", counter, sep="")
                    assign(id, formatC(bmodel[[i]], digits=digits, format="E"), envir=env)
                    bmodel[[i]] <- id
                }
            } else {
                bmodel[[i]] <- replaceID(bmodel[[i]], env, digits = digits)
            }
        }
        bmodel
    }
    bmodel <- deparse(replaceID(bmodel, env, digits = digits), control = NULL)
    for(i in ls(env)){
        bmodel <- gsub(paste('"', i, '"', sep=''), get(i, envir=env), bmodel, fixed=TRUE)
    }
    bmodel
}


