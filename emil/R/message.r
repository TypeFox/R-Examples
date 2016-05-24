#' Print a timestamped and indented log message
#'
#' To suppress messages below a given indentation level set the global
#' \code{\link{option}} setting \code{emil_max_indent}, as in the example below.
#'
#' @param indent Indentation level. Messages with \code{indent=0} are
#'   suppressed.
#' @param ... Sent to \code{\link{sprintf}}.
#' @param time Whether or not to print timestamp.
#' @param domain See \code{\link{message}}.
#' @param appendLF Whether to finish the message with a linebreak or not.
#' @examples
#' equipment <- c("flashlight", "snacks", "pick")
#' {
#'     log_message(1, "Begin descent")
#'     log_message(2, "Oh no, forgot the %s!", sample(equipment, 1))
#'     log_message(2, "Hello? Can you throw it down to me?", time=FALSE)
#'     log_message(1, "Aw shucks, I'm coming back up.")
#' }
#'
#' for(verbose in c(TRUE, FALSE)){
#'     cat("It's", verbose, "\n")
#'     for(i in 0:3)
#'         log_message(indent(verbose, i), "Down")
#' }
#'
#' options(emil_max_indent = 2)
#' for(i in 1:3)
#'     log_message(i, "Down")
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
log_message <- function(indent=1, ..., time=TRUE, domain="R-emil", appendLF=TRUE){
    if(indent > 0){
        for(msg in sprintf(...)){
        message(
            # timestamp
            if(time) format(Sys.time(), "%d %b %H:%M") else
                     paste(rep(" ", 12), collapse=""),
            # indent
            rep("  ", indent),
            # message
            msg,
            # linebreak
            domain = domain,
            appendLF = appendLF)
        }
    }
    if(indent > 4 && is.null(getOption("emil_max_indent"))){
        notify_once("deep_indent", "Use `options(emil_max_indent = 4)` to hide excessive logging.",
                    fun=warning)
    }
}

#' Increase indentation
#'
#' @param base Base indentation level of the function printing the message.
#' @param indent Extra indentation of this message.
#' @return An integer that can be used to specify the indentation level of
#'   messages printed with \code{\link{log_message}}.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
indent <- function(base, indent){
    if(base > 0 && base + indent <= getOption("emil_max_indent", Inf)){
        base + indent
    } else {
        0
    }
}

#' Print a warning message if not printed earlier
#' 
#' To avoid flooding the user with identical warning messages, this function
#' keeps track of which have already been shown.
#' 
#' @param id Warning message id. This is used internally to refer to the
#'   message.
#' @param ... Sent to \code{\link{warning}}.
#' @param if_top_level If \code{TRUE} the notifications will only be reset if
#'   \code{reset_notification} was called from a top-level function call.
#'   This behaviour prevents the notifications from being reset multiple times
#'   during nested calls to functions such as \code{\link{fit}} and
#'   \code{\link{evaluate}}.
#' @param fun Function to display the notification with. Typical choices are
#'  \code{\link{message}} or \code{\link{warning}}.
#' @author Christofer \enc{Bäcklin}{Backlin}
#' @export
notify_once <- function(id, ..., fun=log_message){
    if(!id %in% getOption("emil_notification")){
        match.fun(fun)(...)
        options(emil_notification = c(getOption("emil_notification"), id))
    }
}
#' @rdname notify_once
#' @export
reset_notification <- function(id, if_top_level=TRUE){
    if(!if_top_level || identical(parent.frame(), globalenv())){
        if(missing(id) || is.null(id)){
            options(emil_notification = NULL)
        } else {
            options(emil_notification = setdiff(getOption("emil_notification"), id))
        }
    }
}

