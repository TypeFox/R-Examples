## Option managment shamelessly taken from the lattice package.
.SendmailREnv <- new.env(parent=emptyenv())
.SendmailREnv$options <- list()

.update_list <- function (x, val) {
  if (is.null(x)) 
    x <- list()
  modifyList(x, val)
}

##' Specify global sendmail options so that subsequent calls to
##' \code{sendmail()} do not have to set them in the \code{control}
##' argument.  
##'
##' List of options:
##' \itemize{
##' \item{smtpServer}{SMTP server to contact. This can either be the
##'   mail server responsible for the destination addresses domain or a
##'   smarthost provided by your ISP or institution. SMTP AUTH is
##'   currently unsupported.}
##' \item{smtpPort}{SMTP port to use. Usually 25 but some institutions
##'   require the use of the submission service (port 587).}
##' \item{verbose}{Show detailed information about message
##'   submission. Useful for debugging.}
##' }
##'
##' @param ... Any options can be defined, using \code{name=value} or
##' by passing a list of such tagged values.  However, only the ones
##' below are used in base sendmailR.
##' @return For \code{sendmail_options()}, a list of all set options
##' sorted by name. For \code{sendmail_options(name)}, a list of length
##' one containing the set value, or 'NULL' if it is unset.  For uses
##' setting one or more options, a list with the previous values of
##' the options changed (returned invisibly).  
##' 
##' @title Set package specific options.
##' @export
##' @author Olaf Mersmann \email{olafm@@datensplitter.net}
sendmail_options <- function(...) {
  new <- list(...)
  if (is.null(names(new)) && length(new) == 1 && is.list(new[[1]])) 
    new <- new[[1]]
  old <- .SendmailREnv$options
  if (length(new) == 0) 
    return(old)
  nm <- names(new)
  if (is.null(nm)) 
    return(old[unlist(new)])
  isNamed <- nm != ""
  if (any(!isNamed)) 
    nm[!isNamed] <- unlist(new[!isNamed])
  retVal <- old[nm]
  names(retVal) <- nm
  nm <- nm[isNamed]
  .SendmailREnv$options <- .update_list(old, new[nm])
  invisible(retVal)
}

##' @export
##' @rdname sendmail_options
sendmailOptions <- function(...) {
  .Deprecated("sendmail_options")
  sendmail_options(...)
}
