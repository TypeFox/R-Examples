##---------------------------------------------------
## Modified version sendmailR
##
## sendmailR: http://cran.r-project.org/web/packages/sendmailR/index.html
## Modified date: 2013.08.12
##
##
##---------------------------------------------------

.mime_part_finalizer <- function(x) {
  if (!is.null(x$file))
    file.remove(x$file)
}

.mime_part <- function(headers, file=NULL, text=NULL) {
  if (!is.null(file) && !is.null(text))
    stop("Can only provide file or text for mime part.")

  e <- environment()
  reg.finalizer(e, .mime_part_finalizer, onexit=TRUE)
  class(e) <- "mime_part"
  e
}

.write_mime_part <- function(mp, con=stdout()) {
  writeLines(paste(names(mp$headers), unlist(mp$headers), sep=": "),
             con, sep="\r\n")
  writeLines("", con, sep="\r\n")
  if (is.null(mp$file))
    writeLines(mp$text, con)
  else
    writeLines(readLines(mp$file), con, sep="\r\n")
}

.file_attachment <- function(fn, name,
                             type="application/octet-stream",
                             disposition="attachment") {
  if (missing(name))
    name <- basename(fn)

  text <- base64encode(fn, linewidth=72, newline="\n")
  headers <- list("Content-ID"=sprintf("<%s>",name),
                  "Content-Type"=type,
                  "Content-Disposition"=sprintf("%s; filename=%s",
                    disposition, name),
                  "Content-Transfer-Encoding"="base64")
  
  .mime_part(headers=headers, text=text)
}

.plot_attachment <- function(plt, name=deparse(substitute(plt)), device, ...) {
  fn <- tempfile()
  device(file=fn, ...)
  print(plt)
  dev.off()
  ## FIXME: Guess content type from device!
  res <- .file_attachment(fn, name, type="application/pdf")
  file.remove(fn)
  res
}

##' Create a MIME part
##'
##' @param x Object to include
##' @param name Name of mime part. Usually the filename of the
##'   attachment as displayed by the e-mail client.
##' @param ... Possible further arguments for \code{mime_part}
##'   implementations.
##' @return An S3 \code{mime_part} object.
##' @export
mime_part <- function(x, name, ...)
  UseMethod("mime_part", x)

##' Default MIME part method
##'
##' Creates a string representation of the object \code{x} using
##' \code{dput}. This representation is then turned into a file
##' attachment.
##'
##' @param x R object
##' @param name Filename used for attachment (sans the .R extension)
##' @param ... Ignored.
##' @return An S3 \code{mime_part} object.
##'
##' @method mime_part default
##' @S3method mime_part default
mime_part.default <- function(x, name, ...) {
  str <- dput(x)
  .mime_part(headers=list(
               "Content-Type"="text/plain",
               "Content-Disposition"=sprintf("attachment; file=%s.R", name)),
             text=str)
}

##' Creates a MIME part from a trellis plot object
##'
##' Writes a PDF file of the plot defined by \code{x} and turns this
##' PDF file into a file attachment.
##'
##' @param x A \code{trellis} (lattice) object
##' @param name Name of attachment (sans .pdf extension).
##' @param device Graphics device used to render the plot. Defaults to
##'   \code{pdf}.
##' @param ... Ignored.
##' @return An S3 \code{mime_part} object.
##' 
##' @method mime_part trellis
##' @S3method mime_part trellis
mime_part.trellis <- function(x, name=deparse(substitute(x)), device=pdf, ...)
  .plot_attachment(x, name=name, device=device, ...)

##' Creates a MIME part from a ggplot2 plot object
##'
##' Writes a PDF file of the plot defined by \code{x} and turns this
##' PDF file into a file attachment.
##'
##' @param x A \code{ggplot} object
##' @param name Name of attachment (sans .pdf extension).
##' @param device Graphics device used to render the plot. Defaults to
##'   \code{pdf}.
##' @param ... Ignored.
##' @return An S3 \code{mime_part} object.
##' 
##' @method mime_part ggplot
##' @S3method mime_part ggplot
mime_part.ggplot <- function(x, name=deparse(substitute(x)), device=pdf, ...)
  .plot_attachment(x, name=name, device=device, ...)

##' Create a MIME part from a matrix.
##'
##' @param x Matrix
##' @param name Basename of file attachment that is generated.
##' @param ... Ignored.
##' @return An S3 \code{mime_part} object
##' 
##' @method mime_part matrix
##' @S3method mime_part matrix
mime_part.matrix <- function(x, name=deparse(substitute(x)), ...) {
  f <- tempfile()
  on.exit(file.remove(f))
  write.table(x, file=f, ...)
  .file_attachment(f, name=sprintf("%s.txt", name), type="text/plain")
}

##' Create a MIME part from a \code{data.frame}.
##' 
##' @param x A \code{data.frame}.
##' @param name Basename of file attachment that is generated.
##' @param ... Ignored.
##' @return An S3 \code{mime_part} object.
##'
##' @method mime_part data.frame
##' @S3method mime_part data.frame
mime_part.data.frame <- function(x, name=deparse(substitute(x)), ...) {
  f <- tempfile()
  on.exit(file.remove(f))
  write.table(x, file=f, ...)
  .file_attachment(f, name=sprintf("%s.txt", name), type="text/plain")
}

##' Create a MIME part from a character string. If the string matches
##' a filename, a MIME part containing that file is returned instead.
##' 
##' @param x Character string, possibly a filename.
##' @param name Name of attachment.
##' @param ... Ignored.
##' @return An S3 \code{mime_part} object.
##' 
##' @method mime_part character
##' @S3method mime_part character
mime_part.character <- function(x, name, ...) {
  if (length(x) == 1 && file.exists(x)) {
    .file_attachment(x, name, ...)
  } else {
     .mime_part(headers=list(
                 "Content-Type"="text/plain",
                 "Content-Disposition"="inline"),
               text=paste(x, collapse="\r\n"))
  }
}


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


##
## sendmailR.r - send email from within R
##
## Author:
##  Olaf Mersmann (OME) <olafm@datensplitter.net>
##

.rfc2822_date <- function(time=Sys.time()) {
  lc <- Sys.getlocale("LC_TIME")
  on.exit(Sys.setlocale("LC_TIME", lc))
  Sys.setlocale("LC_TIME", "C")
  strftime(time, format="%a, %d %b %Y %H:%M:%S -0000",
           tz="UTC", use.tz=TRUE)
}

.write_mail <- function(headers, msg, sock) {
  if (!is.list(msg))
    msg <- list(msg)
  content.type <- ifelse("Content-Type" %in% names(headers), headers$`Content-Type`, "text/plain; format=flowed")

  ## Generate MIME headers:
  boundary <- paste(packBits(sample(0:1, 256, TRUE)), collapse="")
  headers$`MIME-Version` <- "1.0"
  headers$`Content-Type` <- sprintf("multipart/mixed; boundary=\"%s\"", boundary)

  writeLines(paste(names(headers),
                   unlist(headers), sep=": "),
             sock, sep="\r\n")
  writeLines("", sock, sep="\r\n")

  writeLines("This is a message with multiple parts in MIME format.", sock, sep="\r\n")

  for (part in msg) {
    writeLines(sprintf("--%s", boundary), sock, sep="\r\n")
    if (inherits(part, "mime_part"))
      .write_mime_part(part, sock)
    else if (is.character(part)) { ## Legacy support for plain old string
      ## writeLines(sprintf("--%s", boundary), sock, sep="\r\n")
      writeLines(paste("Content-Type: ",content.type,"\r\n",sep=""), sock, sep="\r\n")
      writeLines(part, sock, sep="\r\n")
    }
  }
  writeLines(sprintf("--%s--", boundary), sock, sep="\r\n")
}

.smtp_submit_mail <- function(server, port, headers, msg, verbose=FALSE) {
  stopifnot(is.character(headers$From), is.character(headers$To))
  
  wait_for <- function(lcode) {
    done <- FALSE
    while (!done) {
      line <- readLines(con=sock, n=1)
      if (verbose)
        message("<< ", line)
      code <- substring(line, 1, 3)
      msg <- substring(line, 5)
      if (code == lcode) {
        done <- TRUE
      } else {
        if (code >= 500 & code <= 599)
          stop("SMTP Error: ", msg)
        else
          message("Unknown SMTP code: ", code)
      }
      
    }
    return(list(code=code, msg=msg))
  }

  send_command <- function(cmd, code) {
    if (verbose)
      message(">> ", cmd)
    writeLines(cmd, sock, sep="\r\n")
    wait_for(code)
  }

  nodename <- Sys.info()[4]
  sock <- socketConnection(host=server,
                           port=port,
                           blocking=TRUE)
  if (!isOpen(sock))
    stop(sprintf("Could not connect to smtp server '%s' on port '%i'.",
                 server, port))
  on.exit(close(sock))
  ## << 220 <hostname> ESMTP
  wait_for(220)
  ## >> HELO localhost
  ## << 250 mail.statistik.uni-dortmund.de
  send_command(paste("HELO ", nodename), 250)
  ## >> MAIL FROM: <foo@bah.com>
  ## << 250 2.1.0 Ok
  send_command(paste("MAIL FROM: ", headers$From), 250)
  ## >> RCPT TO: <bah@baz.org>
  ## << 250 2.1.5 Ok
  to.list <- unlist(strsplit(headers$To,','))
  for (to in to.list){
      send_command(paste("RCPT TO: ", to), 250)
  }
  ## >> DATA
  ## << 354 blah fu
  send_command("DATA", 354)
  ## >> <actual message + headers + .>
  if (verbose)
    message(">> <message data>")

  .write_mail(headers, msg, sock)
  
  writeLines(".", sock, sep="\r\n")
  
  wait_for(250)
  ## << 250 2.0.0 Ok: queued as XXXXXXXX
  ## >> QUIT
  ## << 221 2.0.0 Bye
  send_command("QUIT", 221)
}

##' Simplistic sendmail utility for R. Uses SMTP to submit a message
##' to a local SMTP server.
##'
##' @title Send mail from within R
##'
##' @param from From whom the mail message is (RFC2822 style address).
##' @param to Recipient of the message (valid RFC2822 style address).
##' @param subject Subject line of message.
##' @param msg Body text of message or a list containing
##'   \code{\link{mime_part}} objects.
##' @param \dots ...
##' @param headers Any other headers to include.
##' @param control List of SMTP server settings. Valid values are the
##'   possible options for \code{\link{sendmail_options}}.
##'
##' @seealso \code{\link{mime_part}} for a way to add attachments.
##' @keywords utilities
##' 
##' @examples
##' \dontrun{
##' from <- sprintf("<sendmailR@@\\%s>", Sys.info()[4])
##' to <- "<olafm@@datensplitter.net>"
##' subject <- "Hello from R"
##' body <- list("It works!", mime_part(iris))
##' sendmail(from, to, subject, body,
##'          control=list(smtpServer="ASPMX.L.GOOGLE.COM"))
##' }
##'
##' @export
sendmail <- function(from, to, subject, msg, ...,
                     headers=list(),
                     control=list()) {
  ## Argument checks:
  stopifnot(is.list(headers), is.list(control))
  if (length(from) != 1)
    stop("'from' must be a single address.")
  
  if (length(to) != 1)
    stop("'to' must be a single address.")
  
  get_value <- function(n, default="") {
    if (n %in% names(control)) {
      return(control[[n]])
    } else if (n %in% names(.SendmailREnv$options)) {
      return(.SendmailREnv$options[[n]])
    } else {
      return(default)      
    }
  }
  
  headers$From <- from
  headers$To <- to
  headers$Subject <- subject

  ## Add Date header if not explicitly set. This fixes the annoyance,
  ## that apparently Thunderbird does not sort mails correctly if they
  ## do not have a Date header.
  if (is.null(headers$Date))
    headers$Date <- .rfc2822_date()
  
  transport <- get_value("transport", "smtp")
  verbose <- get_value("verbose", FALSE)
  if (transport == "smtp") {
    server <- get_value("smtpServer", "localhost")
    port <- get_value("smtpPort", 25)
    
    .smtp_submit_mail(server, port, headers, msg, verbose)
  } else if (transport == "debug") {
    .write_mail(headers, msg, stdout())
  }
}
