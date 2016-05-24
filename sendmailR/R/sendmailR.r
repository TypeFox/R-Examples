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

.get_recipients <- function(headers) {
  res <- headers$To
  if (!is.null(headers$Cc)) {
    res <- c(res, headers$Cc)
  }
  if (!is.null(headers$Bcc)) {
    res <- c(res, headers$Bcc)
  }
  res
}

.write_mail <- function(headers, msg, sock) {
  if (!is.list(msg))
    msg <- list(msg)

  ## Generate MIME headers:
  boundary <- paste(packBits(sample(0:1, 256, TRUE)), collapse="")
  headers$`MIME-Version` <- "1.0"
  headers$`Content-Type` <- sprintf("multipart/mixed; boundary=\"%s\"", boundary)

  headers$To <- paste(headers$To, collapse=", ")
  if (!is.null(headers$Cc))
    headers$Cc <- paste(headers$Cc, collapse=", ")
  ## Do not include BCC recipients in headers, after all, it is a
  ## _blind_ carbon-copy.
  headers$Bcc <- NULL
  
  writeLines(paste(names(headers), unlist(headers), sep=": "),
             sock, sep="\r\n")
  writeLines("", sock, sep="\r\n")

  writeLines("This is a message with multiple parts in MIME format.", sock, sep="\r\n")

  for (part in msg) {
    writeLines(sprintf("--%s", boundary), sock, sep="\r\n")
    if (inherits(part, "mime_part"))
      .write_mime_part(part, sock)
    else if (is.character(part)) { ## Legacy support for plain old string
      ## writeLines(sprintf("--%s", boundary), sock, sep="\r\n")
      writeLines("Content-Type: text/plain; format=flowed\r\n", sock, sep="\r\n")
      writeLines(part, sock, sep="\r\n")
    }
  }
  writeLines(sprintf("--%s--", boundary), sock, sep="\r\n")
}

.smtp_submit_mail <- function(server, port, headers, msg, verbose=FALSE) {
  stopifnot(is.character(headers$From))
  
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
  recipients <- .get_recipients(headers)
  sapply(recipients, function(x) send_command(paste("RCPT TO: ", x), 250))
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
##' @param to Recipient of the message (vector of valid RFC2822 style addresses).
##' @param cc Carbon-copy recipients (vector of valid RFC2822 style addresses).
##' @param bcc Blind carbon-copy recipients (vector of valid RFC2822 style addresses).
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
sendmail <- function(from, to, subject, msg, cc, bcc, ...,
                     headers=list(),
                     control=list()) {
  ## Argument checks:
  stopifnot(is.list(headers), is.list(control))
  if (length(from) != 1)
    stop("'from' must be a single address.")
  
  if (length(to) < 1)
    stop("'to' must contain at least one address.")
  
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
  if (!missing(cc)) 
    headers$Cc <- cc
  if (!missing(bcc))
    headers$Bcc <- bcc
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
    message("Recipients: ", paste(.get_recipients(headers), collapse=", "))
    .write_mail(headers, msg, stderr())
  }
}
