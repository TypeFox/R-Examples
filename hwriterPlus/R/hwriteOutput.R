hwriteOutput <- function(output, page = NULL, fontSize = "10pt", ...,
                         link = NULL, name = NULL, br = NULL, div = NULL) {

  ## Purpose: Addition to hwriter to allow easy printing of R output
  ## ----------------------------------------------------------------------
  ## Arguments: Same as for hwriteString, plus fontSize argument
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date: 24 Oct 2010, 12:23

  ## default arguments
  if (is.null(br)) br <- FALSE
  if (is.null(div)) div <- FALSE
  args <- list(...)

  ## Add line endings
  output <- paste(output, collapse = "\n")
  ### Add tags for preformatting
  output <- paste("<pre style = 'font-size:", fontSize, "'>",
                  output, "</pre>", sep = "")

  hwrite(output, page = page, ..., link = link, name = name,
         br = br, div = div)

}

hwriteScript <- function(scriptFile, page = NULL, fontSize = "10pt",
                         trim = TRUE, ...,
                         link = NULL, name = NULL, br = NULL, div = NULL) {

  ## Purpose: Addition to hwriter to allow easy printing of R session
  ## ----------------------------------------------------------------------
  ## Arguments: Similar to hwriteString and hwriteOutput
  ##            scriptFile is file containing text created by script
  ##            trim trims 1 and last 2 lines of text created by script
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date: 22 August, 2012

  ## Default arguments
  if (is.null(br)) br <- FALSE
  if (is.null(div)) div <- FALSE
  args <- list(...)

  ## Read script output from file
  session <- readLines(scriptFile)

  ## Trim off first and last 2 lines of text created by script
  if (trim) {
      nSession <- length(session)
      session <- session[-c(1,nSession - 1,nSession)]
  }

  ## Add line endings
  session <- paste(session, collapse = "\n")
  ### Add tags for preformatting
  session <- paste("<pre style = 'font-size:", fontSize, "'>",
                  session, "</pre>", sep = "")

  hwrite(session, page = page, ..., link = link, name = name,
         br = br, div = div)

}
