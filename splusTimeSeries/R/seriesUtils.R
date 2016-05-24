par.uin <- function() 
  par("pin") / {u <- par("usr"); c(diff(u[1:2]), diff(u[3:4]))}

canonical.name <-
function(name, choices, ignore.case = TRUE,
         error.name = deparse(substitute(name))[1])
{
  name <- as.character(name)
  ch.names <- names(choices)
  choices <- as.character(choices)
  if(ignore.case)
    f <- casefold
  else f <- function(x)
    x
  i <- charmatch(f(name), f(choices))
  msg <- ""
  if(any(is.na(i)))
    msg <- paste("Value \"", name, "\" for ", error.name, 
                 " doesn't match any allowed choice. ", sep = "")
  else if(any(i == 0))
    msg <- paste("Ambiguous value \"", name, "\" for ", error.name,
                 ". ", sep = "")
  if(nchar(msg)) {
    opts <- paste(choices, collapse = ", ")
    msg <- paste(msg, "Choices are ", opts, ".", sep = "")
    call <- list(as.name("stop"), msg)
    mode(call) <- "call"
    eval(call, envir = sys.parent())
  }
  else {
    if(is.null(ch.names))
      choices[i]
    else {
      ch.names[ch.names == ""] <- choices[ch.names == ""]
      ch.names[i]
    }
  }
}

strip.blanks <-
function(str)
{
  ## strip leading and trailing spaces and tabs, but not interior ones
  strip.blanks.character <- function(str)
    {
      str <- sub("^[ \t]+", "", str)
      str <- sub("[ \t]+$", "", str)
      str
    }
  if(is.recursive(str))
    sapply(str, strip.blanks.character)
  else strip.blanks.character(str)
}
