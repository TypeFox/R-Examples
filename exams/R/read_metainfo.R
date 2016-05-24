extract_environment <- function(x, env, value = TRUE, markup = c("latex", "markdown"))
{
  markup <- match.arg(markup)
  if(markup == "latex") {
    b <- grep(paste0("\\\\(begin|end)\\{", env, "\\}"), x)
    if(length(b) == 0L) return(NULL)
    if(length(b)!= 2L) stop("no unique begin/end pair for", sQuote(env), "found")
    if(value) return(x[(b[1L] + 1L):(b[2L] - 1L)]) else return(b)
  } else {
    ## get all sections and subsections
    seclines <- grep("^====", x)
    sublines <- grep("^----", x)
    alllines <- sort(c(seclines, sublines))
    ## match environment names
    x[alllines - 1L] <- tolower(x[alllines - 1L])
    x[alllines - 1L] <- gsub("-", "", x[alllines - 1L], fixed = TRUE)
    x[alllines - 1L] <- gsub("questionlist", "answerlist", x[alllines - 1L], fixed = TRUE)
    x[alllines - 1L] <- gsub("solutionlist", "answerlist", x[alllines - 1L], fixed = TRUE)
    ## find desired environment
    wi <- which(env == x[alllines - 1L])
    if(length(wi) < 1L) return(NULL)

    ## begin/end
    b <- alllines[wi] - 1L
    e <- if(substr(x[b + 1L], 1L, 1L) == "=") {
      seclines[seclines > b + 1L]
    } else {
      alllines[alllines > b + 1L]
    }
    e <- if(length(e) > 0L) min(e) - 3L else length(x)
    if(value) return(x[(b + 2L):e]) else return(c(b, e))
  }
}

extract_command <- function(x, command, type = c("character", "logical", "numeric"), markup = c("latex", "markdown"))
{
  ## return type and markup type
  type <- match.arg(type)
  markup <- match.arg(markup)

  ## find command line
  command <- if(markup == "latex") paste0("\\", command) else paste0(command, ":")
  rval <- x[grep(command, x, fixed = TRUE)]
  if(length(rval) < 1L) {
      if(type == "logical") return(FALSE) else return(NULL)
  }
  if(length(rval) > 1L) {
    warning("command", sQuote(command), "occurs more than once, last instance used")
    rval <- tail(rval, 1L)
  }
  
  if(markup == "latex") {
    ## strip off everything in brackets
    ## omit everything before \command{
    rval <- strsplit(rval, paste(command, "{", sep = ""), fixed = TRUE)[[1L]][2L]
    ## omit everything after last }
    rval <- gsub("[^\\}]+$", "", rval)
    ## get everthing within brackets
    rval <- gsub("{", "", strsplit(rval, "}")[[1L]], fixed = TRUE)
    ## split further with respect to other symbols (currently only |)
    rval <- unlist(strsplit(rval, "|", fixed = TRUE))
  } else {
    ## strip off command
    rval <- gsub(command, "", rval, fixed = TRUE)
    ## omit leading and trailing white space
    rval <- gsub("^[ \t]+", "", rval)
    rval <- gsub("[ \t]+$", "", rval)
    ## split further with respect to other symbols (currently only |)
    rval <- unlist(strsplit(rval, "|", fixed = TRUE))
  }

  ## convert to return type
  do.call(paste("as", type, sep = "."), list(rval))
}

extract_extra <- function(x, markup = c("latex", "markdown"))
{
  ## markup type
  markup <- match.arg(markup)

  ## search for extra commands
  comm0 <- if(markup == "latex") "\\exextra[" else "exextra["
  comm <- x[grep(comm0, x, fixed = TRUE)]  
  if(length(comm) < 1L) return(list())
  
  ## extract command and type
  comm <- sapply(strsplit(comm, "comm0", fixed = TRUE), "[", 2L)
  comm <- sapply(strsplit(comm, "]", fixed = TRUE), "[", 1L)
  nam <- strsplit(comm, ",", fixed = TRUE)
  typ <- sapply(nam, function(z) if(length(z) > 1L) z[2L] else "character")
  nam <- sapply(nam, "[", 1L)
  
  ## call extract_command
  rval <- lapply(seq_along(comm), function(i) extract_command(x,
    command = paste0("exextra[", comm[i], "]"), type = typ[i], markup = markup))

  names(rval) <- nam
  return(rval)
}

extract_items <- function(x, markup = c("latex", "markdown"))
{
  ## markup type
  markup <- match.arg(markup)

  ## map markdown to tex
  if(markup == "markdown") x <- gsub("^\\* ", "\\\\item ", x)
    
  ## make sure we get items on multiple lines right
  x <- paste(x, collapse = " ")
  x <- gsub("^ *\\\\item *", "", x)
  x <- strsplit(x, " *\\\\item")[[1L]]
  x <- gsub("^ ", "", x)
  gsub(" +$", "", x)
}

read_metainfo <- function(file)
{
  ## read file
  x <- readLines(file)
  markup <- switch(tools::file_ext(file),
    "tex" = "latex",
    "md" = "markdown"
  )

  ## Description ###################################
  extype <- match.arg(extract_command(x, "extype", markup = markup), ## exercise type: schoice, mchoice, num, string, or cloze
    c("schoice", "mchoice", "num", "string", "cloze"))  
  exname <- extract_command(x, "exname", markup = markup)            ## short name/description, only to be used for printing within R
  extitle <- extract_command(x, "extitle", markup = markup)          ## pretty longer title
  exsection <- extract_command(x, "exsection", markup = markup)      ## sections for groups of exercises, use slashes for subsections (like URL)
  exversion <- extract_command(x, "exversion", markup = markup)      ## version of exercise

  ## Question & Solution ###########################
  exsolution <- extract_command(x, "exsolution", markup = markup)    ## solution, valid values depend on extype
  extol <- extract_command(x, "extol", "numeric", markup = markup)   ## optional tolerance limit for numeric solutions
  exclozetype <- extract_command(x, "exclozetype", markup = markup)  ## type of individual cloze solutions

  ## E-Learning & Exam ###################################
  expoints  <- extract_command(x, "expoints",  "numeric", markup = markup) ## default points
  extime    <- extract_command(x, "extime",    "numeric", markup = markup) ## default time in seconds
  exshuffle <- extract_command(x, "exshuffle", "logical", markup = markup) ## shuffle schoice/mchoice answers?
  exsingle  <- extract_command(x, "exsingle",  "logical", markup = markup) ## use radio buttons?
  exmaxchars  <- extract_command(x, "exmaxchars", markup = markup)         ## maximum number of characters in string answers
  exabstention <- extract_command(x, "exabstention", markup = markup)      ## string for abstention in schoice/mchoice answers

  ## User-Defined ###################################
  exextra <- extract_extra(x, markup = markup)

  ## process valid solution types (in for loop for each cloze element)
  slength <- length(exsolution)
  if(slength < 1L) stop("no exsolution specified")
  exsolution <- switch(extype,
    "schoice" = string2mchoice(exsolution, single = TRUE),
    "mchoice" = string2mchoice(exsolution),
    "num" = as.numeric(exsolution),
    "string" = exsolution,
    "cloze" = {
      if(is.null(exclozetype)) {
        warning("no exclozetype specified, taken to be string")
	exclozetype <- "string"
      }
      if(length(exclozetype) > 1L & length(exclozetype) != slength)
        warning("length of exclozetype does not match length of \\exsolution{}")
      exclozetype <- rep(exclozetype, length.out = slength)
      exsolution <- as.list(exsolution)
      for(i in 1L:slength) exsolution[[i]] <- switch(match.arg(exclozetype[i], c("schoice", "mchoice", "num", "string", "verbatim")),
        "schoice" = string2mchoice(exsolution[[i]], single = TRUE),
        "mchoice" = string2mchoice(exsolution[[i]]),
        "num" = as.numeric(exsolution[[i]]),
        "string" = exsolution[[i]],
        "verbatim" = exsolution[[i]])
      exsolution
    })
  slength <- length(exsolution)

  ## lower/upper tolerance value
  if(is.null(extol)) extol <- 0
  extol <- rep(extol, length.out = slength)

  ## compute "nice" string for printing solution in R
  string <- switch(extype,
    "schoice" = paste(exname, ": ", which(exsolution), sep = ""),                                                      ## FIXME: currently fixed
    "mchoice" = paste(exname, ": ", paste(if(any(exsolution)) which(exsolution) else "-", collapse = ", "), sep = ""), ## FIXME: currently fixed
    "num" = if(max(extol) <= 0) {
      paste(exname, ": ", exsolution, sep = "")
    } else {
      if(slength == 1L) {
        paste(exname, ": ", exsolution, " (", exsolution - extol, "--", exsolution + extol, ")", sep = "")
      } else {
	paste(exname, ": [", exsolution[1L], ", ", exsolution[2L], "] ([", exsolution[1L] - extol[1L], "--", exsolution[1L] + extol[1L], ", ",
	  exsolution[2L] - extol[2L], "--", exsolution[2L] + extol[2L], "])", sep = "")
      }
    },
    "string" = paste(exname, ": ", paste(exsolution, collapse = "\n"), sep = ""),
    "cloze" = paste(exname, ": ", paste(sapply(exsolution, paste, collapse = ", "), collapse = " | "), sep = "")
  )

  ## points should be a vector for cloze
  if(!is.null(expoints) & extype == "cloze") {
    expoints <- rep(expoints, length.out = slength)
  }

  ## possible char setting options
  if(!is.null(exmaxchars)) {
    exmaxchars <- rep(exmaxchars, length.out = slength)
    exmaxchars <- lapply(exmaxchars, function(x) {
      x <- gsub("\\s", ",", x)
      x <- strsplit(x, ",")[[1]]
      x <- x[x != ""]
      if(any(x == "NA"))
        x[x == "NA"] <- NA
      mode(x) <- "integer"
      x
    })
    if(slength < 2)
      exmaxchars <- exmaxchars[[1]]
  }

  ## return everything (backward compatible with earlier versions)
  rval <- list(
    file = tools::file_path_sans_ext(file),
    markup = markup,
    type = extype,
    name = exname,
    title = extitle,
    section = exsection,
    version = exversion,
    solution = exsolution,
    tolerance = extol,
    clozetype = exclozetype,
    points = expoints,
    time = extime,
    shuffle = exshuffle,
    single = exsingle,
    length = slength,
    string = string,
    maxchars = exmaxchars,
    abstention = exabstention
  )
  rval <- c(rval, exextra)
  return(rval)
}
