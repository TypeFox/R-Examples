## Unfortunately, FAME does not provide any way at all to discover
## how many attributes an object has, or what their names are. You have to
## know what you're looking for.  

fameAttribute <- function(attribute, fname, db){
  path <- getFamePath(db)
  if(!fameRunning()) fameStart()
  
  openCommand <- paste('open <access read> "', path, '" as blah', sep = "")
  attrCommand <- paste("display ", attribute, "(", fname, ")", sep = "")
  
  if(fameCommand(openCommand) != 0) stop("couldn't open db")
  on.exit(fameCommand("close blah"))
  
  strings <- fameCommand(attrCommand, silent = TRUE, capture = TRUE)
  status <- attr(strings, "status")
  if(status != 0) stop(paste("bad status", status, "returned by", attrCommand))
  retval <- tolower(strings[strings != ""])
  rvname <- retval[1]
  retval <- retval[-1]
  if(length(retval) > 1) 
    retval <- list(retval)
  names(retval) <- rvname
  retval
}

fameSetAttribute <- function(attribute, value, fname, db){
  closeBlah <- function(){
    ## an internal function to make sure the BLAH channel gets closed
    ## close channel BLAH if it was left open
    zz <- fameCommand("type @open.db", capture = TRUE)
    if(tolower(as.vector(zz)) == "blah") fameCommand("close blah")
  }

  path <- getFamePath(db)
  if(!fameRunning()) fameStart()
  
  openCommand <- paste('open <access shared> "', path, '" as blah', sep = "")
  closeBlah()
  if(fameCommand(openCommand) != 0) stop("couldn't open db")
  on.exit(closeBlah())

  if(is.numeric(value)){
    oldOpt <- options(scipen = 15)
    on.exit(options(oldOpt))
  } else {
    value <- dQuote(value)
  }
  attrCommand <- paste("attribute ", attribute, "(", fname, ") = ", value, sep = "")
  
  strings <- fameCommand(attrCommand, silent = TRUE, capture = TRUE)
  closeBlah()
  status <- attr(strings, "status")
  if(status != 0) stop(paste("bad status", status, "returned by", attrCommand))
  retval <- tolower(strings[strings != ""])
  retval
}


fameUpdated <- function(fname, db){
  attString <- fameAttribute("updated", fname, db)
  attString <- gsub("\\..*", "", attString)
  if(grepl("-", attString)) fmt <- "%d-%b-%y"
  else                      fmt <- "%Y%m%d:%H:%M:%S"
  strptime(attString, format = fmt)
}

fameCreated <- function(fname, db){
  attString <- fameAttribute("created", fname, db)
  attString <- gsub("\\..*", "", attString)
  if(grepl("-", attString)) fmt <- "%d-%b-%y"
  else                      fmt <- "%Y%m%d:%H:%M:%S"
  strptime(attString, format = fmt)
}

fameDbStatus <- function(db){
  path <- getFamePath(db)
  if(!fameRunning()) fameStart()
  
  openCommand <- paste("open <access read> \"", path, "\" as blah", sep = "")
  if(fameCommand(openCommand) != 0) stop("couldn't open db")
  on.exit(fameCommand("close blah"))
  
  strings <- fameCommand("dbstatus blah", capture = T)
  status <- attr(strings, "status")
  attr(strings, "status") <- NULL
  if(status != 0) 
    stop(paste("bad status", status, "returned by dbstatus blah"))
  else
    strings
}

fameUserDefinedAttributes <- function(db){
  strings <- tolower(fameDbStatus(db))
  i <- grep("user defined attribute names", strings)
  if(length(i) == 0) return(NULL)
  strings <- c(strings[-(1:(i+1))], "")
  ii <- grep("attribute_names:", strings, ignore.case = TRUE)
  if(length(ii) == 0) return(NULL)
  else atts <- list()

  blank <- seq(along = strings)[strings == ""]
  for(i in seq(along = ii)){
    start <- ii[i]
    type <- gsub("_attribute_names:", "", strings[start])
    nameString <- gsub(" *", "", paste(strings[(start+1):blank[i]], collapse = ""))
    atts[[type]] <- strsplit(nameString, ",")[[1]]
  }
  atts
}

fameAddAttribute <- function(name,
                             type = c("string", "date", "boolean",
                               "precision", "numeric", "namelist"),
                             db){
  type <- match.arg(type)
  attList <- list()

  currentAtts <- fameUserDefinedAttributes(db)
  ctypes <- names(currentAtts)
  if(length(currentAtts) > 0){
    for(ctype in ctypes){
      cAtts <- currentAtts[[ctype]]
      if(name %in% cAtts){
        if(ctype == type){
          cat(paste("Db already has attribute", name), "\n")
          return(NULL)
        }
        else stop(paste("Db already has attribute", name,
                        "but of type", ctype, "not", type))
      }
    }
    if(type %in% ctypes)
      attList[[type]] <- c(currentAtts[[type]], name)
    else
      attList[[type]] <- name
  }
  else attList[[type]] <- name
  fameSetDbAttributes(attList, db)
}

fameSetDbAttributes <- function(attributeList, db){
  atts <- attributeList
  if(!is.list(atts)) stop("attributeList is not a list")
  fameAttTypes <- c("string", "date", "boolean", "precision", "numeric", "namelist")
  if(!all(names(atts) %in% fameAttTypes)) stop("attributeList has bad names")
  
  path <- getFamePath(db)
  if(!fameRunning()) fameStart()
  
  openCommand <- paste('open <access shared> "', path, '" as blah; overwrite on', sep = "")
  if(fameCommand(openCommand) != 0) stop("couldn't open db")
  on.exit(fameCommand("overwrite off; close blah"))
  outstrings <- character(0)

  for(nm in names(atts)){
    cmd <- paste(nm, "_attribute_names = {",
                 paste(atts[[nm]], collapse = ", "),
                 "}", sep = "")
    strings <- fameCommand(cmd, silent = TRUE, capture = TRUE)
    status <- attr(strings, "status")
    outstrings <- c(outstrings, tolower(strings))
    attr(outstrings, "status") <- status
    if(status != 0)
      stop(paste(paste("bad status", status, "returned by", cmd),
                 outstrings, sep = "\n"))
  }
  invisible(outstrings)
}

