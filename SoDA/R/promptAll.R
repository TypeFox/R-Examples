mergeSections <- function(input, keepFirst = 1, keepLast = 1, sep = "") {
   item <- input[[1]]
   start <- item[seq(length=keepFirst)]
   end <- item[seq(length=keepLast, to=length(item))]
   middle <- character(); all <- c(start, end)
  for(i in seq(along=input)) {
    item <- input[[i]]
    item <- item[is.na(match(item, all))]
    middle <- c(middle, sep, item)
    all <- c(all, item)
  }
   c(start, middle, end)
}

.fdummy <- function()NULL
.fdoctype <- elNamed(prompt(.fdummy, NA), "doctype")
.isFunctionDoc <- function(object)
  identical(elNamed(object, "doctype"), .fdoctype)

.functionDocEl <- function(x, name) {
    if(.isFunctionDoc(x))
      elNamed(x, name)
    else
      character()
  }

promptAll <- function(objects, name, filename = paste(objects[[1]],
  "Rd", sep="."), where = topenv(parent.frame()), ...) {
  if(length(objects) == 0)
    stop("No objects to document!")
  allDocs <- as.list(objects)
  for(i in seq(along = objects)) {
    WHAT <- as.name(objects[[i]])
    allDocs[[i]] <- eval(substitute(prompt(WHAT, NA, ...)), where)
  }
  doc <- allDocs[[1]]
  aliases <- lapply(allDocs, elNamed, name = "aliases")
  elNamed(doc, "aliases") <- mergeSections(aliases, 0, 1, character())
  usage <- lapply(allDocs, .functionDocEl, name = "usage")
  elNamed(doc, "usage") <- mergeSections(usage, 1, 2, "")
  arguments <- lapply(allDocs, .functionDocEl, name = "arguments")
  elNamed(doc, "arguments") <- mergeSections(arguments, 1, 1, "")
  if(!missing(name)) {
    elNamed(doc, "name") <- paste("\\name{", name, "}", sep="")
    if(missing(filename))
      filename <- paste(name, "Rd", sep=".")
  }
    if(!is.na(filename)) {
      cat(unlist(doc), file = filename, sep = "\n")
      message("Wrote documentation to \"", filename, "\"")
      invisible(filename)
    }
    else
      doc
}

  
  
