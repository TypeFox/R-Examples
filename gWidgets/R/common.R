##  Some functions that seem to be useful. These are not exported,
## hence repeated in gWidgetsR****. This should likely be changed, but
## for now it isn't

gwCat <- function(...) {
  doCat <- getOption("gWidgetsDebug")
  if(!is.null(doCat)) message(...)
}

Paste = function(..., sep="", collapse="") {
  x = unlist(list(...))
  x = x[!is.na(x)]
  x = x[x != "NA"]
  paste(x, sep=sep, collapse=collapse)
}


PasteWithComma = function(...) {
  args = unlist(list(...))
  args = args[!is.na(args)]
  paste(args, sep="", collapse=", ")
}

quoteIfNeeded = function(str) {
  if(length(grep('^\\".*\\"$', str, perl=TRUE)) > 0 ||
     length(grep("^\\'.*\\'$", str, perl=TRUE)) > 0 )
    return(str)
  else
    return(paste('"',str,'"',sep="",collapse=""))
}

## ReadParseEvaL -- saves typing
rpel = function(string, envir=.GlobalEnv) {
  eval(parse(text=string), envir=envir)
}

## return type
str1 <- function(obj) {
  md <- mode(obj)
  lg <- length(obj)
  objdim <- dim(obj)
  if (length(objdim) == 0) 
    dim.field <- paste("length:", lg)
  else {
    dim.field <- "dim:"
    for (i in 1:length(objdim)) dim.field <- paste(dim.field, 
                                                   objdim[i])
    if (is.matrix(obj)) 
      md <- "matrix"
  }
  obj.class <- oldClass(obj)
  if (!is.null(obj.class)) {
    md <- obj.class[1]
    if (inherits(obj, "factor")) 
      dim.field <- paste("levels:", length(levels(obj)))
  }
  list( type = md, dim.field = dim.field)
}


## untaint a variable name so that $ can be used
untaintName = function(objName) {
  if (length(grep(" |\\+|\\-|\\*|\\/\\(|\\[|\\:",objName)) > 0) {
    objName=Paste("\"",objName,"\"")
  }
  return(objName)
}

## as name says
stripWhiteSpace = function(str) {
  sub('[[:space:]]+$', '', str) ## from ?gsub
  sub('^[[:space:]]+', '', str) ## from ?gsub
  return(str)
}

### We should make this configurable. Likely using options
### Use this for filtering
.datasets = c(
  "numeric","logical","factor","character",
  "data.frame","matrix","list",
  "table","xtabs",
  "nfnGroupedData","nffGroupedData","nmGroupedData"
  )
.models = c("lm","glm","lqs","aov","anova",
    "lme","lmList","gls",
  "ar","arma","arima0","fGARCH","fAPARCH"
    )
.ts = c("ts", "mts", "timeSeries", "its", "zoo")
.functions=c("function")
.plots = c("recordedplot")

knownTypes = list(
  "data sets and models"=c(.datasets, .models, .ts),
  "data sets"= c(.datasets,.ts),
  "model objects" = .models,
  "time series objects" = .ts,
  "functions"=.functions,
  "saved plots" = .plots,
  "all" = NULL
  )


## get does not work with name$component, this gets around that
## returns NULL if not available
getObjectFromString = function(STRING="", envir=.GlobalEnv) {
  tmp = try(get(STRING, envir), silent = TRUE)
  if(!inherits(tmp, "try-error")) return(tmp)
  
  tmp = try(rpel(STRING,envir), silent=TRUE)
  if(!inherits(tmp, "try-error"))  return(tmp)

  ## out of chances
  return(NULL)
}


## Find the name of the object by pasting toghther the pieces
## better to do name$name, but value may be a numeric
makeObjectName = function(root,value) {
  if(is.null(root)) return(untaintName(value))

  ## now decide between $ and [[]]
  if(value == make.names(value)) {
    return(Paste(root,"$",untaintName(value)))
  } else {
    return(Paste(root,"[['",value,"']]"))
  }
}



## a function to get objects and their types
## filter is a vector of classes
getObjectsWithType = function(root=NULL, filter = NULL, envir=.GlobalEnv) {

  if(is.null(root)) {
    objects = ls(envir=envir)
  } else {
    string = Paste("with(",root,",ls())")
    objects = try(rpel(string,envir=envir), silent=TRUE)
  }
  ## objects is character vector of components of root.
  badnames = grep("[[<-]|\\*",objects)
  if(length(badnames) > 0)
    objects = objects[-badnames]

  objectsWithRoot = sapply(objects,function(i) makeObjectName(root,i))

  
  type = sapply(objectsWithRoot, function(i) {
    string = Paste("str2(",i,")")
    rpel(string, envir=envir)
  })

  objects = data.frame(Name=objects,Type=type,stringsAsFactors=FALSE)

  ## filter
  if(!is.null(filter))
    objects = objects[type %in% filter,]

  return(objects)
  
  
}


## get the names of the object, if available (datastores)
getNamesofObject = function(string="") {
  ## if empty string, get variables in .GlobalEnv
  if(string == "") {
    ## return objects of certain type
    objects = getObjectsWithType(root=NULL, filter=knownTypes[['data sets']])
    return(unlist(objects$Name))
  } 
  obj = getObjectFromString(string)
  if(!is.null(obj)) {
    if(is.list(obj)) {
      return(names(obj))
    } else if(is.matrix(obj)) {
      return(colnames(obj))
    } else{
      return(NULL)
    }
  } else {
    return(NULL)
  }
}



## used to check output 
is.empty = function(obj) {
  
  if(is.null(obj) || is.na(obj) || (is.character(obj) && obj == "")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


## fix font mess up DEPRECATED!!
.fixFontMessUp = function(val) {
  if(is.vector(val)) {
    tmp = val
    val = list()
    for(i in names(tmp)) val[[i]] = tmp[i]
  }

  weights = c("normal","oblique", "italic")
  styles = c("ultra-light","light","normal","bold","ultra-bold", "heavy")

  a = val$weight; b = val$style
  if((!is.null(a) && a %in% styles) || (!is.null(b) &&b %in% weights)) {
    tmp = a
    val$weight <- b
    val$style <- a
  }
  return(val)
}

## what type of object is thixs and a size
str2 <- function(obj) {
  md <- mode(obj)
  if (is.matrix(obj))  md <- "matrix"
  obj.class <- oldClass(obj)
  if (!is.null(obj.class)) {
    md <- obj.class[1]
  }
  return(md)
}

## blurb about installation -- put in so can be updated easily
## look to file to update
installing_gWidgets_toolkits <- function() {
  file <- system.file("install/Installing_gWidgets_Toolkits.txt", package="gWidgets")
  tmp <- readLines(file)
  for(i in tmp) cat(i,"\n")

}

##' return x unless null then give default
getWithDefault <- function(x, default)
  ifelse(is.null(x) || is.na(x) || x == "", default, x)


## ID is used by the ANY widgets
n=0;assignInNamespace("n",0,"gWidgets")
getNewID = function() {                 # get new one, incremented
  n = getFromNamespace("n",ns="gWidgets")
  assignInNamespace("n",n+1,ns="gWidgets")
  return(n+1)
}
