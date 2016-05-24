## Common functions


#Paste = function(x,...) paste(x,...,sep="",collapse="")
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
stripWhiteSpace = function(str) {
  sub('[[:space:]]+$', '', str) ## from ?gsub
  sub('^[[:space:]]+', '', str) ## from ?gsub
  return(str)
}

toupperFirst <- function(str="") {
  if(nchar(str) == 0)
    return(str)
  out <- toupper(substr(str, 0, 1))
  if(nchar(str) > 1) {
    out <- paste(out, substr(str,2, nchar(str)), sep="")
  }
  return(out)
}
           
quoteIfNeeded = function(str) {
  if(length(grep('^\\".*\\"$', str, perl=TRUE)) > 0 ||
     length(grep("^\\'.*\\'$", str, perl=TRUE)) > 0 )
    return(str)
  else
    return(paste('"',str,'"',sep="",collapse=""))
}


showErrorMessage = function() {       # leave here for scoping on command
  message = Paste("Error:",
    "\n\t",geterrmessage())
  gmessage(message,icon="error")
  stop()
}
  

## Push and Pop -- for convenience
Push = function(v,d) c(v,d)
Pop = function(v) ifelse(length(v) > 1, v[-length(v)], NA)


### is functions
is.RGtkObject = function(obj) {
  is(obj,"RGtkObject") 
}

is.guiWidget = function(obj) {
  is(obj,"guiWidget")
}
is.gWidget = function(obj) {
  is(obj,"gWidgetRGtk")
}
is.gWindow = function(obj) {
  is(obj,"gWindowRGtk")
}
is.gComponent = function(obj) {
  is(obj,"gComponentRGtk")
}
is.gContainer = function(obj) {
  is(obj,"gContainer")
}

is.gImage = function(obj) {
  is(obj,"gImageRGtk")
}
is.gLabel = function(obj) {
  is(obj,"gLabelRGtk") 
}

is.gMenu = function(obj) {
  is(obj,"gMenuRGtk") 
}
is.gEditDataFrame=function(obj) {
  stop("deprecated, use is.gGrid")
}
is.gGrid = function(obj) {
  is(obj,"gGridRGtk")
}

is.invalid = function(obj) {
  widget = getWidget(obj)
  parent = widget$GetParentWindow()
  ifelse(inherits(parent,"<invalid>"), TRUE, FALSE)
}
## used to check output 
is.empty = function(obj) {
  if(is.null(obj) || is.na(obj) || obj == "") {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


## for showing only possible values
is.dataframelike = function(obj) {
  if(is.data.frame(obj) || is.matrix(obj) ||
     is.numeric(obj) || is.logical(obj) ||
     is.factor(obj)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

## check if a gtkTreeViewCOlumn, make no GTK language
is.gdataframecolumn = function(obj) {
  if(class(obj)[1] == "GtkTreeViewColumn")
    return(TRUE)
  else
    return(FALSE)
}

## Function to convert back and forth between R classes and GObject classes


### these are used by gvarbrowser
## This is from browseEnv in base
## what type of object is thixs and a size
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

## untaint a variable name so that $ can be used
untaintName = function(objName) {
  if (length(grep(" |\\+|\\-|\\*|\\/\\(|\\[|\\:",objName)) > 0) {
    objName=Paste("\"",objName,"\"")
  }
  return(objName)
}

## try to stip off data frame stuff in fron to DND target
findDataParent = function(x) {
  child = sub(".*]]","",x)
  child = sub(".*\\$","",child)
  parent = sub(Paste(child,"$"),"",x)
  parent = sub("\\$$","",parent)
  return(list(child=child,parent=parent))
}


## basically repeat findDataParent until no parent
findRootObject = function(x) {
  x = sub("\\[\\[.*","",x)
  x = sub("\\$.*","", x)
  return(x)
}

## ReadParseEvaL -- saves typing
rpel = function(string, envir=.GlobalEnv) {
  eval(parse(text=string), envir=envir)
}


## get does not work with name$component, this gets around that
## returns NULL if not available
getObjectFromString = function(string="", envir=.GlobalEnv) {
  tmp = gtktry(get(string, envir), silent = TRUE)
  if(!inherits(tmp, "try-error")) return(tmp)
  
  tmp = gtktry(rpel(string,envir), silent=TRUE)
  if(!inherits(tmp, "try-error"))  return(tmp)

  ## out of chances
  return(NULL)
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

## a function to get objects and their types
## filter is a vector of classes
getObjectsWithType = function(root=NULL, filter = NULL, envir=.GlobalEnv) {

  if(is.null(root)) {
    objects = ls(envir=envir)
  } else {
    string = Paste("with(",root,",ls())")
    objects = gtktry(rpel(string,envir=envir), silent=TRUE)
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

  objects = data.frame(Name=I(objects),Type=I(type))

  ## filter
  if(!is.null(filter))
    objects = objects[type %in% filter,]

  return(objects)
  
  
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



######
## send a file to csv mode for editing
"browseDataAsCSV" <-
  function(x) {

    x = try(as.data.frame(x))
    if(inherits(x,"try-error")) {
      stop("Can not coerce data into a data frame")
    }

    tmpfile = paste(tempfile(),".csv",sep="",collapse="")
    write.csv(x,file=tmpfile)
    browseURL(paste("file://",tmpfile,sep="",collapse=""))

  }

## help out with gtree
byReturnVector = function(df, FUN,...) {
  tmp = by(df, factor(1:nrow(df)), FUN)
  sapply(tmp, function(x) x)
}

hack.as.data.frame = function(items) {
  ## check rectangular, or coerce to rectangules
  if(!(is.data.frame(items) || is.matrix(items) || is.vector(items))) {
    warning("Needs rectangular data, either a vector, matrix or data.frame")
    return(NA)
  }
  
  ## coerce to data frame
  if(is.vector(items)) {
    itemsName = deparse(substitute(items))
    items = data.frame(I(items))
    names(items) = itemsName
  }
  if(is.matrix(items)) {
    items = hack.as.data.frame.matrix(items) # fun in common.R
  }

  ## make a data frame (CO2)
  items <- as.data.frame(items)
  
  return(items)
}

## no easy way to not convert character vectors. This is a hack.
hack.as.data.frame.matrix = 
  function (x, row.names = NULL, optional = FALSE) 
  {
    d <- dim(x)
    nrows <- d[1]
    ir <- seq(length = nrows)
        ncols <- d[2]
    ic <- seq(length = ncols)
    dn <- dimnames(x)
    if (missing(row.names)) 
      row.names <- dn[[1]]
    collabs <- dn[[2]]
    if (any(empty <- nchar(collabs) == 0)) 
      collabs[empty] <- paste("V", ic, sep = "")[empty]
    value <- vector("list", ncols)
    if (mode(x) == "character") {
      for (i in ic) value[[i]] <- as.character(x[, i])
    }
    else {
      for (i in ic) value[[i]] <- as.vector(x[, i])
    }
    if (length(row.names) != nrows) 
      row.names <- if (optional) 
        character(nrows)
      else as.character(ir)
    if (length(collabs) == ncols) 
      names(value) <- collabs
    else if (!optional) 
      names(value) <- paste("V", ic, sep = "")
        attr(value, "row.names") <- row.names
    class(value) <- "data.frame"
    value
  }

##################################################
## timestamp function for objects made with pmg
## Modified from R mailing list, value is comment. Need <- to act in
## OO manner comment needs to be a character vector. If a list were
## okay (say serialize()) then this could be different
"Timestamp<-" <- function(obj,value) {
  currentStamp = Timestamp(obj)
  currentStamp = c(currentStamp, timestamp=as.character(Sys.time()),comment=value)
  comment(obj) <- currentStamp
  return(obj)
}

Timestamp = function(obj,k=1) {
  currentComment= comment(obj)
  allStamps =comment(obj)[names(comment(obj)) %in% "timestamp"]
  n = length(allStamps)
  if(n > 0)
    return(allStamps[(max(1,n+1-k)):n])
  else
    return(NA)
}

#################################################
## functions to deal with icons
## class to icon translation -- return stock name
## with prefix
stockIconFromClass = function(theClass=NULL) {
  default = "symbol_star"
  
  if(is.null(theClass) ||
     is.na(theClass) ||
     length(theClass) == 0
     )
    return(NA)


  theClass = theClass[1]

  if(theClass %in% .models)
    return(getstockiconname("lines"))
  if(theClass %in% .ts)
    return(getstockiconname("ts"))
  if(theClass %in% .functions)
    return(getstockiconname("function"))

  ret = switch(theClass,
    "numeric"= "numeric",
    "integer"= "numeric",
    "logical" = "logical",
    "character"="select-font",
    "matrix" = "matrix",
    "data.frame" = "dataframe",
    "list" = "dataframe",
    "complex"="numeric",
    "factor"="factor",
    "recordedplot" = "plot",
    NA)
  
  return(getstockiconname(ret))
}

stockIconFromObject = function(obj)
  stockIconFromClass(class(obj)[1])

## get with default
getWithDefault <- function(x, default) {
  if(is.null(x) || (is.atomic(x) && length(x) ==1 && is.na(x)))
    default
  else
    x
}
