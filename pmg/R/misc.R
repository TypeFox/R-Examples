## some miscellaneous functions
## Use this to filter by type
## knownTypes in common
### Use this for filtering by (gvarbrowser, gvarbrowsertree)
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
  "data sets"= .datasets,
  "model objects" = .models,
  "time series objects" = .ts,
  "functions"=.functions,
  "saved plots" = .plots,
  "all" = NULL
  )

## list of some type
lsType = function(type, envir=.GlobalEnv) {
  x = with(.GlobalEnv, sapply(ls(), function(i) class(get(i))))
  objects = names(x)[sapply(x, function(i) any(i %in% type))]
  return(objects)
}
lsDatasets = function(envir=.GlobalEnv)  lsType(.datasets, envir)
lsModels = function(envir=.GlobalEnv)  lsType(.models, envir)
lsTs = function(envir=.GlobalEnv)  lsType(.ts, envir)
lsFunctions = function(envir=.GlobalEnv)  lsType(.functions, envir)

###  These should be in gWidgets or gWidgetsRGtk2, but arent
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




##################################################
##
## make a fancy summary function for showing on double click
## in varbrowser
## make generic, but not needed
pmgSummary = function(obj,...) UseMethod("pmgSummary")
pmgSummary.default = function(obj, ...) {
  ## what is object?
  objName = deparse(substitute(obj))

  if(is.character(obj) && length(obj) == 1) {
    ## assume it is a string containing object
    objName = obj
    obj = svalue(obj)

  }

  
  group = ggroup(horizontal = FALSE,...)

  ## should I export this function?
  icon = stockIconFromClass(class(obj))
  add(group, gimage(icon, dirname="stock", size="DIALOG"))
  table = glayout(adjust="left")
  add(group, table)

  table[1,1] = glabel("<b>Name:</b>", markup=TRUE)
  table[1,2] = glabel(objName)

  table[2,1] = glabel("<b>Kind:</b> ", markup=TRUE)
  table[2,2] = glabel(paste(class(obj),sep="",collapse=", "))


  table[3,1] = glabel("<b>Size:</b>",markup=TRUE)
  if(!is.function(obj)) {
    theSize = str1(obj)$dim.field
    table[3,2] = glabel(theSize)
  } else {
    table[3,2] = glabel("NA")
  }

  stamp = Timestamp(obj)
  if(!is.na(stamp)) {
    table[4,1] = glabel("<b>Last modified:</b>", markup=TRUE)
    table[4,2] = glabel(format(as.Date(stamp), "%B %d, %Y"))
  }

  table[5,1] = glabel("<b>Preview:</b>", markup=TRUE)
  theValue = capture.output(eval(obj))
  if(length(theValue) > 10)
    theValue = c(theValue[1:10],"... 8< snipped >8 ...")
  theHead = gtext(font.attr=c(style="monospace"))
  add(theHead,theValue)
  enabled(theHead) <- FALSE
  add(group, theHead, expand=TRUE)
    
  visible(table) <- TRUE

  return(group)
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
  while(!is.RGtkObject(obj))
    obj = obj@block
  ifelse("<invalid>" %in% class(obj), TRUE, FALSE)
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
  ## is this making windows bug out?
  if(class(obj)[1] == "GtkTreeViewColumn")
    return(TRUE)
  else
    return(FALSE)
}

## Function to convert back and forth between R classes and GObject classes
RtoGObjectConversion = function(obj) {
  if("gComponent" %in% class(obj)) return("GObject")
  if(is.list(obj)) return("GObject")
  
  Klasse = class(obj)[1]                # silly name?
  switch(Klasse,
         "integer"="gint",
         "numeric"="gdouble",
         "gtk"="GObject",
         "logical" = "gboolean",
         "gchararray"
         )
}


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


## get does not work with name$component, this gets around that
## returns NULL if not available
getObjectFromString = function(string="", envir=.GlobalEnv) {
  tmp = try(get(string, envir), silent = TRUE)
  if(!inherits(tmp, "try-error")) return(tmp)
  
  tmp = try(rpel(string,envir), silent=TRUE)
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

Paste = function(..., sep="", collapse="") {
  x = unlist(list(...))
  x = x[!is.na(x)]
  x = x[x != "NA"]
  paste(x, sep=sep, collapse=collapse)
}

stripWhiteSpace = function(str) {
  sub('[[:space:]]+$', '', str) ## from ?gsub
  sub('^[[:space:]]+', '', str) ## from ?gsub
  return(str)
}
## ReadParseEvaL -- saves typing
rpel = function(string, envir=.GlobalEnv) {
  eval(parse(text=string), envir=envir)
}



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



##################################################
## define skewness and kurtosis
skewness = function(x, na.rm=TRUE,...) UseMethod("skewness")
### FROM http://finzi.psych.upenn.edu/R/Rhelp02a/archive/44065.html
skewness.factor <- function(x, na.rm=TRUE, ...) NA
skewness.character <- skewness.factor
skewness.list = function(x, na.rm=TRUE, ...) sapply(x,skewness)
skewness.data.frame = function(x, na.rm=TRUE, ...) sapply(x,skewness)
skewness.default =  function(x, na.rm=TRUE, ...)  {
  ## Remove NAs:
  if (na.rm) x = x[!is.na(x)]

  ## Warnings:
  if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("argument is not numeric or logical: returning NA")
    return(as.numeric(NA))}
  
  
  ## Skewness:
  n = length(x)
  if (is.integer(x)) x = as.numeric(x)
  skewness = sum((x-mean(x))^3/sqrt(var(x))^3)/length(x)
  
  ## Return Value:
  skewness
}
kurtosis <- function(x, na.rm=TRUE, ...) UseMethod("kurtosis")
kurtosis.list <- function(x, na.rm=TRUE, ...) sapply(x, kurtosis) # lazy == na.rm?
kurtosis.factor <- function(x, na.rm=TRUE, ...) return(NA)
kurtosis.character <- kurtosis.factor 
kurtosis.data.frame = function(x, na.rm=TRUE, ...) sapply(x, kurtosis)
kurtosis.default = function(x, na.rm=TRUE, ...) {
  ## Remove NAs:
  if (na.rm) x = x[!is.na(x)]
  
  n = length(x)
  if (is.integer(x)) x = as.numeric(x)
  kurtosis = sum((x-mean(x))^4/var(x)^2)/length(x) - 3
                                        
  kurtosis
}

