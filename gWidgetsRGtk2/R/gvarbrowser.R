## Use this to filter by type
## knownTypes in common
### Use this for filtering by (gvarbrowser, gvarbrowsertree)

## This is *ugly* -- how to get a reasonable set of values here?
.datasets = c(
  "numeric","logical","factor","character","integer",
  "data.frame","matrix","list",
  "table","xtabs",
  "nfnGroupedData","nffGroupedData","nmGroupedData",
  "POSIXct","POSIXlt","POSIXt"
  )
.models = c("lm","glm","lqs","aov","anova",
  "lme","lmList","gls",
  "ar","arma","arima0","fGARCH","fAPARCH"
    )
.ts = c("ts", "mts", "timeSeries", "its", "zoo","xts")
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

##' function to capture summary (in string) of object
ourStr <- function(x) UseMethod("ourStr")
ourStr.default <- function(x) ""        # say nothing if nothing good to say -- thanks Mom.
ourStr.character <- ourStr.logical <- ourStr.numeric <- function(x) sprintf("length %s", length(x))
ourStr.matrix <- function(x) sprintf("%s by %s", nrow(x), ncol(x))
ourStr.data.frame <- function(x) sprintf("%s variables, %s observations", length(x), nrow(x))
ourStr.list <- function(x) sprintf("%s components", length(x))
ourStr.lm <- function(x) deparse(x$call)






##' Make offspring data frame
offspring = function(path=c(), data=NULL) {
  emptyDf <- data.frame(names="",hasSubTree=FALSE,type="", summary="", stringsAsFactors=FALSE)

  
  if(!is.null(data) && is.function(data)) data <- data()

  ## data is knownClass value. This checks through inheritance but still the
  ## question of what classes to show is hardcoded -- eh
  .inClass <- function(x,data) {
    if(is.null(data)) return(TRUE)
    any(sapply(1:length(data), function(i) {
      out <- is(x,data[i])
      if(inherits(out,"try-error")) return(FALSE)
      return(out)
    }))
  }
  
  if(length(path) == 0) {
    fullx <- x <- ls(envir=.GlobalEnv)
  } else {
    string <- paste(path,collapse="$")
    obj <- getObjectFromString(string)

    x <- with(obj, ls())
    fullx <- paste(string,x,sep="$")
  }

  if(length(x) == 0) {
    return(emptyDf)
  }

  objType <- newNames <- objSummary <- character(0)
  hasTree <- logical(0);
    
  for(i in seq_along(x)) {
    y <-  getObjectFromString(fullx[i])
    if(.inClass(y,data)) {
      j <- length(objType)+ 1
      objType[j] <- str2(y)
      hasTree[j] <- hasSubTree(y)
      newNames[j] <- x[i]
      objSummary[j] <- ourStr(y)
    }
  }
  
  if(length(objType) == 0) {
    return(emptyDf)
  }

  allValues <-  data.frame(names=I(newNames),
                           hasSubTree=hasTree,
                           type=I(objType),
                           summary=I(objSummary), stringsAsFactors=FALSE)
  ## Thanks Stephanie
  if(!is.null(data)) { 
    return(allValues[allValues$type %in% data, ,drop=FALSE]) 
  } else { 
    return(allValues) 
  } 
  
  
}

hasSubTree = function(x) {
  tmp  = gtktry(is.list(x) &&
    !is.guiWidget(x) && !is.gWidget(x) && !is.RGtkObject(x) &&
    !is.null(names(x)), silent=TRUE)
  if(!inherits(tmp,"try-error") && tmp)
    return(TRUE)
  else
    return(FALSE)
}


setClass("gVarbrowserRGtk",
         representation(filter="guiComponent"),
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )

##' Toolkit constructor for gvarbrowser widget
##'
##' Call the update method to update the object
##' 
##' 
##' @example
##' Make a popup menu for actions. Issue: works with selection, but selection updated first by left click
##' is needed.
##' library(gWidgets)
##' options(guiToolkit="RGtk2")
##' v <- gvarbrowser(container =gwindow("Object broser"), handler=function(h,...) {
##' varname <- h$obj[]
##' if(length(varname) == 1) {
##' do.call("fix", list(varname))
##' }
##' })
##' Helper function to get object from argument h passed in to menulist
## getObjectFrom_h <- function(h) {
##   varname <- h$action[]                 # note action, not obj
##   obj <- get(varname[1], envir=.GlobalEnv)
##   if(length(varname) > 1)
##     obj <- obj[[varname[-1]]]
##   obj
## }

## ##' a list of gaction items or separators
## ml <- list(
##            summary=gaction("summary...", action=v, handler=function(h,...) {
##              obj <- getObjectFrom_h(h)
##              print(summary(obj))
##            }),
##            plot=gaction("plot...", action=v, handler=function(h,...) {
##              obj <- getObjectFrom_h(h)
##              try(plot(obj), silent=TRUE)
##            }),
##            sep=gseparator(),
##            remove=gaction("remove", action=v, handler=function(h,...) {
##              varname <- h$action[]
##              print(varname)
##              if(gconfirm(sprintf("Really delete %s?", varname[1])))
##                rm(list=varname[1], envir=.GlobalEnv)
##            })
##            )
## add3rdMousePopupmenu(v, menulist=ml)

setMethod(".gvarbrowser",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   handler = NULL,
                   action = "summary",
                   container = NULL,
                   ...) {

            force(toolkit)

            theArgs <- list(...)
            if(!is.null(theArgs$inteval)) theArgs$interval <- theArgs$interval ## typo fix. Remove later
            interval <- ifelse(is.null(theArgs$interval), 2000, theArgs$interval)

            ## fix up known types
            if(!is.null(theArgs$knownTypes))
              knownTypes <- theArgs$knownTypes
            else if(!is.null(getOption("knownTypes"))) {
              knownTypes <- getOption("knownTypes")
            }

            multiple <- getWithDefault(theArgs$multiple, TRUE)
            
            ## fix handler if action is non-null
            if(is.null(handler) && !is.null(action)) {
              handler = function(h, ...) {
                values <- h$obj[]
                value <- paste(values, collapse = "$")
                if (!is.null(action))
                        print(do.call(h$action, list(svalue(value))))
              }
            }

            ## begin
            group <- ggroup(horizontal=FALSE, container=container,...)
            filterGroup <- ggroup(container=group)
            glabel("Filter by:",container=filterGroup)
            filterPopup <- gdroplist(names(knownTypes), container=filterGroup)
            val <- ifelse("data sets" %in% names(knownTypes), "data sets", names(knownTypes)[1])
            svalue(filterPopup) <- val

            
            ## main tree
            tree = gtree(offspring=offspring,
              offspring.data=function() knownTypes[[svalue(filterPopup)]],
              col.types=data.frame(Name="string",Type="string",Summary="String"),
              icon.FUN = function(d,...) {
                ## Allow user to change stock icon
                FUN <- getWithDefault(getOption("gWidgetsStockIconFromClass"), stockIconFromClass)
                byReturnVector(d,function(x) FUN(x[,'type']))
              },
              multiple=multiple,
              container = group, expand=TRUE
              )

            updateGroup <- ggroup(container =group)
            autoUpdate <- gcheckbox("Auto update", checked=TRUE, use.togglebutton=TRUE, container =updateGroup,
                                    handler=function(h,...) {
                                      enabled(refreshButton) <- !svalue(h$obj)
                                    })
            refreshButton <- gimage("refresh", dirname="stock", container =updateGroup, handler=function(h,...) {
              key <- svalue(filterPopup)
              offspring.data <- knownTypes[[key]]
              update(obj, offspring.data)
            })
            enabled(refreshButton) <- FALSE
            tooltip(refreshButton) <- "Click to refresh display"
            visible(updateGroup) <- FALSE
            
            ## update the tree this way
            addhandlerclicked(filterPopup,
                              handler = function(h,...) {
                                key = svalue(filterPopup)
                                offspring.data = knownTypes[[key]]
                                update(h$action, offspring.data)
                              },
                              action=tree)
            
            
            
            
            ## drop handler
            adddropsource(tree,handler=function(h,...) {
              
              values = h$obj[]
              values = sapply(values, untaintName)
              return(paste(values,collapse="$"))
            })
            
            tag(tree,"view")$SetEnableSearch(TRUE)
            tag(tree,"view")$SetHeadersClickable(TRUE)

            obj <- new("gVarbrowserRGtk",block=group, widget=tree, filter=filterPopup, toolkit=toolkit)

            tag(obj, "filterPopup") <- filterPopup
            
            ## ### In place of an idleHandler, we use a taskCallback
            ## ### This is a little fragile, as remove taskCallback can remove
            ## updateCallback <- function(x) {
            ##   function(expr, value, ...) {
                
            ##     if(!isExtant(x)) return(FALSE)      # need widget to be available, otherwise shut off
                
            ##     if(is.call(expr)) {
            ##       FUN <- deparse(expr[[1]])
            ##       if(FUN %in% c("=","<-", "assign", "rm")) {
            ##         update(x)
            ##       }
            ##     }
            ##     return(TRUE)
            ##   }
            ## }
            ## addTaskCallback(updateCallback(obj), name="gvarbrowser")

            
             ## add an idle handler for updating tree every  second (or interval)
             idleHandler <- function(h,...) {

               visible(updateGroup) <- tag(h$action, "logsize") > 2 #  50 or more

               if(!svalue(autoUpdate))
                 return()


               key = svalue(filterPopup)
               offspring.data = knownTypes[[key]]
               update(h$obj, offspring.data)

               ## do we make timeframe longer bigger?
               n <- ceiling(log(1+ length(.GlobalEnv), 7)) 
               if(n != tag(h$action, "logsize")) {
                 tag(h$action, "logsize") <- n
                 idleid <- tag(h$action, "idleid")
                 gSourceRemove(idleid)
                 tag(h$action, "idleid") <-
                   addhandleridle(tree, interval=2^n*1000, handler = idleHandler, action=h$action)
               }
             } 
             idleid <- addhandleridle(tree, interval=interval, handler = idleHandler, action=obj)
             tag(obj, "idleid") <- idleid
             tag(obj, "logsize") <- 1    # ceiling(log(1+ length(.GlobalEnv), 10))
             addhandlerunrealize(tree, handler = function(h,...) {
               idleid <- tag(h$action, "idleid")
               gSourceRemove(idleid)
             },
                                 action=obj)

            
            ## override how we compare items. Default is just by name, here we want
            ## to include class and summary
            tag(tree, "isStillThere") <- function(old, new) {
              if(length(old) && length(new)) {
                identical(any(ind <- (old[1] == new[,1, drop=TRUE])) &&
                          (old[2] %in% new[which(ind),3, drop=TRUE]) &&
                          (old[3] %in% new[which(ind),4, drop=TRUE]), # for wxf
                          TRUE)         # Tom Taverner change
              } else {
                FALSE
              }
            }
            
            if(!is.null(handler)) {
              id = addhandlerdoubleclick(tree,
                handler=handler, action=action)
            }
            
            ## all done
            return(obj)
          })

### methods
## push methods and handlers down to tree in this case
setMethod(".update",
          signature(toolkit="guiWidgetsToolkitRGtk2",object="gVarbrowserRGtk"),
          function(object, toolkit, ...) {
            filterPopup <- tag(object, "filterPopup")
            key <- svalue(filterPopup)
            offspring.data <- knownTypes[[key]]
            update(object@widget, offspring.data)
          })

setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gVarbrowserRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {

            ## check if any selection
            value <- NA
            x <- svalue(obj@widget)
            if(!(is.atomic(x) && length(x) == 1 && is.na(x))) {
              f <- function(x) paste(x, collapse="$")
              values <- obj@widget[]       # from tree
              if(is.list(values))
                value <- sapply(values, f)
              else
                value <- f(values)
            }
            return(value)
          })
setMethod("[",
          signature(x="gVarbrowserRGtk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x,guiToolkit("RGtk2"), i, j, ..., drop=drop)
          })
setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gVarbrowserRGtk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            if(missing(i))
              x@widget[...]
            else
              x@widget[i,...]
          })




