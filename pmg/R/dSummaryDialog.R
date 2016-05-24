
dSummaryDialog = function(container=NULL) {

  defaultMsg = "Drop variable here"
  dynamicWarning = "Editing values doesn't work with dynamic data"
  
  ## default is to have a window
  if(is.null(container))
    container = pmgWC$new("Dynamic Summaries", visible=TRUE)

  gp = ggroup(horizontal = FALSE, container=container, raise.on.dragmotion = TRUE)
  obj = gp
#  class(obj) = c("gDSummary","gComponent","gWidget")


  tag(obj,"dropHandlers") <- list()
  
  tbl = glayout()

  tbl[1,1] = glabel("x:")
  tbl[2,1] = glabel("[y]:")
  tbl[3,1] = glabel("[group by]:")

  
  xVar = glabel(defaultMsg,editable=TRUE)
  font(xVar) <- c(style="bold")
  tag(obj,"xVarData") <- NULL
  tbl[1,2] = xVar

  yVar = glabel(defaultMsg,editable=TRUE)
  font(yVar) <- c(style="bold")
  enabled(yVar)<-FALSE                  # initially not enabled
  tag(obj,"yVarData") <- NULL
  tbl[2,2] = yVar

  groupingVar = glabel(defaultMsg,editable=TRUE)
  font(groupingVar) <- c(style="bold")
  tag(obj,"groupingVarData") <- NULL
  tbl[3,2] = groupingVar

  add(gp,tbl)
  visible(tbl) <- TRUE

  popupGroup = ggroup(container=gp)
  addSpring(popupGroup)
  gbutton("clear",container=popupGroup, handler=function(h,...) clear())
  glabel("Select summary:",container=popupGroup)

  ## changes here need to propogate to indices below
  univariateSummaries = c("summary","length","mean","median","sd","IQR","mad","range","skewness","kurtosis")
  bivariateSummaries = c("cor")
  actionPopup = gdroplist(
#    c("Clear variables",
    c("--- Summary ---",
      univariateSummaries[1:2],
      "--- Center ---",
      univariateSummaries[3:4],
      "--- Spread ---",
      univariateSummaries[5:8],
      "--- Shape ---",
      univariateSummaries[9:length(univariateSummaries)],
      "--- Bivariate ---",
      bivariateSummaries
      ),selected=3, container=popupGroup)

  gseparator(container=gp)
  summaryArea = gtext()
  add(gp,summaryArea, expand=TRUE)
  add(summaryArea,
      c("Add variables above either by clicking on bold-faced values, or dragging and dropping values.",
        "The summary appears in this window.",
        "The 'Clear variables' options resets the values."))

  ## add handlers
  addhandlerchanged(xVar,
                    handler = function(h,...) {
                      cat(dynamicWarning,"\n")
                      ids = tag(obj,"dropHandlers")
                      if(length(ids) > 0) {
                        removehandler(obj,ids)
                        tag(obj,"dropHandlers") <- list()
                      }
                      tag(obj, "xVarData") <- svalue(h$obj)
                      ## put popup on 1
                      ## svalue(tag(obj,"actionPopup"),index=TRUE) <- 1
                      
                      update()
                    })
  adddroptarget(xVar,
                handler=function(h, ...) {
                  tag(obj,"xVarData") <- h$dropdata
                  svalue(xVar) <- id(h$dropdata)
                  ## put popup on 1
                  ## svalue("actionPopup",index=TRUE) <- 1
                  
                  update()

                  ## now bind to be dynamic *if* a treeviewcolumn
                  if(is.gdataframecolumn(h$dropdata)) {
                    view.col = h$dropdata
                    id = addhandlerchanged(view.col,
                      signal = "edited",
                      handler=function(h,...) update()
                      )
                    dropHandlers = tag(obj,"dropHandlers")
                    dropHandlers[[length(dropHandlers)+1]] = list(
                                  view.col = view.col,
                                  id = id
                                  )
                    tag(obj,"dropHandlers") <- dropHandlers
                  }
                })
  ## yvar
  addhandlerchanged(yVar,
                    handler = function(h,...) {
                      cat(dynamicWarning,"\n")
                      ids = tag(obj,"dropHandlers")
                      if(length(ids) > 0) {
                        removehandler(obj,ids)
                        tag(obj,"dropHandlers") <- list()
                      }
                      tag(obj,  "yVarData") <-  svalue(h$obj)
                      ## put popup on 1
                      ## svalue(tag(obj,"actionPopup"),index=TRUE) <- 1
                      
                      update()
                    })
  adddroptarget(yVar,
                handler=function(h, ...) {
                  tag(obj,"yVarData") <- h$dropdata
                  svalue(yVar) <-  id(h$dropdata)
                  ## put popup on 1
                  ## svalue("actionPopup", index=TRUE) <- 1
                  
                  update()

                  ## now bind to be dynamic *if* a treeviewcolumn
                  if(is.gdataframecolumn(h$dropdata)) {
                    view.col = h$dropdata
                    id = addhandlerchanged(view.col,
                      signal = "edited",
                      handler=function(h,...) update()
                      )
                    dropHandlers = tag(obj,"dropHandlers")
                    dropHandlers[[length(dropHandlers)+1]] = list(
                                  view.col = view.col,
                                  id = id
                                  )
                    tag(obj,"dropHandlers") <- dropHandlers
                  }
                })


    addhandlerchanged(groupingVar,
                    handler = function(h,...) {
                      cat(dynamicWarning,"\n")
                      ids = tag(obj,"dropHandlers")
                      if(length(ids) > 0) {
                        removehandler(obj,ids)
                        tag(obj,"dropHandlers") <- list()
                      }
                      tag(obj,  "groupingVarData") <-  svalue(h$obj)
                      ## put popup on 1
                      ## svalue(tag(obj,"actionPopup"),index=TRUE) <- 1
                      
                      update()
                    })
  adddroptarget(groupingVar,
                handler=function(h, ...) {
                  tag(obj,"groupingVarData") <- h$dropdata
                  svalue(groupingVar) <-  id(h$dropdata)
                  ## put popup on 1
                  ## svalue("actionPopup",index=TRUE) <- 1
                  
                  update()

                  ## now bind to be dynamic *if* a treeviewcolumn
                  if(is.gdataframecolumn(h$dropdata)) {
                    view.col = h$dropdata
                    id = addhandlerchanged(view.col,
                      signal = "edited",
                      handler=function(h,...) update()
                      )
                    dropHandlers = tag(obj,"dropHandlers")
                    dropHandlers[[length(dropHandlers)+1]] = list(
                                  view.col = view.col,
                                  id = id
                                  )
                    tag(obj,"dropHandlers") <- dropHandlers
                  }
                })

  addhandlerchanged(actionPopup, handler = function(h,...) {
    actionVal = svalue(actionPopup)
    if(actionVal == "Select summary:") {

    } else if(actionVal == "Clear variables") {
      clear()
      enabled(yVar)<-FALSE
    } else if(length(grep("^---",actionVal))>0) {
      ##cat("No summary selected")
      clear()
      enabled(yVar)<-FALSE
    } else {
      ## gray out y value unless action popup is in bivariate
      if(actionVal %in% bivariateSummaries) {
        enabled(yVar)<-TRUE
      } else {
        enabled(yVar)<-FALSE
      }
      update()                          # call update with new action
    }
  })
  
  ## three actions: clear, update, unrealize
  clear = function() {
    svalue(xVar) <- defaultMsg
    svalue(yVar) <- defaultMsg
    svalue(groupingVar) <- defaultMsg
    dispose(summaryArea)
    
    tag(obj,"xVarData") <- NULL
    tag(obj,"yVarData") <- NULL
    tag(obj,"groupingVarData") <-  NULL

    clearDropHandlers()
  }

  update = function() {
    ## get variables
    xVarData = tag(obj,"xVarData")
    yVarData = tag(obj,"yVarData")
    groupingVarData = tag(obj,"groupingVarData")
    actionVal = svalue(actionPopup)

    ## proceed?
    if(is.null(xVarData)) {
      cat("Need data in the x variable\n")
      return()
    }
    if(actionVal %in%   c("Select summary:","Clear variables") ||
       length(grep("^---", actionVal)) >0 
       ) {
      return()
    }

    ## update summaryArea
    envir = environment()               # use eval(parse... here
    assign(id(xVarData),svalue(xVarData), envir=envir)
    if(actionVal %in% univariateSummaries) { ##is.null(yVarData)
      if(is.null(groupingVarData)) {
        command = Paste(actionVal,"(",id(xVarData),")")
      } else {
        assign(id(groupingVarData),svalue(groupingVarData), envir=envir)
        command = Paste("sapply(split(",id(xVarData),",",
          id(groupingVarData),"),",actionVal,")")
      }
    } else {
      ## yVar is good, how about group
      assign(id(yVarData),svalue(yVarData), envir=envir)
      if(is.null(groupingVarData)) {
        command = Paste(actionVal,"(",
          id(xVarData),",",id(yVarData),
          ")")
      } else {
        ## bivariate summary
        assign(id(groupingVarData),svalue(groupingVarData), envir=envir)
        d = data.frame(svalue(xVarData),svalue(yVarData))
        names(d) = c(id(xVarData),id(yVarData))
        assign("d",d, envir=envir)
        command = Paste("sapply(split(d,", id(groupingVarData),"),",
          "function(x)", actionVal,"(x[,1],x[,2])",
          ")")
      }
    }
    ## cat(options("prompt"), command,"\n")
    out = capture.output(eval(parse(text=command), envir=envir))

    dispose(summaryArea)
    add(summaryArea, command,font.attr=c(style="monospace",color="blue"))
    add(summaryArea, out, font.attr=c(style="monospace"))

  }
  
  ## clear out view.col handlers
  clearDropHandlers = function(...) {
    dropHandlers = tag(obj,"dropHandlers")
    if(length(dropHandlers) > 0) {
      for(i in 1:length(dropHandlers)) {
        removehandler(dropHandlers[[i]]$view.col,dropHandlers[[i]]$id)
      }
    }
  }
  addhandlerunrealize(obj, handler = clearDropHandlers)
  
  return(obj)
}
  
