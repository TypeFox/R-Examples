## interface to variables  used by ggenericwidget
##TODO:

## * responsevar-- droplist (update with data), with dnd, close subset dialog, formula editor dialog
## * formulavar -- droplist updated with data, with dnd, also has button for foumulat editor
## * datavar -- droplist updated with dnd, when updated chnages droplist in responsevar and formulavar
## * subsetvar -- updated with dnd, has subsetEditor
## * formulaEditor -- how to do? Must have following
##   adding binary terms: +, * ....
##   grouping via ()
##   wrapping terms: I(), tsvar(), ...
##   how to integrate lm, glm, sspir usages, nlm?


##################################################
## some specific widgets for handling variables

## call this instead of the others
## gvarariables(

gvariables = function(variableType=NULL,..., toolkit=guiToolkit()) {
  if(variableType %in% c("univariate","univariatetable","bivariate","model","lattice")) {
    tmp = switch(variableType,
      "univariate" = gunivariate(...,toolkit=toolkit),
      "univariatetable" = gunivariatetable(...,toolkit=toolkit),
      "fileurl" = gfileurl(...,toolkit=toolkit),
      "bivariate" = gbivariate(...,toolkit=toolkit),
      "model" = gmodel(...,toolkit=toolkit),
      "lattice" = glattice(...,toolkit=toolkit),
      "lmer" = glmer(...,toolkit=toolkit)
      )
    return(tmp)
  } else {
    ## can't do anything here yet
    ## expects something so we give back a box
    return(ggroup(...,toolkit=toolkit))
  }
}



## For these we return a string ready to go.

## take an idwidget, wrap an argument around it, return as a container
## with value method giving "arg=val"

setClass("gAddargANY",
         representation(argument="character"),
         contains="gComponentANY"
         )

addArg = function(argument,widget,container = NULL, toolkit=guiToolkit()) {
  
  
  delete(container,widget)                 # a hack, take away parent to reparent
  label = glabel(text=Paste(argument,"= "), container=container)
  add(container,widget)

  obj = new("gAddargANY",block=container, widget=widget,
    toolkit=toolkit,
    argument=argument)
  return(obj)
}



setMethod(".svalue",
          signature(toolkit="ANY",obj="gAddargANY"),
          function(obj, toolkit, index=NULL, drop=NULL,  ...) {
            val = svalue(obj@widget)
            arg = obj@argument


            val = stripWhiteSpace(val)

            if(val == "") {
              return(NA)
            } else if(
                      is(obj@widget,"gEditListANY") ||
                      is(obj@widget,"gEditNamedListANY") ||
                      arg == "...") {
              return(val)
            } else {
              if(is.null(drop) || drop==FALSE)
                return(Paste(arg,"=",val))
              else
                return(val)
            }
          })

setReplaceMethod(".svalue",
                 signature(toolkit="ANY",obj="gAddargANY"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   svalue(obj@widget) <- value
                   return(obj)
                 })

## the most basic drop handler 
basicdrophandler = function(h,...) svalue(h$obj) <- h$dropdata
##################################################
### UNIVARIATE
## returns gedit vaRlue
setClass("gUnivariateANY",
         representation(widgets="list"),
         contains="gComponentANY"
         )

gunivariate = function(xlabel="x",container=NULL, ..., toolkit=guiToolkit()) {
  frame = gframe(text = "data", horizontal=TRUE, container=container, expand=TRUE)
  font(frame) <- c(weight="bold", size="small")
  
  xentry = gedit(text="",width=30, container=frame) # was NULL
  xarg = addArg(argument=xlabel, xentry, container=frame)

  obj = new("gUnivariateANY",block=frame, widget=frame,
    toolkit=toolkit, 
    widgets=list(xarg))
  return(obj)
}

setMethod(".svalue",
          signature(toolkit="ANY",obj="gUnivariateANY"),
          function(obj, toolkit, index=NULL, drop=NULL,  ...) {
            svalue(obj@widgets[[1]])
          })

##################################################
#######################################################################
## returns gedit value, perhaps in table()
setClass("gUnivariateTableANY",
         representation(widgets="list"),
         contains="gComponentANY"
         )

gunivariatetable = function(xlabel="x",container=NULL, ..., toolkit=guiToolkit()) {

  frame = gframe(text = "data", horizontal=TRUE, container=container,
    expand=TRUE)
  font(frame) <- c(weight="bold")
  
  xentry = gedit(text="",width=30, container=frame)
  xarg = addArg(argument=xlabel, xentry, container=frame)
  glabel("Tabulate data?",container=frame)
  doTable=gcombobox(c(TRUE,FALSE),container=frame)
  obj = new("gUnivariateTableANY",block=frame, widget=frame, toolkit=toolkit, widgets=list(x=xarg, doTable=doTable))
  return(obj)
}

setMethod(".svalue",
          signature(toolkit="ANY",obj="gUnivariateTableANY"),
          function(obj, toolkit, index=NULL, drop=NULL,  ...) {
            vals = svalue(obj@widgets[[1]], drop=TRUE)
            doTable = as.logical(svalue(obj@widgets[[2]]))
            if(doTable)
              vals = paste("table(",vals,")", sep="", collapse="")
            return(vals)
          })


#################################################
#################################################
## a File browser with a switch for wrapping in url
## returns gedit value, perhaps in table()
setClass("gFileURLANY",
         representation(widgets="list"),
         contains="gComponentANY"
         )

gfileurl = function(xlabel="file",container=NULL, ..., toolkit=guiToolkit()) {

  frame = gframe(text = "file", horizontal=TRUE, container=container,
    expand=TRUE)
  font(frame) <- c(weight="bold")
  
  xentry = gfilebrowse(text="",width=40, container=frame)
  xarg = addArg(argument=xlabel, xentry, container=frame)
  glabel("A url?",container=frame)
  doURL=gcombobox(c(FALSE,TRUE),container=frame)
  obj = new("gFileURLANY",block=frame, widget=frame, toolkit=toolkit, widgets=list(x=xarg, doURL=doURL))
  return(obj)
}

setMethod(".svalue",
          signature(toolkit="ANY",obj="gFileURLANY"),
          function(obj, toolkit, index=NULL, drop=NULL,  ...) {
            vals = svalue(obj@widgets[[1]])
            doURL = as.logical(svalue(obj@widgets[[2]])) # TRUE of FALSE
            if(doURL)
              vals = paste("url(",quoteIfNeeded(vals),")", sep="", collapse="")
            return(vals)
          })

##################################################
## bivariate
setClass("gBivariateANY",
         representation(widgets="list"),
         contains="gComponentANY"
         )

gbivariate = function(xlabel = "x", ylabel = "y", container=NULL, ...,
  toolkit=guiToolkit()) {
  
  frame = gframe(text = "data", horizontal=TRUE, container=container,
    expand=TRUE)
  font(frame) <- c(weight="bold")

  xentry = gedit(text="",width=30, container=frame)
  ##  adddroptarget(xentry, handler = basicdrophandler)
  xarg = addArg(argument=xlabel, xentry, container=frame)
  
  yentry = gedit(text="",width=30, container=frame)
  ##  adddroptarget(yentry, handler = basicdrophandler)
  yarg = addArg(argument=ylabel, yentry, container=frame)
  
  obj = new("gBivariateANY",block=frame, widget=frame,
    toolkit=toolkit,  widgets=list(xarg,yarg))
  return(obj)
}

setMethod(".svalue",
          signature(toolkit="ANY",obj="gBivariateANY"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            return(PasteWithComma(
                                  svalue(obj@widgets[[1]]), # xvalue
                                  svalue(obj@widgets[[2]])  # yvalue
                                  )
                   )
          })

setReplaceMethod(".svalue",
                 signature(toolkit="ANY",obj="gBivariateANY"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   value = rep(value,2)[1:2]
                   lapply(1:2, function(i) svalue(obj@widgets[[i]]) <- value[i])
                   return(obj)
                 })


##################################################
## model formula -- this has drag and drop stuff, and formula editing
setClass("gModelANY",
         representation(widgets="list"),
         contains="gComponentANY",
         )

gmodel = function(lattice=FALSE, container=NULL,...,toolkit=guiToolkit()) {
  ## containers
  frame = gframe(text = "data",  horizontal=FALSE, container=container,
    expand=TRUE, anchor=c(-1,1))
  font(frame) <- c(weight="bold", size="small")

  
  ## we have 4 main things: response, predictor(s), data, subset
  
  editPredictorHandler = function(h,...) {
    editFormulaDialog(data=h$action$dataEntry,
                      responsewidget=h$action$responseEntry,
                      predictorwidget=h$action$predictorEntry)
  }
  
  dataDropHandler = function(h,...) {
     val.str = h$dropdata
     ## in gWidgetsRGtk2 dnd is a little messy
#JV     if(getOption("guiToolkit") != "RGtk2")
       svalue(h$obj) <- val.str ## gets on enter handler?
     names = try(getNamesofObject(val.str))
     if(!inherits(names,"try-error") &&
        !is.null(names)) {
       ## set dropdown values:
       responseEntry[] <- names
       predictorEntry[] <- names
       if(lattice)
         conditionEntry[] <- names
       ## what else?
     }
   }
    
  dataEnterHandler = function(h,...) {
     val.str = svalue(h$obj)
     names = getNamesofObject(val.str)
     if(!is.null(names)) {
       ## set dropdown values:
       responseEntry[] <- names
       predictorEntry[] <- names
       if(lattice)
         conditionEntry[] <- names
       ## what else?
     }
   }

  ## subsetDropHandler = basicdrophandler
  editSubsetHandler = function(h,...) {
    editSubsetDialog(data=h$action$dataEntry,
                     widget=h$action$subsetEntry)
  }
  editConditionHandler = function(h,...) {
    editConditionDialog(data=h$action$dataEntry,
                        widget = h$action$conditionEntry)
  }
  
  ## Build up the widgets now
  variableNames = c("",getNamesofObject())
  
  
  tbl = glayout(container=frame, anchor=c(-1,1), width=500,height=100)
  tbl[1,1] <- (responseEntry = gcombobox(variableNames,width=40,editable=TRUE, container =tbl))
  tbl[2,1] <- "response"
  size(responseEntry) <- c(200, -1)
  
  ## add ~
  tbl[1,2] <- " ~ "

  tbl[1,3] <- (predictorEntry = gcombobox(variableNames,width=50, editable=TRUE, container =tbl))
  tbl[2,3] <- "predictor(s)"
  size(predictorEntry) <- c(200, -1)
  
  tbl[1,4] <- (predictorEdit = gbutton("edit",container=tbl))
  
  
  
  
  ## lattice has a conditioning variable
  ## define this even if not needed
  conditionEntry = NULL
  if(lattice) {
      tbl[3,2] <- " | "
      tbl[3,3] <- (conditionEntry = gcombobox(variableNames, editable=TRUE,
                       container=tbl))
      tbl[4,3] <- "conditioning variable(s)"
    
    
    tbl[3,4] <- (conditionEdit = gbutton("edit",container=tbl))
  }
  visible(tbl) <- TRUE



               
  ## visual separation
  gseparator(container=frame)
  
  ## Now for data frame
  tbl = glayout(container=frame, anchor = c(-1,1), width=500,height=100)
  tbl[1,1, anchor=c(1,0)] = "data="
  tbl[1,2, anchor=c(-1,0)] = (dataEntry <-  gedit("", container=tbl))

  addhandlerchanged(dataEntry, handler = dataEnterHandler)
  adddroptarget(dataEntry, handler = dataDropHandler)
  
  ## subset has extra button for editing
  tbl[2,1, anchor=c(1,0)] = "subset="
  tbl[2,2, anchor=c(-1,0)] = (subsetEntry <- gedit("", container =tbl))
  
  tbl[2,3] = (subsetEdit <- gbutton("edit",container=tbl)
              )
  visible(tbl) <- TRUE                # for RGtk2
  
  ## add handlers --after dataEntry is defined
  addHandlerChanged(predictorEdit,
                    handler=editPredictorHandler,
                    action=list(dataEntry=dataEntry,
                      responseEntry=responseEntry,
                      predictorEntry=predictorEntry))
  

  addHandlerChanged(subsetEdit,
                    handler=editSubsetHandler,
                    action=list(
                      dataEntry=dataEntry,
                      subsetEntry=subsetEntry
                      )
                    )
  
  if(lattice) {
    addHandlerChanged(conditionEdit,
                      handler = editConditionHandler,
                      action = list(dataEntry=dataEntry,
                        conditionEntry =conditionEntry))
  }
  
  obj = new("gModelANY",block=frame,widget=frame,
    toolkit=toolkit,
    widgets =
    list(response  = responseEntry,
         predictor = predictorEntry,
         data      = dataEntry,
         subset    = subsetEntry,
         lattice   = conditionEntry
         )
    )
  return(obj)
}

setMethod(".svalue",
          signature(toolkit="ANY",obj="gModelANY"),
          function(obj, toolkit,index=NULL, drop=NULL, ...) {
            lst = lapply(obj@widgets, function(i) {
              if(!is.null(i)) svalue(i)
            })
            names(lst) = c("response","predictor","data","subset","lattice")[1:length(lst)]
            ## return a string of type r ~ p, data=..., subset=...
            str = with(lst, Paste(response," ~ ",predictor))
            if(!is.null(lst$lattice) && lst$lattice != "")
              str = Paste(str, " | ", lst$lattice)
            if(lst$data != "") str = PasteWithComma(str,Paste("data=",lst$data))
            if(lst$subset != "") str = PasteWithComma(str,Paste("subset=",lst$subset))
            
            return(str)
          })

setReplaceMethod(".svalue",
                 signature(toolkit="ANY",obj="gModelANY"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   if(is.list(value)) unlist(value)
                   lapply(1:4, function(i) svalue(obj@widjets[[i]]) <- value[i])
                   return(obj)
                 })


##################################################
## interface for lattice graphics
glattice = function(..., toolkit=guiToolkit()) return(gmodel(lattice = TRUE, ...,toolkit=toolkit))

##################################################
setClass("gLmerANY",
         representation(formulaEntry="guiWidget",
                        dataEntry="guiWidget"),
         contains="gComponentANY")


## interface for linear mixed effects models models
glmer = function(container=NULL, ..., toolkit=guiToolkit()) {
  
  frame = gframe(text = "data",  horizontal=FALSE, container=container)
  font(frame) <- c(weight="bold")
  
  tbl = glayout(container=frame)
  tbl[1,1] = "formula"
  tbl[1,2] = (formulaEntry <- gedit("", width=40, container =tbl))
  tbl[2,1] = "data="
  tbl[2,2] = (dataEntry <- gedit("", container =tbl))
  visible(tbl) <- TRUE
  
  obj = new("gLmerANY",block=frame,widget=frame,toolkit=toolkit,
    dataEntry = dataEntry, formulaEntry = formulaEntry)

  return(obj)
}

setMethod(".svalue",
          signature(toolkit="ANY",obj="gLmerANY"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            formulaValue = svalue(obj@formulaEntry)
            dataValue = svalue(obj@dataEntry)
            
            if (formulaValue == "")
              return("")
            string = formulaValue
            if(dataValue != "")
              string = Paste(string, ", data=",dataValue)
            
            return(string)
          })





## some special edit fields
## this is used to get "..." into generic widget
setClass("gEditListANY",
         contains="gComponentANY")

geditlist = function(...,toolkit=guiToolkit()) {

  edit = gedit(...,toolkit=toolkit)
  obj = new("gEditListANY", block=edit@widget@block, widget=edit@widget@widget, toolkit=toolkit)
  invisible(obj)

}


setClass("gEditNamedListANY",
         contains="gComponentANY")

geditnamedlist = function(..., toolkit=guiToolkit()) {
  edit = gedit(...,toolkit=toolkit)
  obj = new("gEditNamedListANY",
    block=edit@widget@block, widget=edit@widget@widget, toolkit=toolkit)
  invisible(obj)
}


##################################################
##
## two dialogs for editing


## sample edit Formula dialog, Smilar to one in Splus, but layout is
## not as nice, not quite as feature rich

editFormulaDialog = function(
  data=NULL,                            # within these variables
  responsewidget = NULL,                # response value
  predictorwidget = NULL                # predictor value
  ) {
  ## actually get data, not just a string
  if(is(data,"guiWidget") || is(data,"gComponentANY")) {
    data = svalue(data)
    dataName = deparse(substitute(data))
  }
  if(is.character(data) && length(data) == 1) {
    dataName = data 
    data = getObjectFromString(data)
  }

  if(is.na(data) || is.null(data)) {
    warning("Can't find data set")
    return(NA)
  }

  ## coerce data if possible
  if(!is.data.frame(data)) {
    tmp = try(as.data.frame(data), silent=TRUE)
    if(inherits(tmp,"try-error")) {
      warning("gtable shows data frames of vectors")
      return(NA)
    }
    data = tmp
  }
  varNames = names(data)

  ## define key widgets
  
  
  ## Set up the window
  ## main window
  win = gwindow("Edit model formula values")
  
  group = ggroup(horizontal=FALSE, container=win)
  
  datagroup = ggroup(container=group)
  size(datagroup) <- c(300,200)
  glabel("Dataset: ", container=datagroup)
  tmp = glabel(dataName, container=datagroup);
  font(tmp) <- c(weight="bold")
  addSpace(datagroup, 10)

  variables = gtable(varNames,multiple=TRUE, container =datagroup, expand=TRUE)

  gseparator(container =group)
  ## buttons

  
  buttonGroup = ggroup(container=group)
  glabel("Actions:", container=buttonGroup)
  tbl = glayout(container =buttonGroup, expand=TRUE)
  tbl[1,1,anchor=c(-1,0)] = (addresponse <- gbutton("Response",container =tbl))
  tbl[1,2,anchor=c(-1,0)] = (addterm <-  gbutton("+ (main effect)",container =tbl))
  tbl[2,1,anchor=c(-1,0)] = (addwithin <- gbutton(": (interaction)",container =tbl))
  tbl[2,2,anchor=c(-1,0)] = (addinteraction <- gbutton("* (main + interaction)",container =tbl))
  tbl[3,1,anchor=c(-1,0)] = (addsecondpowers <- gbutton("^2 (second-order)",container =tbl))
  tbl[3,2,anchor=c(-1,0)] = (subtractintercept <-  gbutton("remove intercept",container =tbl))
  visible(tbl) <- TRUE

  gseparator(container =group)
  tbl = glayout(container=group)

  response = gedit(svalue(responsewidget),container =tbl)
  predictor = gedit(svalue(predictorwidget), container =tbl)

  tbl[1,1] <- response
  tbl[1,2] <- glabel(" ~ ",container =tbl)
  tbl[1,3:6] <- predictor
  tbl[2,1] <- "response"
  tbl[2,3:6] <- "predictor formula"
  visible(tbl) <- TRUE


  

  buttonbox = ggroup(container=group)
  addSpring(buttonbox)
  okbutton = gbutton("ok", container =buttonbox)
  addSpace(buttonbox,15)  
  clearbutton = gbutton("clear", container =buttonbox)
  cancelbutton = gbutton("cancel", container =buttonbox)

  ## Now add handlers
  addhandlerclicked(addresponse, handler=function(h,...) {
    vals = svalue(variables)
    if(!is.null(vals)) {
      svalue(response) <- vals[1]
    }
  })
  addhandlerclicked(addinteraction,handler=function(h,...) {
    vars = svalue(variables)
    if(!is.null(vars)) {
      oldval = svalue(predictor)
      if(!is.null(oldval) && oldval !="")
        oldval = Paste(oldval, " + ")
      else
        oldval = ""
      svalue(predictor) <- Paste(oldval, paste(vars, sep="", collapse=" * "))
    }
  })
  addhandlerclicked(addterm,handler=function(h,...) {
    vars = svalue(variables)
    
    if(!is.null(vars)) {
      oldval = svalue(predictor)
      if(!is.null(oldval) && oldval !="")
        oldval = Paste(oldval, " + ")
      else
        oldval = ""
      svalue(predictor) <- Paste(oldval, paste(vars, sep="", collapse=" + "))
    }
  })
  addhandlerclicked(addsecondpowers,handler=function(h,...) {
    vars = svalue(variables)
    if(!is.null(vars)) {
      oldval = svalue(predictor)
      if(!is.null(oldval) && oldval !="")
        oldval = Paste(oldval, " + ")
      else
        oldval = ""
      svalue(predictor) <- Paste(oldval,
                                 " (",
                                 paste(vars, sep="", collapse=" + "),
                                 ")^2")
    }
  })
  addhandlerclicked(addwithin,handler=function(h,...) {
    vars = svalue(variables)
    if(!is.null(vars)) {
      oldval = svalue(predictor)
      if(!is.null(oldval) && oldval !="")
        oldval = Paste(oldval, " + ")
      else
        oldval = ""
      svalue(predictor) <- Paste(oldval,
                                 paste(vars, sep="", collapse=":")
                                 )
    }

  })
  addhandlerclicked(subtractintercept,handler=function(h,...) {
    svalue(predictor) <- Paste(svalue(predictor), " -1")
  })
  addhandlerclicked(okbutton, handler = function(h,...) {
    svalue(responsewidget) <- svalue(response)
    svalue(predictorwidget) <- svalue(predictor)
    dispose(win)
  })
  addhandlerclicked(clearbutton, handler=function(h,...) {
    svalue(response) <- ""
    svalue(predictor) <- ""
  })
  addhandlerclicked(cancelbutton,handler=function(h,...) {
    dispose(win)
  })
}

editSubsetDialog = function(
  data=NULL,
  widget = NULL                         # what to write to, start with
  ) {
  
  ## get data values
  if(is(data,"guiComponent") || is(data,"gComponentANY")) {
    data = svalue(data)
  }

  
  if(is.character(data)) {
    if(data == "") {
      warning("A data set needs to be set")
      return()
    }
    dataName = data
    data = svalue(data)
  } else {
    dataName = deparse(substitute(data))
  }

  if(!is.data.frame(data)) {
    tmp = try(as.data.frame(data), silent=TRUE)
    if(inherits(tmp,"try-error")) {
      warning("gtable shows data frames of vectors")
      return(NA)
    }
    data = tmp
  }
  varNames = names(data)
  
  
  ## main widgets
  
                                        #  preview = gtable(head(data))
  
  
  ## layout
  
  win = gwindow("Edit subset value")                    # main window
  group = ggroup(horizontal=FALSE, container=win)
  glabel(Paste("Data set:", dataName,""), container=group)

  tbl = glayout(container=group)
  tbl[1,1] <- (andOrPopup = gcombobox(c("","&","|"), container =tbl))
  tbl[1,2] <- (notPopup  = gcombobox(c("","!"), container =tbl))
  tbl[1,3] <- (var1 = gcombobox(varNames, editable=TRUE, container =tbl))
  tbl[1,4] <- (logicalPopup = gcombobox(c("","<","<=","==","!=",">=",">","%in%"), container =tbl))
  tbl[1,5] <- (var2 = gcombobox(varNames, editable = TRUE, container =tbl))

  tbl[2,1] <- "join with"
  tbl[2,2] <- "negate"
  tbl[2,3] <- ""
  tbl[2,4] <- "compare with"
  tbl[2,5] <- ""
  visible(tbl) <- TRUE
  
  buttonGroup = ggroup(container=group)
  addSpring(buttonGroup)

  addButton = gbutton("add", container =buttonGroup)
  addSpace(buttonGroup,10)
  clearButton =   gbutton("clear",container =buttonGroup)

  outputGroup = ggroup(container=group)
  glabel("Subset by",container=outputGroup)
  output = glabel("", container =group)
  font(output) <- c(weight="bold")  # make bold
  status = glabel("",container =group)
  
  gseparator(container=group)
  spacingGroup = ggroup(container=group)
  glabel("Preview:", container=spacingGroup)
  addSpring(spacingGroup)

  preview = gtable(data, container =group, expand=TRUE)
  enabled(preview) <- FALSE
  
  buttonGroup = ggroup(container=group)
  addSpring(buttonGroup)
  okButton = gbutton("ok", container =buttonGroup)
  addSpace(buttonGroup,10)
  cancelButton = gbutton("cancel", container =buttonGroup)

  ## now give some actions
  addHandlerAction = function(h,...) {
    str = Paste(
      svalue(var1),
      svalue(logicalPopup),
      svalue(var2)
      )
    if(!is.null(svalue(notPopup)) && length(svalue(notPopup)) > 0 && svalue(notPopup) == "!")
      str = Paste("!( ",str," )")
    ## error check
    tmp = try(rpel(Paste("subset(data,subset=",str,")"), envir=environment(NULL)))
    if(! inherits(tmp,"try-error")) {
      combine = svalue(andOrPopup)
      if(is.null(combine)) combine = ""
      previous = svalue(output)
      if(length(previous) == 0 || nchar(previous) == 0) {
        combine=""; previous=""
      }
      if(nchar(previous) > 0 && combine == "") {
        svalue(status) <-"You need an operator to combine terms"
        return(FALSE)
      }
        
      svalue(output) <- Paste(previous," ",combine," ",str)
      svalue(status) <- ""
      ## update preview
      subsetVal = svalue(output)

      tmpVals = rpel(Paste("head(subset(data,subset=",subsetVal,"))"),
        envir=environment(NULL))
      preview[,] = tmpVals
    } else {
      svalue(status) <- Paste(str," failed")
    }
  }

  ## when logical is %in$ and leftside if factor, we can do some  business
  addhandlerchanged(logicalPopup,handler = function(h,...) {
    ## check that all conditions are met
    if(length(svalue(logicalPopup)) > 0 && svalue(logicalPopup) == "%in%") {
      varName = svalue(var1)
      var = getObjectFromString(Paste(dataName,"$",varName))
      if(is.factor(var)) {
        varLevels = levels(var)
        selectVectorEntriesDialog(varLevels, var2)
      }
    }
  }
                    )
  
  addhandlerclicked(clearButton, handler = function(h,...) {
    svalue(output) <- ""
    svalue(status) <- ""
    preview[] <- head(data)
  })
  
  addhandlerclicked(addButton, handler=addHandlerAction)
  addhandlerclicked(okButton, handler = function(h,...) {
    outputValue = svalue(output)
    ## strip leading spaces
    svalue(widget) <- sub("^ +","",outputValue)
    dispose(win)
  })
  addhandlerclicked(cancelButton, handler=function(h,...) {
    dispose(win)
  })
}



## add terms 
editConditionDialog <- function(data,widget) {
    data <- svalue(data)
  if(is.na(data) || is.null(data)) {
    warning("empty dataset")
    return()
  } else if(data == "") {
    varNames = ls(envir=.GlobalEnv)
  } else {
    ## get data values
    if(is(data,"gComponentANY")) {
      data = svalue(data)
    }
    if(is.character(data)) {
      dataName = data
      data = rpel(Paste("get(\"",data,"\")", sep=" "))
    } else {
      dataName = deparse(substitute(data))
    }
    
    if(!is.data.frame(data)) {
      tmp = try(as.data.frame(data), silent=TRUE)
      if(inherits(tmp,"try-error")) {
        warning("gtable shows data frames of vectors")
        return(NA)
      }
      data = tmp
    }
    varNames = names(data)
  }

  win = gwindow("Select variables for conditioning variable",
    visible=TRUE)
  group = ggroup(horizontal=FALSE, container=win)

  vars = gtable(varNames, multiple=TRUE, handler=function(h,...) {
    values = svalue(h$obj)
    if(length(values) > 0) {
      out = paste(values, collapse=" + ")
      svalue(widget) <- out
    }
    dispose(win)
  },
    container =group, expand=TRUE)

  buttonGroup = ggroup(container=group)
  addSpring(buttonGroup)
  gbutton("ok",container=buttonGroup, action=vars,
          handler=function(h,...) {
            values = svalue(h$action)
            if(length(values) > 0) {
              out = paste(values, collapse=" + ")
              svalue(widget) <- out
            }
            dispose(win)
          })
  gbutton("close",container=buttonGroup, handler = function(h,...) dispose(win))
  status = gstatusbar("Select one or more variables to condition by.",
    container=group)

}

## select asks user for which columns to consider
editSelectDialog = function(data,widget) {
  
    ## get data values
  if(is(data,"gComponentANY")) {
    data = svalue(data)
    if(data == "") {
      warning("A data set needs to be set")
      return()
    }
  }
  if(is.character(data)) {
    if(data == "") {
      warning("A data set needs to be set")
      return()
    }
    dataName = data
    data = rpel(Paste("get(\"",data,"\")", sep=" "))
  } else {
    dataName = deparse(substitute(data))
  }

  if(!is.data.frame(data)) {
    tmp = try(as.data.frame(data), silent=TRUE)
    if(inherits(tmp,"try-error")) {
      warning("gtable shows data frames of vectors")
      return(NA)
    }
    data = tmp
  }
  varNames = names(data)

  win <- gwindow("Select variables", visible=TRUE)

  group <- ggroup(horizontal=FALSE,container=win)
  frame <- gframe("Select all the desired variables",
    container=group, expand=TRUE)
  font(frame) <- c(weight="bold")

  varList <- gtable(varNames, multiple=TRUE, container=frame, expand=TRUE)


  buttonGroup = ggroup(container=group)
  addSpring(buttonGroup)
  submitButton = gbutton("submit",container=buttonGroup,
    handler = function(h,...) {
    values = svalue(varList)
    ## wrap it up
    str = Paste("c(",paste(sapply(values,function(str) paste("'",str,"'",sep="")),sep=" ",collapse=", "),")")
    svalue(widget) <-  str
    ## lclose dialog
    dispose(win)
  })
}
  
  
  ## put values into a vector with paste
  selectVectorEntriesDialog <- function(vals, widget) {
    
    
    
    win=gwindow("Select values", visible=T)
    group = ggroup(horizontal = FALSE, container=win)
    varNamesList = gtable(vals, multiple = TRUE, group, expand=TRUE)

    buttonGroup = ggroup(container=group)
    addSpring(buttonGroup)
    okButton = gbutton("ok",container=buttonGroup)
    cancelButton = gbutton("cancel", container=buttonGroup)

    addhandlerclicked(cancelButton, handler = function(h,...) {
      dispose(win)
    })
    addhandlerclicked(okButton, handler = function(h,...) {
      varNames = svalue(varNamesList)
      varNamesAsString = Paste("c(\"",paste(varNames,sep="'",collapse="\",\""),"\")")
      svalue(widget) <- varNamesAsString
      ## tidy up
      dispose(win)
    })
  }
