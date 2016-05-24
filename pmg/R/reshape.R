## gui for reshape package

## startup stuff
## options("guiToolkit"="RGtk2")

## provides two functions: meltGUI and castGUI.
## How to improve these?

pmg.meltGUI = function(container = pmgWC$new("Melt a data frame")) {
    

  ## Helper functions
  idVarsDF = function(df) {
    d = data.frame(
      "variable name"=names(df),
      "variable type"=sapply(df, function(i) class(i)[1]),
      stringsAsFactors=FALSE)
    return(d)
  }

  guessIdVars = function(d) {
   ## d has been through idVarsDF
    which(d[,2] %in% c("factor","integer"))
  }

  getMeltedObject = function() {

    theName = svalue(theDF)
    theVals = try(get(theName),silent=TRUE)
    ## verify data frame is good
    if(inherits(theVals,"try-error") ||
       !(is.data.frame(theVals) || is.matrix(theVals))) {
      gmessage(paste(theName,"does not refer to a data frame or matrix",sep=" ",collapse=" "))
      return()
    }

    varNames = svalue(theIDVars)
    theList = list(data=theVals,
      id.var = varNames,
      variable_name = svalue(theVarName),
      preserve.na = svalue(thePreserveNA)
      )

    newMelt = do.call("melt",theList)

    return(newMelt)
  }
 


  ## GUI
  g = ggroup(horizontal=FALSE, container=container, raise.on.dragmotion = TRUE)
  glabel("Melt a data frame", container=g)

  layout = glayout(container=g)

  ## theDF -- for data frame name
  layout[1,1] = glabel("data frame:")
  theDF = gedit()
  layout[1,2] = theDF

  ## theIDVars -- for selecting id variables
  layout[2,1] = glabel("id.var")
  dummyDF = data.frame("variable Name" = "","variable type"="",stringsAsFactors=FALSE)
  theIDVars = gtable(dummyDF, multiple=TRUE)
  
  size(theIDVars) <- c(300,200)
  layout[2,2] = theIDVars

  ## variable Name
  layout[3,1] = glabel("variable name")
  theVarName = gedit("variable")
  layout[3,2] = theVarName

  ## preserve.na
  layout[4,1] = glabel("preserve.na")
  thePreserveNA = gdroplist(c(TRUE,FALSE))
  layout[4,2] = thePreserveNA

  ## update
  theUpdate = gbutton("update")
  layout[5,2] = theUpdate
  
  visible(layout) <- TRUE

  
  ## preview area
  add(g,gseparator(horizontal=TRUE))
  previewGroup = gexpandgroup("Preview",container=g)
  thePreview  =  glabel("")
  add(previewGroup, thePreview, expand=TRUE) # use delete/add to chnage
  visible(previewGroup) <- TRUE              # open as default
  
  ## saveAs area
  add(g,gseparator(horizontal=TRUE))
  saveAsGroup = ggroup(container=g)
  saveAsButton = gbutton("Save output as:",container=saveAsGroup)
  saveAs = gedit("",container=saveAsGroup)
  enabled(saveAsGroup) <- FALSE
  ## End of layout


  ## data frame
  addhandlerchanged(theDF,handler = function(h,...) {
    theName = svalue(theDF)

    ## trust but verify
    theVals = try(get(theName))
    if(inherits(theVals,"try-error") ||
       !(is.data.frame(theVals) || is.matrix(theVals))) {
      gmessage(paste(theName,"does not refer to a data frame or matrix",sep=" ",collapse=" "))
    } else {

      ## updateIDVars area
      tmp <- idVarsDF(theVals)
      theIDVars[,] <- tmp
      svalue(theIDVars) <- guessIdVars(tmp)
    }
  })

  ## updatebutton
  addhandlerchanged(theUpdate, handler = function(h,...) {
    newMelt = getMeltedObject()

    ## update preview
    delete(previewGroup, thePreview)
    thePreview <<- gtable(head(newMelt, n = 15))
    
    add(previewGroup, thePreview, expand=TRUE)
    enabled(thePreview) <- FALSE
    
    ## make output area visible
    enabled(saveAsGroup) <- TRUE
    
  })

  ## saveAs
  saveHandler = handler = function(h,...) {
    newMelt = getMeltedObject()

    varName = svalue(saveAs)
    ## check
    if(varName == "") {
      gmessage("No variable name specified")
      return()
    }
    if(exists(varName, envir=.GlobalEnv)) {
      val = gconfirm(paste("Overwrite value for",varName,"?",sep=" "))
      if(val == FALSE)
        return()
    }
                         
    assign_global(varName, newMelt)

    enabled(saveAsGroup) <- FALSE
  }

  ## clicking button, or enter after editing variable will do it.
  addhandlerchanged(saveAsButton, handler=saveHandler)
  addhandlerchanged(saveAs, handler=saveHandler)

  invisible()
}


pmg.castGUI = function(container=pmgWC$new("Cast data")) {

  
  g = ggroup(horizontal=FALSE, container=container, raise.on.dragmotion = TRUE)

  theData = gedit("", width=75)

  theVariables = gtable(data.frame(ID.vars="", stringsAsFactors=FALSE))
  adddropsource(theVariables, handler=function(h,...) svalue(theVariables))

  ## formula
  defColFormText = "Drop column variable(s) here"
  defRowFormText = "Drop row variable(s) here"
  colFormula = glabel(defColFormText, editable=TRUE)
  font(colFormula) <- c(style="bold")
  rowFormula = glabel(defRowFormText, editable=TRUE)
  font(rowFormula) <- c(style="bold")  

  aggregateFuns = c("length","mean","median","IQR","sd","range","summary")
  theAggregateFun = gdroplist(aggregateFuns, editable=TRUE)

  ## Should have "TRUE" here as well, but get wierd condense error
  defMarginVals = c("FALSE","TRUE","grand_col","grand_row")
  theMargins = gdroplist(defMarginVals)

  theSubset = gedit("", width=75)

  possDotsVals = c("","na.rm = TRUE")
  theDots = gdroplist(possDotsVals, editable=TRUE, width=75)

  clearFormulaButton = gbutton("clear")
  editSubsetButton = gbutton("edit")
  updateButton = gbutton("update")
  
  ## the layout
  glabel("Cast a melted data set", container=g)
  
  layout = glayout(container=g)

  layout[1,1] = glabel("data:")
  layout[1,2] = theData

  layout[2,1] = glabel("variables:")
  layout[2,2] = theVariables
      
  layout[3,1] = glabel("formula:")
  layout[3,2] = colFormula
  layout[3,3] = clearFormulaButton

  layout[4,2] = glabel(" ~ ")
  layout[5,2] = rowFormula
  
  layout[6,1] = glabel("fun.aggregate:")
  layout[6,2] = theAggregateFun

  layout[7,1] = glabel("margins:")
  layout[7,2] = theMargins
  
  layout[8,1] = glabel("subset:")
  layout[8,2] = theSubset
  layout[8,3] = editSubsetButton

  layout[9,1] = glabel("...")
  layout[9,2] = theDots

  layout[10,2] = updateButton
  
  visible(layout) <- TRUE

  ## preview
  add(g,gseparator(horizontal=TRUE))
  previewGroup = gexpandgroup("Preview",container=g)
  thePreview  =  glabel("")
  add(previewGroup, thePreview, expand=TRUE) # use delete/add to chnage
  visible(previewGroup) <- TRUE              # open as default

  ## saveAs area
  add(g,gseparator(horizontal=TRUE))
  saveAsGroup = ggroup(container=g)
  saveAsButton = gbutton("Save output as:",container=saveAsGroup)
  saveAs = gedit("",container=saveAsGroup)

  ##################################################
  ## helper functions
  getCast = function() {
    ## gather pieces and call cast. Return FALSE if there is an error

    ## get the data set
    theName = svalue(theData)
    theVals = try(get(theName),silent=TRUE)
    if(inherits(theVals,"try-error") || !is.data.frame(theVals)) {
      msg = paste(theName,"does not refer to a data frame or matrix",sep=" ",collapse="")
      return(list(value=msg, flag=FALSE))
    }

    ## get ready
    ## margins
    marVal<-svalue(theMargins)
    if(is.null(marVal)) marVal = FALSE  # no margins
    if(marVal %in% c("FALSE","TRUE")) marVal <- as.logical(marVal)
    
    ## subset is a character, need to get logical
    subsetVal = svalue(theSubset)
    if(subsetVal == "")
      subsetVal = TRUE
    else
      subsetVal = eval(parse(text=subsetVal), envir=theVals)

    ## get the formula
    ## don't do this if not set yet
    if(svalue(colFormula) == defColFormText ||
       svalue(rowFormula) == defRowFormText) {
      cat("Drop more variables into formula\n")
      return(list(value="",flag=NULL))
    }
    
    theFormula = as.formula(paste(
      svalue(colFormula),
      "~",
      svalue(rowFormula),
      sep="",collapse="")
      )

    
    theArgs = list(
      data=theVals,
      formula = theFormula,
      "fun.aggregate" = svalue(theAggregateFun), 
      margins = marVal,
      subset = subsetVal
      )

    ## the dots
    if((theDotsVal <- svalue(theDots)) != "") {
      ## need to split on "="
      tmp = splitAndStrip(theDotsVal,"=")
      ## assign
      theDotsValue = eval(parse(text=tmp[2]),envir=.GlobalEnv)
      theArgs[[tmp[1]]] <- theDotsValue
    }

    
    ## this is so errors do show up
    theCast = do.call("cast",theArgs)
    
#    theCast = try(do.call("cast",theArgs), silent=TRUE)
    if(inherits(theCast,"try-error")) {
      ## error
      return(list(value=theCast,flag=FALSE))
    } else {
      return(list(value=theCast,flag=TRUE))
    }
  }
  
  cleanPreview = function() {
    delete(previewGroup, thePreview)
    thePreview <<- glabel("")
    add(previewGroup, thePreview, expand=TRUE)
##    enabled(thePreview) <- FALSE
  }    

  
  ##################################################
  ## handlers

  addhandlerchanged(theData, handler=function(h,...) {
    theName = svalue(theData)
    theVals = try(get(theName),silent=TRUE)
    ## verify data frame is good

    ## CHECK::HOW TO CHECK IF IS A MELTED OBJECT?
    if(inherits(theVals,"try-error") || !is.data.frame(theVals)) {
      gmessage(paste(theName,"does not refer to a data frame or matrix",sep=" ",collapse=" "))
      return()
    }

    ## otherwise, this is good
    ## add to variable list
    ID.vars = c(".",sort(names(theVals)),"...")
    ID.vars = ID.vars[ID.vars != "value"]

    
    theVariables[,] <- data.frame(ID.vars = ID.vars, stringsAsFactors=FALSE)


    ## clearout values
    svalue(colFormula) <- defColFormText
    svalue(rowFormula) <- defRowFormText
    svalue(theAggregateFun, index=TRUE) <-1
    ## add variables to marings
    theMargins[]<-c(defMarginVals, rev(rev(names(theVals))[-(1:2)]))
    

    svalue(theSubset) <- ""
    cleanPreview()
  })


  formulaDropHandler = function(h,...) {
    curText = svalue(h$obj)
    if(curText == defColFormText || curText == defRowFormText) {
      svalue(h$obj) <- h$dropdata
    } else {
      svalue(h$obj) <- paste(curText,h$dropdata, sep=" + ",collapse="")
    }

    updateHandler(list())
  }

  adddroptarget(colFormula, handler=formulaDropHandler)
  adddroptarget(rowFormula, handler=formulaDropHandler)

  updateHandler = function(h,...) {
    theCast = getCast()

    if(is.null(theCast$flag)) return(FALSE)
    
    if(theCast$flag == FALSE) {
      gmessage(theCast$value)
      return()
    }
    theCast = theCast$value             
    theCast = as.data.frame(theCast)    # chop off class cast_df gives gdf fits

    ## now update preview
    if(is.data.frame(theCast)) {
      delete(previewGroup, thePreview)
      thePreview <<- gtable(head(theCast, n = 15))
      
      add(previewGroup, thePreview, expand=TRUE)
##      enabled(thePreview) <- FALSE
    } else {
      ## something more complicated
      cat("DEBUG: something more complicated\n")
    }
    ## allow saving
    enabled(saveAsGroup) <- TRUE
  }

  ## active, these things
  addhandlerchanged(colFormula, handler = updateHandler) 
  addhandlerchanged(rowFormula, handler = updateHandler)
  addhandlerchanged(theAggregateFun, handler = updateHandler)
  addhandlerchanged(theMargins, handler = updateHandler)
  addhandlerchanged(theSubset, handler = updateHandler)
  addhandlerchanged(theDots, handler = updateHandler) 
  addhandlerchanged(updateButton, handler = updateHandler) 
  
  ## clear
  clearFormulaHandler = function(h,...) {
    svalue(colFormula) <- defColFormText
    svalue(rowFormula) <- defRowFormText

    delete(previewGroup, thePreview)
    thePreview <<- glabel("")
      
    add(previewGroup, thePreview, expand=TRUE)
##    enabled(thePreview) <- FALSE
    
  }
  
  addhandlerchanged(clearFormulaButton, handler=clearFormulaHandler)

  editSubsetHandler = function(h,...) {
    
    ## we need to have theData set properly
    theName = svalue(theData)
    theVals = try(get(theName),silent=TRUE)
    ## verify data frame is good

    ## CHECK::HOW TO CHECK IF IS A MELTED OBJECT?
    if(inherits(theVals,"try-error") || !is.data.frame(theVals)) {
      gmessage("first set a data value")
      return()
    }

    ## this is exported from gWidgets
    editSubsetDialog(theName, widget=theSubset)
    
  }
  addhandlerchanged(editSubsetButton, handler=editSubsetHandler)

  ## Save as
  saveAsHandler = function(h,...) {
    theVals = getCast()
    
    if(theVals$flag == FALSE) {
      gmessage(paste("Can't save, an error",
                     theVals$value,
                     sep="",collapse=""))
      return()
    }

    varName = svalue(saveAs)
    ## check
    if(varName == "") {
      gmessage("No variable name specified")
      return()
    }
    if(exists(varName, envir=.GlobalEnv)) {
      val = gconfirm(paste("Overwrite value for",varName,"?",sep=" "))
      if(val == FALSE)
        return()
    }
    
    assign_global(varName, theVals$value)
    enabled(saveAsGroup) <- FALSE

  }
  addhandlerchanged(saveAsButton, handler=saveAsHandler)
  addhandlerchanged(saveAs, handler=saveAsHandler)

  invisible()
}

### Helpers
splitAndStrip = function(x, pat) {
  tmp = unlist(strsplit(x,pat, perl=TRUE))
  sub('\\s+$', '', tmp, perl = TRUE) # trim white space
  sub('^\\s+', '', tmp, perl = TRUE) 

  return(tmp)
}
