## Simple Dynamic Models dialog
## output has regular model methods defined for it except update()
dModelsDialog = function() {
  ## list of avail. models
  models = list(
    "linear regression"= "gui.lm",
    "robust regresion" = "gui.rlm",
    "-----" = "-----",
    "Analysis of variance"="gui.aov"
    )

  ## main groups
  win = pmgWC$new("Models",visible=TRUE)
  size(win) <- c(600, 400)
  gp = ggroup(horizontal=FALSE, container=win, raise.on.dragmotion = TRUE)
  popupGroup = ggroup(container=gp)
  addSpring(popupGroup)
  testPopup = gdroplist(c("",names(models)), container=popupGroup)
  gseparator(container=gp)
  testWindow = ggroup(container=gp)

  obj = list(ref=testWindow)
  class(obj) = c("gDynamicModelDialog")
  
  dialogList = list()
  dialogList[["FirstOne"]] = glabel("Select a model from popup")
  
  add(testWindow,dialogList[["FirstOne"]], expand=TRUE)
  tag(testWindow, "dialogList") <- dialogList
  tag(testWindow,"currentTest") <- dialogList[["FirstOne"]]
  
  addhandlerchanged(testPopup, handler = function(h,...) {
    popupValue = svalue(testPopup)
    if(!is.empty(popupValue) || popupValue != "-----") {
      delete(testWindow,tag(testWindow,"currentTest"))
      dialogList = tag(testWindow, "dialogList")
      if(is.null(dialogList[[popupValue]])) {
        dialogList[[popupValue]] <- do.call(models[[popupValue]],list())
        tag(testWindow, "dialogList") <-  dialogList
      }
      add(testWindow,dialogList[[popupValue]]$ref, expand=TRUE)
      tag(testWindow,"currentTest") <- dialogList[[popupValue]]$ref
  }
  })

  ## tidy up on uneralize. Need to get afresh
  addhandlerdestroy(win, handler = function(h,...) {
    for(i in dialogList) {
      dropHandlers = tag(i,"dropHandlers")
      if(length(dropHandlers) > 0) {
        for(i in 1:length(dropHandlers)) {
          removehandler(dropHandlers[[i]]$view.col,dropHandlers[[i]]$id)
        }
      }
    }
  })

  return(obj)
}


## The regular R model methods
summary.gDynamicModelDialog = function(object, ...)
  summary(tag(object$ref,"currentTest"), ...)
anova.gDynamicModelDialog = function(object, ...)
  anova(tag(object$ref,"currentTest"), ...)
coefficients.gDynamicModelDialog = function(object,...)
  coefficients(tag(object$ref,"currentTest"), ...)
effects.gDynamicModelDialog  = function(object,...)
  effects(tag(object$ref,"currentTest"), ...)
fitted.values.gDynamicModelDialog = function(object,...)
  fitted.values(tag(object$ref,"currentTest"), ...)
residuals.gDynamicModelDialog = function(object, ...) 
  residuals(tag(object$ref,"currentTest"),  ...)
predict.gDynamicModelDialog = function(object, ...) 
  predict(tag(object$ref,"currentTest"),  ...)
plot.gDynamicModelDialog = function(x, ...)  {
  print(list(...))
  plot(tag(x$ref,"currentTest"),  ...)
}
##################################################
## wrappers
## use lm function
gui.lm = function(container=NULL) {
  actions =list(
    "drop1" = list(
      FUN = "drop1"
      ),
    "plot: Residuals vs Fitted" =list(
      FUN="plot",
      ARGS = list(which=1)
      ),
    "plot: Normal Q-Q"=list(
      FUN="plot",
      ARGS = list(which=2)
      ),
    "plot: Scale-Location"=list(
      FUN="plot",
      ARGS = list(which=3)
      ),
    "plot: Cook's distance"=list(
      FUN="plot",
      ARGS = list(which=4)
      ),
    "plot: Residuals vs Leverage"=list(
      FUN="plot",
      ARGS = list(which=5)
      ),
    "plot: Cook's distance vs Leverage"=list(
      FUN="plot",
      ARGS = list(which=6)
      )
    )
  dynamicModelWidget(FUN="lm", actions=actions, container=container)
}

gui.rlm = function(container=NULL) {
  ## actions inherit from lm. -- sublcass?
  actions =list(
    "plot: Residuals vs Fitted" =list(
      FUN="plot",
      ARGS = list(which=1)
      ),
    "plot: Normal Q-Q"=list(
      FUN="plot",
      ARGS = list(which=2)
      ),
    "plot: Scale-Location"=list(
      FUN="plot",
      ARGS = list(which=3)
      ),
    "plot: Cook's distance"=list(
      FUN="plot",
      ARGS = list(which=4)
      ),
    "plot: Residuals vs Leverage"=list(
      FUN="plot",
      ARGS = list(which=5)
      ),
    "plot: Cook's distance vs Leverage"=list(
      FUN="plot",
      ARGS = list(which=6)
      )
    )
  dynamicModelWidget(FUN="rlm", actions=actions, container=container)
}

## use aov
gui.aov = function(container=NULL) {
  ## actions inherit from lm. -- sublcass?
  actions =list(
    "plot: Residuals vs Fitted" =list(
      FUN="plot",
      ARGS = list(which=1)
      ),
    "plot: Normal Q-Q"=list(
      FUN="plot",
      ARGS = list(which=2)
      ),
    "plot: Scale-Location"=list(
      FUN="plot",
      ARGS = list(which=3)
      ),
    "plot: Cook's distance"=list(
      FUN="plot",
      ARGS = list(which=4)
      ),
    "plot: Residuals vs Leverage"=list(
      FUN="plot",
      ARGS = list(which=5)
      ),
    "plot: Cook's distance vs Leverage"=list(
      FUN="plot",
      ARGS = list(which=6)
      )
    )
 
  dynamicModelWidget(FUN="aov", actions=actions, container=container)
}
##################################################

### workhorse functino
dynamicModelWidget = function(
  FUN = "lm",
  extra.args = NULL,
  actions = c(),                        # for actions window. Called on obj
  container = NULL,
  ...) {

  group = ggroup(horizontal=FALSE, container = container)
  obj = list(ref=group)
  class(obj) = c("gDynamicModel","gComponent","gWidget")

  ## store values
  tag(obj$ref,"FUN") <- FUN
  tag(obj$ref,"extra.args") <- extra.args

  formulaGroup = ggroup(container=group)
  ## key widgets:
  responseVar = glabel("response", container=formulaGroup,editable=TRUE)
  font(responseVar) <-  c(style="bold")
  tag(obj$ref,"responseVar") <- responseVar
  tag(obj$ref, "responseVarData") <- NA

  glabel(" ~ ", container=formulaGroup)
  
  intercept = gdroplist(c("1","-1"), container=formulaGroup)
  tag(obj$ref,"intercept") <- intercept

  predictorVars = glabel("predictor(s)",container=formulaGroup, editable=TRUE)
  font(predictorVars) <- c(style="bold")
  tag(obj$ref,"predictorVars") <- predictorVars
  tag(obj$ref,"predictorVarsData") <- list()

  addSpring(formulaGroup)
  actionPopup = gdroplist(c(
    "Select an action","Clear formula","Save model object",
    names(actions)),
    container = formulaGroup)
  tag(obj$ref,"actionPopup") <- actionPopup
  tag(obj$ref,"actions") <- actions
  
  gseparator(container=group)  
  outputArea = gtext("")
  size(outputArea) <-  c( 300,300)
  add(group, outputArea, expand=TRUE)
  tag(obj$ref,  "outputArea") <- outputArea

  tag(obj$ref, "res") <- NA
  ## store the drop handlers in the main object.
  tag(obj$ref, "dropHandlers")  <- list()

  ## add droptargets, handlers
  addhandlerchanged(responseVar,
                    handler = function(h,...) {
#                      cat("This doesn't work with dynamic data\n")
                      ids = tag(obj$ref,"dropHandlers")
                      if(length(ids) > 0) {
                        removehandler(obj$ref,ids)
                        tag(obj$ref,"dropHandlers") <- list()
                      }
                      
                      tag(obj$ref,  "responseVarData") <- svalue(responseVar)

                      
                      ## put popup on 1
                      svalue(tag(obj$ref,"actionPopup"),index=TRUE) <- 1
                      
                      update(obj)
                    })
  adddroptarget(responseVar,
                handler=function(h, ...) {
                  
                  tag(obj$ref,"responseVarData") <- h$dropdata
                  svalue(tag(obj$ref, "responseVar")) <- id(h$dropdata)
                  ## put popup on 1
                  svalue(tag(obj$ref,"actionPopup"),index=TRUE) <- 1
                  
                  update(obj)

                  ## now bind to be dynamic *if* a treeviewcolumn
                  if(is.gdataframecolumn(h$dropdata)) {
                    view.col = h$dropdata
                    id = addhandlerchanged(view.col,
                      signal = "edited",
                      handler=function(h,...) update(obj)
                      )
                    dropHandlers = tag(obj$ref,"dropHandlers")
                    dropHandlers[[length(dropHandlers)+1]] = list(
                                  view.col = view.col,
                                  id = id
                                  )
                    tag(obj$ref,"dropHandlers") <- dropHandlers
                  }
                })

  addhandlerchanged(intercept, handler=function(h,...) {
    ## put popup on 1
    svalue(tag(obj$ref,"actionPopup"),index=TRUE) <- 1

    update(obj)
  })
  addhandlerchanged(predictorVars,
                    handler = function(h,...) {
                      ## we need to be careful here. This overwrites any drop data.
                      cat("This doesn't work with dynamic data\n")
                      ids = tag(obj$ref,"dropHandlers")
                      if(length(ids) > 0) {
                        removehandler(obj$ref,ids)
                        tag(obj$ref,"dropHandlers") <- list()
                      }

                      vals = svalue(predictorVars)

                      ## clear out any leading or trailing +
                      vals = sub("\\s*[+]\\s*","",vals)
                      if(vals == "")    # in case of cleaning
                        vars = NULL
                      else 
                        vars = sapply(strsplit(vals,"\\+"),stripWhiteSpace)
                      tag(obj$ref,  "predictorVarsData") <- vars

                      ## add to predictorVars
                      if(is.null(vars))
                        svalue(predictorVars) <- "predictor(s)" # leave target
                      else
                        svalue(predictorVars) <- paste("+",paste(vars,collapse=" + "),collapse="  ")
                      ## put popup on 1
                      svalue(tag(obj$ref,"actionPopup"),index=TRUE) <- 1
                      
                      update(obj)
                    })

  adddroptarget(predictorVars,
                handler = function(h,...) {
                  varList = tag(obj$ref,"predictorVarsData")
                  if(!is.list(varList)) 
                    varList = list()
                  n = length(varList)
                  varList[[n+1]] = h$dropdata
                  tag(obj$ref,"predictorVarsData") <- varList
                  predictorVars = tag(obj$ref,"predictorVars")
                  curLabel = svalue(predictorVars)
                  if(curLabel == "predictor(s)")
                    curLabel = ""                   
                  newLabel = Paste(curLabel," + ", id(h$dropdata))
                  svalue(predictorVars) <- newLabel
                  update(obj)

                  ## now bind to be dynamic *if* a treeviewcolumn
                  if(is.gdataframecolumn(h$dropdata)) {
                    view.col = h$dropdata
                    id = addhandlerchanged(view.col,
                      signal = "edited",
                      handler=function(h,...) update(obj)
                      )
                    dropHandlers = tag(obj$ref,"dropHandlers")
                    dropHandlers[[length(dropHandlers)+1]] = list(
                                  view.col = view.col,
                                  id = id
                                  )
                    tag(obj$ref,"dropHandlers") <- dropHandlers
                  }
                })


  actionPopupHandler = handler = function(h,...) {
    if(svalue(h$obj) == "Select an action") {
      ## do nothing
    } else if(svalue(h$obj) == "Clear formula") {
      ## clear out values, labels
      tag(obj$ref,"responseVarData") <- NA
      svalue(responseVar) <-"response"
      
      svalue(intercept,index=TRUE) <- 1
      
      tag(obj$ref,"predictorVarsData") <- list()
      svalue(predictorVars) <- "predictor(s)"
      
      dispose(outputArea)
      ## would like to reset popup, bu tthis causes infinite loop
      ## clear out res
      tag(obj$ref, "res") <- NA
      
    } else if(svalue(h$obj) == "Save model object") {
      ## pop up a dialog to save the model object
      saveHandler = function(h,...) {
        varName = svalue(theName)
        if(!is.empty(varName)) {
          varName = make.names(varName)
          assign_global(varName, tag(obj$ref,"res"))
          dispose(win)
        }
      }

      win = pmgWC$new("Save model object as...", visible=TRUE)
      gp = ggroup(horizontal=FALSE, container=win)
      glabel("Specify a variable name for the object:", container=gp)
      theName = gedit("",container=gp, handler=saveHandler)
      if(length(lsModels()) > 0) theName[]<- lsModels() # add to type ahead
      buttonGroup = ggroup(container=gp)
      addSpring(buttonGroup)
      gbutton("ok", container=buttonGroup, handler =saveHandler)
      gbutton("cancel", container=buttonGroup, handler = function(h,...) dispose(win))
      ## that's all
    } else {
      argList = c(list(obj),actions[[svalue(h$obj)]]$ARGS)
      do.call(actions[[svalue(h$obj)]]$FUN, argList)
      
    }
  }

  ## how to set value of popup back to start without call ing handler?
  addhandlerchanged(actionPopup,handler=actionPopupHandler)
  
  return(obj)
}


## regular model methods
summary.gDynamicModel = function(object, ...)
  summary(tag(object$ref,"res"), ...)
anova.gDynamicModel = function(object, ...)
  anova(tag(object$ref,"res"), ...)
coefficients.gDynamicModel = function(object,...)
  coefficients(tag(object$ref,"res"), ...)
effects.gDynamicModel  = function(object,...)
  effects(tag(object$ref,"res"), ...)
fitted.values.gDynamicModel = function(object,...)
  fitted.values(tag(object$ref,"res"), ...)
residuals.gDynamicModel = function(object, ...) 
  residuals(tag(object$ref,"res"),  ...)
predict.gDynamicModel = function(object, ...) 
  predict(tag(object$ref,"res"),  ...)
plot.gDynamicModel = function(x, ...) 
  plot(tag(x$ref,"res"),  ...)

## This is main function
## call function, update widgets
## this is *not* the update for models traditionally expected
update.gDynamicModel = function(object, ...) {
  obj = object$ref                        # for s3 consistency

  
  FUN = tag(obj,"FUN")
  extra.args = tag(obj,"extra.args")

  responseVarData = tag(obj,"responseVarData")
  intercept = tag(obj,"intercept")
  intercept.val = svalue(intercept) # a string
  predictorVarsData = tag(obj,"predictorVarsData") # a list

  
  ## make formula
  ## how to avoid eval(parse... here?
  if(is.na(responseVarData)) {
    cat("Need a response variable\n")
    return()
  }

  env = environment()
  assign(id(responseVarData),svalue(responseVarData), envir=env)

  if(length(predictorVarsData) > 0) {
    sapply(predictorVarsData, function(i)
           assign(id(i), svalue(i), envir=env) ) 
  }

  formula = Paste(id(responseVarData), " ~ ", intercept.val)
  if(length(predictorVarsData) > 0) {
    formula = Paste(formula," + ",
      paste(sapply(predictorVarsData, function(i) id(i)), collapse=" + "))
  }

  command = Paste(FUN,"(",formula)
  if(!is.null(extra.args))
    command = Paste(command,",", extra.args)
  command = Paste(command, ")")

  
  res = eval(parse(text=command), envir=env)
  tag(obj,"res") <- res

  if(FUN == "aov")
    out = c(capture.output(res),"\n","Summary:","\n", capture.output(summary(res)))
  else
    out = capture.output(summary(res))
  

  outputArea = tag(obj,"outputArea")
  dispose(outputArea)
  for(i in out) {
    if(length(grep(":$",i)) > 0)
      add(outputArea,i,font.attr=c(style="monospace",color="blue"))
    else
      add(outputArea,i,font.attr=c(style="monospace"))
  }
}

