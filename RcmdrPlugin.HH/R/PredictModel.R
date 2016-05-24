"PredictModel" <-
function(){
    env <- environment()
    initializeDialog(title=gettextRcmdr("Enter X data for new observations"))
    confLevelVar <- tclVar(".95")
    confLevelEntry <- tkentry(top, width="6", textvariable=confLevelVar)
    checkBoxes(frame="seFrame", boxes=c("se"), initialValues=rep(0,1),
        labels=gettextRcmdr(c("Standard Error")))
    outerTableFrame <- tkframe(top)
    assign(".tableFrame", tkframe(outerTableFrame), envir=env)
    setUpTable <- function(...) {
        tkdestroy(get(".tableFrame", envir=env))
        assign(".tableFrame", tkframe(outerTableFrame), envir=env)
        nrows <- as.numeric(tclvalue(rowsValue))
        ncols <- length(colsNames)
        make.col.names <- "tklabel(.tableFrame, text='')"
        for (j in seq(length.out=ncols)) {
            col.varname <- colsNames[j] ## paste(".colname.", j, sep="")
            assign(col.varname, tclVar(colsNames[j]), envir=env)
            make.col.names <- paste(make.col.names, ", ", "tkentry(.tableFrame, width='12', textvariable=", 
                    col.varname, ")", sep="")
            }
        if (ncols > 0) {
eval(parse(text=paste("tkgrid(", make.col.names, ")", sep="")), envir=env)
        for (i in 1:nrows) {   
            varname <- paste(".tab.", i, ".1", sep="") 
            assign(varname, tclVar("") , envir=env)
            row.varname <- paste(".rowname.", i, sep="")
            assign(row.varname, tclVar(i), envir=env)
            make.row <- paste("tkentry(.tableFrame, width='12', textvariable=",
                row.varname, ")", sep="")
            make.row <- paste(make.row, ", ", "tkentry(.tableFrame, width='12', textvariable=", 
                varname, ")", sep="")
            if (ncols > 1)
            for (j in 2:ncols) {
                varname <- paste(".tab.", i, ".", j, sep="")
                assign(varname, tclVar(""), envir=env)
                make.row <- paste(make.row, ", ", "tkentry(.tableFrame, width='12', textvariable=", 
                    varname, ")", sep="")
## recover()
                }
            eval(parse(text=paste("tkgrid(", make.row, ")", sep="")), envir=env)
            }
        tkgrid(get(".tableFrame", envir=env), sticky="w")
      }
      }

    .active.model <- ActiveModel()

    colsNames <- if (is.null(.active.model)) NULL else dimnames(attr(get(.active.model)$terms,"factors"))[[2]]

    rowColFrame <- tkframe(top)
    rowsValue <- tclVar("1")
    rowsSlider <- tkscale(rowColFrame, from=1, to=10, showvalue=FALSE, variable=rowsValue,
        resolution=1, orient="horizontal", command=setUpTable)
    rowsShow <- tklabel(rowColFrame, textvariable=rowsValue, width=2, justify="right")

    onOK <- function() {
    if (is.null(.active.model)) {
      errorCondition(recall=PredictModel, message=sprintf(gettextRcmdr("No active model.  Please press Cancel.")))
      return()
    }
        nrows <- as.numeric(tclvalue(rowsValue))
        ncols <- length(colsNames)

    col.names <- dimnames(attr(get(.active.model)$terms,"factors"))[[2]]
    newdata <- get(activeDataSet())[1:nrows, col.names, drop=FALSE]
    for (i in names(newdata)) newdata[[i]][] <- NA
    ## sapply(newdata, levels)
    row.names <- rep("", nrows)
    for (i in 1:nrows) row.names[i] <- 
      eval(parse(text=paste("tclvalue(", paste(".rowname.", i, sep=""),")", sep="")))
        for (i in 1:nrows){
            for (j in 1:ncols){
                varname <- paste(".tab.", i, ".", j, sep="")
                cell.value <-
                  eval(parse(text=paste("tclvalue(", varname,")", sep="")))
                aDSj <- get(activeDataSet())[[col.names[j]]]
                if (is.factor(aDSj) || is.character(aDSj))
                  newdata[i,j] <- cell.value
                else
                  newdata[i,j] <- as.numeric(cell.value)
                }
            }
        newdata <- na.omit(newdata)
        if (prod(dim(newdata)) != nrows*ncols) {
            errorCondition(recall=PredictModel,
                           message=sprintf(gettextRcmdr("Number of valid entries (%d)\nnot equal to number of rows (%d) * number of columns (%d)."),
                             prod(dim(newdata)), nrows, ncols))
            return()
            }
        if (length(unique(row.names)) != nrows) {
            errorCondition(recall=PredictModel, message=gettextRcmdr("Row names are not unique."))
            return()
            }     
        if (length(unique(col.names)) != ncols) {
            errorCondition(recall=PredictModel, message=gettextRcmdr("Column names are not unique."))
            return()
            }     
    closeDialog()

    confLevel <- as.numeric(tclvalue(confLevelVar))
    j <- 1
    command <- paste(".NewData <- data.frame(", col.names[j], "=",
                     deparse(newdata[,j]), sep="")
    if (ncols > 1)
      for (j in 2:ncols) {
        command <- paste(command, ", ", col.names[j], "=",
                         deparse(newdata[,j]), sep="")
      }
    command <- paste(command, ", ", "row.names=",
                     deparse(row.names),")", sep="")
    doItAndPrint(command)
    doItAndPrint(".NewData  # Newdata")

    command <- paste('predict(',
                     .active.model,
                     ', newdata=.NewData, interval="',
                     tclvalue(predictVariable),
                     '", level=', tclvalue(confLevelVar),
                     ', se.fit=', ("1"==tclvalue(seVariable)) ,
                     ')', sep="")
    doItAndPrint(command)
    ## logger("remove(.NewData)") 
    ## remove(.NewData, envir=.GlobalEnv)                                                      
    tkfocus(CommanderWindow())
  }

    OKCancelHelp(helpSubject="predict")
    
    radioButtons(name="predict",
                 buttons=c("none",
                   "confidence",
                   "prediction"),
                 values=c("none", "confidence", "prediction"),
                 initialValue="confidence",
                 labels=gettextRcmdr(c
                   ("point estimate only",
                    "confidence interval for mean",
                    "prediction interval for individual")),
                 title=gettextRcmdr("Prediction Interval"))
    
    if (length(colsNames) > 0) {
    
    tkgrid(tklabel(rowColFrame,
                   text=gettextRcmdr("Number of Rows:")),
           rowsSlider, rowsShow, sticky="w")
    tkgrid(rowColFrame, sticky="w")
    tkgrid(tklabel(top,
                   text=gettextRcmdr("Enter X values:"),
                   fg="blue"), sticky="w")
    tkgrid(outerTableFrame, sticky="w")
    tkgrid(predictFrame, sticky="w")
    tkgrid(tklabel(top, text=gettextRcmdr("Confidence Level")),
           confLevelEntry, sticky="w")
    tkgrid(seFrame, sticky="w")
  }
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=7, columns=2)
    }

## dimnames(attr(RegModel.1$terms,"factors"))[[2]]


## source("~/HH-R.package/RcmdrPlugin.HH/R/PredictModel.R")
