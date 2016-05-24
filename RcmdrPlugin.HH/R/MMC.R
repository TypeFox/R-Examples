MMCmenu <- function() {
    initializeDialog(title=gettextRcmdr("MMC Plot"))
    modelBox <- variableListBox(top, listAOVModels(),
                                title=gettextRcmdr("AOV model (pick one)"))
    checkBoxes(frame="tieFrame", boxes=c("tie"), initialValues=rep(0),
        labels=gettextRcmdr(c("Tiebreaker Plot")))

    onOK <- function() {
      model <- getSelection(modelBox)

      if (length(model) != 1) {
        errorCondition(recall=MMCmenu,
                       message=sprintf(gettextRcmdr("Please select an AOV model.")))
        return()
      }

      if (get(model)$df.residual == 0) {
        errorCondition(recall=MMCmenu,
                       message=sprintf(gettextRcmdr("Model has 0 df for residual.\nPlease pick a different model or cancel the MMC command.")))
        return()
      }

      closeDialog()

      ## doItAndPrint(paste(".activeAOVModel <- '", model, "'", sep=""))
      putRcmdr(".activeAOVModel", model)

      command.factors <-
        paste(model, ".factors <- ", "names(", model, "$xlevels)", sep="")
      eval(parse(text=command.factors))

      command.factors.length <-
        paste("length(", model, ".factors)", sep="")
      focus.and.lmatrows <-  if (eval(parse(text=command.factors.length)) > 1)
        {
        putRcmdr(".MMC2.result", "cancelled") ## to check later
        MMC2menu()
        getRcmdr(".MMC2.result")
      }
      else
        ""
      if (focus.and.lmatrows=="cancelled") return()

      c0a <- paste("right.omd <- .8")
      eval(parse(text=c0a))
      
      command0<- paste("old.omd <- par(omd=c(0, ", right.omd, ", 0,1))", sep="")
      doItAndPrint(command0)

      command1 <- paste(model, ".mmc <- mmc(", model,
                        focus.and.lmatrows, ")", sep="")

      doItAndPrint(command1)

      command2 <- paste(model, ".mmc", sep="")
      doItAndPrint(command2)

      c2a <- paste(command2, "$mca$table[,'estimate']", sep="")
      c2b <- paste("x.left <- -max(abs(",
                   c2a, "))",
                   sep="")
      eval(parse(text=c2b))

      c2c <- paste("x.right <- max(", command2,
                   "$mca$table[,'upper'])",
                   sep="")
      eval(parse(text=c2c))

      eval(parse(text=(paste("x.offset <- (x.left+x.right)/2"))))


      c2e <- paste("y.range <- range(", command2,
                   "$none$table[,'estimate'])",
                   sep="")
      eval(parse(text=c2e))

      eval(parse(text=("ry <- mean(y.range) + (y.range - mean(y.range))*(-2*x.left+x.offset)/(-2*x.left)*1.4")))

      command3 <- paste("plot(", model, ".mmc, x.offset=",
                        x.offset, ", ry=", deparse(ry), ")", sep="")
      doItAndPrint(command3)

      if (tclvalue(tieVariable)=="1") {
        command4 <-
          paste("plotMatchMMC(", model, ".mmc$mca, xlabel.print=FALSE)",
                sep="")
        doItAndPrint(command4)
      }

      command5<- "par(old.omd)"
      doItAndPrint(command5)

      
      logger("## The default placement of MMC plots may not be ideal.
## You may need to adjust the x.offset and/or ry arguments.
## You may need to adjust the par('omd') parameter before the plot.
## See the help file for MMC for examples.")
    }

    OKCancelHelp(helpSubject="MMC")
    tkgrid(getFrame(modelBox), sticky="w")
    tkgrid(tieFrame, sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=3, columns=2)
  }

MMC2menu <- function() {
    model <- getRcmdr(".activeAOVModel")
    initializeDialog(title=gettextRcmdr(paste("MMC Plot for", model)))
    ## checkBoxes(frame="tieFrame", boxes=c("tie"), initialValues=rep(0),
    ##     labels=gettextRcmdr(c("Tiebreaker Plot")))

    xlevels <- eval(parse(text=paste(model, "$xlevels", sep="")))
    factors <- names(xlevels)
    focusBox <-
      variableListBox(top,
                      factors,
                      title=gettextRcmdr("focus factor (pick one)"))

    tmp <- eval(parse(text=paste("glht(get(model), linfct=mcp(",
                        factors[1], "='Tukey'))", sep="")))
    lmat.rownames <- dimnames(tmp$linfct)[[2]]

    checkBoxes(frame="choose.lmat.rowsFrame", boxes=c("lmatrows"),
               initialValues=rep(1),
        labels=gettextRcmdr(c("Accept the default lmat.rows.")))

    lmat.rowsBox <- variableListBox(top,
                                lmat.rownames,
                                title=gettextRcmdr("Or select all the lmat rows corresponding to the focus factor.")
                                    , selectmode="multiple")

    onOK <- function() {

      focus <- getSelection(focusBox)
        if (length(focus) != 1) {
            errorCondition(recall=MMC2menu,
                           message=sprintf(gettextRcmdr("Please select a focus factor.")))
            return()
            }

      if (tclvalue(lmatrowsVariable)=="1") {
        lmat.rows <- grep(paste("^", focus, sep=""),
                          gsub(".*:.*", "", lmat.rownames))

        logger(paste("Default lmat.rows:",
                     paste(lmat.rownames[lmat.rows], collapse=" ")))
        if (length(lmat.rows)==0) {
          logger("We can't verify that the rows you selected are correct.")
          errorCondition(recall=MMC2menu,
                         message=sprintf(gettextRcmdr("Please select the lmat rowfrom the list below.")))
        }
      }
      else {
        lmat.rowselected <- getSelection(lmat.rowsBox)
        lmat.rows <- match(lmat.rowselected, lmat.rownames)

        if (!((length(xlevels[[focus]]) == length(lmat.rows)+1) &&
              all(substring(lmat.rowselected, 1, nchar(focus)) ==  focus) &&
              !any(grep("\\[", sub("\\[", "\\*", lmat.rownames)) %in% lmat.rows) &&
              all(diff(lmat.rows)==1))) {
          errorCondition(recall=MMC2menu,
                         message=sprintf(gettextRcmdr("Please select the lmat rows corresponding to the focus factor.")))
          return()
        }
      }
      
      putRcmdr(".MMC2.result",
               paste(", focus='", focus,
                     "', lmat.rows=", deparse(lmat.rows),
                     sep=""))
      closeDialog()


    }

    OKCancelHelp(helpSubject="MMC")
    tkgrid(getFrame(focusBox), sticky="w")
    tkgrid(choose.lmat.rowsFrame, sticky="w")
    tkgrid(getFrame(lmat.rowsBox), sticky="w")
    tkgrid(buttonsFrame, columnspan=2, sticky="w")
    dialogSuffix(rows=3, columns=2)
}

AOVModelsP <-
  function (n = 1) 
  length(listAOVModels()) >= n


## source("~/HH-R.package/RcmdrPlugin.HH/R/MMC.R")
