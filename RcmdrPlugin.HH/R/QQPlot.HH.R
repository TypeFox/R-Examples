"QQPlot.HH" <-
function()
  ## this function modified by Martin Maechler
  ## this function modified by Richard M. Heiberger to add Shapiro-Wilk test
{
    initializeDialog(title=gettextRcmdr("Quantile-Comparison (QQ) Plot"))
    xBox <- variableListBox(top, Numeric(), title=gettextRcmdr("Variable (pick one)"))
    onOK <- function(){
        x <- getSelection(xBox)
        closeDialog()
       if (0 == length(x)) {
            errorCondition(recall=QQPlot.HH, message=gettextRcmdr("You must select a variable.")) 
            return()
            }
        dist <- tclvalue(distVariable)
        save <- options(warn=-1)
        on.exit(options(save))
        retryMe <- function(msg) {
            Message(message= msg, type="error")
            QQPlot.HH()
        }
        switch(dist,
               "norm" = {
                 args <- 'dist= "norm"'
                 ShapWilk <- as.numeric(tclvalue(ShapWilkVariable))
               },
               "t" =  {
                   df <- tclvalue(tDfVariable)
                   df.num <- as.numeric(df)
                   if (is.na(df.num) || df.num < 1) {
                       retryMe(gettextRcmdr("df for t must be a positive number."))
                       return()
                   }
                   args <- paste('dist="t", df=', df, sep="")
               },
               "chisq" = {
                   df <- tclvalue(chisqDfVariable)
                   df.num <- as.numeric(df)
                   if (is.na(df.num) || df.num < 1) {
                       retryMe(gettextRcmdr("df for chi-square must be a positive number."))
                       return()
                   }
                   args <- paste('dist="chisq", df=', df, sep="")
               },
               "f" = {
                   df1 <- tclvalue(FDf1Variable)
                   df2 <- tclvalue(FDf2Variable)
                   df.num1 <- as.numeric(df1)
                   df.num2 <- as.numeric(df2)
                   if (is.na(df.num1) || df.num1 < 1 ||
                       is.na(df.num2) || df.num2 < 1) {
                       retryMe(gettextRcmdr("numerator and denominator \ndf for F must be positive numbers."))
                       return()                            
                   }
                   args <- paste('dist="f", df1=', df1, ', df2=', df2, sep="")
               },
               ## else -- other `dist' :
           {
               dist <- tclvalue(otherNameVariable)
               params <- tclvalue(otherParamsVariable)
               args <- paste('dist="', dist,'", ', params, sep="")
           }) # end{switch}
        .activeDataSet <- ActiveDataSet()
        labels <-
            if ("1" == tclvalue(identifyVariable))
                paste("rownames(", .activeDataSet, ")", sep="")
            else "FALSE"
        command <- paste("qq.plot", "(", .activeDataSet, "$", x, ", ", args,
                          ", labels=", labels, ")", sep="")
        command2 <- paste("shapiro.test", "(", .activeDataSet, "$", x, ")", sep="")
        doItAndPrint(command)
        if (dist=="norm" && ShapWilk) doItAndPrint(command2)
        activateMenus()
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject="qq.plot")
    distFrame <- tkframe(top)
    distVariable <- tclVar("norm")
    normalButton <- tkradiobutton(distFrame, variable=distVariable, value="norm")
    tButton <- tkradiobutton(distFrame, variable=distVariable, value="t")
    chisqButton <- tkradiobutton(distFrame, variable=distVariable, value="chisq")
    FButton <- tkradiobutton(distFrame, variable=distVariable, value="f")
    otherButton <- tkradiobutton(distFrame, variable=distVariable, value="other")
    checkBoxes(frame="normFrame", boxes=c("ShapWilk"), initialValues=c("1"),
        labels=gettextRcmdr(c("Shapiro-Wilk test of normality")))
    tDfFrame <- tkframe(distFrame)
    tDfVariable <- tclVar("")
    tDfField <- tkentry(tDfFrame, width="6", textvariable=tDfVariable)
    chisqDfFrame <- tkframe(distFrame)
    chisqDfVariable <- tclVar("")
    chisqDfField <- tkentry(chisqDfFrame, width="6", textvariable=chisqDfVariable)
    FDfFrame <- tkframe(distFrame)
    FDf1Variable <- tclVar("")
    FDf1Field <- tkentry(FDfFrame, width="6", textvariable=FDf1Variable)
    FDf2Variable <- tclVar("")
    FDf2Field <- tkentry(FDfFrame, width="6", textvariable=FDf2Variable)
    otherParamsFrame <- tkframe(distFrame)
    otherParamsVariable <- tclVar("")
    otherParamsField <- tkentry(otherParamsFrame, width="30", textvariable=otherParamsVariable)
    otherNameVariable <- tclVar("")
    otherNameField <- tkentry(otherParamsFrame, width="10", textvariable=otherNameVariable)
    identifyVariable <- tclVar("0")
    identifyFrame <- tkframe(top)
    identifyCheckBox <- tkcheckbutton(identifyFrame, variable=identifyVariable)
    tkgrid(getFrame(xBox), sticky="nw")
    tkgrid(tklabel(identifyFrame, text=gettextRcmdr("Identify observations with mouse")),
           identifyCheckBox, sticky="w")
    tkgrid(identifyFrame, sticky="w")
    tkgrid(tklabel(distFrame, text=gettextRcmdr("Distribution"), fg="blue"), columnspan=6, sticky="w")
    tkgrid(tklabel(distFrame, text=gettextRcmdr("Normal")), normalButton, normFrame, sticky="w")
    tkgrid(tklabel(tDfFrame, text=gettextRcmdr("df = ")), tDfField, sticky="w")
    tkgrid(tklabel(distFrame, text="t"), tButton, tDfFrame, sticky="w")
    tkgrid(tklabel(chisqDfFrame, text=gettextRcmdr("df = ")), chisqDfField, sticky="w")
    tkgrid(tklabel(distFrame, text=gettextRcmdr("Chi-square")), chisqButton,
           chisqDfFrame, sticky="w")
    tkgrid(tklabel(FDfFrame, text=gettextRcmdr("Numerator df = ")), FDf1Field,
           tklabel(FDfFrame, text=gettextRcmdr("Denominator df = ")), FDf2Field, sticky="w")
    tkgrid(tklabel(distFrame, text="F"), FButton, FDfFrame, sticky="w")
    tkgrid(tklabel(otherParamsFrame, text=gettextRcmdr("Specify: ")),
           otherNameField, tklabel(otherParamsFrame, text=gettextRcmdr("Parameters: ")),
           otherParamsField, sticky="w")
    tkgrid(tklabel(distFrame, text=gettextRcmdr("Other")), otherButton,
           otherParamsFrame, sticky="w")
    tkgrid(distFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix(rows=5, columns=1)
    }

