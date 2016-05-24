# Last modified Feb 16, 2008

`checkMultiLevelFactors` <-
function (n = 1) 
{
    if (length(MultiLevelFactors()) < n) {
        if (n > 1) 
            Message(message = sprintf(gettextRcmdr("There fewer than %d multi-level factors in the active data set."), 
                n), type = "error")
        else Message(message = gettextRcmdr("There are no multi-level factors in the active data set."), 
            type = "error")
        tkfocus(CommanderWindow())
        FALSE
    }
    else TRUE
}


`listMultiLevelFactors` <-
function (dataSet = ActiveDataSet()) 
{
    factors <- listFactors(dataSet)
    if (length(factors) == 0) 
        return(NULL)
    factors[sapply(factors, function(.x) 2 < length(levels(eval(parse(text = .x), 
        envir = eval(parse(text = dataSet), envir = .GlobalEnv)))))]
}


`MultiLevelFactors` <-
function (names) 
{
    if (missing(names)) 
        getRcmdr("multiLevelFactors")
    else putRcmdr("multiLevelFactors", names)
}


`multiLevelFactorsP` <-
function (n = 1) 
activeDataSetP() && length(listMultiLevelFactors()) >= n
`multiSampleProportionsTest` <-
function () 
{
    require("abind")
    initializeDialog(title = gettextRcmdr("Test equality of several proportions..."))
    .multifactors <- listMultiLevelFactors()
    .twoLevelFactors <- TwoLevelFactors()
    groupsBox <- variableListBox(top, .multifactors, title = gettextRcmdr("Groups (pick one)"))
    xBox <- variableListBox(top, .twoLevelFactors, title = gettextRcmdr("Response Variable (pick one)"))
    onOK <- function() {
        groups <- getSelection(groupsBox)
        if (length(groups) == 0) {
            errorCondition(recall = multiSampleProportionsTest, 
                message = gettextRcmdr("You must select a groups variable."))
            return()
        }
        x <- getSelection(xBox)
        if (length(x) == 0) {
            errorCondition(recall = multiSampleProportionsTest, 
                message = gettextRcmdr("You must select a response variable."))
            return()
        }
        if (x == groups) {
            errorCondition(recall = multiSampleProportionsTest, 
                message = gettextRcmdr("Groups and response variables must be different."))
            return()
        }
        closeDialog()
        command <- paste(".Table <- xtabs(~", groups, "+", x, ", data=", 
            ActiveDataSet(), ")", sep = "")
        logger(paste(".Table <-", command))
        justDoIt(command)
        doItAndPrint("rowPercents(.Table)")
        doItAndPrint(paste("prop.test(.Table)", sep = ""))
        logger("remove(.Table)")
        remove(.Table, envir = .GlobalEnv)
        tkfocus(CommanderWindow())
    }
    OKCancelHelp(helpSubject = "prop.test")
    rightFrame <- tkframe(top)
    tkgrid(getFrame(groupsBox), getFrame(xBox), sticky = "nw")
    tkgrid(buttonsFrame, columnspan = 2, sticky = "w")
    dialogSuffix(rows = 5, columns = 2)
}
