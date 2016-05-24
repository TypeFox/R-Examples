modelFormulaDoE <- defmacro(frame = modelFrame, hasLhs = TRUE, rhschr="", expr={
    currentModel <- FALSE
    # not for all cases useful to have current model
    # therefore local false version
    # rhs can be given to the macro as rhschr
    checkAddOperator <- function(rhs) {
        rhs.chars <- rev(strsplit(rhs, "")[[1]])
        if (length(rhs.chars) < 1)
            return(FALSE)
        check.char <- if ((rhs.chars[1] != " ") || (length(rhs.chars) ==
            1))
            rhs.chars[1]
        else rhs.chars[2]
        !is.element(check.char, c("+", "*", ":", "/", "-", "^",
            "(", "%"))
    }
    .variables <- Variables()
    word <- paste("\\[", gettextRcmdr("factor"), "\\]", sep = "")
    variables <- paste(.variables, ifelse(is.element(.variables,
        Factors()), paste("[", gettextRcmdr("factor"), "]", sep = ""),
        ""))
    xBox <- variableListBox(frame, variables, selectmode = "multiple",
        title = gettextRcmdr("Variables (double-click to formula)"))
    onDoubleClick <- if (!hasLhs) {
        function() {
            var <- getSelection(xBox)
            tkselection.clear(xBox$listbox, "0", "end")
            if (length(grep(word, var)) == 1)
                var <- sub(word, "", var)
            tkfocus(rhsEntry)
            rhs <- tclvalue(rhsVariable)
            rhs.chars <- rev(strsplit(rhs, "")[[1]])
            check.char <- if (length(rhs.chars) > 0) {
                if ((rhs.chars[1] != " ") || (length(rhs.chars) ==
                  1))
                  rhs.chars[1]
                else rhs.chars[2]
            }
            else ""
           tclvalue(rhsVariable) <- if (rhs == "" || is.element(check.char,
                c("+", "*", ":", "/", "-", "^", "(", "%")))
                paste(rhs, var, sep = "") else paste(rhs, "+", var)
            tkicursor(rhsEntry, "end")
            tkxview.moveto(rhsEntry, "1")
        }
    }
    else {
        function() {
            var <- getSelection(xBox)
            which <- tkcurselection(xBox$listbox)
            tkselection.clear(xBox$listbox, "0", "end")
            if (length(grep(word, var)) == 1)
                var <- sub(word, "", var)
            lhs <- tclvalue(lhsVariable)
            if (lhs == "" || tclvalue(tkselection.present(lhsEntry)) ==
                "1") {
                tclvalue(lhsVariable) <- var
                tkselection.clear(lhsEntry)
                tkfocus(rhsEntry)
            }
            else {
                tkfocus(rhsEntry)
                rhs <- tclvalue(rhsVariable)
                rhs.chars <- rev(strsplit(rhs, "")[[1]])
                check.char <- if (length(rhs.chars) > 0) {
                  if ((rhs.chars[1] != " ") || (length(rhs.chars) ==
                    1))
                    rhs.chars[1]
                  else rhs.chars[2]
                }
                else ""
            putRcmdr(tclvalue(rhsVariable)) <- if (rhs == "" || is.element(check.char,
                c("+", "*", ":", "/", "-", "^", "(", "%")))
                paste(rhs, var, sep = "")
                else paste(rhs, "+", var)
            }
            tkicursor(rhsEntry, "end")
            tkxview.moveto(rhsEntry, "1")
        }
    }
    tkbind(xBox$listbox, "<Double-ButtonPress-1>", onDoubleClick)
    onPlus <- function() {
        rhs <- tclvalue(rhsVariable)
        var <- getSelection(xBox)
        tkselection.clear(xBox$listbox, "0", "end")
        if ((check <- !checkAddOperator(rhs)) && length(var) ==
            0)
            return()
        if (length(var) > 1) {
            if (length(grep(word, var)) > 0)
                var <- sub(word, "", var)
            if (length(var) > 1)
                var <- paste(var, collapse = " + ")
        }
        tclvalue(rhsVariable) <- paste(rhs, if (!check)
            " + ", var, sep = "")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
    }
    onTimes <- function() {
        rhs <- tclvalue(rhsVariable)
        var <- getSelection(xBox)
        tkselection.clear(xBox$listbox, "0", "end")
        if ((check <- !checkAddOperator(rhs)) && length(var) ==
            0)
            return()
        if (length(var) > 1) {
            if (length(grep(word, var)) > 0)
                var <- sub(word, "", var)
            var <- trim.blanks(var)
            if (length(var) > 1)
                var <- paste(var, collapse = "*")
            tclvalue(rhsVariable) <- paste(rhs, if (!check)
                " + ", var, sep = "")
        }
        else tclvalue(rhsVariable) <- paste(rhs, if (!check)
            "*", sep = "")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
    }
    onColon <- function() {
        rhs <- tclvalue(rhsVariable)
        var <- getSelection(xBox)
        tkselection.clear(xBox$listbox, "0", "end")
        if ((check <- !checkAddOperator(rhs)) && length(var) ==
            0)
            return()
        if (length(var) > 1) {
            if (length(grep(word, var)) > 0)
                var <- sub(word, "", var)
            var <- trim.blanks(var)
            if (length(var) > 1)
                var <- paste(var, collapse = ":")
            tclvalue(rhsVariable) <- paste(rhs, if (!check)
                " + ", var, sep = "")
        }
        else tclvalue(rhsVariable) <- paste(rhs, if (!check)
            ":", sep = "")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
    }
    onSlash <- function() {
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs))
            return()
        tclvalue(rhsVariable) <- paste(rhs, "/", sep = "")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
    }
    onIn <- function() {
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs))
            return()
        tclvalue(rhsVariable) <- paste(rhs, "%in% ")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
    }
    onMinus <- function() {
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs))
            return()
        tclvalue(rhsVariable) <- paste(rhs, "- ")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
    }
    onPower <- function() {
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs))
            return()
        tclvalue(rhsVariable) <- paste(rhs, "^", sep = "")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
    }
    onLeftParen <- function() {
        tkfocus(rhsEntry)
        rhs <- tclvalue(rhsVariable)
        tclvalue(rhsVariable) <- paste(rhs, "(", sep = "")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
    }
    onRightParen <- function() {
        rhs <- tclvalue(rhsVariable)
        if (!checkAddOperator(rhs))
            return()
        tclvalue(rhsVariable) <- paste(rhs, ")", sep = "")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
    }
    outerOperatorsFrame <- tkframe(frame)
    operatorsFrame <- tkframe(outerOperatorsFrame)
    plusButton <- buttonRcmdr(operatorsFrame, text = "+", width = "3",
        command = onPlus)
    timesButton <- buttonRcmdr(operatorsFrame, text = "*", width = "3",
        command = onTimes)
    colonButton <- buttonRcmdr(operatorsFrame, text = ":", width = "3",
        command = onColon)
    slashButton <- buttonRcmdr(operatorsFrame, text = "/", width = "3",
        command = onSlash)
    inButton <- buttonRcmdr(operatorsFrame, text = "%in%", width = "5",
        command = onIn)
    minusButton <- buttonRcmdr(operatorsFrame, text = "-", width = "3",
        command = onMinus)
    powerButton <- buttonRcmdr(operatorsFrame, text = "^", width = "3",
        command = onPower)
    leftParenButton <- buttonRcmdr(operatorsFrame, text = "(",
        width = "3", command = onLeftParen)
    rightParenButton <- buttonRcmdr(operatorsFrame, text = ")",
        width = "3", command = onRightParen)
    tkgrid(plusButton, timesButton, colonButton, slashButton,
        inButton, minusButton, powerButton, leftParenButton,
        rightParenButton, sticky = "w")
    formulaFrame <- tkframe(frame)
    if (hasLhs) {
        tkgrid(labelRcmdr(outerOperatorsFrame, text = gettextRcmdr("Model Formula:     "),
            fg = "blue"), operatorsFrame)
        lhsVariable <- if (currentModel)
            tclVar(currentFields$lhs)
        else tclVar("")
        rhsVariable <- if (currentModel)
            tclVar(currentFields$rhs)
        else tclVar("")
        rhsEntry <- ttkentry(formulaFrame, width = "50", textvariable = rhsVariable)
        rhsXscroll <- ttkscrollbar(formulaFrame, orient = "horizontal",
            command = function(...) tkxview(rhsEntry, ...))
        tkconfigure(rhsEntry, xscrollcommand = function(...) tkset(rhsXscroll,
            ...))
        lhsEntry <- ttkentry(formulaFrame, width = "10", textvariable = lhsVariable)
        lhsScroll <- ttkscrollbar(formulaFrame, orient = "horizontal",
            command = function(...) tkxview(lhsEntry, ...))
        tkconfigure(lhsEntry, xscrollcommand = function(...) tkset(lhsScroll,
            ...))
        tkgrid(lhsEntry, labelRcmdr(formulaFrame, text = " ~    "),
            rhsEntry, sticky = "w")
        tkgrid(lhsScroll, labelRcmdr(formulaFrame, text = ""),
            rhsXscroll, sticky = "w")
        tkgrid.configure(lhsScroll, sticky = "ew")
    }
    else {
        rhsVariable <- if (currentModel)
            tclVar(currentFields$rhs)
        else tclVar(rhschr)
        rhsEntry <- ttkentry(formulaFrame, width = "50", textvariable = rhsVariable)
        rhsXscroll <- ttkscrollbar(formulaFrame, orient = "horizontal",
            command = function(...) tkxview(rhs, ...))
        tkconfigure(rhsEntry, xscrollcommand = function(...) tkset(rhsXscroll,
            ...))
        tkgrid(labelRcmdr(formulaFrame, text = "   ~ "), rhsEntry,
            sticky = "w")
        tkgrid(labelRcmdr(formulaFrame, text = ""), rhsXscroll,
            sticky = "w")
    }
    tkgrid.configure(rhsXscroll, sticky = "ew")
    tkgrid(operatorsFrame, sticky="w")
    tkgrid(getFrame(xBox), outerOperatorsFrame, sticky="w")
    tkgrid(formulaFrame, sticky="w", columnspan="2")
}
)