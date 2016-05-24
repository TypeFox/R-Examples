rsmFormula <- defmacro(frame=top, hasLhs=TRUE, expr={
    ## prepares formula menu
    .activeDataSet <- ActiveDataSet()
    di <- design.info(get(.activeDataSet))
    fn <- names(di$factor.names)
    bn <- di$block.name
    rn <- di$response.names
    if (getRcmdr("coded")) .codenames <- paste("x", 1:di$nfactors,sep="")
    ## looks at last character
    checkAddOperator <- function(rhs) {
        rhs.chars <- rev(strsplit(rhs, "")[[1]])
        if (length(rhs.chars) < 1)
            return(FALSE)
        check.char <- if ((rhs.chars[1] != " ") || (length(rhs.chars) ==
            1))
            rhs.chars[1]
        else rhs.chars[2]
        !is.element(check.char, c("+", ":", "-", "^","(",","))
    }
    .variables <- c(fn, bn, rn)
    word <- paste("\\[", gettextRcmdr("factor"), "\\]", sep = "")
    variables <- paste(.variables, ifelse(is.element(.variables,
        Factors()), paste("[", gettextRcmdr("factor"), "]", sep = ""),
        ""))
        if (getRcmdr("coded")){ variables[1:di$nfactors] <-
                  paste(.codenames,variables[1:di$nfactors], sep="<")
                  .variables[1:di$nfactors] <- .codenames}
    xBox <- variableListBox(frame, variables, selectmode = "multiple",
        title = gettextRcmdr("Variables (double-click to formula)"))
    onDoubleClick <- if (!hasLhs) {
        function() {
            var <- getSelection(xBox)
            which <- as.numeric(tkcurselection(xBox$listbox))+1
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
                c("+", ":", "-", "^", "(", ",")))
                paste(rhs, paste(.variables[which],collapse="+"), sep = "")
            else paste(rhs, "+", paste(.variables[which],collapse="+"))
            tkicursor(rhsEntry, "end")
            tkxview.moveto(rhsEntry, "1")
        }
    }
    else {
        function() {
            var <- getSelection(xBox)
            which <- as.numeric(tkcurselection(xBox$listbox))+1
            tkselection.clear(xBox$listbox, "0", "end")
            if (length(grep(word, var)) == 1)
                var <- sub(word, "", var)
            lhs <- tclvalue(lhsVariable)
            if (lhs == "")
                tclvalue(lhsVariable) <- var
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
                tclvalue(rhsVariable) <- if (rhs == "" || is.element(check.char,
                  c("+", ":", "-", "^", "(")))
                  paste(rhs, paste(.variables[which],collapse="+"), sep = "")
                else paste(rhs, "+", paste(.variables[which],collapse="+"))
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
    onFO <- function() {
        rhs <- tclvalue(rhsVariable)
        which <- as.numeric(tkcurselection(xBox$listbox))+1
            tkselection.clear(xBox$listbox, "0", "end")
        if (!checkAddOperator(rhs))
            tclvalue(rhsVariable) <- paste(rhs, " FO(", paste(.variables[which],collapse=","),")",sep="")
        else
            tclvalue(rhsVariable) <- paste(rhs, "+ FO(", paste(.variables[which],collapse=","),")",sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
    }
    onTWI <- function() {
        rhs <- tclvalue(rhsVariable)
        which <- as.numeric(tkcurselection(xBox$listbox))+1
            tkselection.clear(xBox$listbox, "0", "end")
        if (!checkAddOperator(rhs))
            tclvalue(rhsVariable) <- paste(rhs, " TWI(", paste(.variables[which],collapse=","),")",sep="")
        else
            tclvalue(rhsVariable) <- paste(rhs, "+ TWI(", paste(.variables[which],collapse=","),")",sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
    }
    onPQ <- function() {
        rhs <- tclvalue(rhsVariable)
        which <- as.numeric(tkcurselection(xBox$listbox))+1
            tkselection.clear(xBox$listbox, "0", "end")
        if (!checkAddOperator(rhs))
            tclvalue(rhsVariable) <- paste(rhs, " PQ(", paste(.variables[which],collapse=","),")",sep="")
        else
            tclvalue(rhsVariable) <- paste(rhs, "+ PQ(", paste(.variables[which],collapse=","),")",sep="")
        tkicursor(rhsEntry, "end")
        tkxview.moveto(rhsEntry, "1")
    }
    onSO <- function() {
        rhs <- tclvalue(rhsVariable)
        which <- as.numeric(tkcurselection(xBox$listbox))+1
            tkselection.clear(xBox$listbox, "0", "end")
        if (!checkAddOperator(rhs))
            tclvalue(rhsVariable) <- paste(rhs, " SO(", paste(.variables[which],collapse=","),")",sep="")
        else
            tclvalue(rhsVariable) <- paste(rhs, "+ SO(", paste(.variables[which],collapse=","),")",sep="")
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
    FOButton <- buttonRcmdr(operatorsFrame, text = "FO: first order",
        command = onFO)
    TWIButton <- buttonRcmdr(operatorsFrame, text = "TWI: 2-factor-interaction",
        command = onTWI)
    PQButton <- buttonRcmdr(operatorsFrame, text = "PQ: pure quadratic",
        command = onPQ)
    SOButton <- buttonRcmdr(operatorsFrame, text = "SO: full second order",
        command = onSO)
    plusButton <- buttonRcmdr(operatorsFrame, text = "+", width = "3",
        command = onPlus)
  # rsm seems to ignore -
  #  minusButton <- buttonRcmdr(operatorsFrame, text = "-", width = "3",
  #      command = onMinus)
    colonButton <- buttonRcmdr(operatorsFrame, text = ":", width = "3",
        command = onColon)
    powerButton <- buttonRcmdr(operatorsFrame, text = "^", width = "3",
        command = onPower)
    leftParenButton <- buttonRcmdr(operatorsFrame, text = "(",
        width = "3", command = onLeftParen)
    rightParenButton <- buttonRcmdr(operatorsFrame, text = ")",
        width = "3", command = onRightParen)
    tkgrid(plusButton, colonButton, powerButton, leftParenButton,
        rightParenButton, sticky = "w")
  #  tkgrid(plusButton, minusButton, colonButton, powerButton, leftParenButton,
  #      rightParenButton, sticky = "w")
    tkgrid(FOButton, TWIButton, PQButton, SOButton, sticky = "w")

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
        rhsEntry <- ttkentry(formulaFrame, width = "80", textvariable = rhsVariable)
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
        else tclVar("")
        rhsEntry <- ttkentry(formulaFrame, width = "80", textvariable = rhsVariable)
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
})
