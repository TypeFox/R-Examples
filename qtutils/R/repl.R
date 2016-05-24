
qselectedText_QTextEdit <- function(x)
{
    ans <- x$textCursor()$selection()$toPlainText()
    if (is.null(ans)) ""
    else ans
}



## not exactly a REPL, but something that looks like it. The idea is
## to have a widget where one can input R commands, and evaluate them
## on demand.  The typical use case would be a limited REPL as used in
## browser() and recover().


## Two editors in a splitter, one for input and one for output.


tryComplete <- function(text, cursor = nchar(text))
{
    utils:::.win32consoleCompletion(linebuffer = text, cursorPosition = cursor,
                                    check.repeat = TRUE, 
                                    minlength = -1)
}


tryParseEval <- function(text, env)
{
    exprs <- try(parse(text = text), silent = TRUE)
    esrc <- lapply(attr(exprs, "srcref"), as.character)
    if (is(exprs, "try-error")) return (exprs)
    ans <- vector(mode = "list", length = length(exprs))
    for (i in seq_along(exprs))
    {
        ein <- esrc[[i]]
        .GlobalEnv$.qexpr <- exprs[[i]] # env$.expr <- exprs[[i]]
        output <-
            capture.output(evis <- try(evalq(withVisible(eval(.GlobalEnv$.qexpr)),
                                             envir = env),
                                        silent = TRUE))
        ##eout <- evis
        ans[[i]] <- list(ein = paste("> ", paste(ein, collapse = "\n+ "), sep = ""),
                         evis = evis,
                         output = output)
    }
    ans
}




##' Creates a Qt Widget that emulates the R REPL.
##'
##' This widget tries to emulate the behaviour of the R command-line
##' interface (Read-Eval-Print-Loop) in a GUI.  The current
##' implementation is essentially a proof-of-concept, and not meant
##' for serious use.
##'
##' Two versions are available.  \code{qrepl} emulates the
##' conventional REPL interface where commands are typed at a command
##' prompt and evaluated when Enter is pressed (except that parse
##' errors are trapped and not evaluated).   \code{qrepl} provides an
##' alternative interface with two command specification modes and a
##' common output area.   Commands can be entered either in input
##' mode, where one or more commands may be typed and then  executed
##' using Ctrl+Enter, or in edit mode, where commands are typed in a
##' file, and selections may be executed similarly.  The latter mode
##' is intended to facilitate the recommended method of working in
##' ESS, where commands are entered and modified within an editor, and
##' executed as necessary.  When working in edit mode, commands
##' executed are automatically added to the editor area, keeping a
##' record of those commands.
##'
##' All code editing interfaces support command completion and code
##' indentation (but the rules are not yet customizable).
##' 
##' @title A Qt based REPL emulator
##' @param env The evaluation environment for the REPL. 
##' @param ... Further arguments, passed on to \code{\link{RCodeEditor}}.
##' @param family Font family to be used.
##' @param pointsize Font pointsize to be used.
##' @param incolor The color used for code that is evaluated.
##' @param outcolor The color used for output.
##' @param msgcolor The color used for messages.
##' @param html.preferred Logical flag indicating whether HTML output is preferred (for table-like objects).
##' @param history Logical flag indicating whether command history
##' should be available.  If enabled, pressing Ctrl+Up or Down arrows
##' allow navigation  through previous commands.   History is
##' currently not retained across invocations.  
##' @param title Character string giving window title.
##' @return A QWidget instance.
##' @author Deepayan Sarkar
qrepl <- function(env = .GlobalEnv,
                  ...,
                  family = "monospace", pointsize = 12,
                  incolor = "red",
                  outcolor = "blue",
                  msgcolor = "black",
                  html.preferred = FALSE,
                  history = TRUE,
                  title = sprintf("Qt REPL %s", capture.output(print(.GlobalEnv))))
{
    html.preferred <- html.preferred && require(xtable)
    font <- qfont(family = family, pointsize = pointsize)
    informat <- Qt$QTextCharFormat()
    informat$setForeground(qbrush(incolor))
    informat$setFont(font)
    outformat <- Qt$QTextCharFormat()
    outformat$setForeground(qbrush(outcolor))
    outformat$setFont(font)
    msgformat <- Qt$QTextCharFormat()
    msgformat$setForeground(qbrush(msgcolor))
    msgformat$setFont(font)

    ## input1: REPL-like mode, type and execute code
    ined1 <- qeditor(rsyntax = TRUE, richtext = FALSE,
                     family = family, pointsize = pointsize, ...)
    ## input2: Editor mode, select and execute code.  input1 text gets appended here
    ined2 <- qeditor(rsyntax = TRUE, richtext = FALSE,
                     family = family, pointsize = pointsize, ...)
    ## output
    outed <- Qt$QTextEdit()
    outed$readOnly <- TRUE
    ## tabbed widget holding the two input editors
    intab <- Qt$QTabWidget()
    intab$addTab(ined1, label = "Input mode")
    intab$addTab(ined2, label = "Edit mode")

    intab$setSizePolicy(Qt$QSizePolicy$Expanding,
                        Qt$QSizePolicy$Preferred)
    outed$setSizePolicy(Qt$QSizePolicy$Expanding,
                        Qt$QSizePolicy$Expanding)
    ## qsetExpanding(intab, vertical = FALSE)
    ## qsetExpanding(outed, vertical = TRUE)
    ## messages
    msg <- Qt$QLabel("")
    msg$wordWrap <- TRUE
    ## container
    container <- Qt$QSplitter(Qt$Qt$Vertical)
    container$addWidget(outed)
    container$addWidget(intab)
    container$addWidget(msg)
    container$setStretchFactor(0L, 10L)
    container$setStretchFactor(1L, 0L)
    container$setStretchFactor(2L, 0L)
    msg$text <- "Type code, press Ctrl+Return to evaluate"
    ined1$setContextMenuPolicy(Qt$Qt$ActionsContextMenu)
    ined2$setContextMenuPolicy(Qt$Qt$ActionsContextMenu)
    outed$setContextMenuPolicy(Qt$Qt$ActionsContextMenu)

    ## History mechanism: Ctrl+Up/Down.  Unused if history=FALSE.
    lhist <- list() # list of commands
    ihist <- 0L     # index
    cur_command <- ""
    if (history)
    {
        addToHistory <- function(s)
        {
            n <- length(lhist)
            if (n == 0 || s != lhist[[n]]) lhist[[n + 1L]] <<- s
        }
        prevHistoryItem <- function()
        {
            n <- length(lhist)
            ihist <<-
                if (ihist == 0) n
                else ihist - 1L 
            if (ihist == 0) cur_command else lhist[[ihist]]
        }
        nextHistoryItem <- function()
        {
            n <- length(lhist)
            ihist <<- ihist + 1L
            if (ihist == n + 1L) ihist <<- 0L
            if (ihist == 0) cur_command else lhist[[ihist]]
        }
        showPrevItem <- function(...)
        {
            ## str(list(lhist = lhist, ihist = ihist, cur_command = cur_command))
            ## if editing something new, make that current before changing
            if (ihist == 0L) cur_command <<- ined1$plainText
            ined1$setPlainText(prevHistoryItem())
        }
        showNextItem <- function(...)
        {
            ## str(list(lhist = lhist, ihist = ihist, cur_command = cur_command))
            ## if editing something new, make that current before changing
            if (ihist == 0L) cur_command <<- ined1$plainText
            ined1$setPlainText(nextHistoryItem())
        }

        prevHistoryAct <- Qt$QAction(text = "Previous command", parent = ined1)
        prevHistoryAct$setShortcut(Qt$QKeySequence("Ctrl+Up"))
        prevHistoryAct$setShortcutContext(Qt$Qt$WidgetShortcut)
        qconnect(prevHistoryAct, signal = "triggered", handler = showPrevItem)
        ined1$addAction(prevHistoryAct)
        
        nextHistoryAct <- Qt$QAction(text = "Next command", parent = ined1)
        nextHistoryAct$setShortcut(Qt$QKeySequence("Ctrl+Down"))
        nextHistoryAct$setShortcutContext(Qt$Qt$WidgetShortcut)
        qconnect(nextHistoryAct, signal = "triggered", handler = showNextItem)
        ined1$addAction(nextHistoryAct)
    }

    ## function to perform code execution
    executeCode <- function(text, mode = c("input", "edit"))
    {
        mode <- match.arg(mode)
        pe <- tryParseEval(text = text, env = env)
        if (is(pe, "try-error"))
        {
            msg$text <- paste(strsplit(as.character(pe), "\n", fixed = TRUE)[[1]], collapse = "\\n")
        }
        else
        {
            msg$text <- ""
            if (mode == "input")
            {
                if (history) addToHistory(text)
                ihist <<- 0L
                ined1$selectAll()
                ined2$appendPlainText(text)
            }
            for (i in seq_along(pe))
            {
                ein <- pe[[i]]$ein
                output <- pe[[i]]$output
                evis <- pe[[i]]$evis
                ## input
                outed$setCurrentCharFormat(informat)
                outed$append(ein)
                ## output
                outed$setCurrentCharFormat(outformat)
                ## any captured output (by product of evaluation)
                if (length(output))
                    outed$append(paste(output, collapse = "\n"))
                ## return value of evaluation (may need to be printed)
                if (inherits(evis, "try-error"))
                {
                    ## outed$moveCursor(Qt$QTextCursor$End)
                    outed$setCurrentCharFormat(msgformat)
                    outed$append(paste(strsplit(as.character(evis), "\n")[[1]],
                                       collapse = "\n")) # remove final newline
                }
                else if (evis$visible)
                {
                    if (html.preferred &&
                        !inherits(try(xtab <- xtable(evis$value), silent = TRUE),
                                  "try-error"))
                    {
                        ## FIXME: need something append-like (add to end)
                        outed$moveCursor(Qt$QTextCursor$End)
                        html.output <- capture.output(print(xtab, type = "html"))
                        outed$insertHtml(paste(html.output, collapse = "\n"))
                    }
                    else
                    {
                        text.output <- capture.output(evis$value)
                        outed$append(paste(text.output, collapse = "\n"))
                    }
                }
            }
            outed$moveCursor(Qt$QTextCursor$End)
        }
    }

    clearAct <- Qt$QAction(text = "Clear contents", parent = outed)
    clearAct$setShortcut(Qt$QKeySequence("Ctrl+L"))
    clearHandler <- function(checked) { outed$setPlainText("") }
    ## clearAct$setShortcutContext(Qt$Qt$WidgetShortcut) ## only triggered when widget has focus
    qconnect(clearAct, signal = "triggered", handler = clearHandler)
    outed$addAction(clearAct)


    ##qaction(desc = "Execute", shortcut = "Ctrl+Return", parent = ined1)
    runAct1 <- Qt$QAction(text = "Execute", parent = ined1)
    runAct1$setShortcut(Qt$QKeySequence("Ctrl+Return"))
    runHandler1 <- function(checked) { executeCode(ined1$plainText, mode = "input") }
    runAct1$setShortcutContext(Qt$Qt$WidgetShortcut) ## only triggered when widget has focus
    qconnect(runAct1, signal = "triggered", handler = runHandler1)
    ined1$addAction(runAct1)
    
    runAct2 <- Qt$QAction(text = "Execute selection", parent = ined2)
    runAct2$setShortcut(Qt$QKeySequence("Ctrl+Return"))
    ## runAct2 <- qaction(desc = "Execute selection", shortcut = "Ctrl+Return", parent = ined2)
    runHandler2 <- function(checked)
    {
        ## execute if selection exists, else select minimal parseable
        ## input starting backwords from current line
        sel <- qselectedText_QTextEdit(ined2)
        if (nzchar(sel))
        {
            executeCode(sel, mode = "edit")
            ined2$moveCursor(Qt$QTextCursor$StartOfLine)
            ined2$moveCursor(Qt$QTextCursor$Down)
        }
        else
        {
            oldCursor <- ined2$textCursor()
            ## oldPos <- ined2$textCursor()$position()
            ined2$moveCursor(Qt$QTextCursor$EndOfLine)
            ined2$moveCursor(Qt$QTextCursor$StartOfLine, Qt$QTextCursor$KeepAnchor)
            parseable <- !is(try(parse(text = qselectedText_QTextEdit(ined2)), silent = TRUE), "try-error")
            reached0 <- ined2$textCursor()$position() == 0L ## qcursorPosition(ined2) == 0L
            while (!parseable && !reached0)
            {
                ined2$moveCursor(Qt$QTextCursor$Up, Qt$QTextCursor$KeepAnchor)
                ## qmoveCursor(ined2, "up", select = TRUE)
                parseable <- !is(try(parse(text = qselectedText_QTextEdit(ined2)),
                                     silent = TRUE), "try-error")
                reached0 <- ined2$textCursor()$position() == 0L
            }
            if (!parseable) ## qsetCursorPosition(ined2, oldPos) ## restore
                ined2$setTextCursor(oldCursor)
        }
    }
    runAct2$setShortcutContext(Qt$Qt$WidgetShortcut) ## only triggered when widget has focus
    qconnect(runAct2, signal = "triggered", handler = runHandler2)
    ined2$addAction(runAct2)
    
    ## ## add action for text-completion (not yet in edit mode)
    ## compAct1 <- Qt$QAction(text = "Complete", parent = ined1)
    ## compAct1$setShortcut(Qt$QKeySequence("Ctrl+I"))
    ## ## qaction(desc = "Complete", shortcut = "Ctrl+I", parent = ined1)
    ## compAct1$setShortcutContext(Qt$Qt$WidgetShortcut)
    ## compHandler1 <- function(checked) {
    ##     comps <- tryComplete(text = ined1$plainText, ined1$textCursor()$position())
    ##     ined1$insertPlainText(comps$addition)
    ##     msg$text <-
    ##         if (nzchar(comps$addition) || any(nzchar(comps$comps))) paste(comps$comps)
    ##         else "No completions."
    ## }
    ## qconnect(compAct1, signal = "triggered", handler = compHandler1)
    ## ined1$addAction(compAct1)

    ## save file in edit mode
    saveFileName <- NULL
    saveAct <- Qt$QAction(text = "Save", parent = ined2)
    saveAct$setShortcut(Qt$QKeySequence("Ctrl+S"))
    saveAct$setShortcutContext(Qt$Qt$WidgetShortcut)
    saveHandler <- function(checked) {
        if (is.null(saveFileName))
        {
            file <- qfile.choose(caption = "Choose output file", filter = "R code (*.R *.r *.S *.s);;All files (*)", allow.new = TRUE)
            if (nzchar(file)) saveFileName <<- file
        }
        if (!is.null(saveFileName)) cat(ined2$plainText, file = saveFileName)
    }
    qconnect(saveAct, signal = "triggered", handler = saveHandler)
    ined2$addAction(saveAct)

    ## load file in edit mode
    loadAct <- Qt$QAction(text = "Append contents from file", parent = ined2)
    loadAct$setShortcut(Qt$QKeySequence("Ctrl+O"))
    loadAct$setShortcutContext(Qt$Qt$WidgetShortcut)
    loadHandler <- function(checked) {
        file <- qfile.choose(caption = "Choose file", filter = "R code (*.R *.r *.S *.s);;All files (*)", allow.new = FALSE)
        if (nzchar(file)) {
            ined2$appendPlainText(paste(readLines(file), collapse = "\n"))
            saveFileName <<- NULL
        }
    }
    qconnect(loadAct, signal = "triggered", handler = loadHandler)
    ined2$addAction(loadAct)

    ## return containing splitter
    container$resize(600, 400)
    container$setWindowTitle(title)
    container
}



## ## text zoom: doesn't wok because font reset during each eval
## zoominAct1 <- qaction(desc = "Increase text size", shortcut = "Ctrl++", parent = ined1)
## zoominAct1$shortcutContext <- 0 ## only triggered when widget has focus
## qconnect(zoominAct1, signal = "triggered", handler = function() ined1$zoomIn())
## qaddAction(ined1, zoominAct1)

## zoomoutAct1 <- qaction(desc = "Decrease text size", shortcut = "Ctrl+-", parent = ined1)
## zoomoutAct1$shortcutContext <- 0 ## only triggered when widget has focus
## qconnect(zoomoutAct1, signal = "triggered", handler = function() ined1$zoomOut())
## qaddAction(ined1, zoomoutAct1)

## zoominAct2 <- qaction(desc = "Increase text size", shortcut = "Ctrl++", parent = ined2)
## zoominAct2$shortcutContext <- 0 ## only triggered when widget has focus
## qconnect(zoominAct2, signal = "triggered", handler = function() ined2$zoomIn())
## qaddAction(ined2, zoominAct2)

## zoomoutAct2 <- qaction(desc = "Decrease text size", shortcut = "Ctrl+-", parent = ined2)
## zoomoutAct2$shortcutContext <- 0 ## only triggered when widget has focus
## qconnect(zoomoutAct2, signal = "triggered", handler = function() ined2$zoomOut())
## qaddAction(ined2, zoomoutAct2)



##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param env 
##' @param ... 
##' @param history 
##' @param eval.on.newline 
##' @param incolor 
##' @param outcolor 
##' @param msgcolor 
##' @return 
##' @author Deepayan Sarkar
qreplu <- function(env = .GlobalEnv,
                   ...,
                   history = TRUE,
                   eval.on.newline = TRUE,
                   incolor = "red",
                   outcolor = "blue",
                   msgcolor = "black")
{
    informat <- Qt$QTextCharFormat()
    informat$setForeground(qbrush(incolor))
    outformat <- Qt$QTextCharFormat()
    outformat$setForeground(qbrush(outcolor))
    msgformat <- Qt$QTextCharFormat()
    msgformat$setForeground(qbrush(msgcolor))
    ## ed: REPL-like mode, type and execute code
    ed <- qeditor(richtext = FALSE, rsyntax = FALSE, ..., comp.tooltip = FALSE)
    
    ## messages
    msg <- Qt$QLabel("")
    msg$wordWrap <- TRUE
    ## container
    container <- Qt$QSplitter(Qt$Qt$Vertical)
    container$addWidget(ed)
    container$addWidget(msg)
    container$setStretchFactor(0L, 10L)
    container$setStretchFactor(1L, 0L)
    msg$text <- sprintf("%s Type code, press Ctrl+Return to evaluate",
                        capture.output(env))

    ed$setContextMenuPolicy(Qt$Qt$ActionsContextMenu)

    ed$appendHtml("<p><i>Welcome to qreplu().</i></p> <br>")

    ## ed$setCurrentFont(font)
    ed$setCurrentCharFormat(informat)
    ed$insertPlainText("\n> ")
    ed$moveCursor(Qt$QTextCursor$End)
    lastPosition <- ed$textCursor()$position()

    evalenv <- env
    setEnv <- function(e)
    {
        evalenv <<- e
        m <- capture.output(print(e))
        ## ed$setTextColor(messagecolor)
        ## ed$insertPlainText(sprintf("\n\nChanged active environment to %s.\n", m))
        ## ed$setCurrentFont(font)
        ## ed$setTextColor(incolor)
        ## ed$insertPlainText("\n> ")
        ## ed$moveCursor(Qt$QTextCursor$End)
        ## lastPosition <<- ed$textCursor()$position()
        msg$text <- sprintf("\n\nChanged active environment to %s.\n", m)
   }

    setReadOnlyMode <- function(undo = FALSE)
    {
        ## ed$setReadOnly(ed$textCursor()$position() < lastPosition)
        if (ed$textCursor()$position() < lastPosition)
        {
            if (undo) ed$undo()
            ic <- ed$textCursor()
            ic$setPosition(lastPosition)
            ed$setTextCursor(ic)
        }
    }
    
    qconnect(ed, "cursorPositionChanged", setReadOnlyMode)
    qconnect(ed, "textChanged", setReadOnlyMode, user.data = TRUE)
    msgCompletions <- function(comps) {
        msg$text <- if (nzchar(comps)) comps else "No completions."
    }
    qconnect(ed, "completionsAvailable", msgCompletions)

    ## History mechanism: Ctrl+Up/Down.  Unused if history=FALSE.
    lhist <- list() # list of commands
    ihist <- 0L     # index
    cur_command <- ""
    if (history)
    {
        addToHistory <- function(s)
        {
            n <- length(lhist)
            if (n == 0 || s != lhist[[n]]) lhist[[n + 1L]] <<- s
        }
        prevHistoryItem <- function()
        {
            n <- length(lhist)
            ihist <<-
                if (ihist == 0) n
                else ihist - 1L 
            if (ihist == 0) cur_command else lhist[[ihist]]
        }
        nextHistoryItem <- function()
        {
            n <- length(lhist)
            ihist <<- ihist + 1L
            if (ihist == n + 1L) ihist <<- 0L
            if (ihist == 0) cur_command else lhist[[ihist]]
        }
        showPrevItem <- function(...)
        {
            ## str(list(lhist = lhist, ihist = ihist, cur_command = cur_command))
            ed$moveCursor(Qt$QTextCursor$End)
            ic <- ed$textCursor()
            ic$setPosition(lastPosition, Qt$QTextCursor$KeepAnchor)
            ed$setTextCursor(ic)
            ## If editing something new, make that current before changing
            if (ihist == 0L) cur_command <<- qselectedText_QTextEdit(ed)
            ic$removeSelectedText()
            ed$insertPlainText(prevHistoryItem())
        }
        showNextItem <- function(...)
        {
            ## str(list(lhist = lhist, ihist = ihist, cur_command = cur_command))
            ed$moveCursor(Qt$QTextCursor$End)
            ic <- ed$textCursor()
            ic$setPosition(lastPosition, Qt$QTextCursor$KeepAnchor)
            ed$setTextCursor(ic)
            ## If editing something new, make that current before changing
            if (ihist == 0L) cur_command <<- qselectedText_QTextEdit(ed)
            ic$removeSelectedText()
            ed$insertPlainText(nextHistoryItem())
        }
        
        prevHistoryAct <- Qt$QAction(text = "Previous command", parent = ed)
        prevHistoryAct$setShortcut(Qt$QKeySequence("Ctrl+Up"))
        prevHistoryAct$setShortcutContext(Qt$Qt$WidgetShortcut)
        qconnect(prevHistoryAct, signal = "triggered", handler = showPrevItem)
        ed$addAction(prevHistoryAct)
        
        nextHistoryAct <- Qt$QAction(text = "Next command", parent = ed)
        nextHistoryAct$setShortcut(Qt$QKeySequence("Ctrl+Down"))
        nextHistoryAct$setShortcutContext(Qt$Qt$WidgetShortcut)
        qconnect(nextHistoryAct, signal = "triggered", handler = showNextItem)
        ed$addAction(nextHistoryAct)
    }
    
    processCodeIfReady <- function()
    {
        ed$moveCursor(Qt$QTextCursor$End)
        ic <- ed$textCursor()
        ic$setPosition(lastPosition, Qt$QTextCursor$KeepAnchor)
        ed$setTextCursor(ic)
        intext <- qselectedText_QTextEdit(ed)
        parseable <- !is(try(parse(text = intext), silent = TRUE), "try-error")
        if (parseable)
        {
            ihist <<- 0L
            ic$removeSelectedText()
            executeCode(intext)
        }
        else 
        {
            msg$text <- "Invalid or incomplete input"
            ed$moveCursor(Qt$QTextCursor$End)
        }
    }

    ## add action to execute code

    runAct <- Qt$QAction(text = "Execute", parent = ed)
    runAct$setShortcut(Qt$QKeySequence("Ctrl+Return"))
    ## runAct$setShortcut(Qt$QKeySequence(Qt$Qt$Key_Return))
    runAct$setShortcutContext(Qt$Qt$WidgetShortcut) ## only triggered when widget has focus
    qconnect(runAct, signal = "triggered", handler = processCodeIfReady)
    ed$addAction(runAct)

    ## Also do this everytime Enter is pressed
    if (eval.on.newline)
        qconnect(ed, signal = "enterPressed", handler = processCodeIfReady)
    
    closeAct <- Qt$QAction(text = "Close REPL", parent = ed)
    closeAct$setShortcut(Qt$QKeySequence("Ctrl+Q"))
    closeAct$setShortcutContext(Qt$Qt$WidgetShortcut) ## only triggered when widget has focus
    qconnect(closeAct, signal = "triggered", handler = container$close)
    ed$addAction(closeAct)
    
    ## function to perform code execution

    executeCode <- function(text)
    {
        if (!nzchar(text)) return (FALSE)
        msg$text <- "Evaluating..."
        msg$update()
        pe <- tryParseEval(text = text, env = env)
        if (is(pe, "try-error"))
        {
            msg$text <- paste(strsplit(as.character(pe), "\n", fixed = TRUE)[[1]], collapse = "\\n")
            return (FALSE)
        }
        else
        {
            if (history) addToHistory(text)
            ## code already deleted from prompt onwards. Now to remove prompt as well.
            
            lastPosition <<- 0 # to allow deleting prompt
            ed$moveCursor(Qt$QTextCursor$Left, Qt$QTextCursor$KeepAnchor)
            ed$moveCursor(Qt$QTextCursor$Left, Qt$QTextCursor$KeepAnchor)
            ed$moveCursor(Qt$QTextCursor$Left, Qt$QTextCursor$KeepAnchor)
            ed$textCursor()$removeSelectedText()

            for (i in seq_along(pe))
            {
                ein <- pe[[i]]$ein
                output <- pe[[i]]$output
                evis <- pe[[i]]$evis
                ## input
                ed$setCurrentCharFormat(informat)
                ed$appendPlainText(ein)
                ## ed$insertPlainText("\n")
                ## output
                ed$setCurrentCharFormat(outformat)
                ## any captured output (by-product of evaluation)
                if (length(output))
                    ed$appendPlainText(paste(output, collapse = "\n"))
                ## return value of evaluation (may need to be printed)
                if (inherits(evis, "try-error"))
                {
                    ## ed$moveCursor(Qt$QTextCursor$End)
                    ed$setCurrentCharFormat(msgformat)
                    ed$appendPlainText(paste(strsplit(as.character(evis), "\n")[[1]],
                                             collapse = "\n")) # remove final newline
                }
                else if (evis$visible)
                {
                    text.output <- capture.output(evis$value)
                    ed$appendPlainText(paste(text.output, collapse = "\n"))
                }
            }
            ## ed$setTextColor(incolor)
            ed$setCurrentCharFormat(informat)
            ed$appendPlainText("> ")
            ed$moveCursor(Qt$QTextCursor$End)
            lastPosition <<- ed$textCursor()$position()
            msg$text <- sprintf("%s Type code, press Ctrl+Return to evaluate",
                                capture.output(evalenv))
        }
    }

    ## return containing splitter
    container$resize(600, 400)
    attr(container, "setEnv") <- setEnv
    container
}


