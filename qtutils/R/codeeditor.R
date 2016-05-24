

## Class that subclasses QPlainTextEdit and overrides TAB key etc

qsetClass("RCodeEditor", Qt$QPlainTextEdit,
          constructor = function(family = "monospace", pointsize = 14,
                                 underscore.assign = FALSE, comp.tooltip = TRUE) {
              this$setFont(qfont(family = family, pointsize = pointsize))
              this$centerOnScroll <- FALSE
              this$setLineWrapMode(Qt$QTextEdit$NoWrap)
              ## this$tabMode <- "complete"
              this$uassign <- underscore.assign
              this$ctooltip <- comp.tooltip
          })



isblank <- function(s) 
{
    # is s composed of only spaces?
    identical(unique(charToRaw(s)), as.raw(32))
}

firstNonBlank <- function(s)
{
    n <- countLeadingSpaces(s)
    if (nchar(s) > n) substring(s, n+1, n+1) else ""
}


findLastUnmatched <- function(s)
{
    ## s is like "(" "[" "]" "(" ")" "{" "(" ")" "{" "}"
    ## want to find last unmatched left-open
    nbrace <- nbracket <- nparen <- 0L
    for (i in rev(seq_along(s)))
    {
        switch(s[i],
               "}" = { nbrace <- nbrace + 1L },
               "]" = { nbracket <- nbracket + 1L },
               ")" = { nparen <- nparen + 1L },
               "{" = { if (nbrace == 0L) return(i) else nbrace <- nbrace - 1L },
               "[" = { if (nbracket == 0L) return(i) else nbracket <- nbracket - 1L },
               "(" = { if (nparen == 0L) return(i) else nparen <- nparen - 1L })
    }
    return (0L)
}

countLeadingSpaces <- function(x) # x is already split
{
    if (length(x) == 1) x <- strsplit(x, "")[[1]]
    spaces <- (x == " ")
    if (length(spaces))
    {
        nsp <- which(!spaces)
        if (length(nsp))
        {
            nsp[1] - 1L
        }
        else length(spaces) # all spaces
    }
    else 0
}

computeTabSpaces <- function(text)
{
    ## want intended number of spaces for last line.
    ## Assume current call goes back to last line with 0 leading spaces
    ll <- strsplit(text, "\n")[[1]]
    spaces <- sapply(strsplit(ll, ""), countLeadingSpaces)
    goBackTo <- tail(which(spaces == 0), 1)
    if (length(goBackTo))
    {
        useText <- if (goBackTo == 1L) ll else tail(ll, -(goBackTo-1L))
        ## cat(useText, sep = "\n")
        ## now count number of {, }, ), (
        useTextCombined <- paste(useText, collapse = "\n")
        occurencePattern <- function(p)
        {
            p <- gregexpr(p, useTextCombined, fixed = TRUE)[[1]]
            if (all(p == -1)) integer(0)
            else as.numeric(p)
        }
        linebreaks <- occurencePattern("\n")
        pkey <- c("(", ")", "[", "]", "{", "}")
        poccur <- sapply(pkey, occurencePattern, simplify = FALSE)
        ## str(poccur)
        ord <- order(unlist(poccur))
        osorted <- rep(pkey, sapply(poccur, length))[ord]
        ## print(osorted)
        lastUnmatched <- findLastUnmatched(osorted)
        if (lastUnmatched == 0L) return(0) # apparently complete
        ## what is it and where?
        posUnmatched <- unlist(poccur)[ord][lastUnmatched]
        switch(osorted[lastUnmatched],
               "(" = ,
               "[" = {
                   ## Indent to that position + 1.
                   ## To find position, see difference from previous newline
                   if (any(linebreaks < posUnmatched)) 
                   {
                       return(posUnmatched -
                              linebreaks[tail(which(linebreaks < posUnmatched), 1)])
                   }
                   else return(posUnmatched) # unmatched bracket in first line
               },
               "{" = {
                   ## Indent to indent of line containing that { + 4
                   ## Find line number by counting preceding newlines
                   linenum <- sum(linebreaks < posUnmatched) + 1L
                   return(countLeadingSpaces(strsplit(useText[linenum], "")[[1]]) + 4)
               })
        1
    }
    else # no line starts at beginning 
        0
}


qsetMethod("currentLine", RCodeEditor,
           function(uptocursor = TRUE, remove = FALSE) {
               cc <- textCursor()
               if (!uptocursor)
                   cc$movePosition(Qt$QTextCursor$EndOfLine, Qt$QTextCursor$MoveAnchor)
               cc$movePosition(Qt$QTextCursor$StartOfLine, Qt$QTextCursor$KeepAnchor)
               ans <- cc$selection()$toPlainText() # can be NULL
               if (remove) cc$removeSelectedText()
               if (is.null(ans)) "" else ans
           })

qsetMethod("currentDocument", RCodeEditor,
           function(uptocursor = TRUE, remove = FALSE) {
               cc <- textCursor()
               if (!uptocursor)
                   cc$movePosition(Qt$QTextCursor$End, Qt$QTextCursor$MoveAnchor)
               cc$movePosition(Qt$QTextCursor$Start, Qt$QTextCursor$KeepAnchor)
               ans <- cc$selection()$toPlainText()
               if (remove) cc$removeSelectedText()
               if (is.null(ans)) "" else ans
           })


qsetMethod("indentCurrentLine", RCodeEditor,
           ## return: TRUE if something happened, FALSE if nothing to do.
           function() {
               ## how many leading blanks in current line?
               cur <- countLeadingSpaces(currentLine(uptocursor = FALSE))
               ## how many should it have? Go to beginning and decide
               cc <- textCursor()
               cc$movePosition(Qt$QTextCursor$StartOfLine, Qt$QTextCursor$MoveAnchor)
               cc$movePosition(Qt$QTextCursor$Start, Qt$QTextCursor$KeepAnchor)
               selText <- cc$selection()$toPlainText()
               wanted <- 
                   if (is.null(selText) || !nzchar(selText)) # first line
                       0L
                   else
                       computeTabSpaces(cc$selection()$toPlainText())
               ## adjust if first non-blank char is } or ) or ]
               adj <- switch(firstNonBlank(currentLine(uptocursor = FALSE)),
                             "}" = 4L,
                             ")" = , "]" = 1L, 0L)
               wanted <- max(0L, wanted - adj)
               ## base::print(c(cur, wanted))
               if (wanted == cur)
               {
                   ## if cursor in initial blank space, move it
                   ## forward. Otherwise nothing to do.
                   cursorcol <- textCursor()$positionInBlock()
                   if (cursorcol >= cur)
                       return(FALSE)
                   else
                   {
                       for (i in seq_len(cur - cursorcol))
                           moveCursor(Qt$QTextCursor$Right)
                       return(TRUE)
                   }
               }
               cc <- textCursor()
               cc$movePosition(Qt$QTextCursor$StartOfLine, Qt$QTextCursor$MoveAnchor)
               for (i in seq_len(cur)) cc$deleteChar()
               for (i in seq_len(wanted)) cc$insertText(" ")
               ## Could have done this instead, but then cursor would not be at indent.
               ## if (wanted > cur)
               ## {
               ##     for (i in seq_len(wanted-cur)) cc$insertText(" ")
               ## }
               ## else if (wanted < cur)
               ## {
               ##     for (i in seq_len(cur-wanted)) cc$deleteChar()
               ## }
               return (TRUE)
           })


qsetMethod("cursorGlobalPosition", RCodeEditor,
           function() {
               cm <- as.matrix(cursorRect())
               globalpos <- mapToGlobal(pos)
               ## qpoint(window()$geometry$x() + cm[2, 1],
               ##        window()$geometry$y() + cm[2, 2])
               qpoint(globalpos$x() + cm[2, 1],
                      globalpos$y() + cm[2, 2])
           })

qsetSignal("completionsAvailable(QString character)", RCodeEditor)

qsetSignal("enterPressed()", RCodeEditor)

qsetMethod("keyPressEvent", RCodeEditor,
           function(e) {
               if (ctooltip && Qt$QToolTip$isVisible()) Qt$QToolTip$hideText()
               if (e$modifiers() == Qt$Qt$ControlModifier)
               {
                   ek <- e$key()
                   if (ek == Qt$Qt$Key_E) {
                       moveCursor(Qt$QTextCursor$EndOfLine)
                       return()
                   }
                   else if (ek == Qt$Qt$Key_A) {
                       moveCursor(Qt$QTextCursor$StartOfLine)
                       return()
                   }
                   else return(super("keyPressEvent", e))
               }
               et <- e$text()
               ## base::print(et)
               if (et == "\r") # Enter/Return
               {
                   insertPlainText("\n")
                   ## spaces <- computeTabSpaces(document()$toPlainText())
                   spaces <- computeTabSpaces(currentDocument(uptocursor = TRUE))
                   insertPlainText(base::paste(rep(" ", spaces), collapse = ""))
                   ## emit signal to indicate new line (REPL may want to execute line)
                   enterPressed()
               }
               else if (et == "}")
               {
                   ## go back by 4 if this is the first character of line
                   cl <- currentLine(uptocursor = TRUE)
                   if (nchar(cl) >= 4 && isblank(cl))
                   {
                       for (i in 1:4) textCursor()$deletePreviousChar()
                   }
                   insertPlainText(et)
               }
               else if (et %in% c(")", "]"))
               {
                   ## go back by 4 if this is the first character of line
                   cl <- currentLine(uptocursor = TRUE)
                   if (nchar(cl) >= 1 && isblank(cl))
                   {
                       textCursor()$deletePreviousChar()
                   }
                   insertPlainText(et)
               }
               else if (et == "\t") 
               {
                   ## TAB mode: complete/indent/indent-then-complete
                   ## we will do only third for now.  Others are easy with flags. 
                   if (!indentCurrentLine()) 
                   {
                       comps <- tryComplete(currentLine(uptocursor = TRUE))
                       if (comps$addition == "" && nzchar(comps$comps)) 
                       {
                           ## indicate available comps$comps
                           if (ctooltip)
                               Qt$QToolTip$showText(cursorGlobalPosition(),
                                                    base::paste(strwrap(comps$comps, 80),
                                                                collapse = "\n"),
                                                    this, this$rect)
                           else # emit signal
                               completionsAvailable(comps$comps)
                       }
                       else insertPlainText(comps$addition)
                   }
               }
               else if (uassign && et == "_")
                   insertPlainText(" <- ")
               else
               {
                   ## base::print(et)
                   super("keyPressEvent", e)
               }
           })



## rm(RCodeEditor)
## (foo <- RCodeEditor())


