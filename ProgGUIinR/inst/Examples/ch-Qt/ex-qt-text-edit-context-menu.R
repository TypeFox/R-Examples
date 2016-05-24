
###################################################
### code chunk number 419: QTextEditWithMenu
###################################################
qsetClass("QTextEditWithCompletions", Qt$QTextEdit)
##
qsetMethod("contextMenuEvent", QTextEditWithCompletions, 
           function(event) 
           {
  menu <- this$createStandardContextMenu()
  if(this$textCursor()$hasSelection()) {
    selection <- this$textCursor()$selectedText()
    completions <- utils:::matchAvailableTopics(selection)
    completions <- setdiff(completions, selection)
    if(length(completions) > 0 && length(completions) < 25) {
      menu$addSeparator()                  # add actions
      sapply(completions, function(completion) {
        action <- Qt$QAction(completion, this)
        qconnect(action, "triggered", function(checked) {
          insertPlainText(completion)
        })
        menu$addAction(action)
      })
    }
  }
  menu$exec(event$globalPos())
})
text_edit <- QTextEditWithCompletions()


###################################################
### code chunk number 420: raise (eval = FALSE)
###################################################
text_edit$show()
text_edit$raise()
