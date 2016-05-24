

###################################################
### code chunk number 112: gtk-widget-entry-validate
###################################################
validatedEntry <- gtkEntry()
gSignalConnect(validatedEntry, "changed", function(entry) {
  text <- entry$getText()
  if (nzchar(gsub("[a-zA-Z]", "", text))) {
    entry$setIconFromStock("primary", "gtk-no")
    entry$setIconTooltipText("primary", 
                                 "Only letters are allowed")
  } else { 
    entry$setIconFromStock("primary", "gtk-yes")
    entry$setIconTooltipText("primary", NULL)
  }
})
validatedEntry$setIconFromStock("primary", "gtk-yes")


###################################################
### code chunk number 113: BasicComponents.Rnw:430-433
###################################################
w <- gtkWindow(show=FALSE)
w$add(validatedEntry)
w$showAll()

