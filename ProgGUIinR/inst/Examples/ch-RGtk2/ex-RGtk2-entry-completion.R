### R code from vignette source 'ex-RGtk2-entry-completion.Rnw'

###################################################
### code chunk number 1: ex-RGtk2-entry-completion.Rnw:2-3
###################################################
require(RGtk2)


###################################################
### code chunk number 2: AppendWords
###################################################
entry <- gtkEntry(); completion <- gtkEntryCompletion()
entry$setCompletion(completion)


###################################################
### code chunk number 3: SetCompletion
###################################################
model <- rGtkDataFrame(state.name)
completion$setModel(model)
completion$setTextColumn(0)
completion['inline-completion'] <- TRUE
completion['popup-single-match'] <- FALSE


###################################################
### code chunk number 4: SetMatchFunc
###################################################
matchAnywhere <- function(completion, key, iter, user.data) {
  model <- completion$getModel()
  row_value <- model$getValue(iter, 0)$value
  key <- completion$getEntry()$getText() # case sensitivity
  grepl(key, row_value)
}
completion$setMatchFunc(matchAnywhere)


###################################################
### code chunk number 5: notShown
###################################################
## Our basic GUI is basic:
w <- gtkWindow(show=FALSE)
w$setTitle("Test of entry with completion")
w$add(entry)
w$showAll()


