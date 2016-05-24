
###################################################
### code chunk number 61: ex-gWidgets-gedit-validation.Rnw:14-18
###################################################
window <- gwindow("Validation example")
lyt <- glayout(cont = window)
lyt[1,1] <- "R expression:"
lyt[1,2] <- (entry <- gedit("", cont = lyt))


###################################################
### code chunk number 62: ex-gWidgets-gedit-validation.Rnw:33-38
###################################################
require(evaluate)
isValid <- function(e) {
  out <- try(evaluate:::evaluate(e), silent=TRUE)
  !(inherits(out, "try-error") ||  is(out[[2]], "error"))
}


###################################################
### code chunk number 63: validate
###################################################
addHandlerChanged(entry, handler = function(h,...) {
  cur_val <- svalue(entry)
  if(isValid(cur_val)) {
    font(entry) <- c(color = "black")
  } else {
    font(entry) <- c(color = "red")
  }
})

