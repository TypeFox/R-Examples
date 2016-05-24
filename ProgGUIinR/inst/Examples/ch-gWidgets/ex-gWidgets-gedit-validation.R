###################################################
### chunk number 1: 
###################################################
## Validation example
require(gWidgets)


###################################################
### chunk number 2: 
###################################################
w <- gwindow("Validation example")
validRegexpr <- "[[:digit:]]{3}-[[:digit:]]{4}"
tbl <- glayout(cont = w)
tbl[1,1] <- "Phone number (XXX-XXXX)"
tbl[1,2] <- (e <- gedit("", cont = tbl))
tbl[2,2] <- (b <- gbutton("submit", cont = tbl, 
                          handler=function(h,...) print("hi")))
## Blur is focus out event
addHandlerBlur(e, handler = function(h,...) {
  curVal <- svalue(h$obj)
  if(grepl(validRegexpr, curVal)) {
    font(h$obj) <- c(color="black")
  } else {
    focus(h$obj) <- TRUE
    font(h$obj) <- c(color="red")
  }
})


