###################################################
### chunk number 1: 
###################################################
#line 3 "ex-gWidgets-calculator.Rnw"
## A calculator layout with gWidgets
library(gWidgets)


###################################################
### chunk number 2: 
###################################################
#line 20 "ex-gWidgets-calculator.Rnw"
buttons <- rbind(c(7:9, "(", ")"),
                 c(4:6, "*", "/"),
                 c(1:3, "+", "-"))
#
w <- gwindow("glayout for a calculator", visible=FALSE)
g <- ggroup(cont=w, expand=TRUE, horizontal=FALSE)
tbl <- glayout(cont=g, spacing=2)
                                        
tbl[1, 1:5, anchor=c(-1,0)] <-          # span 5 columns
  (eqnArea <- gedit("", cont=tbl))
tbl[2, 1:5, anchor=c(1,0)] <- 
  (outputArea <- glabel("", cont=tbl))
#
bList <- list()
for(i in 3:5) {
  for(j in 1:5) {
    val <- buttons[i-2, j]
    tbl[i,j] <- (bList[[val]] <- gbutton(val, cont=tbl))
  }
}
tbl[6,2] <- (bList[["0"]] <- gbutton("0", cont=tbl))
tbl[6,3] <- (bList[["."]] <- gbutton(".", cont=tbl))
tbl[6,4:5] <- (eqButton <- gbutton("=", cont=tbl))
#
visible(w) <- TRUE


###################################################
### chunk number 3: 
###################################################
#line 53 "ex-gWidgets-calculator.Rnw"
addButton <- function(h, ...) {
  curExpr <- svalue(eqnArea)
  newChar <- svalue(h$obj)              # the button's value
  svalue(eqnArea) <- paste(curExpr, newChar, sep="")
  svalue(outputArea) <- ""              # clear label 
}
#out <- sapply(bList, function(i) 
#              addHandlerChanged(i, handler=addButton))
sapply(bList, addHandlerChanged, handler=addButton)


###################################################
### chunk number 4: 
###################################################
#line 67 "ex-gWidgets-calculator.Rnw"
addHandlerClicked(eqButton, handler = function(h,...) {
  curExpr <- svalue(eqnArea)
  out <- try(capture.output(eval(parse(text=curExpr))), 
             silent=TRUE)
  if(inherits(out, "try-error")) {
    galert("There is an error")
  } else {
    svalue(outputArea) <- out
    svalue(eqnArea) <- ""            # restart
  }
})
                  


