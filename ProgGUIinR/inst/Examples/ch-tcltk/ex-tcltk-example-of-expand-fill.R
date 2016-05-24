
###################################################
### code chunk number 70: expandFill
###################################################
window <- tktoplevel()
tkwm.title(window, "Expand/Fill arguments")
frame <- ttkframe(window, padding = c(3,3,12,12))
tkpack(frame, expand = TRUE, fill = "both")
##
pack_btn <- function(txt, ...) 
  tkpack(button <- ttkbutton(frame, text = txt), ...)
##
pack_btn("Top",    side="top",    expand=TRUE, fill="both") 
pack_btn("Bottom", side="bottom", expand=TRUE, fill="both") 
pack_btn("Left",   side="left",   expand=TRUE, fill="both") 
pack_btn("Right",  side="right",  expand=TRUE, fill="both") 


###################################################
### code chunk number 71: BasicContainers.Rnw:537-539
###################################################
children <- as.character(tkwinfo("children", frame))
sapply(children, tkpack.configure, fill = "none")


###################################################
### code chunk number 72: BasicContainers.Rnw:554-555 (eval = FALSE)
###################################################
tkwm.geometry(window, "")




###################################################
### code chunk number 71: BasicContainers.Rnw:537-539
###################################################
children <- as.character(tkwinfo("children", frame))
sapply(children, tkpack.configure, fill = "none")


###################################################
### code chunk number 72: BasicContainers.Rnw:554-555 (eval = FALSE)
###################################################
tkwm.geometry(window, "")


## Change to "x" fill:
children <- as.character(tkwinfo("children", frame))
sapply(children, tkpack.configure, fill = "x")
tkwm.geometry(window, "")

