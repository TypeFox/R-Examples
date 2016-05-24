
###################################################
### code chunk number 69: anchorPoints
###################################################
## Example of anchor arguments
## not shown in text
w <- tktoplevel()
tkwm.title(w, "anchor arguments")
f <- ttkframe(w, padding = c(3,3,12,12)); tkpack(f, expand = TRUE, fill = "both")
tr <- ttkframe(f); mr <- ttkframe(f); br <- ttkframe(f)
sapply(list(tr, mr, br), tkpack, expand = TRUE, fill = "both")
tkpack(ttklabel(tr, text = "nw"), anchor = "nw", expand = TRUE, side = "left")
tkpack(ttklabel(tr, text = "n"),  anchor = "n",  expand = TRUE, side = "left")
tkpack(ttklabel(tr, text = "ne"), anchor = "ne", expand = TRUE, side = "left")
##
tkpack(ttklabel(mr, text = "w"),       anchor = "w",       expand = TRUE, side = "left")
tkpack(ttklabel(mr, text = "center"),  anchor = "center",  expand = TRUE, side = "left")
tkpack(ttklabel(mr, text = "e"),       anchor = "e",       expand = TRUE, side = "left")
##
tkpack(ttklabel(br, text = "sw"), anchor = "sw", expand = TRUE,  side = "left")
tkpack(ttklabel(br, text = "s"),  anchor = "s",  expand = TRUE,  side = "left")
tkpack(ttklabel(br, text = "se"), anchor = "se", expand = TRUE,  side = "left")

