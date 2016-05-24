

###################################################
### code chunk number 66: padding-pady-ipady
###################################################
## Code to make padding/pady/ipady figure
## padding
w <- tktoplevel();
tkwm.title(w, "Padding/pady/ipady")
f1 <- ttkframe(w, padding = c(3,3,12,12))  # Some breathing room
tkpack(f1, expand = TRUE, fill = "both")     # For window expansion


f <- ttkframe(f1,padding = 20, border = 1);
tkpack(f, side = "left", expand = TRUE, fill = "x")
tkpack(ttkbutton(f, text = "padding"))
tkpack(ttkbutton(f, text = "padding"))

tkpack(ttkseparator(f1, orient = "vertical"), expand = TRUE, fill = "y", side = "left")

##  pady outside border
f <- ttkframe(f1, border = 1)
tkpack(f, side = "left", expand = TRUE, fill = "both")
tkpack(ttkbutton(f, text = "pady"), pady = 10)
tkpack(ttkbutton(f, text = "pady"), pady = 10)

tkpack(ttkseparator(f1, orient = "vertical"), expand = TRUE, fill = "y", side = "left")
## ipady within border
f <- ttkframe(f1, border = 1)
tkpack(f, side = "left", expand = TRUE, fill = "both")
tkpack(ttkbutton(f, text = "ipady"), ipady = 10)
tkpack(ttkbutton(f, text = "ipady"), ipady = 10)
