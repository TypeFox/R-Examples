### tcltk2 examples

### autoscroll
tclRequire("autoscroll")
tt <- tktoplevel()
scrl <- tkscrollbar(tt, orient = "v", command = function(...) tkyview(txt, ...))
txt <- tktext(tt, highlightthickness = 0, yscrollcommand = function(...) tkset(scrl, ...))
tkpack(scrl, side = "right", fill = "y")
tkpack(txt, side = "left", fill = "both", expand = 1)
tcl("::autoscroll::autoscroll", scrl)


### combobox: to eliminate!

### choosefont
### TODO


### ctext
### TODO


### cursor
### TODO


### swaplist
tclRequire("swaplist")
tt <- tktoplevel()
opts <- tclVar()
sl <- tcl("swaplist::swaplist", tt, opts, 1:9, c(1, 3, 5))
cat("You choose:", tclvalue(opts), "\n")
