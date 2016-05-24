# Test of widget 3.0

# The calendar widget (impossible to load)
tclRequire("widget::calendar")
tt <- tktoplevel()
db <- tkwidget(tt, "widget::calendar")
tkpack(db, fill = "both", expand = 1)


# The menuentry widget (OK)
tclRequire("widget::menuentry")
tt <- tktoplevel()
me <- tkwidget(tt, "widget::menuentry")
mnu <- tkmenu(me, tearoff = 0)
foo <- tclVar()
tkadd(mnu, "radiobutton", label = "Name", variable = foo, value = "name")
tkadd(mnu, "radiobutton", label = "Abstract", variable = foo, value = "abstract")
tkadd(mnu, "separator")
tkadd(mnu, "radiobutton", label = "Name and Abstract", variable = foo, value = "name abstract")
tkconfigure(me, menu = mnu)
tkpack(me, fill = "x", expand = 1, padx = 4, pady = 4)
# Change the selection in the menu, and then:
tclvalue(foo)


# The panelframe (does not work!?)
tclRequire("widget::panelframe")
tt <- tktoplevel()
pf <- tkwidget(tt, "widget::panelframe", text = "My Panel")
sf <- tkframe(pf, padx = 4, pady = 4)
txt <- tktext(sf)
tkpack(txt, fill = "both", expand = 1)
tcl(pf, "setwidget", sf)
tkpack(sf, fill = "both", expand =  1, padx = 4, pady = 4)

# The superframe (what is this?)
tclRequire("widget::superframe")
tt <- tktoplevel()
sf <- tkwidget(tt, "widget::superframe", text = "Superframe:")
tkpack(sf)
tkpack(tk2button(sf, text = "A button"))


