#######  read in data  from Reventador Volcano, Ecuador, and plot
######  using PICK.GEN
options(demo.ask=FALSE)

data(KH)
STDLAB = c("DONE",  "zoom in", "zoom out", "refresh", "restore",
 "XTR", "SPEC", "SGRAM" ,"WLET", "FILT",  "Pinfo")

######## swig is a generic Picking program
######## Use left Mouse click to select traces and windows.
######## Use the buttons to operate on the selected traces/windows
########  To zoom, click twice with the left mouse on the traces, then right mouse zooms
########  Right Mouse with no left mouse clicks is equivalent to "Done"

xx =  swig( KH, sel=which(KH$COMPS == "V"), STDLAB = STDLAB)
