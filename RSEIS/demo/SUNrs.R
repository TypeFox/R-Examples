#######  read in data sunspot  and plot
######  using swig
options(demo.ask=FALSE)

data(sunspots)

###  prepare the time series for input to RSEIS

####  step 1:  This puts it into the format used by Rsac
ES = prep1wig(wig=sunspots, dt=1/12, sta="STA", comp="CMP", units="UNITS"    )

####  this puts it into the format used by RSEIS:
####  step 2:
EH=prepSEIS(ES)


########  pop up the signals using swig
STDLAB = c("DONE",  "zoom in", "zoom out", "refresh", "restore",
 "XTR", "SPEC", "SGRAM" ,"WLET", "FILT",  "Pinfo")
######## 
######## 
######## 
######## swig is a generic Picking program
######## Use left Mouse click to select traces and windows.
######## Use the buttons to operate on the selected traces/windows
########  To zoom, click twice with the left mouse on the traces, then right mouse zooms
########  Right Mouse with no left mouse clicks is equivalent to "Done"
######## 

xx =  swig( EH, STDLAB = STDLAB)

