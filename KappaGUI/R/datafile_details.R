datafile_details <-
function() {
# Details de mise en forme sont donnes ici

myenv = new.env() # nouvel environnement qui inclura la variable globale tcltk

start <- tktoplevel()
tkwm.title(start,"About your data file")

# Checkboxes
frameOverall <- tkframe(start)

frame_header <- tkframe(frameOverall)
but_header <- tkcheckbutton(frame_header)
val_header <- tclVar("1")
tkconfigure(but_header,variable=val_header)
tkgrid(tklabel(frame_header,text="the first row contain the names of the variables?"),but_header)

frame_rownames <- tkframe(frameOverall)
but_rownames <- tkcheckbutton(frame_rownames)
val_rownames <- tclVar("1")
tkconfigure(but_rownames,variable=val_rownames)
tkgrid(tklabel(frame_rownames,text="the first column contain the names of the individuals?"),but_rownames)

##
# Affichage
##
tkgrid(tklabel(start,text=paste("In your data file, does...", sep=""),font="bold"))
tkgrid(frame_header)
tkgrid(frame_rownames)

##
# Button OK :
##
PressedOK <- function(){
#caract_files <<- c(tclvalue(val_header), tclvalue(val_rownames)) # old !
assign("caract_files", c(tclvalue(val_header), tclvalue(val_rownames)), envir=myenv)
tkgrab.release(start)
tkdestroy(start)
}


OK.but <- tkbutton(start,text="OK",command=PressedOK,relief="groove",borderwidth=4,width=11,bg="green4",fg="white")
tkgrid(frameOverall)
tkgrid(OK.but)
tkfocus(start)
tkwait.window(OK.but)

resul=mget("caract_files", envir=myenv)$caract_files
return(as.logical(as.numeric(resul)))
}
