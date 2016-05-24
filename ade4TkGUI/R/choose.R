################################
#
# misc user dialog functions
#
################################

################################
# Function to choose the number of axes
################################
"chooseaxes" <- function(eig, rank1)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Number of axes")
	done <- tclVar(0)

	nfvar <- tclVar()
	nf <- 0

	"getnf" <- function()
	{
		if (tclvalue(nfvar) != "") {
			nf  <<- parse(text=tclvalue(nfvar))[[1]]
		} else nf <<- 2
		tkdestroy(tf)
		return(nf)
	}

	plotEig(eig[1:rank1], col.plot="grey", paxes.draw=TRUE, paxes.x.draw=FALSE)

	frame1 <- tkframe(tf, relief="groove", borderwidth=2)	
	nf.entry <- tkentry(frame1, textvariable=nfvar)
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function() getnf())
	tkgrid(tklabel(frame1,text="Select the number of axes : "), nf.entry)
	tkgrid(submit.but)
	tkpack(frame1, fill = "x")
	
	#tkconfigure(nf.entry, state="disabled")

	tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
	tkbind(tf, "<KeyPress-Return>", function() getnf())
	tkfocus(nf.entry)
	
	tkwait.variable(done)
	if(tclvalue(done) == "2") return(nf)

	tkdestroy(tf)
}

################################
# Function to choose the dataframe : builds a listbox containing the dataframes
# that are in the global environment and allows the user to choose one
################################
"choosedf" <- function(df.entry, dfnr.label, dfnc.label)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose dataframe")
	done <- tclVar(0)
	
	vnr <- NULL
	vnc <- NULL
	numi <- 1

	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")

	obj <- ls(globalenv())
#
# For all objects in the global environment, check to see if it is a dataframe
# or a list. If it is a data frame, insert it in the listbox, and if it is a list,
# check its elements.
#
	flb <- function(x1) {
		xobj <- get(x1, envir=globalenv())
		if (is.data.frame(xobj)) {
			tkinsert(tlb, "end", x1)
			cbind(nrow(xobj),ncol(xobj))
		} else if (is.list(xobj)) {
			if (length(names(xobj)) != 0) {
				fn1 <- function(x) {
					sobjn <- paste(x1,"$",x,sep="")
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					if (is.data.frame(sobj)) {
						tkinsert(tlb, "end", sobjn)
					}
				}
				sapply(names(xobj), fn1)
				fn2 <- function(x) {
					sobjn <- paste(x1,"$",x,sep="")
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					if (is.data.frame(sobj)) {
						cbind(nrow(sobj), ncol(sobj))
					}
				}
				res <- sapply(names(xobj), fn2)
				return(res)		
			}
		}
	}
	v <- unlist(lapply(obj, flb))
	if (length(v) > 0) {
		vnr <- v[seq(from=1,to=length(v),by=2)]
		vnc <- v[seq(from=2,to=length(v),by=2)]
	}
	
	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
#
# Get the number of the element choosed by the user
#
	numc <- tclvalue(tkcurselection(tlb))
	numi <- as.integer(numc)+1
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
#
# Put the name of the object in the dataframe text entry
#
	tkdelete(df.entry, 0, "end")
	tkinsert(df.entry, "end", choix)
#
# Put the row and column numbers of the dataframe in the corresponding labels
#
	tkconfigure(dfnr.label, text=as.character(vnr[numi]))
	tkconfigure(dfnc.label, text=as.character(vnc[numi]))

	tkdestroy(tf)
}

################################
# Function to choose the spatial object (contour or area) : builds a listbox containing the contours
# that are in the global environment and allows the user to choose one
################################
"choosesp" <- function(sp.entry)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose sp")
	done <- tclVar(0)
	
	vnr <- NULL
	vnc <- NULL
	numi <- 1

	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...) tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function() tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function() tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")

	obj <- ls(globalenv())
#
# For all objects in the global environment, check to see if it is a dataframe
# or a list. If it is a data frame, insert it in the listbox, and if it is a list,
# check its elements.
#
	flb <- function(x1) {
		xobj <- get(x1, envir=globalenv())
		if (is.data.frame(xobj) && ((ncol(xobj) == 3) | (ncol(xobj) == 4))) {
			tkinsert(tlb, "end", x1)
			cbind(nrow(xobj),ncol(xobj))
		} else if (is.list(xobj)) {
			if (length(names(xobj)) != 0) {
				fn1 <- function(x) {
					sobjn <- paste(x1,"$",x,sep="")
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					if (is.data.frame(sobj) && ((ncol(sobj) == 3) | (ncol(sobj) == 4))) {
						tkinsert(tlb, "end", sobjn)
					}
				}
				sapply(names(xobj), fn1)
				fn2 <- function(x) {
					sobjn <- paste(x1,"$",x,sep="")
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					if (is.data.frame(sobj) && ((ncol(sobj) == 3) | (ncol(sobj) == 4))) {
						cbind(nrow(sobj), ncol(sobj))
					}
				}
				res <- sapply(names(xobj), fn2)
				return(res)		
			}
		}
	}
	v <- unlist(lapply(obj, flb))
  
	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
#
# Get the number of the element choosed by the user
#
	numc <- tclvalue(tkcurselection(tlb))
	numi <- as.integer(numc)+1
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
#
# Put the name of the object in the dataframe text entry
#
	tkdelete(sp.entry, 0, "end")
	tkinsert(sp.entry, "end", choix)

	tkdestroy(tf)
}

################################
# Function to choose the pixmap : builds a listbox containing the pixmaps
# that are in the global environment and allows the user to choose one
################################
"choosepm" <- function(pm.entry)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose pixmap")
	done <- tclVar(0)
	
	vnr <- NULL
	vnc <- NULL
	numi <- 1

	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")

	obj <- ls(globalenv())
#
# For all objects in the global environment, check to see if it is a dataframe
# or a list. If it is a data frame, insert it in the listbox, and if it is a list,
# check its elements.
#
	flb <- function(x1) {
		xobj <- get(x1, envir=globalenv())
		if (class(xobj) == "pixmapIndexed") {
			tkinsert(tlb, "end", x1)
			cbind(nrow(xobj),ncol(xobj))
		} else if (is.list(xobj)) {
			if (length(names(xobj)) != 0) {
				fn1 <- function(x) {
					sobjn <- paste(x1,"$",x,sep="")
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					if (class(xobj) == "pixmapIndexed") {
						tkinsert(tlb, "end", sobjn)
					}
				}
				sapply(names(xobj), fn1)
				fn2 <- function(x) {
					sobjn <- paste(x1,"$",x,sep="")
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					if (class(xobj) == "pixmapIndexed") {
						cbind(nrow(sobj), ncol(sobj))
					}
				}
				res <- sapply(names(xobj), fn2)
				return(res)		
			}
		}
	}
	v <- unlist(lapply(obj, flb))

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
#
# Get the number of the element choosed by the user
#
	numc <- tclvalue(tkcurselection(tlb))
	numi <- as.integer(numc)+1
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
#
# Put the name of the object in the dataframe text entry
#
	tkdelete(pm.entry, 0, "end")
	tkinsert(pm.entry, "end", choix)
#
# Put the row and column numbers of the dataframe in the corresponding labels
#
	#tkconfigure(dfnr.label, text=as.character(vnr[numi]))
	#tkconfigure(dfnc.label, text=as.character(vnc[numi]))

	tkdestroy(tf)
}


################################
# Function to choose a distance matrix : builds a listbox containing the dist mat
# that are in the global environment and allows the user to choose one
################################
"choosedist" <- function(df.entry, dfnr.label, dfnc.label)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose distance matrix")
	done <- tclVar(0)
	
	vnr <- NULL
	vnc <- NULL
	numi <- 1

	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")
	
	obj <- ls(globalenv())
#
# For all objects in the global environment, check to see if it is a dist mat
# or a list. If it is a dist mat, insert it in the listbox, and if it is a list,
# check its elements.
#
	flb <- function(x1) {
		xobj <- get(x1, envir=globalenv())
		if (inherits(xobj, "dist")) {
			tkinsert(tlb, "end", x1)
			cbind(attributes(xobj)$Size,attributes(xobj)$Size)
		} else if (is.list(xobj)) {
			if (length(names(xobj)) != 0) {
				fn1 <- function(x) {
					sobjn <- paste(x1,"$",x,sep="")
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					if (inherits(sobj, "dist")) {
						tkinsert(tlb, "end", sobjn)
					}
				}
				sapply(names(xobj), fn1)
				fn2 <- function(x) {
					sobjn <- paste(x1,"$",x,sep="")
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					if (inherits(sobj, "dist")) {
						cbind(attributes(sobj)$Size,attributes(sobj)$Size)
					}
				}
				res <- sapply(names(xobj), fn2)
				return(res)		
			}
		}
	}

	v <- unlist(lapply(obj, flb))
	vnr <- v[seq(from=1,to=length(v),by=2)]
	vnc <- v[seq(from=2,to=length(v),by=2)]

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	
	if(tclvalue(done)=="2") return(0)
#
# Get the number of the element choosed by the user
#
	numc <- tclvalue(tkcurselection(tlb))
	numi <- as.integer(numc)+1
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
#
# Put the name of the object in the dist mat text entry
#
	tkdelete(df.entry, 0, "end")
	tkinsert(df.entry, "end", choix)
#
# Put the row and column numbers of the dist mat in the corresponding labels
#
	tkconfigure(dfnr.label, text=as.character(vnr[numi]))
	tkconfigure(dfnc.label, text=as.character(vnc[numi]))

	tkdestroy(tf)
}

################################
# Function to choose a neighbouring : builds a listbox containing the neighbouring
# that are in the global environment and allows the user to choose one
################################
"chooseneig" <- function(neig.entry)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose neighbouring")
	done <- tclVar(0)
	
	vnr <- NULL
	vnc <- NULL
	numi <- 1

	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")
	
	obj <- ls(globalenv())
#
# For all objects in the global environment, check to see if it is a neighbouring
# or a list. If it is a dist mat, insert it in the listbox, and if it is a list,
# check its elements.
#
	flb <- function(x1) {
		xobj <- get(x1, envir=globalenv())
		if (inherits(xobj, "neig")) {
			tkinsert(tlb, "end", x1)
			cbind(attributes(xobj)$Size,attributes(xobj)$Size)
		} else if (is.list(xobj)) {
			if (length(names(xobj)) != 0) {
				fn1 <- function(x) {
					sobjn <- paste(x1,"$",x,sep="")
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					if (inherits(sobj, "neig")) {
						tkinsert(tlb, "end", sobjn)
					}
				}
				sapply(names(xobj), fn1)
				fn2 <- function(x) {
					sobjn <- paste(x1,"$",x,sep="")
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					if (inherits(sobj, "neig")) {
						cbind(attributes(sobj)$Size,attributes(sobj)$Size)
					}
				}
				res <- sapply(names(xobj), fn2)
				return(res)		
			}
		}
	}

	v <- unlist(lapply(obj, flb))

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	
	if(tclvalue(done)=="2") return(0)
#
# Get the number of the element choosed by the user
#
	numc <- tclvalue(tkcurselection(tlb))
	numi <- as.integer(numc)+1
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
#
# Put the name of the object in the neig text entry
#
	tkdelete(neig.entry, 0, "end")
	tkinsert(neig.entry, "end", choix)

	tkdestroy(tf)
}

################################
# Function to choose a vector of numerics : builds a listbox containing the numeric vectors
# that are in the global environment and allows the user to choose one
################################
"chooseval" <- function(tt, dfnr.label, rw.entry, urwvar)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose values")
	
	done <- tclVar(0)
	
	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")
	
	nbr <- tclvalue(tkcget(dfnr.label, "-text"))
	obj <- ls(globalenv())
#
# For all objects in the global environment, check to see if it is a vector
# or a list. If it is a dist mat, insert it in the listbox, and if it is a list,
# check its elements.
#
	flb <- function(x1) {
		xobj <- get(x1, envir=globalenv())
		if (is.vector(xobj) && identical(class(xobj), "numeric")) {
			if (length(xobj) == nbr) tkinsert(tlb, "end", x1)
		} else if (is.list(xobj)) {
			if (length(names(xobj)) != 0) {
				fn1 <- function(x) {
					sobjn <- paste(x1,"$",x,sep="")
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					if (is.vector(sobj) && identical(class(sobj), "numeric")) {
						if (length(sobj) == nbr) tkinsert(tlb, "end", sobjn)
					} else if (is.data.frame(sobj)) {
						if (nrow(sobj) == nbr)
							tkinsert(tlb, "end", sobjn)
					}
				}
				sapply(names(xobj), fn1)
			}
		}
	}
	v <- unlist(lapply(obj, flb))

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
	
	numc <- tclvalue(tkcurselection(tlb))
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
			
	tkdelete(rw.entry, 0, "end")
	tkinsert(rw.entry, "end", choix)
	
	tkdestroy(tf)
}

################################
# Function to choose a vector of labels : builds a listbox containing the label vectors
# that are in the global environment and allows the user to choose one
################################
"chooselab" <- function(tt, dfnr.label, rw.entry, urwvar)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose labels")
	
	done <- tclVar(0)
	
	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")
	
	nbr <- tclvalue(tkcget(dfnr.label, "-text"))
	obj <- ls(globalenv())
	
#
# For all objects in the global environment, check to see if it is a vector
# or a list. If it is a dist mat, insert it in the listbox, and if it is a list,
# check its elements.
#
	flb <- function(x1) {
		xobj <- get(x1, envir=globalenv())
		if (is.vector(xobj) && identical(class(xobj), "character")) {
			if (length(xobj) == nbr) tkinsert(tlb, "end", x1)
		} else if (is.list(xobj)) {
			if (length(names(xobj)) != 0) {
				fn1 <- function(x) {
					sobjn <- paste(x1,"$",x,sep="")
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					if (is.vector(sobj) && identical(class(sobj), "character")) {
						if (length(sobj) == nbr) tkinsert(tlb, "end", sobjn)
					}
				}
				sapply(names(xobj), fn1)
			}
		}
	}
	v <- unlist(lapply(obj, flb))

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
	
	numc <- tclvalue(tkcurselection(tlb))
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
			
	tkdelete(rw.entry, 0, "end")
	tkinsert(rw.entry, "end", choix)
	
	tkdestroy(tf)
}

################################
# Function to choose row weights
################################
"chooserw" <- function(tt, dfnr.label, rw.entry, urwvar)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose row weights")
	
	done <- tclVar(0)
	
	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")
	
	nbr <- tclvalue(tkcget(dfnr.label, "-text"))
	obj <- ls(globalenv())
	
#
# For all objects in the global environment, check to see if it is a vector
# or a list. If it is a dist mat, insert it in the listbox, and if it is a list,
# check its elements.
#
	flb <- function(x1) {
		xobj <- get(x1, envir=globalenv())
		if (is.vector(xobj) && identical(class(xobj), "numeric")) {
			if (length(xobj) == nbr) tkinsert(tlb, "end", x1)
		} else if (is.list(xobj)) {
			if (length(names(xobj)) != 0) {
				fn1 <- function(x) {
					sobjn <- paste(x1,"$",x,sep="")
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					if (is.vector(sobj) && identical(class(sobj), "numeric")) {
						if (length(sobj) == nbr) tkinsert(tlb, "end", sobjn)
					}
				}
				sapply(names(xobj), fn1)
			}
		}
	}
	v <- unlist(lapply(obj, flb))

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
	
	numc <- tclvalue(tkcurselection(tlb))
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
			
	tkdelete(rw.entry, 0, "end")
	tkinsert(rw.entry, "end", choix)
	
	tclvalue(urwvar)<-"0"
	urwvar <- tclVar(0)
	rwl.cbut <- tkcheckbutton(tt,text="Uniform row weights", variable=urwvar)

	tkdestroy(tf)
}

################################
# Function to choose column weights
################################
"choosecw" <- function(tt, dfnc.label, cw.entry, ucwvar)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose column weights")
	
	done <- tclVar(0)
	
	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")
	
	nbc <- tclvalue(tkcget(dfnc.label, "-text"))
	obj <- ls(globalenv())
	
#
# For all objects in the global environment, check to see if it is a vector
# or a list. If it is a dist mat, insert it in the listbox, and if it is a list,
# check its elements.
#
	flb <- function(x1) {
		xobj <- get(x1, envir=globalenv())
		if (is.vector(xobj) && identical(class(xobj), "numeric")) {
			if (length(xobj) == nbc) tkinsert(tlb, "end", x1)
		} else if (is.list(xobj)) {
			if (length(names(xobj)) != 0) {
				fn1 <- function(x) {
					sobjn <- paste(x1,"$",x,sep="")
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					if (is.vector(sobj) && identical(class(sobj), "numeric")) {
						if (length(sobj) == nbc) tkinsert(tlb, "end", sobjn)
					}
				}
				sapply(names(xobj), fn1)
			}
		}
	}
	v <- unlist(lapply(obj, flb))

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
	
	numc <- tclvalue(tkcurselection(tlb))
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
			
	tkdelete(cw.entry, 0, "end")
	tkinsert(cw.entry, "end", choix)
	
	tclvalue(ucwvar)<-"0"
	tkcheckbutton(tt,text="Uniform row weights", variable=ucwvar)

	tkdestroy(tf)
}

################################
# Function to choose a data set in the ade4 package
################################
"choosepackage" <- function(show, history)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose ade4 data set")
	
	done <- tclVar(0)
	
	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")
	
	lade4 <- data(package="ade4", envir=env_ade4tkgui)
	for (i in seq(1, length(lade4$result)/4)) {
		dsname <- lade4$result[i,3]
		tkinsert(tlb, "end", dsname)
	}

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tlb, "<Button-3>", function() if (tclvalue(tkcurselection(tlb)) != "") print(help(tclvalue(tkget(tlb, tclvalue(tkcurselection(tlb)))))))
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
	
	numc <- tclvalue(tkcurselection(tlb))
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
	data(list=choix, envir=env_ade4tkgui)

	if (show) {
		pr1 <- substr(options("prompt")$prompt, 1,2)
		cat("data(", choix, ")\n", pr1, sep="")
	}

	if (history)
    rewriteHistory(paste("data(", choix, ")", sep=""))
	
	tkdestroy(tf)
	q <- tkmessageBox(icon="info", title="Data set loaded", type="ok", message=paste("The \"",choix,"\" data set has been successfully loaded.", sep=""))
}

################################
# Function to choose a text file
################################
"readtable" <- function(show, history)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose a text file")
	
	done <- tclVar(0)
	typesvar <- tclVar()
	tabnamevar <- tclVar("")
	filenamevar <- tclVar("clipboard")
	fictrt <- tclVar()
	varnames <- tclVar(1)
	sepVar <- tclVar(1)
	otherSepVar <- tclVar("")
	decSepVar <- tclVar(1)
	colClassVar <- tclVar(4)

	"choosefile" <- function()
	{
		tclvalue(typesvar)="{{Text files} {.txt}}"
		fictrt <- tkgetOpenFile(filetypes=tclvalue(typesvar))
		fpath <- tclvalue(fictrt)
		tkdelete(file.entry, 0, "end")
		tkinsert(file.entry, "end", fpath)
	}
	
	"choosefic" <- function(show, history)
	{	
		if (tclvalue(tabnamevar) != "") {
			tabname  <- parse(text=tclvalue(tabnamevar))[[1]]
		} else tabname <- "untitled"

		if (tclvalue(filenamevar) != "") {
			filename  <- tclvalue(filenamevar)
		} else return()
		
		varn <- as.logical(tclObj(varnames))
		sep <- tclvalue(sepVar)
		if (sep == 1) sepch <- ""
		if (sep == 2) sepch <- ","
		if (sep == 3) sepch <- ";"
		if (sep == 4) {
			if (tclvalue(otherSepVar) != "") {
				otherSep <- tclvalue(otherSepVar)
			} else otherSep <- ""
			sepch <- otherSep
		}
		decSep <- tclvalue(decSepVar)
		if (decSep == 1) decsepch <- "."
		if (decSep == 2) decsepch <- ","

		colClass <- tclvalue(colClassVar)
		if (colClass == 1) colClass <- "numeric"
		if (colClass == 2) colClass <- "character"
		if (colClass == 3) colClass <- "factor"
		if (colClass == 4) colClass <- "auto"

		if (colClass == "auto") {
			rdcom <- paste(tabname," <<- read.table(file='", filename, "', header=",varn,
				", sep='",sepch,"', dec='",decsepch,"')", sep="")
		} else {
			rdcom <- paste(tabname," <<- read.table(file='", filename, "', header=",varn,
				", sep='",sepch,"', dec='",decsepch,"', colClasses='",colClass,"')", sep="")
		}
	
		if (show) {
			pr1 <- substr(options("prompt")$prompt, 1,2)
			cat(rdcom, "\n", pr1, sep="")
		}

		if (history) rewriteHistory(rdcom)

		eval(parse(text=rdcom))
		tkdestroy(tf)
		eval(parse(text=paste("edit(", tabname, ")", sep="")))
	}
	
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	labh <- tklabel(frame1, bitmap="questhead")
	tkbind(labh, "<Button-1>", function() print(help("read.table")))
	tkgrid(tklabel(frame1,text="Read text file", font="Times 18", foreground="red"), labh, columnspan=2)
	tkpack(frame1, fill = "x")

	frame2 <- tkframe(tf, relief="groove", borderwidth=2)	
	tab.entry <- tkentry(frame2, textvariable=tabnamevar)
	file.entry <- tkentry(frame2, textvariable=filenamevar)
	choosefile.but <- tkbutton(frame2, text="Set", command=function() choosefile())
	tkgrid(tklabel(frame2,text="Text file to read : "), file.entry, choosefile.but)
	tkgrid(tklabel(frame2,text="Dataframe to receive the data : "), tab.entry)
	varnames.cbut <- tkcheckbutton(frame2,text="Variables names on the first row of data file", variable=varnames)
	tkgrid(varnames.cbut, columnspan=2, sticky="w")
	
	sepFrame <- tkframe(frame2, relief="groove", borderwidth=2)
	sep.entry <- tkentry(sepFrame, textvariable=otherSepVar, width=10)
  tkgrid(tklabel(sepFrame, text="Field separator:", foreground="blue"))
  tkgrid(tkradiobutton(sepFrame, text="Default", value=1, variable=sepVar), sticky="w")
  tkgrid(tkradiobutton(sepFrame, text="Commas", value=2, variable=sepVar), sticky="w")
  tkgrid(tkradiobutton(sepFrame, text="Semicolon", value=3, variable=sepVar), sticky="w")
  tkgrid(tkradiobutton(sepFrame, text="Other", value=4, variable=sepVar), sep.entry, sticky="w")

  decSepFrame <- tkframe(frame2, relief="groove", borderwidth=2)
  tkgrid(tklabel(decSepFrame, text="Decimal separator:", foreground="blue"))
  tkgrid(tkradiobutton(decSepFrame, text="Period [.]", value=1, variable=decSepVar), sticky="w")
  tkgrid(tkradiobutton(decSepFrame, text="Comma [,]", value=2, variable=decSepVar), sticky="w")

  colClassFrame <- tkframe(frame2, relief="groove", borderwidth=2)
  tkgrid(tklabel(colClassFrame, text="Column types:", foreground="blue"))
  tkgrid(tkradiobutton(colClassFrame, text="Default", value=4, variable=colClassVar), sticky="w")
  tkgrid(tkradiobutton(colClassFrame, text="Numeric", value=1, variable=colClassVar), sticky="w")
  tkgrid(tkradiobutton(colClassFrame, text="Character", value=2, variable=colClassVar), sticky="w")
  tkgrid(tkradiobutton(colClassFrame, text="Factor", value=3, variable=colClassVar), sticky="w")
 	tkgrid(sepFrame, decSepFrame, colClassFrame, sticky="n")

	tkpack(frame2, fill = "x")

	ok.but <- tkbutton(tf, text="Submit", command=function() choosefic(show, history))
	cancel.but <- tkbutton(tf, text="Dismiss", command=function() tkdestroy(tf))
	tkpack(cancel.but, ok.but, side="left", fill="x", expand=1)

	tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
	tkbind(tf, "<KeyPress-Return>", function() choosefic(show, history))
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	tkwait.variable(done)
	if(tclvalue(done) == "2")
    return(0)
	tkdestroy(tf)
}

################################
# Function to display a text file
################################
"displaytable" <- function(tabname)
{
  tt <- tktoplevel()
  tkwm.title(tt, "Dataframe display")
  txt <- tktext(tt, bg="white", font="courier")
  scr <- tkscrollbar(tt, repeatinterval=5, command=function(...)tkyview(txt,...))
  
  ## Safest to make sure scr exists before setting yscrollcommand
  tkconfigure(txt, yscrollcommand=function(...)tkset(scr,...))
  tkpack(txt, side="left", fill="both", expand=TRUE)
  tkpack(scr, side="right", fill="y")

	sink(conn <- file("Rlisting001.tmp", open="w"))
	cheval <- paste(tabname,"[1:min(20,nrow(", tabname, ")), 1:min(5,ncol(", tabname, "))]", sep="")
	cat("First rows and columns of dataframe:\n")
	print(eval(parse(text=cheval)))
	sink()
	close(conn)
  chn <- tclopen(file.path("Rlisting001.tmp"))
  tkinsert(txt, "end", tclread(chn))
  tclclose(chn)
  system("rm Rlisting001.tmp")

  tkconfigure(txt, state="disabled")
  tkmark.set(txt,"insert","0.0")
  tkfocus(txt)
}

################################
# dialog box to display a dudi
################################
"dudisp" <- function(show, history)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Dudi display")

	"calldisp" <- function(show, history)
	{
		dialog.dudi.display(show, history, tclvalue(dudivar))

		if (show) {
			pr1 <- substr(options("prompt")$prompt, 1,2)
			cat(paste("dialog.dudi.display(", show, ",", history, ",eval(expression(", dQuote(tclvalue(dudivar)), ")))", sep=""), "\n", pr1, sep="")
		}

		if (history)
      rewriteHistory(paste("dialog.dudi.display(", show, ",", history, ",eval(expression(", dQuote(tclvalue(dudivar)), ")))", sep=""))

		tkdestroy(tf)
	}
	
	done <- tclVar(0)
	dudivar <- tclVar()

	frame1 <- tkframe(tf, relief="groove", borderwidth=2)	
	labh <- tklabel(frame1, bitmap="questhead")
	tkbind(labh, "<Button-1>", function() print(help("dudi")))
	tkgrid(tklabel(frame1,text="Dudi display", font="Times 18", foreground="red"), labh, columnspan=2)
	tkpack(frame1, fill = "x")

	frame2 <- tkframe(tf, relief="groove", borderwidth=2)	
	dudi.entry <- tkentry(frame2, textvariable=dudivar)
	choosedudi.but <- tkbutton(frame2, text="Set", command=function() choosedudi(dudi.entry))
	tkgrid(tklabel(frame2,text="dudi name : "), dudi.entry, choosedudi.but)
	tkpack(frame2, fill = "x")

	cancel.but <- tkbutton(tf, text="Dismiss", command=function() tkdestroy(tf))
	dispdudi.but <- tkbutton(tf, text="Submit", command=function() calldisp(show, history))
	tkpack(cancel.but, dispdudi.but, expand=1, fill="x", side="left")

	tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	tkbind(tf, "<KeyPress-Return>", function() calldisp(show, history))

	tkwait.variable(done)
	if(tclvalue(done) == "2") return(0)
	tkdestroy(tf)
}

################################
# Function to choose the dudi
################################
"choosedudi" <- function(dudi.entry)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose dudi")
	done <- tclVar(0)
	
	numi <- 1

	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")

	obj <- ls(globalenv())
#
# For all objects in the global environment, check to see if it is a dudi
# or a list. If it is a dudi, insert it in the listbox, and if it is a list,
# check its elements.
#
	flb <- function(x1) {
		xobj <- get(x1, envir=globalenv())
		if (is.dudi(xobj) || class(xobj)=="discrimin" || class(xobj)=="dpcoa") {
			tkinsert(tlb, "end", x1)
		}
	}
	v <- unlist(lapply(obj, flb))

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
#
# Get the number of the element choosed by the user
#
	numc <- tclvalue(tkcurselection(tlb))
	numi <- as.integer(numc)+1
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
#
# Put the name of the object in the dataframe text entry
#
	tkdelete(dudi.entry, 0, "end")
	tkinsert(dudi.entry, "end", choix)

	tkdestroy(tf)
}

################################
# Function to choose the dudi
################################
"choosedudirc" <- function(dudi.entry, dfnr.label, dfnc.label)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose dudi")
	done <- tclVar(0)
	
	vnr <- NULL
	vnc <- NULL
	numi <- 1

	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")

	obj <- ls(globalenv())
#
# For all objects in the global environment, check to see if it is a dudi
# or a list. If it is a dudi, insert it in the listbox, and if it is a list,
# check its elements.
#
	flb <- function(x1) {
		xobj <- get(x1, envir=globalenv())
		if (is.dudi(xobj)) {
			tkinsert(tlb, "end", x1)
			nr <- eval(parse(text=paste("nrow(",x1,"$tab",")",sep="")))
			nc <- eval(parse(text=paste("ncol(",x1,"$tab",")",sep="")))
			return(cbind(nr, nc))
		}
	}
	v <- unlist(lapply(obj, flb))
	if (length(v) > 0) {
		vnr <- v[seq(from=1,to=length(v),by=2)]
		vnc <- v[seq(from=2,to=length(v),by=2)]
	}

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
#
# Get the number of the element choosed by the user
#
	numc <- tclvalue(tkcurselection(tlb))
	numi <- as.integer(numc)+1
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
#
# Put the name of the object in the dataframe text entry
#
	tkdelete(dudi.entry, 0, "end")
	tkinsert(dudi.entry, "end", choix)
#
# Put the row and column numbers of the dataframe in the corresponding labels
#
	tkconfigure(dfnr.label, text=as.character(vnr[numi]))
	tkconfigure(dfnc.label, text=as.character(vnc[numi]))

	tkdestroy(tf)
}

################################
# Function to choose the dudi
################################
"chooseduditest" <- function(dudi.entry)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose dudi")
	done <- tclVar(0)
	
	numi <- 1

	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")

	obj <- ls(globalenv())
#
# For all objects in the global environment, check to see if it is a dudi
# or a list. If it is a dudi, insert it in the listbox, and if it is a list,
# check its elements.
#
	flb <- function(x1) {
		xobj <- get(x1, envir=globalenv())
		if (class(xobj)=="between" || class(xobj)=="discrimin" || class(xobj)=="coinertia" || class(xobj)=="cca" || class(xobj)=="pcaiv" || class(xobj)=="pcaivortho") {
			tkinsert(tlb, "end", x1)
		}
	}
	v <- unlist(lapply(obj, flb))

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
#
# Get the number of the element choosed by the user
#
	numc <- tclvalue(tkcurselection(tlb))
	numi <- as.integer(numc)+1
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
#
# Put the name of the object in the dataframe text entry
#
	tkdelete(dudi.entry, 0, "end")
	tkinsert(dudi.entry, "end", choix)

	tkdestroy(tf)
}

################################
# dialog box to launch the explore function
################################
"exploregraph" <- function(show, history) {
	tf <- tktoplevel()
	tkwm.title(tf, "Graph exploration")
	loclist <- get("cmdlist", envir = env_ade4tkgui)
	
	"callexp" <- function()
	{
		appel <- tclvalue(callvar)[[1]]
		plotcmd <- parse(text = loclist[as.numeric(appel) + 1])		
		explorecmd <- parse(text = paste("explore(", plotcmd, ")", sep=""))		
		if (show) {
			pr1 <- substr(options("prompt")$prompt, 1,2)
			cat(paste("explore(", plotcmd, ")", sep=""), "\n", pr1, sep="")
		}
		if (history) rewriteHistory(paste("explore(", plotcmd, ")", sep=""))
		eval.parent(explorecmd)
		tkdestroy(tf)
	}
	
	done <- tclVar(0)
	callvar <- tclVar()

	frame1 <- tkframe(tf, relief="groove", borderwidth=2)	
	tkgrid(tklabel(frame1,text="Graph exploration", font="Times 18", foreground="red"), columnspan=2)
	tkpack(frame1, fill = "x")

	frame2 <- tkframe(tf, relief="groove", borderwidth=2)	
	call.entry <- tkentry(frame2, textvariable=callvar)
	choosegraph.but <- tkbutton(frame2, text="Set", command=function() choosegraph(call.entry))
	tkgrid(tklabel(frame2,text="Graph function : "), call.entry, choosegraph.but)
	tkpack(frame2, fill = "x")

	explore.but <- tkbutton(tf, text="Submit", command=function() callexp())
	cancel.but <- tkbutton(tf, text="Dismiss", command=function() tkdestroy(tf))
	tkpack(cancel.but, explore.but, expand=1, fill="x", side="left")

	tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
	tkbind(tf, "<KeyPress-Return>", function() callexp())
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))

	tkwait.variable(done)
	if(tclvalue(done) == "2") return(0)
	tkdestroy(tf)
}

################################
# Function to choose the graph
################################
"choosegraph" <- function(graph.entry)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose graph")
	done <- tclVar(0)
	
	numi <- 1

	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")

#
# cmdlist contains the list of graphs that were drawn by the user.
#
	cmdlist1 <- get("cmdlist", envir=env_ade4tkgui)
	if (length(cmdlist1) > 1) {
		for (i in 2:length(cmdlist1)) {
			if (is.call(cmdlist1[[i]])) {
				dcall <- cmdlist1[[i]]
				narg <- length(names(dcall))
				paramlst <- encodeString(dcall)[2:narg]
				arglst <- names(dcall)[2:narg]
				call1 <- paste(encodeString(dcall)[1],"(", paste(arglst, paramlst, sep=" = ",collapse=", "), ")", sep="")
				tkinsert(tlb, "end", call1)
			} else {
				tkinsert(tlb, "end", as.character(cmdlist1[[i]]))
			}
		}
	}	
	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
#
# Get the number of the element choosed by the user
#
	numc <- tclvalue(tkcurselection(tlb))
	numi <- as.integer(numc)+1
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
#
# Put the name of the object in the dataframe text entry
#
	tkdelete(graph.entry, 0, "end")
	tkinsert(graph.entry, "end", numi)
	
	tkdestroy(tf)
}

################################
# Function to choose a factor (s.class version)
################################
"choosefac" <- function(fac.entry, dfnr.label)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose factor")
	
	done <- tclVar(0)
	
	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")
	
	nbr <- tclvalue(tkcget(dfnr.label, "-text"))
	obj <- ls(globalenv())	
	
#
# For all objects in the global environment, check to see if it is a vector
# or a list. If it is a dist mat, insert it in the listbox, and if it is a list,
# check its elements.
#

	flb <- function(x1) {
		xobj <- get(x1, envir=globalenv())
		if (is.factor(xobj)) {
			if (length(xobj) == nbr) tkinsert(tlb, "end", x1)
		} else if (is.list(xobj)) {
			if (length(names(xobj)) != 0) {
				fn1 <- function(x) {
					sobjn <- paste(x1,"$", x, sep="")
					options(show.error.messages = FALSE)
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					options(show.error.messages = TRUE)
					if (is.factor(sobj)) {
						if (length(sobj) == nbr) tkinsert(tlb, "end", sobjn)
					} else if (is.list(sobj)) {
						if (length(names(sobj)) != 0) {
							fn2 <- function(x) {
								ssobjn <- paste(sobjn,"$", x, sep="")
								options(show.error.messages = FALSE)
								ssobj <- try(eval(parse(text=ssobjn)), silent=TRUE)
								options(show.error.messages = TRUE)
								if (is.factor(ssobj)) {
									if (length(ssobj) == nbr) tkinsert(tlb, "end", ssobjn)
								}
							}
							sapply(names(sobj), fn2)
						}
					}
				}
				sapply(names(xobj), fn1)
			}
		}
	}
	v <- unlist(lapply(obj, flb))

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
	
	numc <- tclvalue(tkcurselection(tlb))
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
			
	tkdelete(fac.entry, 0, "end")
	tkinsert(fac.entry, "end", choix)
	
	tkdestroy(tf)
}

################################
# Function to choose a factor (between version)
################################
"choosefac2" <- function(fac.entry, dudi.entry)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose factor")
	
	done <- tclVar(0)
	
	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")
	
	obj <- ls(globalenv())
	dudiname <- ""
	dudi <- tclvalue(tkcget(dudi.entry, "-text"))
	if (tclvalue(dudi) != "") dudiname  <- parse(text=tclvalue(dudi))[[1]]
	if (dudiname != "") {
		nbr <- eval(parse(text=paste("nrow(",dudiname,"$tab)",sep="")))
	} else nbr <- 0
	
	flb <- function(x1) {
		xobj <- get(x1, envir=globalenv())
		if (is.factor(xobj)) {
			if (length(xobj) == nbr) tkinsert(tlb, "end", x1)
		} else if (is.list(xobj)) {
			if (length(names(xobj)) != 0) {
				fn1 <- function(x) {
					sobjn <- paste(x1,"$", x, sep="")
					options(show.error.messages = FALSE)
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					options(show.error.messages = TRUE)
					if (is.factor(sobj)) {
						if (length(sobj) == nbr) tkinsert(tlb, "end", sobjn)
					} else if (is.list(sobj)) {
						if (length(names(sobj)) != 0) {
							fn2 <- function(x) {
								ssobjn <- paste(sobjn,"$", x, sep="")
								options(show.error.messages = FALSE)
								ssobj <- try(eval(parse(text=ssobjn)), silent=TRUE)
								options(show.error.messages = TRUE)
								if (is.factor(ssobj)) {
									if (length(ssobj) == nbr) tkinsert(tlb, "end", ssobjn)
								}
							}
							sapply(names(sobj), fn2)
						}
					}
				}
				sapply(names(xobj), fn1)
			}
		}
	}
	v <- unlist(lapply(obj, flb))

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
	
	numc <- tclvalue(tkcurselection(tlb))
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
			
	tkdelete(fac.entry, 0, "end")
	tkinsert(fac.entry, "end", choix)
	
	tkdestroy(tf)
}

################################
# Function to choose a vector of color names
################################
"choosecol" <- function(col.entry, facvar)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose a vector of color names")
	
	done <- tclVar(0)
	
	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")
	
	obj <- ls(globalenv())
	nlev <- length(levels(eval(parse(text=tclvalue(facvar)))))

	for (i in 1:length(obj)) {
		nomobj <- obj[[i]]
		xobj <- get(nomobj, envir=globalenv())
		if (is.vector(xobj) && is.character(xobj)) {
			if (length(xobj) == nlev) tkinsert(tlb, "end", nomobj)
		}
	}

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
	
	numc <- tclvalue(tkcurselection(tlb))
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
			
	tkdelete(col.entry, 0, "end")
	tkinsert(col.entry, "end", choix)
	
	tkdestroy(tf)
}

################################
# Function to choose a vector of weights
################################
"choosewt" <- function(wt.entry, dfnr.label)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose a vector of weights")
	
	done <- tclVar(0)
	
	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")
	
	obj <- ls(globalenv())
	
	for (i in 1:length(obj)) {
		nomobj <- obj[[i]]
		xobj <- get(nomobj, envir=globalenv())
		if (is.vector(xobj) && is.numeric(xobj)) {
			nbr <- tclvalue(tkcget(dfnr.label, "-text"))
			if (length(xobj) == nbr) tkinsert(tlb, "end", nomobj)
		}
	}

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
	
	numc <- tclvalue(tkcurselection(tlb))
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
			
	tkdelete(wt.entry, 0, "end")
	tkinsert(wt.entry, "end", choix)
	
	tkdestroy(tf)
}

################################
# Function to choose hull levels
################################
"choosechull" <- function(chull.entry)
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose a vector hull levels")
	
	done <- tclVar(0)
	
	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")
	
	obj <- ls(globalenv())
	
	for (i in 1:length(obj)) {
		nomobj <- obj[[i]]
		xobj <- get(nomobj, envir=globalenv())
		if (is.vector(xobj) && is.numeric(xobj)) {
			tkinsert(tlb, "end", nomobj)
		}
	}

	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
	
	numc <- tclvalue(tkcurselection(tlb))
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
			
	tkdelete(chull.entry, 0, "end")
	tkinsert(chull.entry, "end", choix)
	
	tkdestroy(tf)
}

################################
# Function to save a graphic in a file
################################
"outgraph" <- function()
{
#
# Main dialog window with title and frames
#
	tf <- tktoplevel()
	tkwm.title(tf,"Save graphic")
#
# Frames
#
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)	
	frame2 <- tkframe(tf, relief="groove", borderwidth=2)	
	frame3 <- tkframe(tf, relief="groove", borderwidth=2)	
    devframe <- tkframe(frame2, relief="groove", borderwidth=2)
#
# Tcl/Tk variables
#
	done <- tclVar(0)
	formatvar <- tclVar(1)
	widthvar <- tclVar(6)
	heightvar <- tclVar(6)
#
# Save function
#
	"savefic" <- function(formatvar, widthvar, heightvar)
	{
		outform <- tclvalue(formatvar)
		width <- as.numeric(tclvalue(widthvar))
		height <- as.numeric(tclvalue(heightvar))
		odev <- dev.cur()
		if (outform == 1) { # postcript
			filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.ps", defaultextension=".ps",
				title="Save graph...", filetypes="{PostScript {.ps .eps}} {{All Files} {*.*}}"))
			if (filename != "") {
				postscript(file=filename, width=width, height=height)
			}
		} else if (outform == 2) { # pdf
			filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.pdf", defaultextension=".pdf",
				title="Save graph...", filetypes="{PDF {.pdf}} {{All Files} {*.*}}"))
			if (filename != "") {
				pdf(file=filename, width=width, height=height)
			}
		} else if (outform == 3) { # pictex
			filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.tex", defaultextension=".tex",
				title="Save graph...", filetypes="{PicTeX {.tex}} {{All Files} {*.*}}"))
			if (filename != "") {
				pictex(file=filename, width=width, height=height)
			}
		} else if (outform == 4) { # xfig
			filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.fig", defaultextension=".fig",
				title="Save graph...", filetypes="{XFig {.fig}} {{All Files} {*.*}}"))
			if (filename != "") {
				xfig(file=filename, width=width, height=height)
			}
		} else if (outform == 5) { # png
			filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.png", defaultextension=".png",
				title="Save graph...", filetypes="{PNG {.png}} {{All Files} {*.*}}"))
			if (filename != "") {
				png(filename=filename, width=width, height=height)
			}
		} else if (outform == 6) { # jpeg
			filename <- tclvalue(tkgetSaveFile(initialfile="Rplots.jpeg", defaultextension=".jpeg",
				title="Save graph...", filetypes="{JPEG {.jpeg .jpg}} {{All Files} {*.*}}"))
			if (filename != "") {
				jpeg(filename=filename, width=width, height=height)
			}
		}
		ndev <- dev.cur()
		dev.set(odev)
		dev.copy(which=ndev)
		dev.off()
		tkdestroy(tf)
	}
#
# Frames setup
#
	tkgrid(tklabel(tf,text="Save current graphic", font="Times 18"), columnspan=2)

	tkgrid(tklabel(frame2,text="Output format : "), sticky="n")
  tkgrid(tkradiobutton(frame2, text="postscript", value=1, variable=formatvar), sticky="w")
  tkgrid(tkradiobutton(frame2, text="pdf", value=2, variable=formatvar), sticky="w")
  tkgrid(tkradiobutton(frame2, text="pictex", value=3, variable=formatvar), sticky="w")
  tkgrid(tkradiobutton(frame2, text="xfig", value=4, variable=formatvar), sticky="w")
  tkgrid(tkradiobutton(frame2, text="png", value=5, variable=formatvar), sticky="w")
  tkgrid(tkradiobutton(frame2, text="jpeg", value=6, variable=formatvar), sticky="w")
	tkgrid(frame2, rowspan=2, sticky="n")
    
	tkgrid(tklabel(frame3,text="Output size : "))
	width.entry <- tkentry(frame3, textvariable=widthvar, width=10)
	height.entry <- tkentry(frame3, textvariable=heightvar, width=10)
	tkgrid(tklabel(frame3,text="Width : "), width.entry)
	tkgrid(tklabel(frame3,text="Height : "), height.entry)
	tkgrid(frame3, column=1, row=1, sticky="n")

	save.but <- tkbutton(frame1, text="Save", command=function() savefic(formatvar, widthvar, heightvar))
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function() tkdestroy(tf))
	tkgrid(save.but, cancel.but)
	tkgrid(frame1, column=1, row=2, sticky="n")
	
	tkbind(tf, "<Destroy>", function() tclvalue(done) <- 2)
	tkbind(tf, "<KeyPress-Return>", function() savefic(formatvar, widthvar, heightvar))
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	tkwait.variable(done)
	if(tclvalue(done) == "2")
    return(0)
	tkdestroy(tf)
}

################################
# Function to create a new graphic window
################################
"newGr" <- function()
{
	if(.Platform$OS.type == "unix") {
		if (capabilities("aqua")) {
			# quartz()
			dev.new()
		} else if (capabilities("X11")) {
			# x11()
			dev.new()
		}
	} else {
		# windows()
		dev.new()
	}
	assign("winlist", get("winlist", envir=env_ade4tkgui)+1, envir=env_ade4tkgui)
}

################################
# Function to reset the list of graphics
################################
"resetgraph" <- function()
{
	if (exists("cmdlist")) rm("cmdlist", envir=env_ade4tkgui)
	assign("cmdlist", "cmdlist", envir=env_ade4tkgui)
}

################################
# Function to choose the graphic Window
################################
"selectGr" <- function()
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose window")
	done <- tclVar(0)
	
	tkpack(tklabel(tf,text="Activate a graphic window", font="Times 18"))

	numi <- 1

	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")

#
# dev.list contains the list of graphic windows.
#
	devlist <- dev.list()
	if (!is.null(devlist)) {
		for (i in 1:length(devlist)) {
			tkinsert(tlb, "end", paste("Device",devlist[i]))
			if (devlist[i] == dev.cur()) {
				tkitemconfigure(tlb, i-1, foreground="red", selectforeground="red")
			}
		}
	}
		
	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
#
# Get the number of the element choosed by the user
#
	numc <- tclvalue(tkcurselection(tlb))
	numi <- as.integer(numc)+1
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
#
# sets the chosen device as active
#
	dev.set(devlist[numi])
	
	tkdestroy(tf)
}

################################
# Function to ask for save before quit
################################
"askQuit" <- function()
{
	q <- tkmessageBox(message="Do you want to save before quitting?",icon="question",type="yesnocancel",default="yes")
	if (tclvalue(q) != "cancel") q(tclvalue(q))
}

################################
# Function to edit a dataframe
################################
"editdf" <- function()
{
	tf <- tktoplevel()
	tkwm.title(tf,"Choose dataframe")
	done <- tclVar(0)
	
	tlb <- tklistbox(tf)
	scr <- tkscrollbar(tf, repeatinterval=5, command=function(...)tkyview(tlb,...))
	tkconfigure(tlb, yscrollcommand=function(...)tkset(scr,...))
	frame1 <- tkframe(tf, relief="groove", borderwidth=2)
	cancel.but <- tkbutton(frame1, text="Dismiss", command=function()tkdestroy(tf))
	submit.but <- tkbutton(frame1, text="Choose", default="active", command=function()tclvalue(done)<-1)
	tkpack(cancel.but, submit.but, side="left")
	tkpack(frame1, side="bottom")
	tkpack(tlb, side="left", fill="both", expand=TRUE)
	tkpack(scr, side="right", fill="y")

	obj <- ls(globalenv())
#
# For all objects in the global environment, check to see if it is a dataframe
# or a list. If it is a data frame, insert it in the listbox, and if it is a list,
# check its elements.
#
	flb <- function(x1) {
		xobj <- get(x1, envir=globalenv())
		if (is.data.frame(xobj)) {
			tkinsert(tlb, "end", x1)
			cbind(nrow(xobj),ncol(xobj))
		} else if (is.list(xobj)) {
			if (length(names(xobj)) != 0) {
				fn1 <- function(x) {
					sobjn <- paste(x1,"$",x,sep="")
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					if (is.data.frame(sobj)) {
						tkinsert(tlb, "end", sobjn)
					}
				}
				sapply(names(xobj), fn1)
				fn2 <- function(x) {
					sobjn <- paste(x1,"$",x,sep="")
					sobj <- try(eval(parse(text=sobjn)), silent=TRUE)
					if (is.data.frame(sobj)) {
						cbind(nrow(sobj), ncol(sobj))
					}
				}
				res <- sapply(names(xobj), fn2)
				return(res)		
			}
		}
	}
	v <- unlist(lapply(obj, flb))
	if (length(v) > 0) {
		vnr <- v[seq(from=1,to=length(v),by=2)]
		vnc <- v[seq(from=2,to=length(v),by=2)]
	}
	
	tkbind(tlb, "<Double-ButtonPress-1>", function() tclvalue(done)<-1)
	tkbind(tf, "<Destroy>", function() tclvalue(done)<-2)
	tkbind(tf, "<KeyPress-Return>", function() tclvalue(done)<-1)
	tkbind(tf, "<KeyPress-Escape>", function() tkdestroy(tf))
	
	tkwait.variable(done)
	if(tclvalue(done)=="2") return(0)
#
# Get the number of the element choosed by the user
#
	numc <- tclvalue(tkcurselection(tlb))
	numi <- as.integer(numc)+1
	
	if(numc == "") {
		tkdestroy(tf)
		return(0)
	}
	
	choix <- tclvalue(tkget(tlb, numc))
	eval(parse(text=paste(choix, " <<- edit(", choix, ")", sep="")))
	tkdestroy(tf)
}

################################
# Function to modify the history (thanks to Thomas Lumley :)
# http://finzi.psych.upenn.edu/R/Rhelp02a/archive/41067.html
# [R] can one evaluate an expression in a comment? (or insert resultsinto history?) from Thomas Lumley on 2004-11-08
################################
"rewriteHistory" <- function(command)
{
	file1 <- tempfile("Rrawhist") 
	on.exit(unlink(file1)) 
	savehistory(file1) 
	conn<-file(file1,open="a") 
	writeLines(command, con=conn) 
	close(conn) 
	loadhistory(file1) 
}



################################
# Functions to test if one or two values are empty and to return a default value
# Two useful functions for programming optimization
################################
".test1value" <- function(val, default) {
  if(val!="")
    return(parse(text=val)[[1]])
  else
    return(default) 
}

".test2values" <- function(val1, val2, default) {
  if((val1!="") & (val2!=""))
    return(c(eval(parse(text=val1)), eval(parse(text=val2))))
  else
    return(default)
}
