### tk2edit.R - A data frame editor using TkTable
### Copyright (c), Jeffrey J. Hallman
### Licensed under LGPL 3 or above
###
### Changes:
### - 2007-01-01: fisrt version (for tcltk2_1.0-0)
###
### To do:
### - Rework all this!!!

###dim.tclArray <- function (ta)
###{
###	nms <- grep(",", names(ta), value = TRUE)
###	if (length(nms) == 0) return(c(0, 0))
###	c(max(as.numeric(gsub(",.*", "", nms))),
###    	max(as.numeric(gsub(".*,", "", nms)))) + 1
###}

tk2edit <- function (x, title = "Matrix Editor", header = NULL,
maxHeight = 600, maxWidth = 800, fontsize = 9, ...)
{
	if (!is.tk())
		stop("Package Tk is required but not loaded")
	if (!inherits(tclRequire("Tktable", warn = FALSE), "tclObj"))
		stop("Tcl package 'Tktable' must be installed first")
	.Tcl(paste("option add *Table.font {courier", fontsize, "bold}"))
	old <- options(scipen = 7)
	on.exit(options(old))
  
	makeCharMat <- function (x) {
    	## Make sure it's a character matrix
    	mat <- matrix(unlist(x), nrow = nrow(as.matrix(x)))
    	dm <- dim(mat)
    
    	## Check for row and column names
    	hasRownames <- length(rn <- rownames(x)) > 0
    	hasColnames <- length(cn <- colnames(x)) > 0
    	## Fake row and column names if they aren't there
    	if (!hasRownames) rn <- paste("[", 1:nrow(x), ",]", sep = "")
    	if(!hasColnames) cn <- paste("[,", 1:ncol(x), "]", sep = "")
    
    	## Format the columns 
    	mat[] <- apply(unclass(mat), 2, format, justify = "right")
    	mat <- rbind(cn, mat)
    	mat <- cbind(c("", rn), mat)
    	mat
  	}
  
  	fillTclArrayFromCharMat <- function (ta, cm) {
    	## cm[,1] contains column names, while cm[1,] has rownames
    	## cm[1,1] is ignored
    	for (j in 2:ncol(cm)) ta[[0, j - 1]] <- as.tclObj(cm[1, j], drop = TRUE)
    	for(i in 2:nrow(cm))
      		for(j in 1:ncol(cm))
        		ta[[i - 1, j - 1]] <- as.tclObj(cm[i, j], drop = TRUE)
  	}
  
  	tA <- tclArray()
  	cmat <- makeCharMat(x)
  	fillTclArrayFromCharMat(tA, cmat)

  	tt <- tktoplevel()
  	tkwm.title(tt, title)

  	colwidths <- apply(cmat, 2, function(x) max(nchar(x)) + 1 )
  	nTableCols <- ncol(cmat)
  	if((moreWidth <- 60 - sum(colwidths)) > 0) {
    	addEach <- moreWidth %/% length(colwidths)
    	if(addEach < 5) colwidths <- colwidths + addEach + 1
    	else nTableCols <- nTableCols + ceiling(moreWidth / 10)
  	}

  	tktable <- tkwidget(tt, "table", variable = tA,
    	rows = nrow(cmat), cols = nTableCols,
        titlerows = 1, titlecols = 1, selecttitle = 1,
        anchor = "e", multiline = 0, selectmode = "extended",
        rowseparator = dQuote("\n"), colseparator = dQuote("\t"),
        background = "white", maxheight = maxHeight, maxwidth = maxWidth,
        xscrollcommand = function(...) tkset(xscr, ...),
        yscrollcommand = function(...) tkset(yscr, ...))
	xscr <-tkscrollbar(tt, orient = "horizontal",
        command = function(...) tkxview(tktable, ...))
	yscr <- tkscrollbar(tt, command = function(...) tkyview(tktable, ...))

  	## Set column widths
  	for(i in 1:ncol(cmat)) tcl(tktable, "width", i - 1, colwidths[i])
  
  	## Rebind the Backspace key, which somehow gets messed up
  	string <- "bind Table <BackSpace> {
		set ::tk::table::Priv(junk) [%W icursor]
    	if {[string compare {} $::tk::table::Priv(junk)] && $::tk::table::Priv(junk)} {
		%W delete active [expr {$::tk::table::Priv(junk)-1}]
    	}}"
	.Tcl(string)
  
	## Internal functions for buttons
	activeRow <- function () as.numeric(tkindex(tktable, "active", "row"))
	
	activeCol <- function () as.numeric(tkindex(tktable, "active", "col"))
	
	undoEdits <- function () {
    	ta <- tclArray()
    	fillTclArrayFromCharMat(ta, cmat)
    	assign("tA", ta, inherits = TRUE)
    	tkconfigure(tktable, variable = tA)
  	}
  	
	finish <- function () tkdestroy(tt)
  	
	cancel <- function () {
    	undoEdits()
    	tkdestroy(tt)
  	}
  	
	insertRow <- function () {
    	row <- activeRow()
    	col <- activeCol()
    	tkinsert(tktable, "rows", row, 1)
    	newCell <- paste(row + 1, col, sep = ",")
    	tkactivate(tktable, newCell)
    	tksee(tktable, newCell)
  	}
  	
	insertCol <- function () {
    	row <- activeRow()
    	col <- activeCol()
    	tkinsert(tktable, "cols", col, 1)
    	newCell <- paste(row, col + 1, sep = ",")
    	tkactivate(tktable, newCell)
    	tksee(tktable, newCell)
  	}
  	
	deleteRow <- function () {
    	if ((row <- activeRow()) != 0)
      		tkdelete(tktable, "rows", row, 1)
  	}
  	
	deleteCol <- function () {
    	if ((col <- activeCol()) != 0)
      		tkdelete(tktable, "cols", col, 1)
  	}
  	
	copyRow <- function () {
    	src <- activeRow()
    	if (src != 0){
      		insertRow()
      		dst <- activeRow()
      		for (j in 0:(ncol(tA) - 1)) tA[[dst, j]] <- tA[[src, j]]
    	}
  	}
  	
	copyCol <- function () {
    	src <- activeCol()
    	if (src != 0) {
      		insertCol()
      		dst <- activeCol()
      		for (i in 0:(nrow(tA) - 1)) tA[[i, dst]] <- tA[[i,src]]
    	}
  	}

  	finishButton <- tkbutton(tt, text = "Finish", command = finish)
  	cancelButton <- tkbutton(tt, text = "Cancel", command = cancel)
  	undoEditsButton <- tkbutton(tt, text = "Undo Edits", command = undoEdits)
  	insertRowButton <- tkbutton(tt, text = "Insert Row", command = insertRow)
  	copyRowButton <- tkbutton(tt, text = "Copy Row", command = copyRow)
  	deleteRowButton <- tkbutton(tt, text = "Delete Row", command = deleteRow)
  	insertColButton <- tkbutton(tt, text = "Insert Col", command = insertCol)
  	copyColButton <- tkbutton(tt, text = "Copy Col", command = copyCol)
  	deleteColButton <- tkbutton(tt, text = "Delete Col", command = deleteCol)

  	## Layout
  	if (length(header) > 0) {
    	for (label in header)
      		tkgrid(tklabel(tt, text = label), columnspan = 7, sticky = "nw")
  	}
  	tkgrid(tktable, yscr, columnspan = 8)
  	tkgrid.configure(tktable, sticky = "news")
  	tkgrid.configure(yscr, sticky = "nsw")
  	tkgrid(xscr, sticky = "new", columnspan = 8)
  	tkgrid(insertRowButton, copyRowButton, deleteRowButton, sticky = "news")
  	tkgrid(insertColButton, copyColButton, deleteColButton,
    	"x", cancelButton, undoEditsButton, finishButton, sticky = "news")
  	tkgrid.columnconfigure(tt, 3, weight = 1)
  	tkgrid.rowconfigure(tt, length(header), weight = 1)
  	tkactivate(tktable, "0,0")
  	tktag.configure(tktable, "active", background = "lightyellow2")
  	tktag.configure(tktable, "title", state = "normal")

  	tkgrab.set(tt)
  	tkfocus(tt)
  	tkwait.window(tt)

  	outMat <- matrix("", nrow = nrow(tA), ncol = ncol(tA))
                   
  	for (i in 1:nrow(outMat)) {
    	for (j in 1:ncol(outMat)) {
      		val <- tA[[i - 1,j - 1]]
      		if (is.null(val)) val <- ""
      		else val <- tclvalue(val)
      		outMat[i,j] <- val
    	}
    }

  	## Recover row and column names
  	rn <- outMat[, 1][-1]
  	cn <- outMat[1, ][-1]
  	outMat <- outMat[-1, -1, drop = FALSE]

  	## Ignore bad and/or NA row and column names
  	badRownames <- c(grep("\\[.*\\]", rn), (1:length(rn))[is.na(rn)])
  	if (length(badRownames) != length(rn)) {
    	rn[badRownames] <- ""
    	rownames(outMat) <- rn
  	}
  	badColnames <- c(grep("\\[.*\\]", cn), (1:length(cn))[is.na(cn)])
  	if (length(badColnames) != length(cn)) {
    	cn[badColnames] <- ""
    	colnames(outMat) <- cn
  	}
  	mode(outMat) <- mode(x)
  	Sys.sleep(0.1)
  	return(outMat)
}
