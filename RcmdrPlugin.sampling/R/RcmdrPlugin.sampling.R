.packageName<-"RcmdrPlugin.EHESsampling"

# The following function written by J. Fox (with contributions from Richard Heiberger) 
# is included to cause the package to load the Rcmdr if it is not already loaded 

.onAttach <- function(libname, pkgname){
    if (!interactive()) return()
    Rcmdr <- options()$Rcmdr
    plugins <- Rcmdr$plugins
    if (!pkgname %in% plugins)  {
        Rcmdr$plugins <- c(plugins, pkgname)
        options(Rcmdr=Rcmdr)
         if("package:Rcmdr" %in% search()) {
            if(!getRcmdr("autoRestart")) {
                closeCommander(ask=FALSE, ask.save=TRUE)
                Commander()
            }
        }
        else {
        Commander()
        } 
    }  
}


# The following was created 2012-08-28 by J. Fox

if (getRversion() >= '2.15.1') globalVariables(c('top', 'buttonsFrame')) #Need to add more in here


### Internal functions for sampling plugin
### Written by Susie Jentoft and Johan Heldal
### last edited 14.11.2012


mk.calc<-function(x,y,z) {
	mk.opt<-x/y
	mk<-ifelse(mk.opt<2,2,mk.opt)
	mk<-ifelse(z==1,1,mk) 
	mk<-ifelse(z<mk,floor(mk),mk)
	mk<-special.round(mk)
	return(mk)
	}

special.round <-function (x){
	uk<-x-floor(x)
	diff<-sum(x)-floor(sum(x))
	number<-ifelse(diff<0.5,floor(sum(uk)),ceiling(sum(uk)))
	order.x<-order(-uk)
	x<-x[order(-uk)]
	x<-c(ceiling(x[if(number==0) 0 else seq(number)]), 
	floor(x[(number+1):length(x)])) 
	return(x[order(order.x)])
	}

# General Disabling function
disable.fn<-function(x,y,z){ 
	if (tclvalue(x)=="1") {
		tkconfigure(y, state="disabled")
		tkconfigure(z, state="disabled")
		}
	if (tclvalue(x)=="0") {
		tkconfigure(y, state="normal")
		tkconfigure(z, state="normal")
		}
	}

disable.fn2<-function(x,y) {
  if (tclvalue(x)=="1") setstate(y, state="disabled")
  if (tclvalue(x)=="0") setstate(y, state="normal")
                }
                

### New variable list box with title label which allows for disabling state
variableListBox2 <- function (parentWindow, variableList = Variables(), bg = "white", 
    selectmode = "single", export = "FALSE", initialSelection = NULL, 
    listHeight = getRcmdr("variable.list.height"), title) 
{
    if (length(variableList) == 1 && is.null(initialSelection)) 
        initialSelection <- 0
    frame <- tkframe(parentWindow)
    minmax <- getRcmdr("variable.list.width")
    listbox <- tklistbox(frame, height = listHeight, 
        selectmode = selectmode, background = bg, exportselection = export, 
        width = min(max(minmax[1], nchar(variableList)), minmax[2]))
    label <- tklabel(frame, text = title, fg = "blue")
    scrollbar <- ttkscrollbar(frame, command = function(...) tkyview(listbox, 
        ...))
    tkconfigure(listbox, yscrollcommand = function(...) tkset(scrollbar, 
        ...))
    for (var in variableList) tkinsert(listbox, "end", var)
    if (is.numeric(initialSelection)) 
        for (sel in initialSelection) tkselection.set(listbox, 
            sel)
    firstChar <- tolower(substr(variableList, 1, 1))
    len <- length(variableList)
    onLetter <- function(letter) {
        letter <- tolower(letter)
        current <- 1 + round(as.numeric(unlist(strsplit(tclvalue(tkyview(listbox)), 
            " "))[1]) * len)
        mat <- match(letter, firstChar[-(1:current)])
        if (is.na(mat)) 
            return()
        tkyview.scroll(listbox, mat, "units")
    }
    onA <- function() onLetter("a")
    onB <- function() onLetter("b")
    onC <- function() onLetter("c")
    onD <- function() onLetter("d")
    onE <- function() onLetter("e")
    onF <- function() onLetter("f")
    onG <- function() onLetter("g")
    onH <- function() onLetter("h")
    onI <- function() onLetter("i")
    onJ <- function() onLetter("j")
    onK <- function() onLetter("k")
    onL <- function() onLetter("l")
    onM <- function() onLetter("m")
    onN <- function() onLetter("n")
    onO <- function() onLetter("o")
    onP <- function() onLetter("p")
    onQ <- function() onLetter("q")
    onR <- function() onLetter("r")
    onS <- function() onLetter("s")
    onT <- function() onLetter("t")
    onU <- function() onLetter("u")
    onV <- function() onLetter("v")
    onW <- function() onLetter("w")
    onX <- function() onLetter("x")
    onY <- function() onLetter("y")
    onZ <- function() onLetter("z")
    for (letter in c(letters, LETTERS)) {
        tkbind(listbox, paste("<", letter, ">", sep = ""), get(paste("on", 
            toupper(letter), sep = "")))
    }
    onClick <- function() tkfocus(listbox)
    toggleSelection <- function() {
        active <- tclvalue(tkindex(listbox, "active"))
        selected <- tclvalue(tkcurselection(listbox))
        if (selected == active) 
            tkselection.clear(listbox, "active")
        else tkselection.set(listbox, "active")
    }
    tkbind(listbox, "<ButtonPress-1>", onClick)
    if (selectmode == "single") 
        tkbind(listbox, "<Control-ButtonPress-1>", toggleSelection)
    if (title!="")  tkgrid(label, columnspan = 2, sticky = "w")
    tkgrid(listbox, scrollbar, sticky = "nw")
    tkgrid.configure(scrollbar, sticky = "wns")
    tkgrid.configure(listbox, sticky = "ew")
    result <- list(frame = frame, listbox = listbox,  label=label, scrollbar = scrollbar, 
        selectmode = selectmode, varlist = variableList)
    class(result) <- "listbox"
    result
}


###
setstate<-function(boxie, state){
	tkconfigure(boxie[[2]], state=state)
	tkconfigure(boxie[[3]], state=state)
	}

###
setstate2<-function(x, y) for (i in 1:length(x)) tkconfigure(x[[i]], state=y)

### New radio buttons with an added command for clicking
radioButtons2 <- defmacro(window=top, name, buttons, values=NULL, initialValue=..values[1], labels, 
	title="", title.color="blue", right.buttons=FALSE, click.command=function(){},
	expr={
		..values <- if (is.null(values)) buttons else values
		..frame <- paste(name, "Frame", sep="")
		assign(..frame, tkframe(window))
		..variable <- paste(name, "Variable", sep="")
		assign(..variable, tclVar(initialValue))
		if(title != ""){
			..title <- paste(name, "Title", sep="")
			assign(..title, labelRcmdr(eval(parse(text=..frame)), text=title, foreground=title.color))
			tkgrid(get(..title), columnspan=2, sticky="w")
		}
		for (i in 1:length(buttons)) {
			..button <- paste(buttons[i], "Button", sep="")
			..buttonlabel<-paste(buttons[i], "Label", sep="")
			assign(..buttonlabel, labelRcmdr(eval(parse(text=..frame)), text=labels[i], justify="left"))
			assign(..button,
				tkradiobutton(eval(parse(text=..frame)), variable=eval(parse(text=..variable)), value=..values[i]))
			if (right.buttons) tkgrid(eval(parse(text=..buttonlabel)), eval(parse(text=..button)), sticky="w")
			else  tkgrid(eval(parse(text=..button)), eval(parse(text=..buttonlabel)), sticky="w")
			tkconfigure(eval(parse(text=..button)), command=click.command)
		}
	}
) 

### add a blank row in a window
tkblank<-function(window) tkgrid(tklabel(window, text=""))


### check that variables have been selected
check.fn2<-function(v){
	mes<-paste("Please select a ", v, " variable", sep="")
 	tkmessageBox(message=mes)
	return()	
	}

### new macro to include tailored onhelp function
OKCancelHelp2<-defmacro(window = top, onHelp, expr = {
    buttonsFrame <- tkframe(window, borderwidth = 5)
    OKbutton <- buttonRcmdr(buttonsFrame, text = gettextRcmdr("OK"), 
        foreground = "darkgreen", width = "12", command = onOK, 
        default = "active", borderwidth = 3)
    onCancel <- function() {
        if (GrabFocus()) tkgrab.release(window)
        tkdestroy(window)
        tkfocus(CommanderWindow())
    }
    cancelButton <- buttonRcmdr(buttonsFrame, text = gettextRcmdr("Cancel"), 
        foreground = "red", width = "12", command = onCancel, 
        borderwidth = 3)
    helpButton <- buttonRcmdr(buttonsFrame, text = gettextRcmdr("Help"), 
            width = "12", command = onHelp, borderwidth = 3)
    tkgrid(OKbutton, labelRcmdr(buttonsFrame, text = "  "), cancelButton, 
        labelRcmdr(buttonsFrame, text = "            "), helpButton, sticky = "w")
})
