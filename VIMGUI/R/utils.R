# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

## widgets
getEscapeChars <- function() {
  c("", "_", "__", "-", " ", "/", "\\", ",")
}
# checkbox
checkboxes <- function(parent, boxes, initial, labels, 
        title, side = c("left","right"), states) {
    n <- length(boxes)
    if(missing(initial)) initial <- rep(FALSE, n)
    side <- match.arg(side)
    if(missing(states)) states <- rep("normal", n)
    frame <- tkframe(parent)
    variables <- lapply(initial, tclVar)
    if(!missing(title)) {
        tkgrid(tklabel(frame, text=title, fg="blue"), columnspan=2, sticky="w")
    }
    cb <- labs <- vector(n, mode="list")
    for(i in 1:n) {
        cb[[i]] <- tkcheckbutton(frame, 
            variable=variables[[i]], state=states[i])
        labs[[i]] <- tklabel(frame, text=labels[i], state=states[i])
        if(side == "left") tkgrid(cb[[i]], labs[[i]])
        else tkgrid(labs[[i]], cb[[i]])
        tkgrid.configure(labs[[i]], sticky="w")
    }
    names(cb) <- boxes
    names(labs) <- boxes
    result <- list(frame=frame, variables=variables, boxes=cb, labels=labs)
    class(result) <- "checkboxes"
    result
}

# ok and cancel buttons
okCancel <- function(window, onOK, parent = NULL) {
    frame <- tkframe(window)
    okButton <- tkbutton(frame, command=onOK, 
        default="active", fg="darkgreen", text="OK", width=10)
    onCancel <- function() closeDialog(window, parent=parent)
    cancelButton <- tkbutton(frame, command=onCancel, 
        fg="red", text="Cancel", width=10)
    tkgrid(okButton, cancelButton, padx=5)
    return <- list(frame=frame, okButton=okButton, cancelButton=cancelButton)
    class(return) <- "okCancel"
    return
}

# listbox
listbox <- function(parent, variables = getVars(), initial = NULL, 
        title, background = "white", exportselection=FALSE, 
        height = 6, selectmode = "single", state = "normal") {
    frame <- tkframe(parent)
    listbox <- tklistbox(frame, background=background, 
        exportselection=exportselection, height=height, 
        selectmode=selectmode, state=state, width=max(20, nchar(variables)))
    scrollbar <- tkscrollbar(frame, repeatinterval=5, 
        command=function(...) tkyview(listbox, ...))
    tkconfigure(listbox, yscrollcommand=function(...) tkset(scrollbar, ...))
    for(var in variables) tkinsert(listbox, "end", var)
    n <- length(variables)
    if(is.character(initial)) initial <- match(initial, variables, 0) - 1
    if(is.numeric(initial)) {
        initial <- initial[!is.na(initial) & 0 <= initial & initial < n]
        if(length(initial)) {
            for(i in initial) tkselection.set(listbox, i)
            #if(initial[1] > height-1) tkyview.moveto(listbox, initial[1]/n)
            tksee(listbox, initial[1])
        }
    }
    # commands for scrolling to a certain letter
    firstChars <- tolower(substr(variables, 1, 1))
    onLetter <- function(letter) {
        letter <- tolower(letter)
        ypos <- as.numeric(unlist(strsplit(tclvalue(tkyview(listbox))," "))[1])
        current <- 1 + round(ypos*n)
        num <- match(letter, firstChars[-(1:current)])
        if(is.na(num) && current > 1) {
            num <- -match(letter, firstChars[(current-1):1])
        }
        if(!is.na(num)) tkyview.scroll(listbox, num, "units")
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
    for(letter in c(letters, LETTERS)) {
        tkbind(listbox, paste("<", letter, ">", sep=""),
            get(paste("on", toupper(letter), sep="")))
    }
    if(!missing(title)) {
        tkgrid(tklabel(frame, text=title, fg="blue"), columnspan=2, sticky="w")
    }
    tkgrid(listbox, scrollbar)
    tkgrid.configure(scrollbar, sticky="nws")
    tkgrid.configure(listbox, sticky="ew")
    result <- list(frame=frame, listbox=listbox, scrollbar=scrollbar)
    class(result) <- "listbox"
    result
}

# radiobuttons
radiobuttons <- function(parent, buttons, values, initial, labels, 
        title, side = c("left","right"), states) {
    n <- length(buttons)
    if(missing(values)) values <- buttons
    if(missing(initial)) initial <- ""
    side <- match.arg(side)
    if(missing(states)) states <- rep("normal", n)
    frame <- tkframe(parent)
    variable <- tclVar(initial)
    if(!missing(title)) {
        tkgrid(tklabel(frame, text=title, fg="blue"), columnspan=2, sticky="w")
    }
    rb <- labs <- vector(n, mode="list")
    for(i in 1:n) {
        rb[[i]] <- tkradiobutton(frame, value=values[i], 
            variable=variable, state=states[i])
        labs[[i]] <- tklabel(frame, text=labels[i], state=states[i])
        if(side == "left") tkgrid(rb[[i]], labs[[i]])
        else tkgrid(labs[[i]], rb[[i]])
        tkgrid.configure(labs[[i]], sticky="w")
    }
    names(rb) <- buttons
    names(labs) <- buttons
    result <- list(frame=frame, variable=variable, buttons=rb, labels=labs)
    class(result) <- "radiobuttons"
    result
}

# ---------------------------------------

# get plot annotation for variables
getLabel <- function(v) {
  sc <- switch(getVm("scaling"), none="", classical="scaled", 
      MCD="robustly scaled", robust="robustly scaled")
  if(nchar(sc)) paste(v, " (", sc, ")", sep="")
  else v
}
defaultNames <- function(p) paste("Var", 1:p, sep="")
## utility functions for widgets

# bind
bind <- function(object, fun, ...) UseMethod("bind")

bind.listbox <- function(object, fun, ...) {
    state <- tclvalue(tkcget(object$listbox, "-state"))
    if(state == "disabled") fun <- function() NULL
    tkbind(object$listbox, "<<ListboxSelect>>", fun)
}

bind.radiobuttons <- function(object, fun, which, ...) {
    buttons <-  object$buttons
    if(!missing(which)) buttons <- buttons[which]
    for(b in buttons) {
        state <- tclvalue(tkcget(b, "-state"))
        if(state == "disabled") fun <- function() NULL
        tkbind(b, "<ButtonRelease>", fun)
    }
}

# close dialog
closeDialog <- function(window, parent = NULL) {
    tkdestroy(window)  # close window
    if(!is.null(parent)) tkfocus(parent)
}

# deselect all
deselectAll <- function(object) UseMethod("deselectAll")

deselectAll.listbox <- function(object) {
    n <- getSize(object)
    if(n > 0) tkselection.clear(object$listbox, 0, n-1)
}

deselectAll.radiobuttons <- function(object) {
    for(b in object$buttons) tkdeselect(b)
}

# empty listbox
empty <- function(object) UseMethod("empty")

empty.listbox <- function(object) {
    n <- getSize(object)
    if(n > 0) tkdelete(object$listbox, 0, n-1)
}

# get selection
getSelection <- function(object, ...) UseMethod("getSelection")

getSelection.listbox <- function(object, variables = getVars(), ...) {
    variables[as.numeric(tkcurselection(object$listbox)) + 1]
}

getSelection.checkboxes <- function(object, which, ...) {
    as.logical(as.numeric(sapply(object$variables, tclvalue)))[which]
}

getSelection.radiobuttons <- function(object, ...) {
    tclvalue(object$variable)
}

# get size
getSize <- function(object) UseMethod("getSize")

getSize.listbox <- function(object) {
    as.numeric(tclvalue(tksize(object$listbox)))
}

getSize.radiobuttons <- function(object) {
    length(object$buttons)
}

# initialize dialog
initializeDialog <- function(title, offset = c(25,10)) {
    window <- tktoplevel(borderwidth=5)  # create window
    if(missing(title)) title <- " "  # default title is empty
    tkwm.title(window, title)  # set title
    if(existsVm(".ttM")) {
        position <- getPosition()  # get position of VIM GUI
        if(!is.null(position) && all(position >= 0)) {
            position <- paste("+", paste(offset+position, collapse="+"), sep="")
            tkwm.geometry(window, position)
        }
    }
    tkwm.resizable(window, 0, 0)  # window not resizable
    tkfocus(window)  # set focus on the new window
    window
}

# insert variables
insert <- function(object, variables) UseMethod("insert")

insert.listbox <- function(object, variables = getVars()) {
    for(var in variables) tkinsert(object$listbox, "end", var)
}

## pack
#pack <- function(object) {
#    tkpack(object$frame, expand=TRUE, fill="x", padx=3, pady=3, side="left")
#}

# get window position of VIM GUI
getPosition <- function () {
    position <- try({
            ID <- getVm(".ttM")$ID
            as.numeric(c(tclvalue(.Tcl(paste("winfo rootx", ID))),
                    tclvalue(.Tcl(paste("winfo rooty", ID)))))
        }, 
        silent=TRUE)
    if(class(position) == "try-error") NULL
    else position
}

# select all
selectAll <- function(object) UseMethod("selectAll")

selectAll.listbox <- function(object) {
    n <- getSize(object)
    if(n > 0) tkselection.set(object$listbox, 0, n-1)
}

# set state
setState <- function(object, state = "normal", ...) UseMethod("setState")

setState.listbox <- function(object, state = "normal", ...) {
    tkconfigure(object$listbox, state=state)
}

setState.radiobuttons <- function(object, state = "normal", which, ...) {
#    n <- getSize(object)
#    for(i in 1:n) {
#        tkconfigure(object$buttons[[i]], state=state)
#        tkconfigure(object$labels[[i]], state=state)
#    }
    if(missing(which)) which <- 1:getSize(object)
    for(i in which) {
        curState <- tclvalue(tkcget(object$buttons[[i]], "-state"))
        if(!(curState %in% c("active", state))) {
            tkconfigure(object$buttons[[i]], state=state)
            tkconfigure(object$labels[[i]], state=state)
        }
    }
}

setState.okCancel <- function(object, state = "normal", ...) {
    tkconfigure(object$okButton, state=state)
}
# set or get active data set
ActiveDataSet <- function(name) {
  if(missing(name)) getVm("activeDataSet")
  else putVm("activeDataSet", name)
}
# ---------------------------------------

## various checks

# check if active data set is selected
checkActiveData <- function() nchar(ActiveDataSet())

# check if variables are selected
checkVars <- function() {
  checkActiveData() && length(getVm("vars")) > 0
}

# check if highlight variables are selected
checkHighlight <- function() {
  checkActiveData() && length(getVm("highlight")) > 0
}

# check requirements for univariate plots
checkUnivar <- function() {
  checkActiveData() &&  length(getVm("vars")) == 1
}

# check requirements for bivariate plots
checkBivar <- function() {
  checkActiveData() && length(getVm("vars")) == 2
}

# check requirements for multivariate plots
checkMultivar <- function() {
  checkActiveData() && length(getVm("vars")) > 1
}

# check requirements for map of missings
checkMap <- function() {
  checkVars() && nchar(getVm("map")) && all(nchar(getVm("coords")))
}

# check requirements for growing dot map with missings
checkGrowdot <- function() {
  checkActiveData() && length(getVm("vars")) == 1 && 
      nchar(getVm("map")) && all(nchar(getVm("coords")))
}

# check requirements for colored map
checkColormap <- function() {
  checkActiveData() && length(getVm("vars")) == 1 && 
      nchar(getVm("map")) && nchar(getVm("region"))
}

## get state for menu items and dialog elements

checkActiveDataS <- function() {
  if(checkActiveData()) "normal" else "disabled"
}
checkVarsS <- function() {
  if(checkVars()) "normal" else "disabled"
}
checkHighlightS <- function() {
  if(checkHighlight()) "normal" else "disabled"
}
checkUnivarS <- function() {
  if(checkUnivar()) "normal" else "disabled"
}
checkBivarS <- function() {
  if(checkBivar()) "normal" else "disabled"
}
checkMultivarS <- function() {
  if(checkMultivar()) "normal" else "disabled"
}
checkMapS <- function() {
  if(checkMap()) "normal" else "disabled"
}
checkGrowdotS <- function() {
  if(checkGrowdot()) "normal" else "disabled"
}
checkColormapS <- function() {
  if(checkColormap()) "normal" else "disabled"
}

# get state for selection method for highlight variables
getSelectionS <- function() {
  if(checkActiveData() && length(getVm("highlight")) > 1) "normal"
  else "disabled"
}

# ---------------------------------------
# get all data.frames
getDataSets <- function(objects, envir = .GlobalEnv) {
  if(missing(objects)) objects <- ls(envir = envir)
  if(length(objects) == 0) return(objects)
  fun <- function(x) {
    if(exists(x, envir = envir)) is.data.frame(get(x, envir = envir))
    else FALSE
  }
  names(which(sapply(objects, fun)))
}

# get variables of a data set specified by name
getVars <- function(name, envir = .GlobalEnv) {
  if(missing(name)) name <- ActiveDataSet()
  if(nchar(name)) names(get(name, envir=envir))
  else character()
}

# save preferences
savePreferences <- function() {
  prefs <- list(col=getVm("col"), alpha=getVm("alpha"), tkr=getVm("tkr"))
  save(prefs, file=".vmGUIprefs.RData")
}

# load preferences
loadPreferences <- function() {
  if(file.exists(".vmGUIprefs.RData")) {
    load(".vmGUIprefs.RData", envir=vmGUIenv())
  }
}


# get vector indicating missings
isNA <- function(x, selection = c("any","all")) {
  selection <- match.arg(selection)
  if(is.null(dim(x))) is.na(x)
  else if(ncol(x) == 1) as.vector(is.na(x))
  else apply(x, 1, function(x) eval(call(selection, is.na(x))))
}

# returns a vector indicating the imputed missings of the current varibale, a vector indicating if the current variable is imputed
# and a vector indicating if there are imputed missings in the other variables
isImp <- function(x, pos, delimiter, imp_var, selection = c("none","any","all")) {
  selection <- match.arg(selection)
  # character vector for possible prefixes for the delimiter
  escape <- getEscapeChars()
  if(is.null(dim(x)) || is.null(dim(imp_var))) list(misspos = imp_var, impp = TRUE, missh = rep(FALSE, NROW(x)))
  else {
    # does the current Variable have imputed missings
    # search escape-vector for possible prefixes
    for(i in 1:length(escape)) {
      indexp <- colnames(imp_var) %in% paste(colnames(x)[pos],delimiter,sep=escape[i])
      # end loop if a match is found
      if(any(indexp))	break
    }
    if(any(indexp)) {
      misspos <- imp_var[,indexp]
      impp <- TRUE
      imp_var <- imp_var[,!indexp, drop = FALSE]
    } else {
      misspos <- rep(FALSE, nrow(x))
      impp <- FALSE
    }
    
    # are there other Variables with missing-indices in the dataset
    # search escape-vector for possible prefixes
    for(i in 1:length(escape)) {
      indexh <- (paste(colnames(x),delimiter,sep=escape[i])) %in% colnames(imp_var)
      # end loop if a match is found
      if(any(indexh)) {
        escape <- escape[i]
        break
      }
    }
    if(any(indexh)) {
      index <- which(indexh)
      tmp <- matrix(nrow = nrow(x), ncol = length(index))
      
      for (i in 1:length(index)) {
        tmp[,i] <- imp_var[,paste(colnames(x)[index[i]],delimiter,sep=escape)]					
      }
      
      if(length(index) > 1 && selection != "none") {
        missh <- apply(tmp, 1, function(tmp) eval(call(selection, tmp)))
      } else {
        missh <- tmp
        colnames(missh) <- colnames(x[,indexh])
      }
    } else {
      missh <- rep(FALSE, nrow(x))
    }
    
    list(misspos = misspos, impp = impp ,missh = missh)
  }
}

# print out which variables are highlighted
highlightInfo <- function(highlight, selection = c("any","all"), imputed = FALSE) {
  if(!imputed) label <- "missings"
  else label <- "imputed missings"
  
  if(length(highlight) == 0) cat(paste("No ", label, " highlighted.\n",sep=""))
  else if(length(highlight) == 1) {
    cat(paste("Highlighted ", label, " in variable ", highlight, ".\n", sep="'"))
  } else {
    selection <- match.arg(selection)
    hlout <- paste(highlight, collapse="', '")
    cat(paste("Highlighted ", label," in ", selection, 
            " of the variables '", hlout, "'.\n", sep=""))
  }
  invisible()
}

# count infinite values
countInf <- function(x) length(which(is.infinite(x)))

# count missings
countNA <- function(x) length(which(is.na(x)))

# count imputed missings
countImp <- function(x, delimiter, imp_var) {
  # character vector for possible prefixes for the delimiter
  escape <- getEscapeChars()
  # search escape-vector for possible prefixes
  for(i in 1:length(escape)) {
    indexh <- (paste(colnames(x),delimiter,sep=escape[i])) %in% colnames(imp_var)
    # end loop if a match is found
    if(any(indexh)) {
      escape <- escape[i]
      break
    }
  }
  tmp <-integer(ncol(x))
  names(tmp) <- colnames(x)
  for ( i in 1:ncol(x)) {
    tmp[i] <- ifelse(indexh[i],length(which(imp_var[,paste(colnames(x)[i],delimiter,sep=escape)])),0)
  }
  tmp
}


# test means for boxplot with missings
testMeans <- function(x, pos = 1, selection = c("any","all")) {
  selection <- match.arg(selection)
  ind <- isNA(x[, -pos], selection)
  x1 <- x[ind, pos]
  x2 <- x[!ind, pos]
  if(length(which(!is.na(x1))) > 1 && length(which(!is.na(x2)))) {
    list(ind=ind, p.v=t.test(x1, x2)$p.v)
  } else list(ind=ind, p.v=NA)
}