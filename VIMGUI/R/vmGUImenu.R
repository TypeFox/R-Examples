# ----------------------------------------------------------
# Authors: Andreas Alfons, Bernd Prantner and Matthias Templ
#          Vienna University of Technology
# ----------------------------------------------------------



#' GUI for Visualization and Imputation of Missing Values
#' 
#' Graphical user interface for visualization and imputation of missing values.
#' 
#' The \emph{Data} menu allows to select a data set from the workspace or load
#' data into the workspace from \code{RData} files.  Furthermore, it can be
#' used to transform variables, which are then appended to the data set in use.
#' Commonly used transformations in official statistics are available, e.g.,
#' the Box-Cox transformation and the log-transformation as an important
#' special case of the Box-Cox transformation.  In addition, several other
#' transformations that are frequently used for compositional data are
#' implemented. Background maps and coordinates for spatial data can be
#' selected in the data menu as well.
#' 
#' After a data set was chosen, variables can be selected in the main menu,
#' along with a method for scaling.  An important feature is that the variables
#' will be used in the same order as they were selected, which is especially
#' useful for parallel coordinate plots.  Variables for highlighting are
#' distinguished from the plot variables and can be selected separately.  For
#' more than one variable chosen for highlighting, it is possible to select
#' whether observations with missing values in any or in all of these variables
#' should be highlighted.
#' 
#' A plot method can be selected from the \emph{Visualization} menu. Note that
#' plots that are not applicable to the selected variables are disabled, for
#' example, if only one plot variable is selected, multivariate plots cannot be
#' chosen.
#' 
#' The \emph{Imputation} menu offers robust imputation methods to impute
#' variables of the data set.
#' 
#' The \emph{Diagnostics} menu is similar to the \emph{Visualization} menu, but
#' is designed to verify the results after the imputation of missing values.
#' 
#' Last, but not least, the \emph{Options} menu allows to set the colors, alpha
#' channel and the delimiter for imputed variables to be used in the plots.  In
#' addition, it contains an option to embed multivariate plots in \code{Tcl/Tk}
#' windows. This is useful if the number of observations and/or variables is
#' large, because scrollbars allow to move from one part of the plot to
#' another.
#' 
#' Internal information regarding the VIM GUI is stored in the environment
#' \code{vmGUIenvir}.
#' 
#' @aliases vmGUImenu vmGUIenvir
#' @author Andreas Alfons, based on an initial design by Matthias Templ,
#' modifications by Bernd Prantner
#' @references M. Templ, A. Alfons, P. Filzmoser (2012) Exploring incomplete
#' data using visualization tools.  \emph{Journal of Advances in Data Analysis
#' and Classification}, Online first. DOI: 10.1007/s11634-011-0102-y.
#' @keywords multivariate hplot
#' @export vmGUImenu
vmGUImenu <- function() {
    ## VisualizationMenunitializations
    # retrieve values of internal variables if they exist
    activeDataSet <- if(existsVm("activeDataSet")) ActiveDataSet() else ""
    vars <- if(existsVm("vars")) getVm("vars") else character()
	imp_vars <- if(existsVm("imp_vars")) getVm("imp_vars") else character()
    scaling <- if(existsVm("scaling")) getVm("scaling") else "none"
    highlight <- if(existsVm("highlight")) getVm("highlight") else character()
    selection <- if(existsVm("selection")) getVm("selection") else "any"
    tvars <- if(existsVm("tvars")) getVm("tvars") else character()
    transform <- if(existsVm("transform")) getVm("transform") else ""
    alrVar <- if(existsVm("alrVar")) getVm("alrVar") else ""
    map <- if(existsVm("map")) getVm("map") else ""
    coords <- if(existsVm("coords")) getVm("coords") else rep("", 2)
    region <- if(existsVm("region")) getVm("region") else ""
    loadPreferences()  # load preferences from file
    col <- if(existsVm("col")) getVm("col") 
        else c("skyblue","red","skyblue4","red4","orange","orange4")
    alpha <- if(existsVm("alpha")) getVm("alpha") else 0.6
    tkr <- if(existsVm("tkr")) getVm("tkr") else FALSE
	delimiter <- if(existsVm("delimiter")) getVm("delimiter") else "_imp"
    # save initial values in environment
    ActiveDataSet(activeDataSet)
    putVm("vars", vars)
	putVm("imp_vars", imp_vars)
    putVm("scaling", scaling)
    putVm("highlight", highlight)
    putVm("selection", selection)
    putVm("tvars", tvars)
    putVm("transform", transform)
    putVm("alrVar", alrVar)
    putVm("map", map)
    putVm("coords", coords)
    putVm("region", region)
    putVm("col", col)
    putVm("alpha", alpha)
    putVm("tkr", tkr)
	putVm("delimiter", delimiter)
    
    ## function to set state of menu items ("normal" or "disabled")
    activateMenus <- function() {
        ## top menu
        state <- checkVarsS()
        for(i in 2:4) 
            .Tcl(paste(topMenu$ID, "entryconfigure", i, "-state", state))
        ## data menu
        .Tcl(paste(DataMenu$ID, "entryconfigure", 2, 
                "-state", checkActiveDataS()))
        .Tcl(paste(DataMenu$ID, "entryconfigure", 3, 
                "-state", checkActiveDataS()))
        ## visualization and imputation menu
        # uni- and bivariate plots
        state <- checkUnivarS()
        for(i in 1:6) {
            .Tcl(paste(VisualizationMenu$ID, "entryconfigure", i, 
                    "-state", state))
			.Tcl(paste(DiagnosticsMenu$ID, "entryconfigure", i, 
					"-state", state))
		}
        state <- checkBivarS()
        for(i in 7:9) {
            .Tcl(paste(VisualizationMenu$ID, "entryconfigure", i, 
                    "-state", state))
			.Tcl(paste(DiagnosticsMenu$ID, "entryconfigure", i, 
				"-state", state))
		}
        # multivariate plots
        state <- checkMultivarS()
        for(i in 10:14)  {
            .Tcl(paste(VisualizationMenu$ID, "entryconfigure", i, 
                    "-state", state))
			.Tcl(paste(DiagnosticsMenu$ID, "entryconfigure", i, 
					"-state", state))
		}
        # maps
		state <- checkMapS()
        .Tcl(paste(VisualizationMenu$ID, "entryconfigure", 15, 
                "-state", state))
		.Tcl(paste(DiagnosticsMenu$ID, "entryconfigure", 15, 
				"-state", state))
		
		state <- checkGrowdotS()
        .Tcl(paste(VisualizationMenu$ID, "entryconfigure", 16, 
                "-state", state))
		.Tcl(paste(DiagnosticsMenu$ID, "entryconfigure", 16, 
				"-state", state))
		
		state <- checkColormapS()
        .Tcl(paste(VisualizationMenu$ID, "entryconfigure", 17, 
                "-state", checkColormapS()))
		.Tcl(paste(DiagnosticsMenu$ID, "entryconfigure", 17, 
				"-state", state))
    }
    
    ## functions bound to top menu buttons
    # close GUI
    # ----------
    # comment during development
    Quit <- function() {
        yesno <- tkmessageBox(message=paste("Close GUI for visualization", 
                "and imputation of missing values?"), icon="question", 
            type="yesno", default="no", parent=ttM, title="Quit")
        if(tclvalue(yesno) == "yes") {
            savePreferences()
            closeDialog(ttM)
            rmVm(".ttM")
        }
    }
    # ----------
    # during development replaced by following function that quits immediately 
    # (to save time)
#    Quit <- function() {
#        savePreferences()
#        closeDialog(ttM)
#        rmVm(".ttM")
#    }
    # ----------
    # not implemented yet
    NotImplemented <- function() {
        tkmessageBox(message="Not implemented yet.", icon="error", parent=ttM)
    }
    
    ## functions bound to data menu buttons
    # select active data set
    vmGUIdata <- function() {
        # start dialog
        ttD <- initializeDialog("Select Data")
        # initializations
        activeDataSet <- ActiveDataSet()
        putVm(".activeDataSet", activeDataSet)
        dataSets <- getDataSets()
        # listbox
        dataFrame <- tkwidget(ttD, "labelframe", 
            text="Select Data", fg="blue")
        dataBox <- listbox(dataFrame, variables=dataSets, 
            initial=activeDataSet, height=6)
        setData <- function() {
            selDataSet <- getSelection(dataBox, variables=dataSets)
            putVm(".activeDataSet", selDataSet)
            activateOK()
        }
        bind(dataBox, setData)
        # ok and cancel buttons
        onOK <- function() {
            selDataSet <- getVm(".activeDataSet")
            if(selDataSet != activeDataSet) {
				# get imputed variables
				delimiter <- getVm("delimiter")
				vars <- colnames(get(selDataSet, envir=.GlobalEnv))
				imp_vars <- grep(delimiter, colnames(get(selDataSet, envir=.GlobalEnv)), value = TRUE)
				vars <- setdiff(vars, imp_vars)
				ActiveDataSet(selDataSet)
                putVm("vars", character())
                putVm("highlight", character())
				putVm("imp_vars", imp_vars)
                putVm("coords", rep("", 2))
                putVm("region", "")
                activateElements()
                activateMenus()
                empty(varsBox)
                deselectAll(varsButtons)
				insert(varsBox, vars)
                empty(highlightBox)
                deselectAll(highlightButtons)
                insert(highlightBox, vars)
            }
            closeDialog(ttD, parent=ttM)
        }
        buttons <- okCancel(ttD, onOK, parent=ttM)
        activateOK <- function() {
            okS <- if(nchar(getVm(".activeDataSet"))) "normal" else "disabled" 
            setState(buttons, okS)
        }
        # display dialog elements
        tkpack(dataBox$frame, expand=TRUE, 
            fill="x", padx=3, pady=3, side="left")
        tkgrid(dataFrame, padx=10, pady=5, sticky="news")
            tkgrid(buttons$frame)
        activateOK()
    }
    # load R data
    LoadRData <- function() {
        fileName <- tclvalue(tkgetOpenFile(parent=ttM, title="Load R Data", 
                filetypes=paste("{{R Data Files} {.RData .Rdata .rdata .rda}}", 
                    "{{All files} *}")))
        if(nchar(fileName)) {
            nam <- try(load(fileName, envir=.GlobalEnv))
            if(!inherits(class(nam), "try-error")) { 
                dataSet <- getDataSets(nam)
                if(length(dataSet) == 0) {
                    msg <- "File '%s' does not contain a data set.\n"
                    tkmessageBox(message=gettextf(msg, fileName), 
                        icon="error", parent=ttM)
                    
                } else if(length(dataSet) == 1) {
                    if(dataSet != ActiveDataSet()) {
						# get imputed variables
						delimiter <- getVm("delimiter")
						vars <- colnames(get(dataSet, envir=.GlobalEnv))
						imp_vars <- grep(delimiter, colnames(get(dataSet, envir=.GlobalEnv)), value = TRUE)
						vars <- setdiff(vars, imp_vars)
						ActiveDataSet(dataSet)
                        putVm("vars", character())
                        putVm("highlight", character())
                        putVm("coords", rep("", 2))
                        putVm("region", "")
                        activateElements()
                        activateMenus()
                        empty(varsBox)
                        deselectAll(varsButtons)
                        insert(varsBox, vars)
                        empty(highlightBox)
                        deselectAll(highlightButtons)
                        insert(highlightBox, vars)
                    }
                } else {
                    msg <- paste("File '%s' contains more than one data set.", 
                        "Please select one using '%s'.\n")
                    tkmessageBox(message=gettextf(msg, fileName, "Select Data"),
                        icon="warning", parent=ttM)
                }
            }
        }
    }
    # transform variables
    vmGUItransform <- function() {
        # start dialog
        ttT <- initializeDialog("Transform Variables")
        # initializations
        allVars <- getVars()
        tvars <- getVm("tvars")
        transform <- getVm("transform")
        alrVar <- getVm("alrVar")
        putVm(".tvars", tvars)
        putVm(".transform", transform)
        putVm(".alrVar", alrVar)
        # frames
        topFrame <- tkframe(ttT)
        tvarsFrame <- tkwidget(topFrame, "labelframe", 
            text="Select Variables", fg="blue")
        transformFrame <- tkwidget(topFrame, "labelframe", 
            text="Select Transformation", fg="blue")
        ## left side
        # listbox for variables
        tvarsBox <- listbox(tvarsFrame, initial=tvars, 
            height=12, selectmode="extended")
        setTvars <- function() {
            tvars <- getSelection(tvarsBox)
            putVm(".tvars", tvars)
            deselectAll(tvarsButtons)
            activateElements()
        }
        bind(tvarsBox, setTvars)
        # select and deselect all variables radiobuttons
        tvarsButtons <- radiobuttons(tvarsFrame, 
            buttons=c("tvarsSelectAll","tvarsDeselectAll"), 
            labels=c("Select all","Deselect all"))
        tvarsSelectAll <- function() {
            putVm(".tvars", allVars)
            selectAll(tvarsBox)
            activateElements()
        }
        tvarsDeselectAll <- function() {
            putVm(".tvars", character())
            deselectAll(tvarsBox)
            activateElements()
        }
        bind(tvarsButtons, tvarsSelectAll, "tvarsSelectAll")
        bind(tvarsButtons, tvarsDeselectAll, "tvarsDeselectAll")
        tkpack(tvarsBox$frame, tvarsButtons$frame, 
            expand=TRUE, fill="x", padx=3, pady=3, side="left")
        ## right side
        # radiobuttons for transformation
        transformButtons <- radiobuttons(transformFrame, 
            buttons=c("minus","reciprocal","logarithm",
                "exponential","boxcox","clr","ilr","alr"), 
            labels=c("Minus","Reciprocal","Logarithm","Exponential",
                "Box-Cox","Centered logratio","Isometric logratio",
                "Additive logratio with ratio variable"), 
            initial=transform)
        setTransform <- function() {
            transform <- getSelection(transformButtons)
            putVm(".transform", transform)
            activateElements()
        }
        bind(transformButtons, setTransform)
        # ratio variable for logratio transformation
        alrVariable <- tclVar(alrVar)
        setAlrVar <- function() {
            tkfocus(ttT)
            putVm(".alrVar", tclvalue(alrVariable))
            activateElements()
        }
        alrComboBox <- tkwidget(transformButtons$frame, "ComboBox", 
            "-modifycmd", setAlrVar, editable=FALSE, values=allVars, 
            width=15, textvariable=alrVariable)
        tkgrid(alrComboBox, row=8, column=1, sticky="w")
        tkpack(transformButtons$frame, expand=TRUE, 
            fill="x", padx=3, pady=3, side="left")
        # ok and cancel buttons
        onOK <- function() {
            tvars <- getVm(".tvars")
            putVm("tvars", tvars)
            transform <- getVm(".transform")
            putVm("transform", transform)
            activeData <- get(ActiveDataSet(), envir=.GlobalEnv)
            if(transform == "alr") {
                alrVar <- getVm(".alrVar")
                putVm("alrVar", alrVar)
                x <- subset(activeData, select=union(tvars, alrVar))
                tx <- prepare(x, transformation=transform, alrVar=alrVar)
                colnames(tx) <- paste("alr(", colnames(tx), ")", sep="")
            } else {
                x <- subset(activeData, select=tvars)
                tx <- prepare(x, transformation=transform)
            }
            if(transform == "minus") colnames(tx) <- paste("-", tvars, sep="")
            else if(transform == "reciprocal") 
                colnames(tx) <- paste("1/", tvars, sep="")
            else if(!(transform %in% c("ilr","alr"))){
                prefix <- switch(transform, 
                    logarithm="log", exponential="exp", 
                    boxcox="boxcox", clr="clr")
                colnames(tx) <- paste(prefix, "(", tvars, ")", sep="")
            }
            assign(ActiveDataSet(), cbind(activeData, tx), envir=vmGUIenv())
			getVm("ActiveDataSet()")
            insert(varsBox, colnames(tx))
            insert(highlightBox, colnames(tx))
            closeDialog(ttT, parent=ttM)
        }
        # ok and cancel buttons
        buttons <- okCancel(ttT, onOK, parent=ttM)
        activateElements <- function() {
            tvars <- getVm(".tvars")
            transform <- getVm(".transform")
            alrVar <- getVm(".alrVar")
            checkTvars <- length(tvars)
            transformS <- if(checkTvars) "normal" else "disabled"
            logratioS <- if(length(tvars) > 1) "normal" else "disabled"
            setState(transformButtons, transformS, 
                which=setdiff(names(transformButtons$buttons), c("clr","ilr")))
            setState(transformButtons, logratioS, which=c("clr","ilr"))
            checkAlr <- checkTvars && transform == "alr"
            alrS <- if(checkAlr) "normal" else "disabled"
            tkconfigure(alrComboBox, state=alrS)
            checkOK <- 
                if(transform %in% c("clr","ilr")) length(tvars) > 1
                else if(transform == "alr") 
                    nchar(alrVar) && length(union(tvars, alrVar)) > 1
                else length(tvars) && nchar(transform)
            okS <- if(checkOK) "normal" else "disabled"
            setState(buttons, okS)
        }
        # display dialog elements
        tkpack(tvarsFrame, fill="y", padx=10, pady=5, side="left")
        tkpack(transformFrame, fill="y", padx=10, pady=5, side="right")
        tkpack(topFrame, side="top")
        tkpack(buttons$frame, side="bottom")
        activateElements()
    }
    # background map
    vmGUImap <- function() {
        # start dialog
        ttMap <- initializeDialog("Background Map")
        # initializations
        map <- getVm("map")
        coords <- getVm("coords")
        region <- getVm("region")
        putVm(".map", map)
        putVm(".coords", coords)
        putVm(".region", region)
        topFrame <- tkframe(ttMap)
        ## left side
        # background map
        mapFrame <- tkwidget(topFrame, "labelframe", 
            text="Select Background Map", fg="blue")
        maps <- ls(envir=.GlobalEnv)
        mapBox <- listbox(mapFrame, variables=maps, initial=map, height=6)
        setMap <- function() {
            map <- getSelection(mapBox, variables=maps)
            putVm(".map", map)
            activateOK()
        }
        bind(mapBox, setMap)
        tkpack(mapBox$frame, expand=TRUE, fill="x", padx=3, pady=3, side="left")
        ## right side
        toprightFrame <- tkframe(topFrame)
        # coordinates
        coordsFrame <- tkwidget(toprightFrame, "labelframe", 
            text="Select Coordinates", fg="blue")
        xyFrame <- tkframe(coordsFrame)
        vars <- c("", getVars())
        xVariable <- tclVar(coords[1])
        yVariable <- tclVar(coords[2])
        setCoords <- function() {
            tkfocus(ttMap)
            coords <- c(tclvalue(xVariable), tclvalue(yVariable))
            putVm(".coords", coords)
            activateOK()
        }
        xComboBox <- tkwidget(xyFrame, "ComboBox", 
            "-modifycmd", setCoords, editable=FALSE, values=vars, 
            width=15, textvariable=xVariable)
        yComboBox <- tkwidget(xyFrame, "ComboBox", 
            "-modifycmd", setCoords, editable=FALSE, values=vars, 
            width=15, textvariable=yVariable)
        tkgrid(tklabel(xyFrame, text="x-Coordinate: "), xComboBox, sticky="w")
        tkgrid(tklabel(xyFrame, text="y-Coordinate: "), yComboBox, sticky="w")
        tkgrid(xyFrame, padx=3, pady=3)
        # region
        regionFrame <- tkwidget(toprightFrame, "labelframe", 
            text="Set Region Variable", fg="blue")
        rFrame <- tkframe(regionFrame)
        regionVariable <- tclVar(region)
        setRegion <- function() {
            tkfocus(ttMap)
            region <- tclvalue(regionVariable)
            putVm(".region", region)
            activateOK()
        }
        regionComboBox <- tkwidget(rFrame, "ComboBox", 
            "-modifycmd", setRegion, editable=FALSE, values=vars, 
            width=15, textvariable=regionVariable)
        tkpack(tklabel(rFrame, text="Region: "), side="left")
        tkpack(regionComboBox, side="right")
        tkpack(rFrame, fill="x", padx=3, pady=3)
        # ok and cancel buttons
        onOK <- function() {
            putVm("map", getVm(".map"))
            putVm("coords", getVm(".coords"))
            putVm("region", getVm(".region"))
            activateMenus()
            closeDialog(ttMap, parent=ttM)
        }
        buttons <- okCancel(ttMap, onOK, parent=ttM)
        activateOK <- function() {
            coords <- getVm(".coords")
            ncCoords <- nchar(coords)
            ok <- nchar(getVm(".map")) && 
                if(any(ncCoords)) all(ncCoords) && coords[1] != coords[2]
                else nchar(getVm(".region"))
            okS <- if(ok) "normal" else "disabled"
            setState(buttons, okS)
        }
        # display dialog elements
        tkpack(mapFrame, fill="y", padx=10, pady=5, side="left")
        tkpack(coordsFrame, fill="x", padx=10, pady=5, side="top")
        tkpack(regionFrame, fill="x", padx=10, pady=5, side="bottom")
        tkpack(toprightFrame, fill="y", side="right")
        tkpack(topFrame, side="top")
        tkpack(buttons$frame, side="bottom")
        activateOK()
    }
    knnGUI <- function(){
      onOKkNN <- function(){
        data <- get(ActiveDataSet()) 
        vars <- getVm("vars")
        distance <- getVm("distance")
        k <- getVm("k")
        assign(paste(ActiveDataSet(),"_IMPUTED",sep=""),
        kNN(data,variable=vars,dist_var=distance,k=k,
			imp_suffix=getVm("delimiter")), envir=vmGUIenv()) #new
		getVm(paste(ActiveDataSet(),"_IMPUTED",sep=""))  #new 
        cmd <- paste(paste(ActiveDataSet(),"_IMPUTED",sep=""),
                " <- kNN(",ActiveDataSet(),",variable=c(",paste("\"",vars,"\"",sep="",collapse=","),
                "),dist_var=c(",paste("\"",distance,"\"",sep="",collapse=","),"),
				k=",k,",imp_suffix=\"",getVm("delimiter"),"\")", sep="")
        cat(cmd,"\n")
		closeDialog(ttIMP)
      }
      ttIMP <- initializeDialog("kNN - Imputation")
      ## frames
      varsFrame <- tkwidget(ttIMP, "labelframe", 
          text="Select Variables to impute", fg="blue")
      distanceFrame <- tkwidget(ttIMP, "labelframe", 
          text="Select Variables to compute distances", fg="blue")
      kFrame <- tkwidget(ttIMP, "labelframe", 
          text="Select the number of nearest neighbours", fg="blue")
      ## dialog elements
      # listbox for plot variables
      varsBox1 <- listbox(varsFrame, initial=vars, 
          height=6, selectmode="extended")
      distanceBox1 <- listbox(distanceFrame,  
          initial=vars, height=6, selectmode="extended")
      kBox1 <- listbox(kFrame,variables=1:10,  
          initial=2, height=6, selectmode="single")
      # select and deselect all variables radiobuttons
      # selection method for highlight variables
      tkpack(varsBox1$frame,  
          expand=TRUE, fill="x", padx=3, pady=3, side="left")
      tkpack(distanceBox1$frame, 
          expand=TRUE, fill="x", padx=3, pady=3, side="left")
      tkpack(kBox1$frame, 
          expand=TRUE, fill="x", padx=3, pady=3, side="left")
      
      tkgrid(varsFrame,  padx=10, pady=5, sticky="news")
      tkgrid(distanceFrame,  padx=10, pady=5, sticky="news")   
      tkgrid(kFrame,  padx=10, pady=5, sticky="news")
      buttons <- okCancel(ttIMP, onOKkNN, parent=ttM)
      tkgrid(buttons$frame)
      setVars1 <- function() {
        vars <- getSelection(varsBox1)
        putVm("vars", vars)
      }
      setdistance <- function() {
        distance <- getSelection(distanceBox1)
        putVm("distance", distance)
      }
      setk <- function() {
        k <- getSelection(kBox1,variables=1:10)
        putVm("k", k)
      }
      bind(varsBox1, setVars1)
      bind(kBox1, setk)
      bind(distanceBox1, setdistance)
      setVars1();setk();setdistance();
    
    }
    hotdeckGUI  <- function(){
      getVarsF <- function(){
        vv <- vector()
        for(v in getVars()){
          x <- is.numeric(subset(get(ActiveDataSet(), envir=.GlobalEnv),select=v))
          if(x)
            vv <- c(vv,v)
        }
        vv
      }
      getVarsN <- function(){
        vv <- vector()
        for(v in getVars()){
          x <- is.numeric(subset(get(ActiveDataSet(), envir=.GlobalEnv),select=v))
          if(!x)
            vv <- c(vv,v)
        }
        vv
      }
      onOKhotdeck <- function(){
        data <- get(ActiveDataSet()) 
        vars <- getVm("vars")
        sort <- getVm("sort")
        domain <- getVm("domain")
        if(length(domain)<1)
          domain <- NULL
        if(length(sort)<1)
          sort <- NULL
        if(length(vars)<1){
          tkmessageBox(message=c("Select at least one variable to impute!"),
              icon="warning", parent=ttM)
        }else{
			  assign(paste(ActiveDataSet(),"_IMPUTED",sep=""),
			  hotdeck(data,variable=vars,ord_var=sort,domain_var=domain,
					  imp_suffix=getVm("delimiter")), envir=vmGUIenv()) #new
	  getVm(paste(ActiveDataSet(),"_IMPUTED",sep=""))  #new 
	  
          if(is.null(sort))
            ord_var <- "NULL"
          else
            ord_var <- paste("c(",paste("\"",sort,"\"",sep="",collapse=","),")",sep="")
          if(is.null(domain))
            domain_var <- "NULL"
          else
            domain_var <- paste("c(",paste("\"",domain,"\"",sep="",collapse=","),")",sep="")
          cmd <- paste(paste(ActiveDataSet(),"_IMPUTED",sep=""),
                " <- hotdeck(",ActiveDataSet(),",variable=c(",paste("\"",vars,"\"",sep="",collapse=","),
                "),ord_var=",ord_var,",domain_var=",domain_var,",imp_suffix=\"",getVm("delimiter"),"\")",sep="")
          cat(cmd,"\n")
        }
		closeDialog(ttIMP)
      }
      ttIMP <- initializeDialog("Hotdeck - Imputation")
      ## frames
      varsFrame <- tkwidget(ttIMP, "labelframe", 
          text="Select Variables to impute", fg="blue")
      sortFrame <- tkwidget(ttIMP, "labelframe", 
          text="Select Variables to sort", fg="blue")
      domainFrame <- tkwidget(ttIMP, "labelframe", 
          text="Select Variables to build domains", fg="blue")
      ## dialog elements
      # listbox for plot variables
      varsBox1 <- listbox(varsFrame, initial=vars, 
          height=6, selectmode="extended")
      sortBox1 <- listbox(sortFrame,variables=getVarsN(),  
          initial=NULL, height=6, selectmode="extended")
      domainBox1 <- listbox(domainFrame,variables=getVarsF(),  
          initial=NULL, height=6, selectmode="extended")
      # select and deselect all variables radiobuttons
      # selection method for highlight variables
      tkpack(varsBox1$frame,  
          expand=TRUE, fill="x", padx=3, pady=3, side="left")
      tkpack(sortBox1$frame, 
          expand=TRUE, fill="x", padx=3, pady=3, side="left")
      tkpack(domainBox1$frame, 
          expand=TRUE, fill="x", padx=3, pady=3, side="left")
      
      tkgrid(varsFrame,  padx=10, pady=5, sticky="news")
      tkgrid(sortFrame,  padx=10, pady=5, sticky="news")   
      tkgrid(domainFrame,  padx=10, pady=5, sticky="news")
      buttons <- okCancel(ttIMP, onOKhotdeck, parent=ttM)
      tkgrid(buttons$frame)
      setVars1 <- function() {
        vars <- getSelection(varsBox1)
        putVm("vars", vars)
      }
      setsort <- function() {
        sort <- getSelection(sortBox1,variables=getVarsN())
        putVm("sort", sort)
      }
      setdomain <- function() {
        domain <- getSelection(domainBox1,variables=getVarsF())
        putVm("domain", domain)
      }
      bind(varsBox1, setVars1)
      bind(domainBox1, setdomain)
      bind(sortBox1, setsort)
      setVars1();setdomain();setsort();
    }
    irmiGUI  <- function(){
      getVarsN <- function(){
        vv <- vector()
        for(v in getVars()){
          x <- is.numeric(subset(get(ActiveDataSet(), envir=.GlobalEnv),select=v))
          if(!x)
            vv <- c(vv,v)
        }
        vv
      }
      onOKirmi <- function(){
        x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), 
            select=getVm("vars"))
        vars=getVm("vars")
        mixed=getVm("mixed")
        robust=getVm("robust")
		assign(paste(ActiveDataSet(),"_IMPUTED",sep=""),
				irmi(x,mixed=mixed,robust=robust), envir=vmGUIenv()) #new
		getVm(paste(ActiveDataSet(),"_IMPUTED",sep=""))  #new 
		
        if(length(mixed)>1)
          mixed <- paste("c(",paste("\"",mixed,"\"",sep="",collapse=","),")",sep="")
        else mixed <- "NULL"
        cmd <- 
            paste(paste(ActiveDataSet(),"_IMPUTED",sep=""),
                " <- irmi(",ActiveDataSet(),"[,c(",
                paste("\"",vars,"\"",sep="",collapse=","),")],mixed=",
                mixed,",robust=",robust,")"
                ,sep="")
        cat(cmd,"\n")
		closeDialog(ttIMP)
      }
      ttIMP <- initializeDialog("IRMI")
      ## frames
      varsFrame <- tkwidget(ttIMP, "labelframe", 
          text="Select Variables to impute", fg="blue")
      mixedFrame <- tkwidget(ttIMP, "labelframe", 
          text="Select semi-continous variables", fg="blue")
      robustFrame <- tkwidget(ttIMP, "labelframe", 
          text="robust/non-robust", fg="blue")
      # listbox for plot variables
      varsBox1 <- listbox(varsFrame, initial=vars, 
          height=6, selectmode="extended")
      mixedBox1 <- listbox(mixedFrame,variables=getVarsN(),  
          initial=NULL, height=6, selectmode="extended")
      robustButtons <- radiobuttons(robustFrame, 
          buttons=c("robust","classical"), 
          labels=c("robust","classical"),initial="classical")
     
      # select and deselect all variables radiobuttons
      # selection method for highlight variables
      tkpack(varsBox1$frame,  
          expand=TRUE, fill="x", padx=3, pady=3, side="left")
      tkpack(mixedBox1$frame, 
          expand=TRUE, fill="x", padx=3, pady=3, side="left")
      tkpack(robustButtons$frame, 
          expand=TRUE, fill="x", padx=3, pady=3, side="left")
      tkgrid(varsFrame, mixedFrame,robustFrame, padx=10, pady=5, sticky="news")
      buttons <- okCancel(ttIMP, onOKirmi, parent=ttM)
      tkgrid(buttons$frame)
      setVars1 <- function() {
        vars <- getSelection(varsBox1)
        putVm("vars", vars)
      }
      setmixed <- function() {
        mixed <- getSelection(mixedBox1,variables=getVarsN())
        putVm("mixed", mixed)
      }
      setrobust <- function(){
        robust <- tclvalue(robustButtons$variable)=="robust"
        putVm("robust", robust)
      }
      bind(varsBox1, setVars1)
      bind(mixedBox1, setmixed)
      bind(robustButtons,setrobust)
      setVars1();setmixed(); setrobust();    
    }
    ## functions bound to visualization and diagnostic menu buttons
    # aggregate missings
    Aggr <- function() {
        x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), 
            select=getVm("vars"))
        res <- aggr(x, plot=FALSE)
        cat("\n\nPrint method:\n")
        print(res)
        cat("\n\nSummary method:\n")
        print(summary(res))
        cat("\n")
        if(getVm("tkr")) {
            TKRaggr(x, col=getVm("col")[1:2], numbers=TRUE, prop=c(TRUE, FALSE))
        } else {
            dev.new()
            plot(res, col=getVm("col")[1:2], numbers=TRUE, prop=c(TRUE, FALSE))
        }
    }
	# aggregate missings and imputed missings
	AggrImp <- function() {
		vars <- c(getVm("vars"), getVm("imp_vars"))
		x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), 
				select=vars)
		res <- aggr(x, delimiter=getVm("delimiter"), plot=FALSE)
		cat("\n\nPrint method:\n")
		print(res)
		cat("\n\nSummary method:\n")
		print(summary(res))
		cat("\n")
		if(getVm("tkr")) {
			TKRaggr(x, delimiter=getVm("delimiter"), col=getVm("col")[c(1,2,5)], numbers=TRUE, prop=c(TRUE, FALSE))
		} else {
			dev.new()
			plot(res, col=getVm("col")[c(1,2,5)], numbers=TRUE, prop=c(TRUE, FALSE))
		}
	}
	
    # histogram with missings
    HistMiss <- function() {
        vars <- union(getVm("vars"), getVm("highlight"))
        x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
        x[,1] <- prepare(x[,1], scaling=getVm("scaling"))
        dev.new()
        histMiss(x, selection=getVm("selection"), 
            col=getVm("col"), xlab=getLabel(vars[1]))
        #leg <- paste(c("observed in","missing in"), vars[2])
        #legend("topleft", legend=leg, pch=15, col=col, pt.cex=1.5, bty="n")
    }
	# histogram with imputed missings
	HistImp <- function() {
		vars <- union(getVm("vars"), getVm("highlight"))
		vars <- c(vars,getVm("imp_vars"))
		x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
		x[,1] <- prepare(x[,1], scaling=getVm("scaling"))
		dev.new()
		histMiss(x, delimiter=getVm("delimiter"), selection=getVm("selection"), 
				col=getVm("col"), xlab=getLabel(vars[1]))
		#leg <- paste(c("observed in","missing in"), vars[2])
		#legend("topleft", legend=leg, pch=15, col=col, pt.cex=1.5, bty="n")
	}
	
    # spinogram with missings
    SpinogramMiss <- function() {
        vars <- union(getVm("vars"), getVm("highlight"))
        x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
        x[,1] <- prepare(x[,1], scaling=getVm("scaling"))
        dev.new()
        spineMiss(x, selection=getVm("selection"), 
            col=getVm("col"), xlab=getLabel(vars[1]))
    }
	# spinogram with imputed missings
	SpinogramImp <- function() {
		vars <- union(getVm("vars"), getVm("highlight"))
		vars <- c(vars,getVm("imp_vars"))
		x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
		x[,1] <- prepare(x[,1], scaling=getVm("scaling"))
		dev.new()
		spineMiss(x, delimiter=getVm("delimiter"), selection=getVm("selection"), 
				col=getVm("col"), xlab=getLabel(vars[1]))
	}
	
    # barplot with missings
    BarMiss <- function() {
        vars <- union(getVm("vars"), getVm("highlight"))
        x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
        dev.new()
        barMiss(x, selection=getVm("selection"), col=getVm("col"), xlab=vars[1])
    }
	# barplot with imputed missings
	BarImp <- function() {
		vars <- union(getVm("vars"), getVm("highlight"))
		vars <- c(vars,getVm("imp_vars"))
		x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
		dev.new()
		barMiss(x, delimiter=getVm("delimiter"), selection=getVm("selection"),
				col=getVm("col"), xlab=vars[1])
	}
	
    # spine plot with missings
    SpineplotMiss <- function() {
        vars <- union(getVm("vars"), getVm("highlight"))
        x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
        dev.new()
        spineMiss(x, selection=getVm("selection"), 
            col=getVm("col"), xlab=vars[1])
    }
	# spine plot with imputed missings
	SpineplotImp <- function() {
		vars <- union(getVm("vars"), getVm("highlight"))
		vars <- c(vars,getVm("imp_vars"))
		x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
		dev.new()
		spineMiss(x, delimiter=getVm("delimiter"), selection=getVm("selection"), 
				col=getVm("col"), xlab=vars[1])
	}
	
    # boxplot with missings
    BoxMiss <- function() {
        vars <- union(getVm("vars"), getVm("highlight"))
        x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
        x[,1] <- prepare(x[,1], scaling=getVm("scaling"))
#        if(length(vars) == 2) {
#            boxnames <- paste(c("obs. in", "miss. in"), vars[2])
#        } else boxnames <- rep("", 2)
        t <- try(testMeans(x, selection=getVm("selection")), silent=TRUE)
        if(class(t) != "try-error") {
            d <- unlist(options("digits"))
            cat(paste("\np.value: ", round(t$p.v, digits=d), "\n", sep=""))
        }
        dev.new()
#        boxplot(x[!t$ind, 1], x[t$ind, 1], names=boxnames, 
#            col=getVm("col")[1:2], ylab=getLabel(vars[1]))
        pbox(x, selection=getVm("selection"), 
            col=getVm("col")[c(1,2,4)], ylab=getLabel(vars[1]))
    }
	# boxplot with ipmuted missings
	BoxImp <- function() {
		vars <- union(getVm("vars"), getVm("highlight"))
		vars <- c(vars,getVm("imp_vars"))
		x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
		x[,1] <- prepare(x[,1], scaling=getVm("scaling"))
		t <- try(testMeans(x, selection=getVm("selection")), silent=TRUE)
		if(class(t) != "try-error") {
			d <- unlist(options("digits"))
			cat(paste("\np.value: ", round(t$p.v, digits=d), "\n", sep=""))
		}
		dev.new()
		pbox(x, delimiter=getVm("delimiter"), selection=getVm("selection"), 
				col=getVm("col")[c(1,2,4,5,6)], ylab=getLabel(vars[1]))
	}
	
    # parallel boxplots
    Pbox <- function() {
        vars <- union(getVm("vars"), getVm("highlight"))
        x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
        x[,1] <- prepare(x[,1], scaling=getVm("scaling"))
        if(getVm("tkr")) {
            TKRpbox(x, col=getVm("col")[c(1,2,4)], ylab=getLabel(vars[1]))
        } else {
            dev.new()
            pbox(x, col=getVm("col")[c(1,2,4)], ylab=getLabel(vars[1]))
        }
    }
	# parallel boxplots with imputed missings
	PboxImp <- function() {
		vars <- union(getVm("vars"), getVm("highlight"))
		vars <- c(vars,getVm("imp_vars"))
		x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
		x[,1] <- prepare(x[,1], scaling=getVm("scaling"))
		if(getVm("tkr")) {
			TKRpbox(x, delimiter=getVm("delimiter"), 
					col=getVm("col")[c(1,2,4,5,6)], ylab=getLabel(vars[1]))
		} else {
			dev.new()
			pbox(x, delimiter=getVm("delimiter"), 
					col=getVm("col")[c(1,2,4,5,6)], ylab=getLabel(vars[1]))
		}
	}
	
    # marginplot
    Marginplot <- function() {
        vars <- getVm("vars")
        x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
        x <- prepare(x, scaling=getVm("scaling"))
        labs <- getLabel(vars)
        dev.new()
        marginplot(x, col=getVm("col")[c(1,2,4)], 
            alpha=getVm("alpha"), xlab=labs[1], ylab=labs[2])
    }
	# marginplot with imputed missings
	MarginplotImp <- function() {
		vars <- c(getVm("vars"), getVm("imp_vars"))
		x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
		x <- prepare(x, scaling=getVm("scaling"))
		labs <- getLabel(vars)
		dev.new()
		marginplot(x, delimiter=getVm("delimiter"),
				col=getVm("col")[c(1,2,4,5,6)], alpha=getVm("alpha"),
				xlab=labs[1], ylab=labs[2])
	}
	
    # scatterplot with missings
    ScattMiss <- function() {
        vars <- getVm("vars")
        x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
        x <- prepare(x, scaling=getVm("scaling"))
        labs <- getLabel(vars)
        dev.new()
        scattMiss(x, col=getVm("col")[1:2], 
            alpha=getVm("alpha"), xlab=labs[1], ylab=labs[2])
    }
	# scatterplot with imputed missings
	ScattImp <- function() {
		vars <- c(getVm("vars"), getVm("imp_vars"))
		x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
		x <- prepare(x, scaling=getVm("scaling"))
		labs <- getLabel(vars)
		dev.new()
		scattMiss(x, delimiter=getVm("delimiter"), col=getVm("col")[c(1,2,5)], 
				alpha=getVm("alpha"), xlab=labs[1], ylab=labs[2])
	}
	
    # bivariate jitter plot
    ScattJitt <- function() {
        vars <- getVm("vars")
        x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
        dev.new()
        scattJitt(x, col=getVm("col")[c(1,2,4)], xlab=vars[1], ylab=vars[2])
    }
	# bivariate jitter plot with imputed missings
	ScattJittImp <- function() {
		vars <- c(getVm("vars"), getVm("imp_vars"))
		x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), select=vars)
		dev.new()
		scattJitt(x, delimiter=getVm("delimiter"),
				col=getVm("col")[c(1,2,4,5,6)], xlab=vars[1], ylab=vars[2])
	}
	
    # marginplot matrix
    Marginmatrix <- function() {
        x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), 
            select=getVm("vars"))
        x <- prepare(x, scaling=getVm("scaling"))
        if(getVm("tkr")) {
            TKRmarginmatrix(x, col=getVm("col")[c(1,2,4)], alpha=getVm("alpha"))
        } else {
            dev.new()
            marginmatrix(x, col=getVm("col")[c(1,2,4)], alpha=getVm("alpha"))
        }
    }
	# marginplot matrix with imputed missings
	MarginmatrixImp <- function() {
		vars <- c(getVm("vars"), getVm("imp_vars"))
		x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), 
				select=vars)
		x <- prepare(x, scaling=getVm("scaling"))
		if(getVm("tkr")) {
			TKRmarginmatrix(x, delimiter=getVm("delimiter"),
					col=getVm("col")[c(1,2,4,5,6)], alpha=getVm("alpha"))
		} else {
			dev.new()
			marginmatrix(x, delimiter=getVm("delimiter"),
					col=getVm("col")[c(1,2,4,5,6)], alpha=getVm("alpha"))
		}
	}
	
    # scatterplot matrix with missings
    ScattmatrixMiss <- function() {
        vars <- getVm("vars")
        highlight <- getVm("highlight")
        x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), 
            select=union(vars, highlight))
        x[,vars] <- prepare(x[,vars], scaling=getVm("scaling"))
        if(getVm("tkr")) {
            TKRscattmatrixMiss(x, highlight=highlight, 
                selection=getVm("selection"), 
                plotvars=vars, col=getVm("col")[1:2], 
                alpha=getVm("alpha"))
        } else {
            dev.new()
            scattmatrixMiss(x, highlight=highlight, 
                selection=getVm("selection"), 
                plotvars=vars, col=getVm("col")[1:2], 
                alpha=getVm("alpha"))
        }
    }
	# scatterplot matrix with imputed missings
	ScattmatrixImp <- function() {
		vars <- getVm("vars")
		highlight <- getVm("highlight")
		x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), 
				select=c(union(vars, highlight), getVm("imp_vars")))
		x[,vars] <- prepare(x[,vars], scaling=getVm("scaling"))
		if(getVm("tkr")) {
			TKRscattmatrixMiss(x, delimiter=getVm("delimiter"),
					highlight=highlight, selection=getVm("selection"), 
					plotvars=vars, col=getVm("col")[c(1,2,5)], 
					alpha=getVm("alpha"))
		} else {
			dev.new()
			scattmatrixMiss(x, delimiter=getVm("delimiter"),
					highlight=highlight, selection=getVm("selection"), 
					plotvars=vars, col=getVm("col")[c(1,2,5)], 
					alpha=getVm("alpha"))
		}
	}
	
    # parallel coordinate plot with missings
    ParcoordMiss <- function() {
        vars <- getVm("vars")
        highlight <- getVm("highlight")
        x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), 
            select=union(vars, highlight))
        x[,vars] <- prepare(x[,vars], scaling=getVm("scaling"))
        if(getVm("tkr")) {
            TKRparcoordMiss(x, highlight=highlight, 
                selection=getVm("selection"), plotvars=vars, col=getVm("col"), 
                alpha=getVm("alpha"))
        }else {
            dev.new()
            parcoordMiss(x, highlight=highlight, selection=getVm("selection"), 
                plotvars=vars, col=getVm("col"), alpha=getVm("alpha"))
        }
    }
	# parallel coordinate plot with imputed missings
	ParcoordImp <- function() {
		vars <- getVm("vars")
		highlight <- getVm("highlight")
		x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), 
				select=c(union(vars, highlight), getVm("imp_vars")))
		x[,vars] <- prepare(x[,vars], scaling=getVm("scaling"))
		if(getVm("tkr")) {
			TKRparcoordMiss(x, delimiter=getVm("delimiter"), highlight=highlight, 
					selection=getVm("selection"), plotvars=vars, col=getVm("col"), 
					alpha=getVm("alpha"))
		}else {
			dev.new()
			parcoordMiss(x, delimiter=getVm("delimiter"), highlight=highlight, selection=getVm("selection"), 
					plotvars=vars, col=getVm("col"), alpha=getVm("alpha"))
		}
	}
	
    # matrix plot
    Matrixplot <- function() {
        x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), 
            select=getVm("vars"))
        x <- prepare(x, scaling=getVm("scaling"))
        if(getVm("tkr")) TKRmatrixplot(x, col=getVm("col")[2])
        else {
            dev.new()
            matrixplot(x, col=getVm("col")[2])
        }
    }
	# matrix plot with imputed missings
	MatrixplotImp <- function() {
		vars <- c(getVm("vars"), getVm("imp_vars"))
		x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), 
				select=vars)
		x <- prepare(x, scaling=getVm("scaling"))
		if(getVm("tkr")) TKRmatrixplot(x, delimiter=getVm("delimiter"), col=getVm("col")[c(2,5)])
		else {
			dev.new()
			matrixplot(x, delimiter=getVm("delimiter"), col=getVm("col")[c(2,5)])
		}
	}
	
    # mosaic plot
    MosaicMiss <- function() {
        vars <- getVm("vars")
        highlight <- getVm("highlight")
        x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), 
            select=union(vars, highlight))
        dev.new()
        mosaicMiss(x, highlight=highlight, selection=getVm("selection"), 
            plotvars=vars, col=getVm("col")[1:2], miss.labels=FALSE)
    }
	# mosaic plot with imputed Missings
	MosaicImp <- function() {
		vars <- getVm("vars")
		highlight <- getVm("highlight")
		x <- subset(get(ActiveDataSet(), envir=.GlobalEnv), 
				select=c(union(vars, highlight), getVm("imp_vars")))
		dev.new()
		mosaicMiss(x, delimiter=getVm("delimiter"), highlight=highlight,
				selection=getVm("selection"), plotvars=vars, 
				col=getVm("col")[c(1,2,5)], miss.labels=FALSE)
	}
	
    # map of missings
    MapMiss <- function() {
        activeData <- get(ActiveDataSet(), envir=.GlobalEnv)
        x <- subset(activeData, select=getVm("vars"))
        coords <- subset(activeData, select=getVm("coords"))
#        main <- paste("Selected variables:\n", paste(colnames(x), 
#                collapse=", "))
        dev.new()
        mapMiss(x, coords, map=getVm("map"), selection=getVm("selection"), 
            col=getVm("col")[1:2], alpha=getVm("alpha"), legend=TRUE, 
            cex.main=1)
    }
	# map of imputed missings
	MapImp <- function() {
		vars <- c(getVm("vars"), getVm("imp_vars"))
		activeData <- get(ActiveDataSet(), envir=.GlobalEnv)
		x <- subset(activeData, select=vars)
		coords <- subset(activeData, select=getVm("coords"))
#        main <- paste("Selected variables:\n", paste(colnames(x), 
#                collapse=", "))
		dev.new()
		mapMiss(x, coords, map=getVm("map"), delimiter=getVm("delimiter"),
				selection=getVm("selection"), col=getVm("col")[c(1,2,5)],
				alpha=getVm("alpha"), legend=TRUE, cex.main=1)
	}
	
    # growing dot map with missings
    GrowdotMiss <- function() {
        activeData <- get(ActiveDataSet(), envir=.GlobalEnv)
        x <- subset(activeData, select=union(getVm("vars"), getVm("highlight")))
        x[,1] <- prepare(x[,1], scaling=getVm("scaling"))
        coords <- subset(activeData, select=getVm("coords"))
        alpha <- getVm("alpha")
        border <- if(alpha == 1) "white" else "transparent"
        dev.new()
        growdotMiss(x, coords, map=getVm("map"), selection=getVm("selection"), 
            col=getVm("col"), alpha=alpha, border=border, legend=TRUE)
    }
	# growing dot map with imputed missings
	GrowdotImp <- function() {
		vars <- c(getVm("vars"), getVm("imp_vars"))
		activeData <- get(ActiveDataSet(), envir=.GlobalEnv)
		x <- subset(activeData, select=union(vars, getVm("highlight")))
		x[,1] <- prepare(x[,1], scaling=getVm("scaling"))
		coords <- subset(activeData, select=getVm("coords"))
		alpha <- getVm("alpha")
		border <- if(alpha == 1) "white" else "transparent"
		dev.new()
		growdotMiss(x, coords, map=getVm("map"), delimiter=getVm("delimiter"),
				selection=getVm("selection"), col=getVm("col"), alpha=alpha,
				border=border, legend=TRUE)
	}
	
    # map with colored regions
    ColormapMiss <- function() {
        activeData <- get(ActiveDataSet(), envir=.GlobalEnv)
        x <- activeData[, getVm("vars")]
        region <- activeData[, getVm("region")]
        dev.new()
        colormapMiss(x, region, map=getVm("map"), col=getVm("col")[2])
    }
	# map with colored regions
	ColormapImp <- function() {
		activeData <- get(ActiveDataSet(), envir=.GlobalEnv)
		var <- getVm("vars")
		x <- activeData[, var]
		region <- activeData[, getVm("region")]
		# get imputation-index for the variable (if exists)
		delimiter <- getVm("delimiter")
		imp_var <- grep(delimiter, colnames(activeData), value = TRUE)
		imp_var <- imp_var[imp_var %in% paste(var,delimiter, sep="")]
		if(length(imp_var) != 0) imp_var <- activeData[, imp_var]
		else imp_var <- NULL
		dev.new()
		colormapMiss(x, region, map=getVm("map"), imp_index = imp_var,
				col=getVm("col")[c(2,5)])
	}
    
    ## functions bound to options menu
    # preferences (colors and alpha value)
    vmGUIpreferences <- function() {
        # start dialog
        ttP <- initializeDialog("Preferences")
        # combo boxes for plot colors
        colFrame <- tkwidget(ttP, "labelframe", 
            text="Select Plot Colors", fg="blue")
        colsFrame <- tkframe(colFrame)
        col <- getVm("col")
        col1Variable <- tclVar(col[1])
        col2Variable <- tclVar(col[2])
        col3Variable <- tclVar(col[3])
        col4Variable <- tclVar(col[4])
		col5Variable <- tclVar(col[5])
		col6Variable <- tclVar(col[6])
        setCol <- function() tkfocus(ttP)
        col1ComboBox <- tkwidget(colsFrame, "ComboBox", "-modifycmd", setCol, 
            values=colors(), width=12, textvariable=col1Variable)
        col2ComboBox <- tkwidget(colsFrame, "ComboBox", "-modifycmd", setCol, 
            values=colors(), width=12, textvariable=col2Variable)
        col3ComboBox <- tkwidget(colsFrame, "ComboBox", "-modifycmd", setCol, 
            values=colors(), width=12, textvariable=col3Variable)
        col4ComboBox <- tkwidget(colsFrame, "ComboBox", "-modifycmd", setCol, 
            values=colors(), width=12, textvariable=col4Variable)
		col5ComboBox <- tkwidget(colsFrame, "ComboBox", "-modifycmd", setCol, 
			values=colors(), width=12, textvariable=col5Variable)
		col6ComboBox <- tkwidget(colsFrame, "ComboBox", "-modifycmd", setCol, 
			values=colors(), width=12, textvariable=col6Variable)
        tkgrid(tklabel(colsFrame,text="Color 1: "), col1ComboBox)
        tkgrid(tklabel(colsFrame,text="Color 2: "), col2ComboBox)
        tkgrid(tklabel(colsFrame,text="Color 3: "), col3ComboBox)
        tkgrid(tklabel(colsFrame,text="Color 4: "), col4ComboBox)
		tkgrid(tklabel(colsFrame,text="Color 5: "), col5ComboBox)
		tkgrid(tklabel(colsFrame,text="Color 6: "), col6ComboBox)
        tkgrid(colsFrame, padx=3, pady=3)
        # scale for alpha value
        alphaFrame <- tkwidget(ttP, "labelframe", 
            text="Set Alpha Value", fg="blue")
        alphaVariable <- tclVar(getVm("alpha"))
        alphaScale <- tkscale(alphaFrame, orient="horizontal", from=0, to=1, 
            resolution=1/255, length=256, variable=alphaVariable)
        tkgrid(alphaScale, padx=3, pady=3)
        # misc
        miscFrame <- tkwidget(ttP, "labelframe", 
            text="Miscellaneous", fg="blue")
        tkrBox <- checkboxes(miscFrame, boxes="tkr", initial=getVm("tkr"), 
            labels="Embed multivariate plots in Tcl/Tk")
        tkgrid(tkrBox$frame, padx=3, pady=3)
		# Delimiter for imputed Variables
		delimiterFrame <- tkwidget(ttP, "labelframe",
			text = "Set imputation-delimiter", fg="blue")
		delimiter <- getVm("delimiter")
		delimiterEntry <- tkwidget(delimiterFrame, "Entry", width=12, text=delimiter, textvariable=delimiter)
		tkgrid(delimiterEntry, padx=3, pady=3)
		
        # onOK function
        onOK <- function() {
            col <- c(tclvalue(col1Variable), tclvalue(col2Variable), 
                tclvalue(col3Variable), tclvalue(col4Variable), tclvalue(col5Variable),
				tclvalue(col6Variable))
            alpha <- as.numeric(tclvalue(alphaVariable))
            tkr <- getSelection(tkrBox)
			delimiter <- tclvalue(delimiter)
            putVm("col", col)
            putVm("alpha", alpha)
            putVm("tkr", tkr)
			putVm("delimiter", delimiter)
#            activateMenus()
            closeDialog(ttP, parent=ttM)
        }
        # ok and cancel buttons
        buttons <- okCancel(ttP, onOK, parent=ttM)
        # display dialog elements
        tkgrid(colFrame, padx=10, pady=5, sticky="news")
        tkgrid(alphaFrame, padx=10, pady=5, sticky="news")
        tkgrid(miscFrame, padx=10, pady=5, sticky="news")
		tkgrid(delimiterFrame, padx=10, pady=5, sticky="news")
        tkgrid(buttons$frame)
    }
    
    ## initialize dialog and add menu
    ttM <- initializeDialog("Visualization and Imputation of Missing Values")
    tkwm.protocol(ttM, "WM_DELETE_WINDOW", Quit)
    putVm(".ttM", ttM)
    topMenu <- tkmenu(ttM)
    tkconfigure(ttM, menu=topMenu)
    DataMenu <- tkmenu(topMenu, tearoff=FALSE)
    VisualizationMenu <- tkmenu(topMenu, tearoff=FALSE)
    ImputationMenu <- tkmenu(topMenu, tearoff=FALSE)
    OptionsMenu <- tkmenu(topMenu, tearoff=FALSE)
	DiagnosticsMenu <- tkmenu(topMenu, tearoff=FALSE)
	# TO DO:
    # ----------
    #SimulationMenu <- tkmenu(topMenu, tearoff=FALSE)
    # ----------
    tkadd(topMenu, "cascade", label="Data", menu=DataMenu)
    tkadd(topMenu, "cascade", label="Visualization", menu=VisualizationMenu)
    # use following commands until other menus are implemented
    tkadd(topMenu, "cascade", label="Imputation", menu=ImputationMenu)
	#tkadd(topMenu, "command", label="Diagnostics", command=NotImplemented)
	tkadd(topMenu, "cascade", label="Diagnostics", menu=DiagnosticsMenu)
    #tkadd(topMenu, "command", label="Simulation", command=NotImplemented)
    # ----------
    tkadd(topMenu, "cascade", label="Options", menu=OptionsMenu)
    tkadd(topMenu, "command", label="Quit", command=Quit)
    #DataMenu
    tkadd(DataMenu, "command", label="Select Data", command=vmGUIdata)
    tkadd(DataMenu, "command", label="Load R Data", command=LoadRData)
    tkadd(DataMenu, "command", label="Transform Variables", 
        command=vmGUItransform)
    tkadd(DataMenu, "command", label="Background Map", command=vmGUImap)
    ########################################################################
    #VisualizationMenu
    tkadd(VisualizationMenu, "command", 
        label="Aggregate Missings", command=Aggr)
    tkadd(VisualizationMenu, "command", 
        label="Histogram with Missings", command=HistMiss)
    tkadd(VisualizationMenu, "command", 
        label="Spinogram with Missings", command=SpinogramMiss)
    tkadd(VisualizationMenu, "command", 
        label="Barplot with Missings", command=BarMiss)
    tkadd(VisualizationMenu, "command", 
        label="Spine Plot with Missings", command=SpineplotMiss)
    tkadd(VisualizationMenu, "command", 
        label="Boxplot with Missings", command=BoxMiss)
    tkadd(VisualizationMenu, "command", 
        label="Parallel Boxplots", command=Pbox)
    tkadd(VisualizationMenu, "command", 
        label="Marginplot", command=Marginplot)
    tkadd(VisualizationMenu, "command", 
        label="Scatterplot with Missings", command=ScattMiss)
    tkadd(VisualizationMenu, "command", 
        label="Bivariate Jitter Plot", command=ScattJitt)
    tkadd(VisualizationMenu, "command", 
        label="Marginplot Matrix", 
        command=Marginmatrix)
    tkadd(VisualizationMenu, "command", 
        label="Scatterplot Matrix with Missings", 
        command=ScattmatrixMiss)
    tkadd(VisualizationMenu, "command", 
        label="Parallel Coordinate Plot with Missings", 
        command=ParcoordMiss)
    tkadd(VisualizationMenu, "command", 
        label="Matrix Plot", command=Matrixplot)
    tkadd(VisualizationMenu, "command", 
        label="Mosaic Plot with Missings", command=MosaicMiss)
    tkadd(VisualizationMenu, "command", 
        label="Map of Missings", command=MapMiss)
    tkadd(VisualizationMenu, "command", 
        label="Growing Dot Map with Missings", command=GrowdotMiss)
    tkadd(VisualizationMenu, "command", 
        label="Map with Colored Regions", command=ColormapMiss)
    ########################################################################
    #ImputationMenu
    tkadd(ImputationMenu, "command", 
      label="k Nearest Neighbour", command=knnGUI)
    tkadd(ImputationMenu, "command", 
      label="Hot Deck", command=hotdeckGUI)
    tkadd(ImputationMenu, "command", 
      label="IRMI", command=irmiGUI)
	########################################################################
	#DiagnosticsMenu
	tkadd(DiagnosticsMenu, "command", 
			label="Aggregate Missings and imputed Missings", command=AggrImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Histogram with imputed Missings", command=HistImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Spinogram with imputed Missings", command=SpinogramImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Barplot with imputed Missings", command=BarImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Spine Plot with imputed Missings", command=SpineplotImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Boxplot with imputed Missings", command=BoxImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Parallel Boxplots", command=PboxImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Marginplot", command=MarginplotImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Scatterplot with imputed Missings", command=ScattImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Bivariate Jitter Plot", command=ScattJittImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Marginplot Matrix", 
			command=MarginmatrixImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Scatterplot Matrix with imputed Missings", 
			command=ScattmatrixImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Parallel Coordinate Plot with imputed Missings", 
			command=ParcoordImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Matrix Plot", command=MatrixplotImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Mosaic Plot with imputed Missings", command=MosaicImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Map of imputed Missings", command=MapImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Growing Dot Map with imputed Missings", command=GrowdotImp)
	tkadd(DiagnosticsMenu, "command", 
			label="Map with Colored Regions", command=ColormapImp)
	########################################################################
    #OptionsMenu
    tkadd(OptionsMenu, "command", label="Preferences", command=vmGUIpreferences)
    activateMenus()
    
    ## frames
    varsFrame <- tkwidget(ttM, "labelframe", 
        text="Select Variables", fg="blue")
    highlightFrame <- tkwidget(ttM, "labelframe", 
        text="Highlight Variables in Plots", fg="blue")
    scalingFrame <- tkwidget(ttM, "labelframe", 
        text="Scaling", fg="blue")
    selectionFrame <- tkwidget(ttM, "labelframe", 
        text="Selection for Highlighting", fg="blue")
    
    ## dialog elements
    # listbox for plot variables
    varsBox <- listbox(varsFrame, initial=vars, 
        height=6, selectmode="extended")
    #putVm(".varsBox", varsBox)
    # select and deselect all variables radiobuttons
    varsButtons <- radiobuttons(varsFrame, 
        buttons=c("varsSelectAll","varsDeselectAll"), 
        labels=c("Select all","Deselect all"))
    # scaling
    scalingButtons <- radiobuttons(scalingFrame, 
        buttons=c("none","classical","MCD","robust"), 
        labels=c("None","Classical",
            "Robust (MCD)","Robust (median, MAD)"), 
        initial=scaling)
    # listbox for highlight variables
    highlightBox <- listbox(highlightFrame,  
        initial=highlight, height=6, selectmode="extended")
    # select and deselect all variables radiobuttons
    highlightButtons <- radiobuttons(highlightFrame, 
        buttons=c("highlightSelectAll","highlightDeselectAll"), 
        labels=c("Select all","Deselect all"))
    # selection method for highlight variables
    selectionButtons <- radiobuttons(selectionFrame, buttons=c("any","all"), 
        labels=c("any","all"), initial=selection)
    
    ## functions bound to dialog elements
    # variables
    setVars <- function() {
        oldVars <- getVm("vars")
        selVars <- getSelection(varsBox)
        newVars <- setdiff(selVars, oldVars)
        vars <- c(intersect(oldVars, selVars), newVars)
        putVm("vars", vars)
        deselectAll(varsButtons)
        activateMenus()
        activateElements()
    }
    varsSelectAll <- function() {
        oldVars <- getVm("vars")
        allVars <- getVars()
        newVars <- setdiff(allVars, oldVars)
        vars <- c(oldVars, newVars)
        putVm("vars", vars)
        selectAll(varsBox)
        activateMenus()
        activateElements()
    }
    varsDeselectAll <- function() {
        vars <- character()
        putVm("vars", vars)
        deselectAll(varsBox)
        activateMenus()
        activateElements()
    }
    # scaling
    setScaling <- function() {
        scaling <- getSelection(scalingButtons)
        putVm("scaling", scaling)
    }
    # highlight variables
    setHighlight <- function() {
        oldHighlight <- getVm("highlight")
        selHighlight <- getSelection(highlightBox)
        newHighlight <- setdiff(selHighlight, oldHighlight)
        highlight <- c(intersect(oldHighlight, selHighlight), newHighlight)
        putVm("highlight", highlight)
        deselectAll(highlightButtons)
        activateMenus()
        activateElements()
    }
    highlightSelectAll <- function() {
        oldHighlight <- getVm("highlight")
        allVars <- getVars()
        newHighlight <- setdiff(allVars, oldHighlight)
        highlight <- c(oldHighlight, newHighlight)
        putVm("highlight", highlight)
        selectAll(highlightBox)
        activateMenus()
        activateElements()
    }
    highlightDeselectAll <- function() {
        highlight <- character()
        putVm("highlight", highlight)
        deselectAll(highlightBox)
        activateMenus()
        activateElements()
    }
    # selection method for highlight variables
    setSelection <- function() {
        selection <- getSelection(selectionButtons)
        putVm("selection", selection)
    }
    
    ## function to set states of dialog elements
    activateElements <- function() {
        boxesS <- checkActiveDataS()  # state for listboxes
        # variables
        setState(varsBox, boxesS)
        bind(varsBox, setVars)
        setState(varsButtons, boxesS)
        bind(varsButtons, varsSelectAll, "varsSelectAll")
        bind(varsButtons, varsDeselectAll, "varsDeselectAll")
        # scaling
        setState(scalingButtons, checkVarsS())
        bind(scalingButtons, setScaling)
        # highlight variables
        setState(highlightBox, boxesS)
        bind(highlightBox, setHighlight)
        setState(highlightButtons, boxesS)
        bind(highlightButtons, highlightSelectAll, "highlightSelectAll")
        bind(highlightButtons, highlightDeselectAll, "highlightDeselectAll")
        # selection method for highlight variables
        setState(selectionButtons, getSelectionS())
        bind(selectionButtons, setSelection)
    }
    
    ## display dialog elements
    activateElements()
    tkpack(varsBox$frame, varsButtons$frame, 
        expand=TRUE, fill="x", padx=3, pady=3, side="left")
    tkpack(scalingButtons$frame, expand=TRUE, 
        fill="x", padx=3, pady=3, side="left")
    tkgrid(varsFrame, scalingFrame, padx=10, pady=5, sticky="news")
    tkpack(highlightBox$frame, highlightButtons$frame, 
        expand=TRUE, fill="x", padx=3, pady=3, side="left")
    tkpack(selectionButtons$frame, expand=TRUE, 
        fill="x", padx=3, pady=3, side="left")
    tkgrid(highlightFrame, selectionFrame, padx=10, pady=5, sticky="news")
}
