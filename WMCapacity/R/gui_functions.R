
.womNotebookPages <- function(name)
{
COLUMNS <- c(intro=0,data=1,models=2,analysis=3,diagnostics=4,results=5,save=6)
as.integer(COLUMNS[name])
}

womExtractModel <- function(name=1)
{
	return(wommbatAnalysis$Models[[name]])
}

.womCharVectorToAlternatingCol <- function(v)
{
	cols = c("#FFFFFF","#AAAAAA")
	r = rle(v)
	x = r$values
	r$values = (0:(length(x)-1))%%2 + 1
	v = inverse.rle(r)
	cols[v]
}
.womChooseAccCol <- function(accRate)
{
	# good, medium, bad
	cols = c("#339900","#AA7700","#CC0000")
	cols[ (accRate<0.5 | accRate>0.9) + (accRate<0.6 | accRate>0.8) + 1 ]
}

.womFileSafeString <- function(string)
{
	gsub("[^\\w-\\s]","_",string,perl=TRUE)
}


theWidget <- function(name) {
	return(StateEnv$GUI$getObject(name))
}

freezeGUI <- function(echo.to.log=T) {
	StateEnv$win$setSensitive(F)
	StateEnv$win$getWindow()$setCursor(gdkCursorNew("watch"))
	#StateEnv$echo.to.log <- echo.to.log
	#setStatusBar("")
}

thawGUI <- function() {
	StateEnv$win$setSensitive(T)
	StateEnv$win$getWindow()$setCursor(NULL)
	#StateEnv$echo.to.log <- T # default
}

getpackagefile <- function(filename) {
	## Try firstly to load from the installed wommbat package
	## Otherwise, look locally.
	myPath <- system.file("etc", filename, package = "WMCapacity")
	if (identical(myPath, "")) 
		myPath <- file.path("WMCapacity", "WMCapacity", "inst", 
			"etc", filename)
	if (!file.exists(myPath)) stop("could not find file ", filename)
	myPath
}

.womSetStatusBarText <- function(text1=NULL,context="General")
{
	status = theWidget("statusbar1")
	context.id=gtkStatusbarGetContextId(status,context)
	if(!is.null(text1)){
		gtkStatusbarPop(status, context.id)
		gtkStatusbarPush(status, context.id, text1)
	}
}

.present_main_window_after_destroy<-function(widget)
{
	StateEnv$win$present()
}

.womSaneNum <- function(x,prec=0)
{
	if(is.null(x)) return(NULL)
	as.character(round(x,prec))
}


.womSetInitialSensitivity<-function()
{
	.womActiveColumnSelection(FALSE)
	.womActiveModelTab(FALSE)
	.womActiveAnalysisTab(FALSE)
	.womActiveDiagnosticsTab(FALSE)
	.womActiveResultsTab(FALSE)
	.womActiveSaveTab(FALSE)

}

.womGetDataFrameList<-function(envir=globalenv())
{
	n =  sapply(ls(envir=envir),function(x){
			eval(parse(text=paste("is.data.frame(",x,")",sep="")),envir=globalenv())
				}
	)
	n=n[n]
	data.frame(Name=names(n))
}


.womSelectDataFrame <- function(envir = .GlobalEnv) 
{
	listOfObjects <- .womGetDataFrameList(envir=envir)
	if(dim(listOfObjects)[1]<1){
		gWidgets::gmessage(paste("There are no data frames in the R global environment."), title="Data error",
		icon="error",toolkit=gWidgets::guiToolkit("RGtk2"))
		return(0)

	}
	toplevel = gWidgets::gwindow("Select R data.frame...")
	gWidgets::gtable(listOfObjects,
		container = toplevel,
		action=NULL,
		handler = function(h,...) {
			.womSetDataForColumnSelection(eval(parse(text=as.character(gWidgets::svalue(h$obj)))))
			gWidgets::dispose(toplevel)
		},toplevel=toplevel)
}

.womGetDataFrame <- function(name)
{
	if(length(name)>1) return(NULL)
	z = eval(parse(text=name),envir=globalenv())
	if(is.data.frame(z)){
		return(z)
	}else{
		return(NULL)
	}
}

.womIntegerMatrix<-function(m){
  dims=dim(m)
  m2=as.integer(as.matrix(m))
  if(is.null(dims)){
    dim(m2)=c(length(m),1)
  }else{  
    dim(m2)=dims
  }
  return(m2)
}


.notebook_changed_page <- function(notebook,page,page_num)
{
	workingDir = getwd()
	theWidget("saveWorkingDirectoryEntry")$setText(workingDir)
}

.womTreeModelGetNthCol<-function(treemodel,n=0)
{
	totalIters = gtkTreeModelIterNChildren(treemodel, NULL)
	if(totalIters==0) return(NULL)


	retval = 1:totalIters * NA
	
	for(i in 1:totalIters)
	{
		iter = gtkTreeModelGetIterFromString(treemodel, as.character(i-1))$iter
		retval[i] = gtkTreeModelGetValue(treemodel,iter,as.integer(n))$value
	}
	
	return(retval)
}

fileChoose <- function(action="cat", text = "Select a file...", type="open", ...) 
{
	gWidgets::gfile(text=text, type=type, ..., action = action, handler =
		function(h,...) { do.call(h$action, list(h$file)) }, toolkit=gWidgets::guiToolkit("RGtk2")
		)
			
}

.clicked_menu_open <- function(menuitem)
{
	gtkNotebookSetCurrentPage(theWidget("notebook1"), .womNotebookPages("data"))
}


.clicked_menu_saveas <- function(menuitem)
{
	gtkNotebookSetCurrentPage(theWidget("notebook1"), .womNotebookPages("save"))
}

.clicked_menu_help <- function(menuitem)
{
	browseURL(url="http://wmcapacity.r-forge.r-project.org/")
}

clearComboModel <- function(combo)
{
	gtkComboBoxSetActive(combo,-1)
	model = gtkComboBoxGetModel(StateEnv$itersCombo)
	Nelements = gtkTreeModelIterNChildren(model)
	for(i in 1:Nelements)
	{
		gtkComboBoxRemoveText(combo, 0)
	}
}


# Next four functions taken from rattle
Rtxt <- function(...)
{
  # Currently, on Windows we are waiting for 2.12.17 of  RGtk2 with
  # rgtk2_bindtextdomain().

#  if (.Platform$OS.type == "windows")
#    paste(...)
#  else
    gettext(paste(...), domain="R-WMCapacity")
}

# This is used to avoid the string being identified as a translation, as in
# RtxtNT(paste(vals ...))

RtxtNT <- Rtxt


packageIsAvailable <- function(pkg, msg=NULL)
{
  if (pkg %notin% rownames(installed.packages()))
  {
    return(FALSE)
  }
  else
    return(TRUE)
}


"%notin%" <- function(x,y) ! x %in% y
