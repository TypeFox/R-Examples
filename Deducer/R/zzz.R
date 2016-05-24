
.registerDialog <- function(name,generator){
	if(bindingIsLocked(".dialogGenerators", asNamespace("Deducer")))
		unlockBinding(".dialogGenerators", asNamespace("Deducer"))
	dg <- .dialogGenerators
	dg[[name]] <- generator
	.dialogGenerators <<- dg
}

.getDialog <- function(name,newInstance=FALSE){
	if(bindingIsLocked(".dialogs", asNamespace("Deducer")))
		unlockBinding(".dialogs", asNamespace("Deducer"))
	dialog <- .dialogs[[name]]
	if(is.null(dialog) || newInstance){
		dialog <- .dialogGenerators[[name]]()
		if(!newInstance){
			di <- .dialogs
			di[[name]] <- dialog
			.dialogs <<- di
		}
	}
	dialog
}

.deducerExecute<-function(cmd){
    cmds<-parse(text=cmd)
    for(i in 1:length(cmds)){
        out<-eval(parse(text=paste("capture.output(",as.character(cmds[i]),")")),globalenv())
        for(line in out)
            cat(line,"\n")
    }
}

#.jChooserMacLAF<-function(){
#
#	if(Sys.info()[1]!="Darwin") stop("Silly rabbit, Mac Look And Feel is only for Macs!")
#	
#	.jaddClassPath("./System/Library/Java")
#	.jaddLibrary("quaqua",system.file("java/quaqua","libquaqua.jnilib",package="Deducer"))
#	.jaddLibrary("quaqua64",system.file("java/quaqua","libquaqua64.jnilib",package="Deducer"))
#	.jpackage("Deducer","quaqua.jar")
#	
#	#HashSet <- J("java.util.HashSet")
#	QuaquaManager <- J("ch.randelshofer.quaqua.QuaquaManager")
#	
#	#excludes <- new(HashSet)
#	#excludes$add("TitledBorder")
#	#excludes$add("Button")
#	#QuaquaManager$setExcludedUIs(excludes)
#	
#	#J("java.lang.System")$setProperty("Quaqua.visualMargin","0,0,0,0")
#	
#	J("javax.swing.UIManager")$setLookAndFeel(QuaquaManager$getLookAndFeel())
#}

.assign.classnames <- function(){
	DeducerMain <<- J("org.rosuda.deducer.Deducer")

	SimpleRDialog <<- J("org.rosuda.deducer.widgets.SimpleRDialog") 
	SimpleRSubDialog <<- J("org.rosuda.deducer.widgets.SimpleRSubDialog")
	RDialog <<- J("org.rosuda.deducer.widgets.RDialog") 

	VariableSelectorWidget <<- J("org.rosuda.deducer.widgets.VariableSelectorWidget")
	VariableListWidget <<- J("org.rosuda.deducer.widgets.VariableListWidget")
	CheckBoxesWidget <<- J("org.rosuda.deducer.widgets.CheckBoxesWidget")
	ButtonGroupWidget <<- J("org.rosuda.deducer.widgets.ButtonGroupWidget")
	SingleVariableWidget <<- J("org.rosuda.deducer.widgets.SingleVariableWidget")
	SliderWidget <<- J("org.rosuda.deducer.widgets.SliderWidget")
	TextAreaWidget <<- J("org.rosuda.deducer.widgets.TextAreaWidget")
	ComboBoxWidget <<- J("org.rosuda.deducer.widgets.ComboBoxWidget")
	ListWidget <<- J("org.rosuda.deducer.widgets.ListWidget")
	TextFieldWidget <<- J("org.rosuda.deducer.widgets.TextFieldWidget")
	ObjectChooserWidget <<- J("org.rosuda.deducer.widgets.ObjectChooserWidget")
	
	RDialogMonitor <<- J("org.rosuda.deducer.widgets.RDialogMonitor")
	AddRemoveButtons <<- J("org.rosuda.deducer.widgets.AddRemoveButtons")
	
	JLabel <<- J("javax.swing.JLabel")
}

DeducerMain <- NULL

SimpleRDialog <- NULL
SimpleRSubDialog <- NULL
RDialog <- NULL

VariableSelectorWidget <- NULL
VariableListWidget <- NULL
CheckBoxesWidget <- NULL
ButtonGroupWidget <- NULL
SingleVariableWidget <- NULL
SliderWidget <- NULL
TextAreaWidget <- NULL
ComboBoxWidget <- NULL
ListWidget <- NULL
TextFieldWidget <- NULL
ObjectChooserWidget <- NULL
	
RDialogMonitor <- NULL
AddRemoveButtons <- NULL
	
JLabel <- NULL

.jgr <- NULL
.deducer <- NULL
.deducer.loaded <- NULL
.dialogs <- list()
.dialogGenerators <- NULL
.windowsGUI <- NULL


.onAttach <- function(libname, pkgname){
	if(!is.null(.startupMsgs))
		packageStartupMessage(.startupMsgs)
}

.onLoad <- function(libname, pkgname) { 
	
	.imports <- parent.env(topenv())
	.i.par <- parent.env(.imports)
	
	#handle messages on .onAttach
	ipe <- new.env(parent=.i.par)
	attr(ipe, "name") <- "volatiles:Deducer"
	parent.env(.imports) <- ipe
	ipe$.startupMsgs <- NULL
	startupMessage <- function(x){
		if(is.null(ipe$.startupMsgs))
			ipe$.startupMsgs <- x
		else
			ipe$.startupMsgs <- paste(ipe$.startupMsgs,x,sep="\n")
	}
	
	
	.jgr <<- try(any(.jcall("java/lang/System","S","getProperty","main.class")=="org.rosuda.JGR.JGR"),silent=TRUE)
	if(class(.jgr) %in% "try-error"){
		startupMessage("\nNote: Problem initiating rJava. The Java GUI will not be available.\n")
		.deducer <<- .jnull()
		.deducer.loaded <<- TRUE
		.jgr <<-FALSE
		.windowsGUI <<- exists("winMenuAdd")
		return(TRUE)
	}
	

	if(!is.null(getOption("DeducerNoGUI")) && getOption("DeducerNoGUI")){
		startupMessage("\nLoading Deducer without GUI. If you wish the GUI to load, run options(DeducerNoGUI=FALSE) before loading Deducer\n")
		.deducer <<- .jnull()
		.deducer.loaded <<- TRUE
		.windowsGUI <<- exists("winMenuAdd")
		return(TRUE)
	}

	if (nzchar(Sys.getenv("NOAWT")) && .jgr!=TRUE) {
		startupMessage("\nNOTE: Environmental variable NOAWT set. Loading Deducer without GUI.\n")
		.deducer <<- .jnull()
		.deducer.loaded <<- TRUE
		.jgr <<-FALSE
		.windowsGUI <<- exists("winMenuAdd")
		return(TRUE)	
	}
	
	jriOK <- try(.jinit(),silent=TRUE)
	if(class(jriOK) %in% "try-error"){
		startupMessage("\nNote: Problem initiating JRI. The Java GUI will not be available.\n")
		#startupMessage(jriOK)
		.deducer <<- .jnull()
		.deducer.loaded <<- TRUE
		.jgr <<-FALSE
		.windowsGUI <<- exists("winMenuAdd")
		return(TRUE)
	}
	
	jriOK <- try(.jengine(TRUE),silent=TRUE)
	if(class(jriOK) %in% "try-error"){
		startupMessage("\nNote: Problem initiating JRI. Make sure you built R with shared libraries The Java GUI will not be available.\n")
		#startupMessage(jriOK)
		.deducer <<- .jnull()
		.deducer.loaded <<- TRUE
		.jgr <<-FALSE
		.windowsGUI <<- exists("winMenuAdd")
		return(TRUE)
	}
	.jpackage(pkgname,"deducer.jar",lib.loc=libname) 
	.jpackage("JavaGD")
	.jpackage("JGR")
	
	.deducer<-try(.jnew("org/rosuda/deducer/Deducer",.jgr),silent=TRUE)
	if(class(deducer) %in% "try-error"){
		startupMessage("Note: Unable to start Deducer's Java class. The Java GUI will not be available.")
		#startupMessage(deducer)
		.deducer <<- .jnull()
		.deducer.loaded <<- TRUE
		.windowsGUI <<- exists("winMenuAdd")
		return(TRUE)
	}
	.deducer <<- .deducer
	.deducer.loaded <<- TRUE
	.windowsGUI <<- FALSE
	

	
	if(exists("winMenuAdd")){
		temp<-try(utils::winMenuAdd("Deducer"),silent=TRUE)
		if(class(temp)!="try-error"){
			utils::winMenuAddItem("Deducer", "Open Data", "deducer('Open Data')")
			utils::winMenuAddItem("Deducer", "Save Data", "deducer('Save Data')")
			utils::winMenuAddItem("Deducer", "Data viewer", "deducer('Data viewer')")
			utils::winMenuAddItem("Deducer", "Deducer Online Help", "deducer('Deducer Online Help')")
			utils::winMenuAdd("Data")
			utils::winMenuAddItem("Data", "Edit Factor", "deducer('Edit Factor')")
			utils::winMenuAddItem("Data", "Recode Variables", "deducer('Recode Variables')")
			utils::winMenuAddItem("Data", "Transform", "deducer('Transform')")
			utils::winMenuAddItem("Data", "Reset Row Names", "deducer('Reset Row Names')")
			utils::winMenuAddItem("Data", "Sort", "deducer('Sort')")
			utils::winMenuAddItem("Data", "Transpose", "deducer('Transpose')")			
			utils::winMenuAddItem("Data", "Merge", "deducer('Merge')")
			utils::winMenuAddItem("Data", "Subset", "deducer('Subset')")
			utils::winMenuAdd("Analysis")
			utils::winMenuAddItem("Analysis", "Frequencies", "deducer('Frequencies')")
			utils::winMenuAddItem("Analysis", "Descriptives", "deducer('Descriptives')")
			utils::winMenuAddItem("Analysis", "Contingency Tables", "deducer('Contingency Tables')")
			utils::winMenuAddItem("Analysis", "One Sample Test", "deducer('One Sample Test')")
			utils::winMenuAddItem("Analysis", "Two Sample Test", "deducer('Two Sample Test')")
			utils::winMenuAddItem("Analysis", "K-Sample Test", "deducer('K-Sample Test')")
			utils::winMenuAddItem("Analysis", "Correlation", "deducer('Correlation')")
			utils::winMenuAddItem("Analysis", "Linear Model", "deducer('Linear Model')")
			utils::winMenuAddItem("Analysis", "Logistic Model", "deducer('Logistic Model')")
			utils::winMenuAddItem("Analysis", "Generalized Linear Model", "deducer('Generalized Linear Model')")
			utils::winMenuAdd("Plots")
			utils::winMenuAddItem("Plots", "Open plot", "deducer('Open plot')")
			utils::winMenuAddItem("Plots", "Plot builder", "deducer('Plot builder')")			
			startupMessage("\n\nDeducer has been loaded from within the Windows Rgui. 
				For the best experience, you are encouraged to use JGR console which can be
				downloaded from: http://jgr.markushelbig.org/\n")
			.windowsGUI <<- TRUE
		}else if(!.jgr || !.jcall("org/rosuda/deducer/Deducer", "Z", "isJGR") )
			startupMessage("\n\nNote Non-JGR console detected:\n\tDeducer is best used from within JGR (http://jgr.markushelbig.org/).
\tTo Bring up GUI dialogs, type deducer().\n")
	}else if(!.jgr || !.jcall("org/rosuda/deducer/Deducer", "Z", "isJGR") )
		startupMessage("\n\nNote Non-JGR console detected:\n\tDeducer is best used from within JGR (http://jgr.markushelbig.org/).
\tTo Bring up GUI dialogs, type deducer().\n")

	##sort of
	deducer.addMenu("File")
	deducer.addMenu("Data")
	deducer.addMenu("Analysis")
	deducer.addMenu("Plots")
	deducer.addMenu("Help")
	populate.items<-function(ch,cmds,men){
		for(i in 1:length(ch)){
			if(.jgr)
				cmd<-paste(".jcall('org/rosuda/deducer/Deducer',,'runCmdThreaded','" , cmds[i],"')",sep="")
			else
				cmd<-paste(".jcall('org/rosuda/deducer/Deducer',,'runCmd','" , cmds[i],"',TRUE)",sep="")
			deducer.addMenuItem(ch[i],,cmd,men)
		}
	}

	file.choices<-c("Open Data", "Save Data","Data viewer")
	file.cmds<-c("Open Data Set", "Save Data Set","table")
	populate.items(file.choices,file.cmds,"File")

	
	analysis.choices<-c("Frequencies","Descriptives",
				"Contingency Tables","One Sample Test","Paired Test","Two Sample Test",
				"K-Sample Test", "Correlation", "Linear Model", 
				"Logistic Model","Generalized Linear Model")
	analysis.cmds<-c("frequency","descriptives","contingency","onesample","pairedtest",
			"two sample","ksample","corr","linear", "logistic","glm")
	populate.items(analysis.choices,analysis.cmds,"Analysis")
	
	data.choices<-c("Edit Factor","Recode Variables","Transform",
			"Reset Row Names","Sort","Transpose","Merge",
			"Subset")
	data.cmds<-c("factor","recode","transform","reset rows","sort","trans","merge","subset")
	populate.items(data.choices,data.cmds,"Data")
	
	plot.choices<-c("Open plot","Plot builder")
	plot.cmds<-c("Open plot", "plotbuilder")
	populate.items(plot.choices,plot.cmds,"Plots")
	
	help.choices<-c("Deducer Online Help")
	help.cmds<-c("dhelp")
	populate.items(help.choices,help.cmds,"Help")
	
	.dialogGenerators <<- list()
	.dialogs <<- list()


	.assign.classnames()

	.registerDialog("Interactive Scatter Plot", .makeIplotDialog)
	.registerDialog("Interactive Histogram", .makeIhistDialog)
	.registerDialog("Interactive Bar Plot", .makeIbarDialog)
	.registerDialog("Interactive Box Plot Long", .makeIboxLongDialog)
	.registerDialog("Interactive Box Plot Wide", .makeIboxWideDialog)
	.registerDialog("Interactive Mosaic Plot", .makeImosaicDialog)
	.registerDialog("Interactive Parallel Coordinate Plot", .makeIpcpDialog)
	.registerDialog("Paired Test", .makePairedTestDialog)

#	if(J("org.rosuda.deducer.toolkit.DeducerPrefs")$USEQUAQUACHOOSER && Sys.info()[1]=="Darwin")
#		.jChooserMacLAF()

} 
.onUnload <- function(libpath){
	#de <- as.environment(match("package:Deducer", search()))
	if(exists(".deducer") && exists(".deducer.loaded") && .deducer.loaded){
		try(.jcall(.deducer,"V","detach",silent=TRUE),silent=TRUE)
		#rm(".deducer",pos=de)
		#rm(".deducer.loaded",pos=de)
	}
}

deducer<-(function(){
	lastItem<-NULL
	deducer<-function(cmd=NULL){
		if(!.jfield("org/rosuda/deducer/Deducer","Z","started")){
			.jcall(.deducer,,"startNoJGR")
		}
		oldDevice<-getOption("device")
		options(device="JavaGD")
		
		menus<- deducer.getMenus()
		menuNames <- names(menus)
		
		if(!is.null(cmd)){
			for(m in menus){
				for(item in m$items){
					if(item$name == cmd){
						if(item$silent==TRUE){
							eval(parse(text=item$command),globalenv())
						}else{
							.deducerExecute(item$command)
						}
						lastItem<<-item
						options(device=oldDevice)
						return(invisible(NULL))
					}
				}
			}
			stop("Command not present in menus")
		}
		menuChoices<-menuNames
		if(!is.null(lastItem)){
			menuChoices <- c(menuChoices,paste("Last:",lastItem$name))
		}
		menuChoices <- c(menuChoices,"exit")
		m<-menu(menuChoices,,"Deducer")
		if(m == length(menuChoices)){
			options(device=oldDevice)
			return(invisible(NULL))
		}
		if(m == length(menuChoices)-1 && !is.null(lastItem)){
			item<-lastItem
			if(item$silent==TRUE){
				eval(parse(text=item$command),globalenv())
			}else{
				.deducerExecute(item$command)
			}
			options(device=oldDevice)
			return(invisible(NULL))
		}
		men <- menus[[m]]$items
		
		itemChoices<- c(sapply(men,function(x)x$name),"back")
		
		i <- menu(itemChoices,,names(menus)[m])
		if(i == length(itemChoices) || i<=0){
			options(device=oldDevice)
			deducer()
			return(invisible(NULL))
		}
		item <- men[[i]]
		lastItem<<-item
		if(item$silent==TRUE){
			eval(parse(text=item$command),globalenv())
		}else{
			.deducerExecute(item$command)
		}
		options(device=oldDevice)
		return(invisible(NULL))
	}
	deducer
})()


data.viewer<-function(){
	deducer("Data viewer")
}









