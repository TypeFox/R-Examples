# Some Rcmdr dialogs for the steepness package

.onAttach <- function(libname, pkgname){
  if (!interactive()) return()
  Rcmdr <- options()$Rcmdr
  plugins <- Rcmdr$plugins
  if (!pkgname %in% plugins) {
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

if (getRversion() >= '2.15.1') globalVariables(c('checkBoxFrame2', 'notebook',
                                                 'dataTab', 'optionsTab',
                                                 'name.optionsVariable', 'namefile', 'X',
                                                 'Rcmdr.steeptest', 'methodVariable',
                                                 'doitAndPrint', 'checkboxframe2', 
                                                 'methodFrame', 'buttonsFrame', 
                                                 'DijVariable','DSVariable',
                                                 'NormDSVariable', 'ResultsVariable',
                                                 'checkBoxFrame', 'PijVariable'))

Rcmdr.steeptestDij <- function(){
  defecto <- list(Rand.inicial="10000",name.options.inicial="0",Dij.inicial="1",DS.inicial="1",
                  NormDS.inicial="1",Results.inicial="1",tab.inicial=0)
  dialog.valores <- getDialog("Rcmdr.steeptestDij",defecto) 
  initializeDialog(title=gettextRcmdr("Steepness Test"),
                   use.tabs=TRUE,tabs=c('dataTab','optionsTab'))  
onOK <- function(){
        tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
        name.labels <- as.numeric(tclvalue(name.optionsVariable))
        command <- "namefile <- tclvalue(tkgetOpenFile(filetypes='{{Text files} {.txt}} 
                {{Data files} {.dat}} {{All files} *}'))"
        justDoIt(command)
        if (namefile == "") return();
        if (name.labels == 0) {
          justDoIt("temp <- scan(namefile)")
          justDoIt("individuals <- sqrt(length(temp))")
          justDoIt("X <- matrix(temp,nrow=individuals,byrow=T)")
          }
        if (name.labels == 1) {
          justDoIt("X1 <- read.table(namefile)")
          justDoIt("names <- rownames(X1)")
          justDoIt("rownames(X1) <- NULL")
          justDoIt("colnames(X1) <- NULL")
          justDoIt("X <- as.matrix(X1)")
        }
        if (is.numeric(X) == FALSE){
           errorCondition(recall=Rcmdr.steeptestDij, message="Invalid Data Type: Original sociomatrix must be numeric.")
            return()
            }

        command <- paste("replications <- as.numeric(",tclvalue(RandVariable),")", sep="")
        justDoIt(command)
        if ( (is.na(replications)) | (replications < 1) | (replications > 1000000) ) {
            errorCondition(recall=Rcmdr.steeptestDij, message="The number of randomizations must be between 1 and 1000000.")
            return()
            }
        
        if (name.labels == 0) {doItAndPrint("test <- steeptest(X,replications,method='Dij')")}
        if (name.labels == 1) {doItAndPrint("test <- steeptest(X,replications,names,method='Dij')")}
        doItAndPrint("test$Stp")
        doItAndPrint("steep.right.pvalue <- (sum(test$Stp <= test$Stpsim)+1)/(test$rep+1)")
      	doItAndPrint("steep.right.pvalue")
	      doItAndPrint("steep.left.pvalue <- (sum(test$Stp >= test$Stpsim)+1)/(test$rep+1)")
        doItAndPrint("steep.left.pvalue")
        doItAndPrint("test$interc")   
        
        if ((as.numeric(tclvalue(DijVariable))+as.numeric(tclvalue(DSVariable))+
              as.numeric(tclvalue(NormDSVariable)))>0){
          doItAndPrint("newX <- getOrderedMatrix(X,names,method='Dij')")} 
        if (tclvalue(DijVariable) == "1") {
          doItAndPrint("getDij(newX$ordered.matrix,names=newX$ordered.names)")
            }
        if (tclvalue(DSVariable) == "1") {
          doItAndPrint("getDS(newX$ordered.matrix,names=newX$ordered.names)")
            }
        if (tclvalue(NormDSVariable) == "1") {
          doItAndPrint("getNormDS(newX$ordered.matrix,names=newX$ordered.names)")
            }
        if (tclvalue(ResultsVariable) == "1") {
          justDoIt("data <- array(dim=c(test$rep,1))
          data[,1] <- test$Stpsim
          colnames(data) <- c('Stpsim')
          Stp_rightpvalue <- (sum(test$Stp <= data[,'Stpsim'])+1)/(test$rep+1)
          Stp_leftpvalue <- (sum(test$Stp >= data[,'Stpsim'])+1)/(test$rep+1)
          results <- array((c(test$Stp, Stp_rightpvalue,Stp_leftpvalue,test$rep,mean(data[,'Stpsim']),var(data[,'Stpsim']),
                        min(data[,'Stpsim']),quantile(data[,'Stpsim'],.25,names=F),quantile(data[,'Stpsim'],.50,names=F),
                        quantile(data[,'Stpsim'],.75,names=F),max(data[,'Stpsim']))),dim=c(11,1))
          dimnames(results) <- list(c('Empirical value', 'Right p-value', 'Left p-value', 'N simulations', 'Mean',
                                'Variance','Minimum', '25th Pctl','50th Pctl', '75th Pctl','Maximum'),'Stp')
          results <-round(as.data.frame(results),round(log(results[4,],10)))")
          doItAndPrint("results")
            }
        putDialog("Rcmdr.steeptestDij",list(Rand.inicial=replications,
                                            Dij.inicial=as.numeric(tclvalue(DijVariable)),
                                            DS.inicial=as.numeric(tclvalue(DSVariable)),
                                            NormDS.inicial=as.numeric(tclvalue(NormDSVariable)),
                                            Results.inicial=as.numeric(tclvalue(ResultsVariable)),
                                            name.options.inicial=as.numeric(tclvalue(name.optionsVariable)),
                                            tab.inicial=tab))                
        
        closeDialog()
        if (name.labels == 0)    
          remove(list=c('namefile','temp','individuals','X','replications','test',
                        'steep.right.pvalue','steep.left.pvalue',
                        'Stp_rightpvalue','Stp_leftpvalue','newX','data','results'),
                 envir=.GlobalEnv)
        if (name.labels == 1)    
          remove(list=c('namefile','X1','names','X','replications','test',
                        'steep.right.pvalue','steep.left.pvalue',
                        'Stp_rightpvalue','Stp_leftpvalue','newX','data','results'),
                 envir=.GlobalEnv)  
        tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject="steeptest",reset="Rcmdr.steeptestDij",apply="Rcmdr.steeptestDij")
checkBoxes(dataTab,frame="checkBoxFrame",boxes=c("name.options"),
           initialValues=c(dialog.valores$name.options.inicial),
           labels=gettextRcmdr(c("            File Includes Row and Column Names")), 
           title = gettextRcmdr("            Original Sociomatrix will be loaded after OK"))

checkBoxes(optionsTab,frame="checkBoxFrame2",boxes=c("Dij","DS","NormDS","Results"),
           initialValues=c(dialog.valores$Dij.inicial,dialog.valores$DS.inicial,
                           dialog.valores$NormDS.inicial,dialog.valores$Results.inicial),
           labels=gettextRcmdr(c("            Dyadic Dominance Indices", "            David's Scores", 
                                 "            Normalized David's Scores", "            Summary Statistics")), 
           title = gettextRcmdr("            Results Options:"))
  RandFrame <- tkframe(optionsTab)
  RandVariable <- tclVar(dialog.valores$Rand.inicial)
  RandField <- ttkentry(RandFrame, width="12", textvariable=RandVariable)
tkgrid(labelRcmdr(RandFrame,text=gettextRcmdr("            Number of Randomizations"),font="RcmdrTitleFont"),
       RandField,sticky="w")
tkgrid(RandFrame, sticky = "w")
tkgrid.configure(RandField, sticky = "e")
tkgrid(checkBoxFrame, sticky="nw")
tkgrid(checkBoxFrame2, sticky="nw")
dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','optionsTab'),
             tab.names=c("Data","Options"))
}

Rcmdr.steeptestPij <- function(){
  defecto <- list(Rand.inicial="10000",name.options.inicial="0",Pij.inicial="1",DS.inicial="1",
                  NormDS.inicial="1",Results.inicial="1",tab.inicial=0)
  dialog.valores <- getDialog("Rcmdr.steeptestPij",defecto) 
  initializeDialog(title=gettextRcmdr("Steepness Test"),
                   use.tabs=TRUE,tabs=c('dataTab','optionsTab'))  
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
    name.labels <- as.numeric(tclvalue(name.optionsVariable))
    command <- "namefile <- tclvalue(tkgetOpenFile(filetypes='{{Text files} {.txt}} 
{{Data files} {.dat}} {{All files} *}'))"
    justDoIt(command)
    if (namefile == "") return();
    if (name.labels == 0) {
      justDoIt("temp <- scan(namefile)")
      justDoIt("individuals <- sqrt(length(temp))")
      justDoIt("X <- matrix(temp,nrow=individuals,byrow=T)")
    }
    if (name.labels == 1) {
      justDoIt("X1 <- read.table(namefile)")
      justDoIt("names <- rownames(X1)")
      justDoIt("rownames(X1) <- NULL")
      justDoIt("colnames(X1) <- NULL")
      justDoIt("X <- as.matrix(X1)")
    }
    if (is.numeric(X) == FALSE){
      errorCondition(recall=Rcmdr.steeptestPij, message="Invalid Data Type: Original sociomatrix must be numeric.")
      return()
    }
    
    command <- paste("replications <- as.numeric(",tclvalue(RandVariable),")", sep="")
    justDoIt(command)
    if ( (is.na(replications)) | (replications < 1) | (replications > 1000000) ) {
      errorCondition(recall=Rcmdr.steeptestPij, message="The number of randomizations must be between 1 and 1000000.")
      return()
    }
    
    if (name.labels == 0) {doItAndPrint("test <- steeptest(X,replications,method='Pij')")}
    if (name.labels == 1) {doItAndPrint("test <- steeptest(X,replications,names,method='Pij')")}
    doItAndPrint("test$Stp")
    doItAndPrint("steep.right.pvalue <- (sum(test$Stp <= test$Stpsim)+1)/(test$rep+1)")
    doItAndPrint("steep.right.pvalue")
    doItAndPrint("steep.left.pvalue <- (sum(test$Stp >= test$Stpsim)+1)/(test$rep+1)")
    doItAndPrint("steep.left.pvalue")
    doItAndPrint("test$interc")   
    
    if ((as.numeric(tclvalue(PijVariable))+as.numeric(tclvalue(DSVariable))+
           as.numeric(tclvalue(NormDSVariable)))>0){
      doItAndPrint("newX <- getOrderedMatrix(X,names,method='Pij')")} 
    if (tclvalue(PijVariable) == "1") {
      doItAndPrint("getPij(newX$ordered.matrix,names=newX$ordered.names)")
    }
    if (tclvalue(DSVariable) == "1") {
      doItAndPrint("getDS(newX$ordered.matrix,names=newX$ordered.names)")
    }
    if (tclvalue(NormDSVariable) == "1") {
      doItAndPrint("getNormDS(newX$ordered.matrix,names=newX$ordered.names)")
    }
    if (tclvalue(ResultsVariable) == "1") {
      justDoIt("data <- array(dim=c(test$rep,1))
               data[,1] <- test$Stpsim
               colnames(data) <- c('Stpsim')
               Stp_rightpvalue <- (sum(test$Stp <= data[,'Stpsim'])+1)/(test$rep+1)
               Stp_leftpvalue <- (sum(test$Stp >= data[,'Stpsim'])+1)/(test$rep+1)
               results <- array((c(test$Stp, Stp_rightpvalue,Stp_leftpvalue,test$rep,mean(data[,'Stpsim']),var(data[,'Stpsim']),
               min(data[,'Stpsim']),quantile(data[,'Stpsim'],.25,names=F),quantile(data[,'Stpsim'],.50,names=F),
               quantile(data[,'Stpsim'],.75,names=F),max(data[,'Stpsim']))),dim=c(11,1))
               dimnames(results) <- list(c('Empirical value', 'Right p-value', 'Left p-value', 'N simulations', 'Mean',
               'Variance','Minimum', '25th Pctl','50th Pctl', '75th Pctl','Maximum'),'Stp')
               results <-round(as.data.frame(results),round(log(results[4,],10)))")
      doItAndPrint("results")
    }
    putDialog("Rcmdr.steeptestPij",list(Rand.inicial=replications,
                                        Pij.inicial=as.numeric(tclvalue(PijVariable)),
                                        DS.inicial=as.numeric(tclvalue(DSVariable)),
                                        NormDS.inicial=as.numeric(tclvalue(NormDSVariable)),
                                        Results.inicial=as.numeric(tclvalue(ResultsVariable)),
                                        name.options.inicial=as.numeric(tclvalue(name.optionsVariable)),
                                        tab.inicial=tab))                
    
    closeDialog()
    if (name.labels == 0)    
    remove(list=c('namefile','temp','individuals','X','replications','test',
                  'steep.right.pvalue','steep.left.pvalue',
                  'Stp_rightpvalue','Stp_leftpvalue','newX','data','results'),
           envir=.GlobalEnv)
    if (name.labels == 1)    
    remove(list=c('namefile','X1','names','X','replications','test',
                  'steep.right.pvalue','steep.left.pvalue',
                  'Stp_rightpvalue','Stp_leftpvalue','newX','data','results'),
           envir=.GlobalEnv)       
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="steeptest",reset="Rcmdr.steeptestPij",apply="Rcmdr.steeptestPij")
  checkBoxes(dataTab,frame="checkBoxFrame",boxes=c("name.options"),
             initialValues=c(dialog.valores$name.options.inicial),
             labels=gettextRcmdr(c("            File Includes Row and Column Names")), 
             title = gettextRcmdr("            Original Sociomatrix will be loaded after OK"))
  
  checkBoxes(optionsTab,frame="checkBoxFrame2",boxes=c("Pij","DS","NormDS","Results"),
             initialValues=c(dialog.valores$Pij.inicial,dialog.valores$DS.inicial,
                             dialog.valores$NormDS.inicial,dialog.valores$Results.inicial),
             labels=gettextRcmdr(c("            Matrix of Pij values", "            David's Scores", 
                                   "            Normalized David's Scores", "            Summary Statistics")), 
             title = gettextRcmdr("            Results Options:"))
  RandFrame <- tkframe(optionsTab)
  RandVariable <- tclVar(dialog.valores$Rand.inicial)
  RandField <- ttkentry(RandFrame, width="12", textvariable=RandVariable)
  tkgrid(labelRcmdr(RandFrame,text=gettextRcmdr("            Number of Randomizations"),font="RcmdrTitleFont"),
         RandField,sticky="w")
  tkgrid(RandFrame, sticky = "w")
  tkgrid.configure(RandField, sticky = "e")
  tkgrid(checkBoxFrame, sticky="nw")
  tkgrid(checkBoxFrame2, sticky="nw")
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','optionsTab'),
               tab.names=c("Data","Options"))
  }


Rcmdr.steepplot <- function(){
  defecto <- list(name.options.inicial="0",method.inicial="Dij",tab.inicial=0)
  dialog.valores <- getDialog("Rcmdr.steepplot",defecto) 
  initializeDialog(title=gettextRcmdr("Steepness Plot"),
                   use.tabs=TRUE,tabs=c('dataTab','optionsTab'))  
  onOK <- function(){
    tab <- if(as.character(tkselect(notebook)) == dataTab$ID) 0 else 1
    name.labels <- as.numeric(tclvalue(name.optionsVariable))
    command <- "namefile <- tclvalue(tkgetOpenFile(filetypes='{{Text files} {.txt}} 
{{Data files} {.dat}} {{All files} *}'))"
    justDoIt(command)
    if (namefile == "") return();
    if (name.labels == 0) {
      justDoIt("temp <- scan(namefile)")
      justDoIt("individuals <- sqrt(length(temp))")
      justDoIt("X <- matrix(temp,nrow=individuals,byrow=T)")
    }
    if (name.labels == 1) {
      justDoIt("X1 <- read.table(namefile)")
      justDoIt("names <- rownames(X1)")
      justDoIt("rownames(X1) <- NULL")
      justDoIt("colnames(X1) <- NULL")
      justDoIt("X <- as.matrix(X1)")
    }
    if (is.numeric(X) == FALSE){
      errorCondition(recall=Rcmdr.steeptestPij, message="Invalid Data Type: Original sociomatrix must be numeric.")
      return()
    }
    
	method.option <- tclvalue(methodVariable)
        tkfocus(CommanderWindow())
	if (method.option == "Dij"){
          if (name.labels == 0) {doItAndPrint("STP<-steeptest(X,rep=1,method='Dij',order=TRUE)")
                                 doitAndPrint("plot(STP)")}
          if (name.labels == 1) {doItAndPrint("STP<-steeptest(X,rep=1,names,method='Dij',order=TRUE)")
                                 doItAndPrint("plot(STP)")}
        }
	if (method.option == "Pij"){
	  if (name.labels == 0) {doItAndPrint("STP<-steeptest(X,rep=1,method='Pij',order=TRUE)")
	                         doitAndPrint("plot(STP)")}
	  if (name.labels == 1) {doItAndPrint("STP<-steeptest(X,rep=1,names,method='Pij',order=TRUE)")
	                         doItAndPrint("plot(STP)")}
        }
	putDialog("Rcmdr.steepplot",list(method.inicial=method.option,
	                                    name.options.inicial=as.numeric(tclvalue(name.optionsVariable)),
	                                    tab.inicial=tab))
	if (name.labels == 0) remove(list=c('namefile','temp','individuals','X','STP'),
	       envir=.GlobalEnv)
	if (name.labels == 1) remove(list=c('namefile','X1','names','X','STP'),envir=.GlobalEnv)
  closeDialog()
	}
  OKCancelHelp(helpSubject="steepplot",reset="Rcmdr.steepplot",apply="Rcmdr.steepplot")
  checkBoxes(dataTab,frame="checkBoxFrame",boxes=c("name.options"),
             initialValues=c(dialog.valores$name.options.inicial),
             labels=gettextRcmdr(c("            File Includes Row and Column Names")), 
             title = gettextRcmdr("            Original Sociomatrix will be loaded after OK"))
  
  radioButtons(optionsTab,name = "method", buttons = c("Dij","Pij"), 
               values = c("Dij","Pij"),
               labels=gettextRcmdr(c("            Steepness plot based on Dij values",
                                     "            Steepness plot based on Pij values")),
               initialValue = dialog.valores$method.inicial,              
               title=gettextRcmdr("            Choose a method for the steepness plot"))
  tkgrid(checkBoxFrame, sticky="nw")
  tkgrid(methodFrame, sticky="w") 
  dialogSuffix(use.tabs=TRUE,grid.buttons=TRUE,tabs=c('dataTab','optionsTab'),
               tab.names=c("Data","Options"))  
}

Rcmdr.help.steepness <- function(){
   doItAndPrint("help(\"steepness\")")
   invisible(NULL)
}

Rcmdr.help.RcmdrPlugin.steepness <- function(){
   doItAndPrint("help(\"RcmdrPlugin.steepness\")")
   invisible(NULL)
}
