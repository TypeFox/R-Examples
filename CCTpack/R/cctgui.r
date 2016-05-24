######################
#Begin CCTpack
######################
.onLoad <- function(libname, pkgname) {
  assign("pkg_globals", new.env(), envir=parent.env(environment()))
}

######################
#Initializations
######################

options(warn=-3)
### Some users have reported memory allocation errors when a high limit is not set here.
suppressMessages(try(memory.limit(10000),silent=TRUE))
suppressMessages(try(memory.limit(20000),silent=TRUE))
options(warn=0)

######################
#Supplementary Functions
######################

Mode <- function(x) {ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
probit <- function(x) {probval <- qnorm(x,0,1); return(probval)}
invprobit <- function(x) {invprobval <- pnorm(x,0,1); return(invprobval)}

######################
#CCTpack GUI, invoked with the command cctgui()
######################
cctgui <- function(){
  ######################
  #Sets Default GUI Variables and Frames
  ######################
  
  guidat <- list();
  guidat$tt <- tktoplevel(bg="#CDEED6")
  guidat$datframe <- tkframe(guidat$tt, borderwidth = 0,bg="#CDEED6")
  guidat$datframe2 <- tkframe(guidat$tt, borderwidth = 0,bg="#CDEED6")
  guidat$datframe3 <- tkframe(guidat$tt, borderwidth = 0,bg="#CDEED6")
  guidat$datframe4 <- tkframe(guidat$tt, borderwidth = 0,bg="#CDEED6")
  guidat$applyframe1 <- tkframe(guidat$tt, borderwidth = 0,bg="#CDEED6")
  guidat$applyframe2 <- tkframe(guidat$tt, borderwidth = 0,bg="#CDEED6")
  guidat$applyframe3 <- tkframe(guidat$tt, borderwidth = 0,bg="#CDEED6")
  guidat$resultsframe1 <- tkframe(guidat$tt, borderwidth = 0,bg="#CDEED6")
  guidat$resultsframe2 <- tkframe(guidat$tt, borderwidth = 0,bg="#CDEED6")
  guidat$resultsframe3 <- tkframe(guidat$tt, borderwidth = 0,bg="#CDEED6")
  guidat$settingsframe <- tkframe(guidat$tt, borderwidth = 0,bg="#CDEED6")
  
  guidat$polywind <- FALSE
  guidat$mval <- FALSE
  guidat$samplesvar <- tclVar("10000") 
  guidat$chainsvar <- tclVar("3")
  guidat$burninvar <- tclVar("2000")
  guidat$thinvar <- tclVar("1")
  guidat$culturesvar <- tclVar("1") 
  guidat$paravar <- tclVar("1")
  guidat$itemdiffvar <- tclVar("0") 
  guidat$polyvar <- tclVar("0")
  guidat$varnametry <- tclVar("cctfit") 
  
  tkwm.title(guidat$tt,"CCT Model Application Software")
  
  ######################
  #The GUI Grid Setup
  ######################
  guidat$samples.entry <- tkentry(guidat$settingsframe, textvariable=guidat$samplesvar,width="6")
  guidat$chains.entry <- tkentry(guidat$settingsframe, textvariable=guidat$chainsvar,width="2")
  guidat$burnin.entry <- tkentry(guidat$settingsframe, textvariable=guidat$burninvar,width="6")
  guidat$thin.entry <- tkentry(guidat$settingsframe, textvariable=guidat$thinvar,width="2")
  
         
  guidat$loaddata.but <- tkbutton(guidat$datframe, text="Load Data", command=loadfilefuncbutton)
  guidat$screeplot.but <- tkbutton(guidat$datframe, text="Scree Plot", command=screeplotfuncbutton)
  guidat$poly.but <- tkcheckbutton(guidat$datframe4,variable=guidat$polyvar,text="Use polychoric correlations",bg="#CDEED6")
  
  guidat$para.but <- tkcheckbutton(guidat$applyframe2,variable=guidat$paravar,text="Parallel run",bg="#CDEED6")
  guidat$item.but <- tkcheckbutton(guidat$applyframe2,variable=guidat$itemdiffvar,text="Item difficulty",bg="#CDEED6")
  guidat$cultdown.but <- tkbutton(guidat$applyframe2, text="<", command=cultdownfuncbutton)
  guidat$cultup.but <- tkbutton(guidat$applyframe2, text=">", command=cultupfuncbutton)
  guidat$applymodel.but <- tkbutton(guidat$applyframe2, text="Apply CCT Model", command=applymodelfuncbutton)
  guidat$cultures.entry <- tkentry(guidat$applyframe2, textvariable=guidat$culturesvar,width="2",disabledforeground="black",disabledbackground="white")
  guidat$varname.entry <- tkentry(guidat$applyframe3, textvariable=guidat$varnametry,width="10")
  
  tkbind(guidat$varname.entry, "<Key>", valid_inputkey) 
  tkbind(guidat$varname.entry, "<FocusOut>", valid_inputfo) 
  tkbind(guidat$varname.entry, "<Return>", valid_inputfo) 

  guidat$plotresults.but <- tkbutton(guidat$resultsframe1, text = "Plot Results", command = plotresultsfuncbutton)
  guidat$doppc.but <- tkbutton(guidat$resultsframe1, text = "Run Checks", command = ppcfuncbutton)
  guidat$exportresults.but <- tkbutton(guidat$resultsframe1, text = "Export Results", command = exportfuncbutton)
  guidat$memb.but <- tkbutton(guidat$resultsframe2, text = "Cluster Membs", command = membfuncbutton)
  guidat$mvest.but <- tkbutton(guidat$resultsframe2, text = "NA Value Est", command = mvestfuncbutton)
  guidat$traceplotdiscrete.but <- tkbutton(guidat$resultsframe2, text = "Traceplot Discrete", command = traceplotdiscretefuncbutton)
  guidat$traceplotall.but <- tkbutton(guidat$resultsframe2, text = "Traceplot All", command = traceplotallfuncbutton)
  guidat$printfit.but <- tkbutton(guidat$resultsframe3, text = "Print Fit", command = printfitfuncbutton)
  guidat$summary.but <- tkbutton(guidat$resultsframe3, text = "Fit Summary", command = summaryfuncbutton)
  guidat$sendconsole.but <- tkbutton(guidat$resultsframe3, text = "Send Console", command = sendconsolefuncbutton)
  
  
  guidat$datafiletxt <- tktext(guidat$tt,bg="white",width=20,height=1)
  guidat$resptxt <- tktext(guidat$tt,bg="white",width=4,height=1)
  guidat$itemtxt <- tktext(guidat$tt,bg="white",width=4,height=1) 
  guidat$dattypetxt <- tktext(guidat$tt,bg="white",width=11,height=1) 
  guidat$modeltxt <- tktext(guidat$tt,bg="white",width=4,height=1)
 
  tkgrid(guidat$datframe,columnspan=3,row=1,column=0)
  tkgrid(guidat$datframe2,columnspan=4,row=2,column=0)
  tkgrid(guidat$datframe3,columnspan=4,row=3,column=0)
  tkgrid(guidat$datframe4,columnspan=5,row=4,column=0)
  tkgrid(guidat$applyframe1,columnspan=10,row=5,column=0)
  tkgrid(guidat$applyframe2,columnspan=10,row=6,column=0)
  tkgrid(guidat$applyframe3,columnspan=10,row=7,column=0)
  tkgrid(guidat$resultsframe1,columnspan=3,row=8)
  tkgrid(guidat$resultsframe2,columnspan=4,row=9)
  tkgrid(guidat$resultsframe3,columnspan=3,row=10) 
  tkgrid(guidat$settingsframe,columnspan=8,row=11)
  
  tkgrid(tklabel(guidat$datframe,text="",bg="#CDEED6"),columnspan=3, pady = 0) 
  tkgrid(tklabel(guidat$datframe,text="Data Input",bg="#CDEED6"),columnspan=3, pady = 2) 
  
  tkgrid(guidat$loaddata.but,guidat$datafiletxt,guidat$screeplot.but,pady= 10, padx= 10)
  tkgrid(tklabel(guidat$datframe2,text="Number of Respondents",bg="#CDEED6"),guidat$resptxt,tklabel(guidat$datframe2,text="Number of Items",bg="#CDEED6"),guidat$itemtxt, padx = 2, pady = 5) 
  tkgrid(tklabel(guidat$datframe3,text="Data Type Detected",bg="#CDEED6"),guidat$dattypetxt,tklabel(guidat$datframe3,text="CCT Model",bg="#CDEED6"),guidat$modeltxt, padx = 2, pady = 5) 
  tkgrid(tklabel(guidat$applyframe1,text="Model Application",bg="#CDEED6"),columnspan=3, pady = 5) 
  tkgrid(guidat$cultdown.but, guidat$cultures.entry, guidat$cultup.but, tklabel(guidat$applyframe2,text="Cultures",bg="#CDEED6"),guidat$item.but,guidat$para.but,guidat$applymodel.but,pady= 10, padx= 2)
  tkgrid(guidat$varname.entry,tklabel(guidat$applyframe3,text="Name of Fit",bg="#CDEED6"),padx= 2)
  
  
  tkconfigure(guidat$screeplot.but, state="disabled") 
  tkgrid(tklabel(guidat$resultsframe1,text="Application Results",bg="#CDEED6"),columnspan=3, pady = 5) 
  tkgrid(guidat$doppc.but,guidat$plotresults.but,guidat$exportresults.but,pady= 5, padx= 10)
  tkgrid(guidat$memb.but,guidat$mvest.but,guidat$traceplotdiscrete.but,guidat$traceplotall.but,pady= 5, padx= 5)
  tkgrid(guidat$printfit.but,guidat$summary.but,guidat$sendconsole.but,pady= 5, padx= 10)
  
  
  tkconfigure(guidat$applymodel.but, state="disabled")
  tkconfigure(guidat$plotresults.but, state="disabled")
  tkconfigure(guidat$doppc.but, state="disabled")
  tkconfigure(guidat$exportresults.but, state="disabled")
  tkconfigure(guidat$printfit.but, state="disabled")
  tkconfigure(guidat$summary.but, state="disabled")
  tkconfigure(guidat$sendconsole.but, state="disabled")
  tkconfigure(guidat$memb.but, state="disabled")
  tkconfigure(guidat$mvest.but, state="disabled")
  tkconfigure(guidat$traceplotdiscrete.but, state="disabled")
  tkconfigure(guidat$traceplotall.but, state="disabled")
  
  tkgrid(tklabel(guidat$settingsframe,text="Sampler Settings (Optional)",bg="#CDEED6"),columnspan=8, pady = 5) 
  tkgrid(tklabel(guidat$settingsframe,text="Samples",bg="#CDEED6"), guidat$samples.entry, tklabel(guidat$settingsframe,text="Chains",bg="#CDEED6"), guidat$chains.entry, tklabel(guidat$settingsframe,text="Burn-in",bg="#CDEED6"), guidat$burnin.entry, tklabel(guidat$settingsframe,text="Thinning",bg="#CDEED6"), guidat$thin.entry, pady= 10, padx= 2)
  
  tkgrid(tklabel(guidat$settingsframe,text="",bg="#CDEED6"),columnspan=8, pady = 0) 
  
  tkconfigure(guidat$cultures.entry, bg="white", state="disabled", background="black")
  tkinsert(guidat$datafiletxt,"end","(load data file)")
  tkconfigure(guidat$datafiletxt, state="disabled")
  
  
  tkconfigure(guidat$resptxt, state="disabled")
  tkconfigure(guidat$itemtxt, state="disabled")
  tkconfigure(guidat$dattypetxt, state="disabled")
  tkconfigure(guidat$modeltxt, state="disabled")
  tkconfigure(guidat$varname.entry, state="disabled")
  
  tkconfigure(guidat$para.but, state="disabled")
  tkconfigure(guidat$item.but, state="disabled")
  tkconfigure(guidat$cultdown.but, state="disabled")
  tkconfigure(guidat$cultup.but, state="disabled")
  
  tkgrid.columnconfigure(guidat$tt,0,weight=1) 
  tkgrid.rowconfigure(guidat$tt,0,weight=1) 
  tkwm.resizable(guidat$tt,0,0)
  
  guidat$guidatnames <- names(guidat)
  list2env(guidat,pkg_globals)
  assign("guidat", guidat, pkg_globals)
  message("\n ...Starting CCT Inference Software\n")
}

######################
#Standalone Functions to use instead of the GUI
######################

######################
#All-in-one function to do each of the GUI tasks, if the gui is not compatible for the user, or not preferred
#- Takes in a data as a respondent by item array or matrix
#- 'clusters' defines # of clusters, and 'itemdiff' if item difficulty is wanted
#- 'jags' defines for the JAGS MCMC, the number of samples, chains, burn-in, thinning
#- 'runchecks' if one wants the posterior predictive checks calculated after inference
#- 'exportfilename' set a name different from "" if one wants to automatically export the results 
#    to the working directory
######################
cctapply <- function(data,clusters=1,itemdiff=FALSE,samples=10000,chains=3,burnin=2000,thinning=1,runchecks=TRUE,exportfilename="",polych=FALSE,parallel=FALSE,seed=NULL,plotr=TRUE){
  if(!is.null(seed)){
    set.seed(1); rm(list=".Random.seed", envir=globalenv())
    set.seed(seed)
  }else{set.seed(1); rm(list=".Random.seed", envir=globalenv())}
  datob <- loadfilefunc(data)
  datob <- screeplotfunc(datob,noplot=TRUE)
  cctfit <- applymodelfunc(datob,clusters=clusters,itemdiff=itemdiff,jags.iter=samples,jags.chains=chains,jags.burnin=burnin,jags.thin=thinning,parallel=parallel)
  if(plotr==TRUE){plotresultsfunc(cctfit)}
  if(runchecks==TRUE){
    cctfit <- ppcfunc(cctfit,polych=polych)
  }
  if(exportfilename!=""){
    exportfunc(cctfit,filename=exportfilename)
  }
  return(cctfit)
}

######################
#Manual function to do the scree plot, equivalent to the 'Scree Plot Button'
#- Takes in a data as a respondent by item array or matrix
######################
cctscree <- function(data,polych=FALSE){  
  datob <- loadfilefunc(data)
  datob <- screeplotfunc(datob,polych=polych)
}

######################
#Manual function to plot the results, equivalent to the 'Plot Results Button'
#- Takes the cctfit object from cctapply() or the 'Apply CCT Model' button
######################
cctresults <- function(cctfit){
  plotresultsfunc(cctfit)
}

######################
#Manual function to plot the posterior predictive check plots, equivalent to the 'Run Checks Button'
#- Takes the cctfit object from cctapply() or the 'Apply CCT Model' button
#- Plots the posterior predictive checks, or calculates them and then plots them (if they have not been calculated yet)
######################
cctppc <- function(cctfit,polych=FALSE){
  
  if(cctfit$whmodel=="LTRM" && cctfit$checksrun==TRUE){
    if(cctfit$polycor!=polych){cctfit$checksrun <- FALSE}
  }
  
  cctfit <- ppcfunc(cctfit,polych=polych)
  return(cctfit)
}

######################
#Manual function to export the results, equivalent to the 'Export Results Button'
#- Takes the cctfit object from cctapply() or the 'Apply CCT Model' button
#- Exports the cctfit object as a .Rdata file, and plots of the scree plot, 
#    posterior results, and posterior predictive checks
######################
cctexport <- function(cctfit,filename="CCTpackdata.Rdata"){
  cctfit <- exportfunc(cctfit,filename=filename)
}

######################
#End of Standalone Functions
######################

#######################
#Background Functions for the GUI and/or Standalone Functions
######################

######################
#Function for the 'Load Data' Button
#- Detects the type of data
#- Detects the number of respondents and number of items
#- Selects the appropriate model for the data
#- Detects the number of missing data points (NA values)
#- Estimates initial estimates for missing data (if any) for the initial scree plot
#- Reports this information back to the GUI or R console
#- Enables the Apply Model Button on the GUI
#- Disables buttons GUI buttons: 'Run Checks', 'Plot Results', 'Export Results' 
#    if they are active from a previous run
######################
loadfilefuncbutton <- function(){
  loadfilefunc(gui=TRUE)
}

loadfilefunc <- function(data=0,gui=FALSE,polych=FALSE){
  if(gui==TRUE){guidat <- get("guidat", pkg_globals)}
  
  datob <- list(); datob$polych <- polych; datob$datfactorsp <- 1; datob$datfactorsc
  
  if(gui==TRUE){
    datob$fileName <- file.path(tclvalue(tkgetOpenFile(filetypes = "{{csv Files} {.csv .txt}}")))
    if (!nchar(datob$fileName)) {
      return()
    }else{
      datob$dat <- as.matrix(read.csv(datob$fileName,header=FALSE))
      
      if(is.list(datob$dat)){
        datob$dat <- matrix(as.numeric(unlist(datob$dat)),dim(datob$dat)[1],dim(datob$dat)[2])
      }
      options(warn=-3)
      if(!is.numeric(datob$dat[1,1])){datob$dat <- matrix(as.numeric(datob$dat),dim(datob$dat)[1],dim(datob$dat)[2])
                                      if(all(is.na(datob$dat[,1]))){datob$dat <- datob$dat[,-1];
                                      }
                                      if(all(is.na(datob$dat[1,]))){datob$dat <- datob$dat[-1,];
                                      }
      }
      
      if(sum(is.na(datob$dat))>0){
        datob$thena <- which(is.na(datob$dat),arr.ind=TRUE)
        datob$thenalist <- which(is.na(datob$dat))
        datob$dat[datob$thena] <- datob$dat[max(which(!is.na(datob$dat)),arr.ind=TRUE)]
        datob$mval <- TRUE
      }else{datob$mval <- FALSE}
      
      options(warn=0)
      if(all(datob$dat[1:dim(datob$dat)[1],1]==1:dim(datob$dat)[1])){datob$dat <- datob$dat[,-1]; if(datob$mval==TRUE){datob$thena[,2] <- datob$thena[,2] - 1}}
      
      setwd(file.path(dirname(datob$fileName)))
      
      tkconfigure(guidat$screeplot.but, state="normal") 
      tkconfigure(guidat$applymodel.but, state="normal") 
    }
    tkconfigure(guidat$datafiletxt, state="normal") 
    tkconfigure(guidat$resptxt, state="normal")
    tkconfigure(guidat$itemtxt, state="normal") 
    tkconfigure(guidat$dattypetxt, state="normal") 
    tkconfigure(guidat$modeltxt, state="normal") 
    tkconfigure(guidat$para.but, state="normal")
    tkconfigure(guidat$item.but, state="normal")
    tkconfigure(guidat$cultdown.but, state="normal")
    tkconfigure(guidat$cultup.but, state="normal")
    
    tkdelete(guidat$datafiletxt,"1.0","800.0")
    tkdelete(guidat$resptxt,"1.0","800.0")
    tkdelete(guidat$itemtxt,"1.0","800.0")
    tkdelete(guidat$dattypetxt,"1.0","800.0")
    tkdelete(guidat$modeltxt,"1.0","800.0")
    
    tkinsert(guidat$datafiletxt,"end",basename(datob$fileName)); datob$nmiss <- 0
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    
    if(!all(is.wholenumber(datob$dat))){ 
      invlogit <- function(x){x <- 1 / (1 + exp(-x)); return(x)}
      logit <- function(x){x <- log(x/(1-x)); return(x)}
      
      tkinsert(guidat$modeltxt,"end","CRM")
      datob$whmodel <- "CRM"
      if(min(datob$dat) >= 0 && max(datob$dat) <= 1){
        datob$datatype <- "Continuous"
        message("\n ...Continuous data detected")
        tkinsert(guidat$dattypetxt,"end","Continuous in [0,1]")
        datob$dat[datob$dat<=0] <- .001; datob$dat[datob$dat>=1] <- .999
        datob$dat <- logit(datob$dat) 
        
        if(datob$mval==TRUE){
          datob$dat[datob$thena] <- NA
          datob$thenalist <- which(is.na(datob$dat)) 
          datob$dat[datob$thena] <- colMeans(datob$dat[,datob$thena[,2]],na.rm=TRUE)
        }
      }else{
        datob$datatype <- "Continuous"
        tkinsert(guidat$dattypetxt,"end","Continuous")
            datob$dat <- invlogit(datob$dat)
            datob$dat[datob$dat<=0] <- .001; datob$dat[datob$dat>=1] <- .999
            datob$dat <- logit(datob$dat) 
            
            if(datob$mval==TRUE){
              datob$dat[datob$thena] <- NA
              datob$thenalist <- which(is.na(datob$dat)) 
              datob$dat[datob$thena] <- colMeans(datob$dat[,datob$thena[,2]],na.rm=TRUE)
            }
            datob$datatype <- "Continuous"
            message("\n ...Continuous data detected")
      }     
    }else{if(min(datob$dat) >= 0 && max(datob$dat) <= 1){  
      tkinsert(guidat$modeltxt,"end","GCM") 
      datob$whmodel <- "GCM"
      datob$datatype <- "Binary"
      tkinsert(guidat$dattypetxt,"end","Binary")
      message("\n ...Binary (dichotomous) data detected")
      if(datob$mval==TRUE){
        datob$dat[datob$thena] <- NA
        datob$thenalist <- which(is.na(datob$dat)) 
        for(i in unique(datob$thena[,2])){
          datob$dat[,i][is.na(datob$dat[,i])] <- Mode(datob$dat[,i][which(!is.na(datob$dat[,i]))]) 
        }
      }
      
    }else{ 
      tkinsert(guidat$modeltxt,"end","LTRM") 
      datob$whmodel <- "LTRM"
      datob$datatype <- "Ordinal"
      tkinsert(guidat$dattypetxt,"end","Ordinal")
      message("\n ...Ordinal (categorical) data detected")
      if(guidat$polywind==FALSE){
        tkgrid(guidat$poly.but,pady= 10, padx= 2)
        guidat$polywind <- TRUE
      }
      
      
      if(datob$mval==TRUE){
        datob$dat[datob$thena] <- NA
        datob$thenalist <- which(is.na(datob$dat)) 
        for(i in unique(datob$thena[,2])){
          datob$dat[,i][is.na(datob$dat[,i])] <- Mode(datob$dat[,i][which(!is.na(datob$dat[,i]))]) 
        }
      }
    }
    
    }
    
    tkinsert(guidat$resptxt,"end",dim(datob$dat)[1])
    tkinsert(guidat$itemtxt,"end",dim(datob$dat)[2])
    
    tkconfigure(guidat$datafiletxt, state="disabled") 
    tkconfigure(guidat$resptxt, state="disabled") 
    tkconfigure(guidat$itemtxt, state="disabled") 
    tkconfigure(guidat$dattypetxt, state="disabled") 
    tkconfigure(guidat$modeltxt, state="disabled") 
    
    datob$datind <- cbind(expand.grid(t(row(datob$dat))),expand.grid(t(col(datob$dat))),expand.grid(t(datob$dat)))
    if(datob$mval==TRUE){
                      datob$nmiss <- dim(datob$thena)[1]
                      datob$datind <- datob$datind[-datob$thenalist,]
                      datob$datna <- datob$dat 
                      datob$datna[datob$thena] <- NA
                      message("\n ...Data has ",datob$nmiss," missing values out of ",length(datob$dat))
    }
    
    if(datob$whmodel=="LTRM"){   
      tkconfigure(guidat$poly.but, state="normal") 
    }else{
      if(guidat$polywind==TRUE){
        tkconfigure(guidat$poly.but, state="disabled") 
      } }  
    
    tkconfigure(guidat$varname.entry, state="normal")
    
    tkconfigure(guidat$plotresults.but, state="disabled") 
    tkconfigure(guidat$doppc.but, state="disabled") 
    tkconfigure(guidat$exportresults.but, state="disabled") 
    tkconfigure(guidat$printfit.but, state="disabled")
    tkconfigure(guidat$summary.but, state="disabled")
    tkconfigure(guidat$sendconsole.but, state="disabled")
    tkconfigure(guidat$memb.but, state="disabled")
    tkconfigure(guidat$mvest.but, state="disabled")
    tkconfigure(guidat$traceplotdiscrete.but, state="disabled")
    tkconfigure(guidat$traceplotall.but, state="disabled")
    
    datob$nobs <- dim(datob$datind)[1]
    
    datob$datobnames <- names(datob)
    #list2env(datob,pkg_globals)
    assign("datob", datob, pkg_globals)
    
    guidat$n <- dim(datob$dat)[1]; guidat$m <- dim(datob$dat)[2]; guidat$mval <- datob$mval
    assign("guidat", guidat, pkg_globals) 
    message("\n ...Data loaded")  
  }
  if(gui==FALSE){
    
    if(is.list(data)){
      datob$dat <- matrix(as.numeric(unlist(data)),dim(data)[1],dim(data)[2])
    }else{
      datob$dat <- data
    }
    
    options(warn=-3)
    if(!is.numeric(datob$dat[1,1])){datob$dat <- matrix(as.numeric(datob$dat),dim(datob$dat)[1],dim(datob$dat)[2])
                                    if(all(is.na(datob$dat[,1]))){datob$dat <- datob$dat[,-1];
                                    }
                                    if(all(is.na(datob$dat[1,]))){datob$dat <- datob$dat[-1,];
                                    }
    }
    
    datob$nmiss <- 0
    if(sum(is.na(datob$dat))>0){
      datob$thena <- which(is.na(datob$dat),arr.ind=TRUE)
      datob$thenalist <- which(is.na(datob$dat))
      datob$dat[datob$thena] <- unlist(datob$dat)[max(which(!is.na(datob$dat)),arr.ind=TRUE)]
      datob$mval <- TRUE
    }else{datob$mval <- FALSE}
    
    options(warn=0)
    if(all(datob$dat[1:dim(datob$dat)[1],1]==1:dim(datob$dat)[1])){datob$dat <- datob$dat[,-1]; if(datob$mval==TRUE){datob$thena[,2] <- datob$thena[,2] - 1}}
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    
    if(!all(is.wholenumber(datob$dat))){ 
      invlogit <- function(x){x <- 1 / (1 + exp(-x)); return(x)}
      logit <- function(x){x <- log(x/(1-x)); return(x)}
      
      datob$whmodel <- "CRM"
      if(min(datob$dat) >= 0 && max(datob$dat) <= 1){
        datob$datatype <- "Continuous"
        message("\n ...Continuous data detected")
        datob$dat[datob$dat<=0] <- .001; datob$dat[datob$dat>=1] <- .999
        #datob$dat <- logit(datob$dat) 
        datob$dat <- logit(datob$dat) 
        
        if(datob$mval==TRUE){
          datob$nmiss <- dim(datob$thena)[1]
          datob$dat[datob$thena] <- NA
          datob$thenalist <- which(is.na(datob$dat)) 
          datob$dat[datob$thena] <- colMeans(datob$dat[,datob$thena[,2]],na.rm=TRUE)
        }
      }else{
        datob$dat <- invlogit(datob$dat)
        datob$dat[datob$dat<=0] <- .001; datob$dat[datob$dat>=1] <- .999
        datob$dat <- logit(datob$dat) 
        
        if(datob$mval==TRUE){
          datob$nmiss <- dim(datob$thena)[1]
          datob$dat[datob$thena] <- NA
          datob$thenalist <- which(is.na(datob$dat)) 
          datob$dat[datob$thena] <- colMeans(datob$dat[,datob$thena[,2]],na.rm=TRUE)
        }
        message("\n ...Continuous data detected")}     
    }else{if(min(datob$dat) >= 0 && max(datob$dat) <= 1){  
      datob$whmodel <- "GCM"
      datob$datatype <- "Binary"
      message("\n ...Binary (dichotomous) data detected")
      if(datob$mval==TRUE){
        datob$nmiss <- dim(datob$thena)[1]
        datob$dat[datob$thena] <- NA
        datob$thenalist <- which(is.na(datob$dat)) 
        for(i in unique(datob$thena[,2])){
          datob$dat[,i][is.na(datob$dat[,i])] <- Mode(datob$dat[,i][which(!is.na(datob$dat[,i]))]) 
        }
      }
    }else{ 
      datob$whmodel <- "LTRM"
      datob$datatype <- "Ordinal"
      message("\n ...Ordinal (categorical) data detected")
      
      if(datob$mval==TRUE){
        datob$nmiss <- dim(datob$thena)[1]
        datob$dat[datob$thena] <- NA
        datob$thenalist <- which(is.na(datob$dat)) 
        for(i in unique(datob$thena[,2])){
          datob$dat[,i][is.na(datob$dat[,i])] <- Mode(datob$dat[,i][which(!is.na(datob$dat[,i]))]) 
        }
      }
    }
    
    }
    
    message("\n ...",dim(datob$dat)[1]," respondents and ",dim(datob$dat)[2]," items")
    
    datob$datind <- cbind(expand.grid(t(row(datob$dat))),expand.grid(t(col(datob$dat))),expand.grid(t(datob$dat)))
    if(datob$mval==TRUE){message("\n ...Data has ",dim(datob$thena)[1]," missing values out of ",length(datob$dat))
                      datob$datind <- datob$datind[-datob$thenalist,]
                      datob$datna <- datob$dat 
                      datob$datna[datob$thena] <- NA
    }
    datob$nobs <- dim(datob$datind)[1]
    return(datob)
  }
  
}

######################
#Function for the 'Scree Plot' Button
#- Uses the data object from loaddatafunc
#- Performs factor analysis on the Pearson correlations of the data 
#    using the fa() function from psych package
#- If polychoric correlations are picked for the LTRM, uses fa() on the polychoric 
#    correlations instead using polychoric() from the polycor package     
#- Creates a plot of 8 eigenvalues with appropriate titles and labels
#- When exporting the results, this function is run and produces .eps and .jpeg's of the plot
######################
screeplotfuncbutton <- function() {
  datob <- get("datob", pkg_globals)
  guidat <- get("guidat", pkg_globals)
  screeplotfunc(datob=datob,saveplots=0,savedir="",polych=as.logical(as.numeric(tclvalue(guidat$polyvar))))
}
screeplotfunc <- function(datob,saveplots=0,savedir="",gui=FALSE,polych=FALSE,noplot=FALSE) {
  options(warn=-3)
  if(datob$whmodel!="LTRM"){
    if(saveplots==0 && noplot==FALSE){
      tmp <- ""
      if(datob$mval==TRUE && datob$whmodel=="GCM"){tmp <- ", missing data handled by mode of respective columns"}
      if(datob$mval==TRUE && datob$whmodel=="CRM"){tmp <- ", missing data handled by mean of respective columns"}
      message("\n ...Producing Scree Plot",tmp)
    }
    if(saveplots==1){jpeg(file.path(gsub(".Rdata","scree.jpg",savedir)),width = 6, height = 6, units = "in", pointsize = 12,quality=100,res=400)}
    if(saveplots==2){postscript(file=file.path(gsub(".Rdata","scree.eps",savedir)), onefile=FALSE, horizontal=FALSE, width = 6, height = 6, paper="special", family="Times")}
    
    if(noplot==FALSE){
      
      if(sum(apply(datob$dat,1,function(x) sd(x,na.rm=TRUE))==0) > 0){
        tmp <- datob$dat
        if(sum(apply(datob$dat,1,function(x) sd(x,na.rm=TRUE))==0)==1){
          tmp[which(apply(tmp,1,function(x) sd(x,na.rm=TRUE))==0),1] <-  min(tmp[which(apply(tmp,1,function(x) sd(x,na.rm=TRUE))==0),1]+.01,.99)
        }else{
          tmp[which(apply(tmp,1,function(x) sd(x,na.rm=TRUE))==0),1] <- sapply(apply(tmp[which(apply(tmp,1,function(x) sd(x,na.rm=TRUE))==0),],1,mean),function(x) min(x+.01,.99))
        }
        par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(1,1))
        suppressMessages(plot(fa(cor(t(tmp)))$values[1:8],las=1,type="b",bg="black",pch=21,xlab="Factor",ylab="Magnitude",main="Scree Plot of Data"))
      }else{
        par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(1,1))
        suppressMessages(plot(fa(cor(t(datob$dat)))$values[1:8],las=1,type="b",bg="black",pch=21,xlab="Factor",ylab="Magnitude",main="Scree Plot of Data"))
      }
      
    }
    
    if(saveplots==1 || saveplots ==2){dev.off()}
  }else{
    if(saveplots==0 && noplot==FALSE){
      tmp <- ""
      if(datob$mval==TRUE && datob$whmodel=="LTRM"){tmp <- ", missing data handled by mode of respective columns"}
      message("\n ...Producing Scree Plot",tmp)
    }
    if(saveplots==1){jpeg(file.path(gsub(".Rdata","scree.jpg",savedir)),width = 6, height = 6, units = "in", pointsize = 12,quality=100,res=400)}
    if(saveplots==2){postscript(file=file.path(gsub(".Rdata","scree.eps",savedir)), onefile=FALSE, horizontal=FALSE, width = 6, height = 6, paper="special", family="Times")}
    
    if(polych==TRUE){
      if(length(datob$datfactorsp)==1){
        message("    note: utilizing polychoric correlations (time intensive)")
        datob$datfactorsp <- suppressMessages(fa(polychoric(t(datob$dat),global=FALSE)$rho)$values[1:8])
      }
      datob$datfactors <- datob$datfactorsp
    }else{
      if(length(datob$datfactorsp)==1){
        datob$datfactorsc <- suppressMessages(fa(cor(t(datob$dat)))$values[1:8])
      }
      datob$datfactors <- datob$datfactorsc
    }
    
    if(noplot==FALSE){
      par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(1,1))
      plot(datob$datfactors,las=1,type="b",bg="black",pch=21,xlab="Factor",ylab="Magnitude",main="Scree Plot of Data")
    }
    
    if(gui==TRUE){
      datob$datobnames <- names(datob)
      #list2env(datob,pkg_globals)
      assign("datob", datob, pkg_globals)
    }
    if(saveplots==1 || saveplots ==2){dev.off()}
  }
  
  if(gui==FALSE){return(datob)}
  options(warn=0)
}

######################
#Functions for the Increasing/Decreasing Cultures Buttons
######################

#For testing
#itemdifffuncbutton <- function(){print(tclvalue(get("guidat", pkg_globals)$itemdiffvar)) }

cultdownfuncbutton <- function(){guidat <- get("guidat", pkg_globals); guidat$culturesvar <- tclVar(as.character(max(1,as.numeric(tclvalue(guidat$culturesvar))-1))); tkconfigure(guidat$cultures.entry, textvariable=guidat$culturesvar); assign("guidat", guidat, pkg_globals)}
cultupfuncbutton <- function(){guidat <- get("guidat", pkg_globals); guidat$culturesvar <- tclVar(as.character(min(floor(guidat$n/4),as.numeric(tclvalue(guidat$culturesvar))+1))); tkconfigure(guidat$cultures.entry, textvariable=guidat$culturesvar); assign("guidat", guidat, pkg_globals)}

#cultdownfuncbutton <- function(){guidat <- get("guidat", pkg_globals); guidat$culturesvar <- tclVar(as.character(max(1,as.numeric(tclvalue(guidat$culturesvar))-1))); tkconfigure(guidat$cultures.entry, textvariable=guidat$culturesvar); assign("guidat", guidat, pkg_globals)}
#cultupfuncbutton <- function(){guidat <- get("guidat", pkg_globals); guidat$culturesvar <- tclVar(as.character(min(floor(dim(datob$dat)[1]/4),as.numeric(tclvalue(guidat$culturesvar))+1))); tkconfigure(guidat$cultures.entry, textvariable=guidat$culturesvar); assign("guidat", guidat, pkg_globals)}

######################
#Functions for Inputting Object Name Entry Field
######################

valid_inputkey <- function() {
  guidat <- get("guidat", pkg_globals);
  val <- gsub(" ","",tclvalue(guidat$varnametry))
  nchars <- nchar(make.names(val))
  if(val != ""){
    if(nchars > 10){
      guidat$varnametry <- tclVar(as.character(substr(make.names(val),start=1,stop=10)))
      tkconfigure(guidat$varname.entry, textvariable=guidat$varnametry)
    }else{
      guidat$varnametry <- tclVar(as.character(make.names(val)))
      tkconfigure(guidat$varname.entry, textvariable=guidat$varnametry)
    }}
  assign("guidat", guidat, pkg_globals)
}

valid_inputfo <- function() {
  guidat <- get("guidat", pkg_globals);
  val <- gsub(" ","",tclvalue(guidat$varnametry))
  nchars <- nchar(make.names(val))
  if(val == ""){
    guidat$varnametry <- tclVar("cctfit")
    tkconfigure(guidat$varname.entry, textvariable=guidat$varnametry)
    assign("guidat", guidat, pkg_globals)
    return() 
  }
  if(nchars>10){
    guidat$varnametry <- tclVar(as.character(substr(make.names(val),start=1,stop=10)))
    tkconfigure(guidat$varname.entry, textvariable=guidat$varnametry)
  }else{
    guidat$varnametry <- tclVar(as.character(make.names(val)))
    tkconfigure(guidat$varname.entry, textvariable=guidat$varnametry)
  }
  assign("guidat", guidat, pkg_globals)
}

######################
#Function for the 'Apply CCT Model' Button
#- Uses the data object from loaddatafunc
#- Applies hierarchical Bayesian inference for the model using JAGS by packages rjags and R2jags
#- Model code for each of the three models are at the bottom of this rcode file
#- Reads in user preferences for model specifications from the GUI
#- Uses these specifications to apply different form(s) of the model
#- During mixture model cases, applies an algorithm that corrects for label-switching/mixing issues
#- Recalculates the statistics for all parameters after correcting for label-switching, as well as
#    the Rhat statistics, and DIC
#- Enables the user to look at traceplots for discrete nodes via dtraceplots() 
#- Provides the model-based clustering by cctfit$respmem   (respondent membership)
#- Provides cctfit$Lik, which is the likelihood of the model evaluated at each sample
#- Enables GUI buttons: 'Run Checks', 'Plot Results', 'Export Results'
######################
applymodelfuncbutton <- function() {
  guidat <- get("guidat", pkg_globals)
  datob <- get("datob", pkg_globals)
  
  guidat$varname <- tclvalue(guidat$varnametry)
  tkconfigure(guidat$varname.entry, state="disabled")
  assign("guidat", guidat, pkg_globals)
  
  assign(guidat$varname,
  applymodelfunc(datob=datob,clusters=as.numeric(tclvalue(guidat$culturesvar)),itemdiff=as.logical(as.numeric(tclvalue(guidat$itemdiffvar))),
                  jags.iter=as.numeric(tclvalue(guidat$samplesvar)),jags.chains= as.numeric(tclvalue(guidat$chainsvar)),jags.burnin=as.numeric(tclvalue(guidat$burninvar)),
                  jags.thin=as.numeric(tclvalue(guidat$thinvar)),parallel= as.logical(as.numeric(tclvalue(guidat$paravar))),gui=TRUE),
  pkg_globals
  )

  #User-friendly (for novices) but not CRAN compliant because it assigns the fit object in the Global environment
  #assign(guidat$varname,
  #applymodelfunc(datob=datob,clusters=as.numeric(tclvalue(guidat$culturesvar)),itemdiff=as.logical(as.numeric(tclvalue(guidat$itemdiffvar))),
  #                jags.iter=as.numeric(tclvalue(guidat$samplesvar)),jags.chains= as.numeric(tclvalue(guidat$chainsvar)),jags.burnin=as.numeric(tclvalue(guidat$burninvar)),
  #                jags.thin=as.numeric(tclvalue(guidat$thinvar)),parallel= as.logical(as.numeric(tclvalue(guidat$paravar))),gui=TRUE),
  #inherits=TRUE
  #)
}

applymodelfunc <- function(datob,clusters=1,itemdiff=FALSE,jags.iter=10000,jags.chains=3,jags.burnin=2000,jags.thin=1,parallel=FALSE,gui=FALSE) {
  
  if(gui==TRUE){guidat <- get("guidat", pkg_globals)}
  ######################
  #Model Codes
  #- Used by the 'Apply CCT Model' button/function
  #- Each model has 2 code versions, 1 without item difficulty, 1 with item difficulty
  ######################
  mcgcm <-
    "model{
  for (l in 1:nobs){
  D[Y[l,1],Y[l,2]] <- (th[Y[l,1]]*(1-lam[Y[l,2]])) / ((th[Y[l,1]]*(1-lam[Y[l,2]]))+(lam[Y[l,2]]*(1-th[Y[l,1]]))) 
  pY[Y[l,1],Y[l,2]] <- (D[Y[l,1],Y[l,2]]*Z[Y[l,2],Om[Y[l,1]]]) +((1-D[Y[l,1],Y[l,2]])*g[Y[l,1]])
  Y[l,3] ~ dbern(pY[Y[l,1],Y[l,2]]) }
  
  for (i in 1:nresp){
  Om[i] ~ dcat(pi) 
  th[i] ~ dbeta(thmu[Om[i]]*thtau[Om[i]],(1-thmu[Om[i]])*thtau[Om[i]])
  g[i] ~ dbeta(gmu[Om[i]]*gtau[Om[i]],(1-gmu[Om[i]])*gtau[Om[i]]) }
  
  for (k in 1:nitem){
  lam[k] <- .5
  for (v in 1:V){
  Z[k,v] ~ dbern(p[v]) }}
  
  #Hyper Parameters
  gsmu <- 10
  gssig <- 10
  dsmu <- 10
  dssig <- 10
  pi[1:V] ~ ddirch(L)
  alpha <- 2
  
  for (v in 1:V){
  p[v] ~ dunif(0,1)
  gmu[v] <- .5
  gtau[v] ~ dgamma(pow(gsmu,2)/pow(gssig,2),gsmu/pow(gssig,2))
  thmu[v] ~ dbeta(alpha,alpha)
  thtau[v] ~ dgamma(pow(dsmu,2)/pow(dssig,2),dsmu/pow(dssig,2))
  L[v] <- 1 }
}"

  mcgcmid <-
    "model{
  for (l in 1:nobs){
  D[Y[l,1],Y[l,2]] <- (th[Y[l,1]]*(1-lam[Y[l,2],Om[Y[l,1]]])) / ((th[Y[l,1]]*(1-lam[Y[l,2],Om[Y[l,1]]]))+(lam[Y[l,2],Om[Y[l,1]]]*(1-th[Y[l,1]])))   
  pY[Y[l,1],Y[l,2]] <- (D[Y[l,1],Y[l,2]]*Z[Y[l,2],Om[Y[l,1]]]) +((1-D[Y[l,1],Y[l,2]])*g[Y[l,1]])
  Y[l,3] ~ dbern(pY[Y[l,1],Y[l,2]]) }  
  
  for (i in 1:nresp){
  Om[i] ~ dcat(pi) 
  th[i] ~ dbeta(thmu[Om[i]]*thtau[Om[i]],(1-thmu[Om[i]])*thtau[Om[i]])
  g[i] ~ dbeta(gmu[Om[i]]*gtau[Om[i]],(1-gmu[Om[i]])*gtau[Om[i]]) }
  
  for (k in 1:nitem){
  for (v in 1:V){
  Z[k,v] ~ dbern(p[v])
  lam[k,v] ~ dbeta(lammu[v]*lamtau[v],(1-lammu[v])*lamtau[v])
  }}
  
  #Hyper Parameters
  lamsmu <- 10
  lamssig <- 10
  gsmu <- 10
  gssig <- 10
  dsmu <- 10
  dssig <- 10
  alpha <- 2
  pi[1:V] ~ ddirch(L)
  
  for (v in 1:V){
  p[v] ~ dunif(0,1)
  lammu[v] <- .5
  lamtau[v] ~ dgamma(pow(lamsmu,2)/pow(lamssig,2),lamsmu/pow(lamssig,2))
  gmu[v] <- .5
  gtau[v] ~ dgamma(pow(gsmu,2)/pow(gssig,2),gsmu/pow(gssig,2))
  thmu[v] ~ dbeta(alpha,alpha)
  thtau[v] ~ dgamma(pow(dsmu,2)/pow(dssig,2),dsmu/pow(dssig,2))
  L[v] <- 1
  }
  }"

  mcltrm <-
    "model{
  for (l in 1:nobs){
  tau[Y[l,1],Y[l,2]] <- pow(E[Y[l,1]],-2)
  pY[Y[l,1],Y[l,2],1] <- pnorm((a[Y[l,1]]*gam[1,Om[Y[l,1]]]) + b[Y[l,1]],T[Y[l,2],Om[Y[l,1]]],tau[Y[l,1],Y[l,2]])
  for (c in 2:(C-1)){pY[Y[l,1],Y[l,2],c] <- pnorm((a[Y[l,1]]*gam[c,Om[Y[l,1]]]) + b[Y[l,1]],T[Y[l,2],Om[Y[l,1]]],tau[Y[l,1],Y[l,2]]) - sum(pY[Y[l,1],Y[l,2],1:(c-1)])}
  pY[Y[l,1],Y[l,2],C] <- (1 - sum(pY[Y[l,1],Y[l,2],1:(C-1)]))
  Y[l,3] ~ dcat(pY[Y[l,1],Y[l,2],1:C])
  }
  
  #Parameters
  for (i in 1:nresp){
  Om[i] ~ dcat(pi) 
  Elog[i] ~ dnorm(Emu[Om[i]],Etau[Om[i]])
  E[i] <- exp(Elog[i])  
  alog[i] ~ dnorm(amu[Om[i]],atau[Om[i]])T(-2.3,2.3)
  a[i] <- exp(alog[i])
  b[i] ~ dnorm(bmu[Om[i]],btau[Om[i]]) }
  
  for (k in 1:nitem){
  for (v in 1:V){
  T[k,v] ~ dnorm(Tmu[v],Ttau[v]) }}
  
  for (v in 1:V){
  gam[1:(C-1),v] <- sort(tgam2[1:(C-1),v])
  for (c in 1:(C-2)){tgam[c,v] ~ dnorm(0,.1)}
  tgam2[1:(C-2),v] <- tgam[1:(C-2),v]
  tgam2[C-1,v] <- -sum(tgam[1:(C-2),v]) }
  
  pi[1:V] ~ ddirch(L)
  
  #Hyperparameters	
  for (v in 1:V){
  L[v] <- 1 
  Tmu[v] ~ dnorm(0,.25)
  Tsig[v] ~ dunif(.25,3)
  Ttau[v] <- pow(Tsig[v],-2)
  Emu[v] ~ dnorm(0,.01)
  Etau[v] ~ dgamma(.01, .01)
  amu[v] <- 0
  atau[v] ~ dgamma(.01, .01)T(.01,)
  bmu[v] <- 0
  btau[v] ~ dgamma(.01, .01)
  }
  }"

  mcltrmid <-
    "model{
  for (l in 1:nobs){
  tau[Y[l,1],Y[l,2]] <- pow(E[Y[l,1]]*exp(lam[Y[l,2],Om[Y[l,1]]]),-2)
  pY[Y[l,1],Y[l,2],1] <- pnorm((a[Y[l,1]]*gam[1,Om[Y[l,1]]]) + b[Y[l,1]],T[Y[l,2],Om[Y[l,1]]],tau[Y[l,1],Y[l,2]])
  for (c in 2:(C-1)){pY[Y[l,1],Y[l,2],c] <- pnorm((a[Y[l,1]]*gam[c,Om[Y[l,1]]]) + b[Y[l,1]],T[Y[l,2],Om[Y[l,1]]],tau[Y[l,1],Y[l,2]]) - sum(pY[Y[l,1],Y[l,2],1:(c-1)])}
  pY[Y[l,1],Y[l,2],C] <- (1 - sum(pY[Y[l,1],Y[l,2],1:(C-1)]))
  Y[l,3] ~ dcat(pY[Y[l,1],Y[l,2],1:C])
  }
  
  #Parameters
  for (i in 1:nresp){
  Om[i] ~ dcat(pi) 
  Elog[i] ~ dnorm(Emu[Om[i]],Etau[Om[i]])
  E[i] <- exp(Elog[i])  
  alog[i] ~ dnorm(amu[Om[i]],atau[Om[i]])T(-2.3,2.3)
  a[i] <- exp(alog[i])
  b[i] ~ dnorm(bmu[Om[i]],btau[Om[i]]) 
  }
  
  for (k in 1:nitem){       
  for (v in 1:V){
  T[k,v] ~ dnorm(Tmu[v],Ttau[v])
  lam[k,v] ~ dnorm(lammu[v],lamtau[v])T(-2.3,2.3)   
  } }
  
  for (v in 1:V){
  gam[1:(C-1),v] <- sort(tgam2[1:(C-1),v])
  for (c in 1:(C-2)){tgam[c,v] ~ dnorm(0,.1)}
  tgam2[1:(C-2),v] <- tgam[1:(C-2),v]
  tgam2[C-1,v] <- -sum(tgam[1:(C-2),v]) 
  }
  
  pi[1:V] ~ ddirch(L)
  
  #Hyperparameters	
  for (v in 1:V){
  L[v] <- 1 
  Tmu[v] ~ dnorm(0,.25)
  Tsig[v] ~ dunif(.25,3)
  Ttau[v] <- pow(Tsig[v],-2)
  Emu[v] ~ dnorm(0,.01)
  Etau[v] ~ dgamma(.01, .01)
  amu[v] <- 0
  atau[v] ~ dgamma(.01, .01)T(.01,)
  bmu[v] <- 0
  btau[v] ~ dgamma(.01, .01)
  lammu[v] <- 0
  lamsig[v] ~ dunif(.25, 2)
  lamtau[v] <- pow(lamsig[v],-2)
  }
  }"

  mccrm <-
    "model{
  for (l in 1:nobs){ 
  Y[l,3] ~ dnorm((a[Y[l,1]]*T[Y[l,2],Om[Y[l,1]]])+b[Y[l,1]],pow(a[Y[l,1]]*E[Y[l,1]],-2))
  }	  
  
  #Parameters
  for (i in 1:nresp){
  Om[i] ~ dcat(pi)
  Elog[i] ~ dnorm(Emu[Om[i]],Etau[Om[i]])
  E[i] <- exp(Elog[i])  
  alog[i] ~ dnorm(amu[Om[i]],atau[Om[i]])T(-2.3,2.3)
  a[i] <- exp(alog[i])
  b[i] ~ dnorm(bmu[Om[i]],btau[Om[i]]) }
  
  for (k in 1:nitem){ 
  for (v in 1:V){
  T[k,v] ~ dnorm(Tmu[v],Ttau[v])
  }}
  
  pi[1:V] ~ ddirch(L)
  
  #Hyperparameters	
  
  for (v in 1:V){
  L[v] <- 1 
  Tmu[v] ~ dnorm(0,0.25)
  Tsig[v] ~ dunif(.25,3)
  Ttau[v] <- pow(Tsig[v],-2)
  Emu[v] ~ dnorm(0,.01)
  Etau[v] ~ dgamma(.01, .01)
  amu[v] <- 0
  atau[v] ~ dgamma(.01, .01)T(.01,)
  bmu[v] <- 0
  btau[v] ~ dgamma(.01, .01)
  }
  }"
 
  mccrmid <-
    "model{
  for (l in 1:nobs){ 
  Y[l,3] ~ dnorm((a[Y[l,1]]*T[Y[l,2],Om[Y[l,1]]])+b[Y[l,1]],pow(a[Y[l,1]]*E[Y[l,1]]*exp(lam[Y[l,2],Om[Y[l,1]]]),-2))
  }    
  
  #Parameters
  for (i in 1:nresp){
  Om[i] ~ dcat(pi) 
  Elog[i] ~ dnorm(Emu[Om[i]],Etau[Om[i]])
  E[i] <- exp(Elog[i])
  alog[i] ~ dnorm(amu[Om[i]],atau[Om[i]])T(-2.3,2.3)
  a[i] <- exp(alog[i])
  b[i] ~ dnorm(bmu[Om[i]],btau[Om[i]]) }
  
  for (k in 1:nitem){ 
  for (v in 1:V){
  T[k,v] ~ dnorm(Tmu[v],Ttau[v])
  lam[k,v] ~ dnorm(lammu[v],lamtau[v])T(-2.3,2.3)
  }}
  
  pi[1:V] ~ ddirch(L)
  
  #Hyperparameters	
  for (v in 1:V){
  L[v] <- 1 
  Tmu[v] ~ dnorm(0,.25)
  Tsig[v] ~ dunif(.25,3)
  Ttau[v] <- pow(Tsig[v],-2)
  Emu[v] ~ dnorm(0,.01)
  Etau[v] ~ dgamma(.01, .01)
  amu[v] <- 0
  atau[v] ~ dgamma(.01, .01)T(.01,)
  bmu[v] <- 0
  btau[v] ~ dgamma(.01, .01)
  lammu[v] <- 0
  lamsig[v] ~ dunif(.25, 2)
  lamtau[v] <- pow(lamsig[v],-2)
  }
  }"

  ## Diffuse Settings
  #  L[v] <- 1 
  #  Tmu[v] ~ dnorm(0,.25)
  #  Ttau[v] ~ dgamma(.01, .01)
  #  Emu[v] ~ dnorm(0,.01)
  #  Etau[v] ~ dgamma(.01, .01)
  #  amu[v] <- 0
  #  atau[v] ~ dgamma(.01, .01)T(.01,)
  #  bmu[v] <- 0
  #  btau[v] ~ dgamma(.01, .01)
  #  lammu[v] <- 0
  #  lamtau[v] ~ dgamma(.01, .01)T(.01,)
  
  ######################
  #Sets up Parameters, Variables, and Data for JAGS for the model selected
  ######################
  if(datob$whmodel=="GCM"){
    Y <- datob$datind; nresp <- dim(datob$dat)[1]; nitem <- dim(datob$dat)[2]; V <- clusters; nobs <- dim(datob$datind)[1]
    jags.data <- list("Y","nresp","nitem","V","nobs")
    
    if(itemdiff==FALSE){
      model.file <- mcgcm
      jags.params <- c("Z","th","g","p","thmu","thtau","gmu","gtau","Om","pi")
      if(clusters>1){
        jags.inits <- function(){ list("Z"=matrix(rbinom(nitem*V,1,.5),nitem,V),"th"= runif(nresp,.2,.8), "g"= runif(nresp,.2,.8),"Om"= sample(1:V,nresp,replace=TRUE) )}
      }else{
        jags.inits <- function(){ list("Z"=matrix(rbinom(nitem*V,1,.5),nitem,V),"th"= runif(nresp,.2,.8), "g"= runif(nresp,.2,.8) )}
      }
      
    }
    if(itemdiff==TRUE){
      model.file <- mcgcmid
      jags.params <- c("Z","th","g","lam","p","thmu","thtau","gmu","gtau","lammu","lamtau","Om","pi")
      if(clusters>1){
        jags.inits <- function(){ list("Z"=matrix(rbinom(nitem*V,1,.5),nitem,V),"th"= runif(nresp,.2,.8), "g"= runif(nresp,.2,.8), "lam"= matrix(runif(nitem*V,.2,.8),nitem,V), "Om"= sample(1:V,nresp,replace=TRUE) )}
      }else{
        jags.inits <- function(){ list("Z"=matrix(rbinom(nitem*V,1,.5),nitem,V),"th"= runif(nresp,.2,.8), "g"= runif(nresp,.2,.8), "lam"= matrix(runif(nitem*V,.2,.8),nitem,V) )}
      }
    }
    if(clusters==1){
      model.file <- gsub(pattern="pi\\[1\\:V\\] ~ ddirch\\(L\\)", "pi <- 1", model.file)
      model.file <- gsub(pattern="Om\\[i\\] ~ dcat\\(pi\\)", "Om\\[i\\] <- 1", model.file)
    }
  }
  
  if(datob$whmodel=="LTRM"){
    Y <- datob$datind; nresp <- dim(datob$dat)[1]; nitem <- dim(datob$dat)[2]; V <- clusters; C <- max(datob$datind[,3]); nobs <- dim(datob$datind)[1]
    jags.data <- list("Y","nresp","nitem","C","V","nobs")
    
    if(itemdiff==FALSE){
      model.file <- mcltrm
      jags.params <- c("T","gam","E","a","b","Tmu","Ttau","Emu","Etau","amu","atau","bmu","btau","Om","pi")
      if(clusters>1){
        jags.inits <- function(){ list("tgam"=matrix(rnorm((C-2)*V,0,1),(C-2),V),"T"=matrix(rnorm(nitem*V,0,1),nitem,V),"Emu"= runif(V,.8,2),"Esig"= runif(V,.4,.6), "asig"= runif(V,.4,.6),"bsig"= runif(V,.4,.6),"Om"= sample(1:V,nresp,replace=TRUE)   )}
      }else{
        jags.inits <- function(){ list("tgam"=matrix(rnorm((C-2)*V,0,1),(C-2),V),"T"=matrix(rnorm(nitem*V,0,1),nitem,V),"Emu"= runif(V,.8,2),"Esig"= runif(V,.4,.6), "asig"= runif(V,.4,.6),"bsig"= runif(V,.4,.6)  )}
      }
      
      if(C==2){
        if(clusters>1){
          jags.inits <- function(){ list("T"=matrix(rnorm(nitem*V,0,1),nitem,V),"Emu"= runif(V,.8,2),"Esig"= runif(V,.4,.6),"bsig"= runif(V,.4,.6),"Om"= sample(1:V,nresp,replace=TRUE)   )}
        }else{
          jags.inits <- function(){ list("T"=matrix(rnorm(nitem*V,0,1),nitem,V),"Emu"= runif(V,.8,2),"Esig"= runif(V,.4,.6),"bsig"= runif(V,.4,.6)   )}
        }
        model.file <- gsub(pattern="a\\[i\\] <- exp\\(alog\\[i\\]\\)", "a\\[i\\] <- 1", model.file)
        model.file <- gsub(pattern="alog\\[i\\] ~ dnorm\\(amu\\[Om\\[i\\]\\],atau\\[Om\\[i\\]\\]\\)T\\(-2.3,2.3\\)", "", model.file)
        model.file <- gsub(pattern="atau\\[v\\] ~ dgamma\\(.01,.01\\)T\\(.01,\\)", "atau\\[v\\] <- 0", model.file)
        model.file <- gsub(pattern="for \\(c in 2\\:\\(C-1\\)\\)", " #for \\(c in 2\\:\\(C-1\\)\\)", model.file)
        model.file <- gsub(pattern="gam\\[1\\:\\(C-1\\),v\\] <- sort\\(tgam2\\[1\\:\\(C-1\\),v\\]\\)", "gam\\[1,v\\] <- 0 }", model.file)
        model.file <- gsub(pattern="for \\(c in 1\\:\\(C-2\\)\\)", "#for \\(c in 1\\:\\(C-2\\)\\)", model.file)
        model.file <- gsub(pattern="\\(1 - sum\\(pY\\[i,k,1\\:\\(C-1\\)\\]\\)\\)", "1 - pY\\[i,k,1\\]", model.file)
        model.file <- gsub(pattern="tgam2\\[1", "#tgam2\\[1", model.file)
        model.file <- gsub(pattern="tgam2\\[C", "#tgam2\\[C", model.file)
      }
    }
    
    if(itemdiff==TRUE){
      model.file <- mcltrmid
      jags.params <- c("T","lam","gam","E","a","b","Tmu","Ttau","Emu","Etau","amu","atau","bmu","btau","lammu","lamtau","Om","pi")
      if(clusters>1){
        jags.inits <- function(){ list("tgam"=matrix(rnorm((C-2)*V,0,1),(C-2),V),"T"=matrix(rnorm(nitem*V,0,1),nitem,V),"Emu"= runif(V,.8,2),"Esig"= runif(V,.4,.6), "asig"= runif(V,.4,.6),"bsig"= runif(V,.4,.6),"lamsig"= runif(V,.25,.30), "Om"= sample(1:V,nresp,replace=TRUE)   )}
      }else{
        jags.inits <- function(){ list("tgam"=matrix(rnorm((C-2)*V,0,1),(C-2),V),"T"=matrix(rnorm(nitem*V,0,1),nitem,V),"Emu"= runif(V,.8,2),"Esig"= runif(V,.4,.6), "asig"= runif(V,.4,.6),"bsig"= runif(V,.4,.6),"lamsig"= runif(V,.25,.30)  )}
      }
      if(C==2){
        if(clusters>1){
          jags.inits <- function(){ list("T"=matrix(rnorm(nitem*V,0,1),nitem,V),"Emu"= runif(V,.8,2),"Esig"= runif(V,.4,.6),"bsig"= runif(V,.4,.6),"lamsig"= runif(V,.25,.30), "Om"= sample(1:V,nresp,replace=TRUE)   )}
        }else{
          jags.inits <- function(){ list("T"=matrix(rnorm(nitem*V,0,1),nitem,V),"Emu"= runif(V,.8,2),"Esig"= runif(V,.4,.6),"bsig"= runif(V,.4,.6),"lamsig"= runif(V,.25,.30)   )}
        }
        model.file <- gsub(pattern="a\\[i\\] <- exp\\(alog\\[i\\]\\)", "a\\[i\\] <- 1", model.file)
        model.file <- gsub(pattern="alog\\[i\\] ~ dnorm\\(amu\\[Om\\[i\\]\\],atau\\[Om\\[i\\]\\]\\)T\\(-2.3,2.3\\)", "", model.file)
        model.file <- gsub(pattern="atau\\[v\\] ~ dgamma\\(.01,.01\\)T\\(.01,\\)", "atau\\[v\\] <- 0", model.file)
        model.file <- gsub(pattern="for \\(c in 2\\:\\(C-1\\)\\)", " #for \\(c in 2\\:\\(C-1\\)\\)", model.file)
        model.file <- gsub(pattern="gam\\[1\\:\\(C-1\\),v\\] <- sort\\(tgam2\\[1\\:\\(C-1\\),v\\]\\)", "gam\\[1,v\\] <- 0 }", model.file)
        model.file <- gsub(pattern="for \\(c in 1\\:\\(C-2\\)\\)", "#for \\(c in 1\\:\\(C-2\\)\\)", model.file)
        model.file <- gsub(pattern="\\(1 - sum\\(pY\\[i,k,1\\:\\(C-1\\)\\]\\)\\)", "1 - pY\\[i,k,1\\]", model.file)
        model.file <- gsub(pattern="tgam2\\[1", "#tgam2\\[1", model.file)
        model.file <- gsub(pattern="tgam2\\[C", "#tgam2\\[C", model.file)
      }
    }
    
    if(clusters==1){ 
      model.file <- gsub(pattern="pi\\[1\\:V\\] ~ ddirch\\(L\\)", "pi <- 1", model.file)
      model.file <- gsub(pattern="Om\\[i\\] ~ dcat\\(pi\\)", "Om\\[i\\] <- 1", model.file)
    }
  }
  
  if(datob$whmodel=="CRM"){
    Y <- datob$datind; nresp <- dim(datob$dat)[1]; nitem <- dim(datob$dat)[2]; V <- clusters; nobs <- dim(datob$datind)[1]
    jags.data <- list("Y","nresp","nitem","V","nobs")
    
    if(itemdiff==FALSE){
      model.file <- mccrm
      jags.params <- c("T","E","a","alog","b","Tmu","Ttau","Emu","Etau","amu","atau","bmu","btau","Om","pi")
      if(clusters>1){
        jags.inits <- function(){ list("T"=matrix(rnorm(nitem*V,0,.5),nitem,V),"Emu"= runif(V,.8,2),"Esig"= runif(V,.5,.8), "asig"= runif(V,.4,.6),"bsig"= runif(V,.4,.6), "Om"= sample(1:V,nresp,replace=TRUE)   )}
      }else{
        jags.inits <- function(){ list("T"=matrix(rnorm(nitem*V,0,1),nitem,V),"Emu"= runif(V,.8,2),"Esig"= runif(V,.5,.8), "asig"= runif(V,.4,.6),"bsig"= runif(V,.4,.6)  )}
      }
    }
    
    if(itemdiff==TRUE){
      model.file <- mccrmid
      jags.params <- c("T","E","a","alog","b","lam","Tmu","Ttau","Emu","Etau","amu","atau","bmu","btau","lammu","lamtau","Om","pi")
      if(clusters>1){
        jags.inits <- function(){ list("T"=matrix(rnorm(nitem*V,0,1),nitem,V),"Emu"= runif(V,.8,2),"Esig"= runif(V,.4,.6), "asig"= runif(V,.4,.6),"bsig"= runif(V,.4,.6),"lamsig"= runif(V,.25,.30), "Om"= sample(1:V,nresp,replace=TRUE)   )}
      }else{
        jags.inits <- function(){ list("T"=matrix(rnorm(nitem*V,0,1),nitem,V),"Emu"= runif(V,.8,2),"Esig"= runif(V,.4,.6), "asig"= runif(V,.4,.6),"bsig"= runif(V,.4,.6),"lamsig"= runif(V,.25,.30)   )}
      }
    }
    if(clusters==1){ 
      model.file <- gsub(pattern="pi\\[1\\:V\\] ~ ddirch\\(L\\)", "pi <- 1", model.file)
      model.file <- gsub(pattern="Om\\[i\\] ~ dcat\\(pi\\)", "Om\\[i\\] <- 1", model.file)
    } 
  }
  
  ######################
  #Runs the Model in JAGS
  #- saves the data object from loadfilefunc within the jags object
  #- saves other useful details used later into the jags object
  ######################
  if(parallel==FALSE){
    cctfit <- jags(data=jags.data, inits=jags.inits, parameters.to.save=jags.params,
                   n.chains=jags.chains, n.iter=(jags.iter+jags.burnin), n.burnin=jags.burnin, 
                   n.thin=jags.thin, model.file=textConnection(model.file)) #textConnection(model.file)
    cctfit$dataind <- datob$datind; cctfit$data <- datob$dat; cctfit$n <- nresp; cctfit$m <- nitem; cctfit$V <- V
    cctfit$mval <- datob$mval; cctfit$itemdiff <- itemdiff; cctfit$checksrun <- FALSE; cctfit$whmodel <- datob$whmodel; cctfit$datob <- datob
    if(cctfit$mval==TRUE){cctfit$datamiss <- datob$datna}
    if(cctfit$whmodel=="LTRM"){cctfit$C <- C;cctfit$polycor <- FALSE}
  }else{
    ### Prepare / Evaluate parameters for Parallel JAGS run
    jags.burnin <- jags.burnin; jags.iter <- jags.iter; jags.chains <- jags.chains; jags.thin <- jags.thin; 
    tmpfn=tempfile()
    tmpcn=file(tmpfn,"w"); cat(model.file,file=tmpcn); close(tmpcn);
    
    message("\n ...Fitting model in parallel \n    note: progress bar currently not available with parallel option")  
    
    cctfit <- jags.parallel(data=jags.data, inits=jags.inits, parameters.to.save=jags.params,
                            n.chains=jags.chains, n.iter=(jags.iter+jags.burnin), n.burnin=jags.burnin, 
                            n.thin=jags.thin, model.file=tmpfn, envir=environment(),jags.seed = abs(.Random.seed[3]))
    cctfit$dataind <- datob$datind; cctfit$data <- datob$dat; cctfit$n <- nresp; cctfit$m <- nitem; cctfit$V <- V
    cctfit$mval <- datob$mval; cctfit$itemdiff <- itemdiff; cctfit$checksrun <- FALSE; cctfit$whmodel <- datob$whmodel; cctfit$datob <- datob
    if(cctfit$mval==TRUE){cctfit$datamiss <- datob$datna}
    if(cctfit$whmodel=="LTRM"){cctfit$C <- C}
  }
  ######################
  #Function used to calculate the Rhats
  ######################
  Rhat1 <- function(mat) {
    m <- ncol(mat)
    n <- nrow(mat)
    b <- apply(mat,2,mean)
    B <- sum((b-mean(mat))^2)*n/(m-1)
    w <- apply(mat,2,var)
    W <- mean(w)
    s2hat <- (n-1)/n*W + B/n
    Vhat <- s2hat + B/m/n 
    covWB <- n /m * (cov(w,b^2)-2*mean(b)*cov(w,b))
    varV <- (n-1)^2 / n^2 * var(w)/m +
      (m+1)^2 / m^2 / n^2 * 2*B^2/(m-1) +
      2 * (m-1)*(n-1)/m/n^2 * covWB
    df <- 2 * Vhat^2 / varV
    R <- sqrt((df+3) * Vhat / (df+1) / W)
    return(R)
  }
  
  Rhat <- function(arr) {
    dm <- dim(arr)
    if (length(dm)==2) return(Rhat1(arr))
    if (dm[2]==1) return(NULL)
    if (dm[3]==1) return(Rhat1(arr[,,1]))
    return(apply(arr,3,Rhat1))
  }
  
  ######################
  #Algorithm that corrects for label switching for the GCM
  ######################
  labelswitchalggcm <- function(cctfit,chnind=0){
    
    cctfit2 <- cctfit
    
    nch <- cctfit2$BUGSoutput$n.chains
    nsamp <- cctfit2$BUGSoutput$n.keep
    
    if(nch!=1){
      ntruths <- cctfit2$V
      truths <- array(NA,c(cctfit2$m,nch,ntruths))
      inds <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Z")]]
      inds <- matrix(inds,cctfit2$m,cctfit2$V)
      for(v in 1:cctfit$V){truths[,,v] <- t(apply(cctfit$BUGSoutput$sims.array[ ,,inds[,v]],c(2,3),mean))}
      
      V <- cctfit2$V
      
      chstart <- 1
      if(length(chnind)==1){ 
        chstart <- 2
        chnind <- array(NA,c(V,nch))
        chnind[1:V,1] <- 1:V
        
        for(v in 1:V){
          for(ch in chstart:nch){
            Tind <- c(1:V)[-chnind[,ch][!is.na(chnind[,ch])]]
            if(length(Tind)==0){Tind <- c(1:V)}
            
            chnind[v,ch] <- which(max(cor(truths[1:cctfit2$m, 1, v],truths[1:cctfit2$m, ch, Tind]))==cor(truths[1:cctfit2$m, 1, v],truths[1:cctfit2$m, ch, ]))
          }}
      }
      
      nsamp <- cctfit$BUGSoutput$n.keep
      
      if(cctfit$itemdiff==TRUE){
        inds2 <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lam")]]
        inds2 <- matrix(inds2,cctfit2$m,cctfit2$V)
        
        inds <- rbind(inds,
                      inds2,
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lammu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lamtau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="thmu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="thtau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="gmu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="gtau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="p")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="pi")]])
      }else{
        inds <- rbind(inds,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="thmu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="thtau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="gmu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="gtau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="p")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="pi")]])
      }
      
      einds <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]]
      tmpeinds <- cctfit$BUGSoutput$sims.array[ ,,einds]
      
      for(v in 1:cctfit$V){
        for(ch in chstart:nch){
          cctfit2$BUGSoutput$sims.array[ ,ch,inds[,v]] <- cctfit$BUGSoutput$sims.array[ ,ch,inds[,chnind[v,ch]]]  #this cctfit2 vs. cctfit difference is intentional
          if(chstart==2){cctfit2$BUGSoutput$sims.array[ ,ch,einds][cctfit$BUGSoutput$sims.array[ ,ch,einds]==chnind[v,ch]] <- chnind[v,1]}
          if(chstart==1){cctfit2$BUGSoutput$sims.array[ ,ch,einds][cctfit$BUGSoutput$sims.array[ ,ch,einds]==v] <- chnind[v,1]}
        }}
      
    }
    
    cctfit2$BUGSoutput$sims.list[["Om"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]]], c(nsamp*nch,cctfit$n))
    cctfit2$BUGSoutput$mean$Om <- apply(cctfit2$BUGSoutput$sims.list[["Om"]],2,mean)
    
    cctfit2$BUGSoutput$sims.matrix <- array(cctfit2$BUGSoutput$sims.array,c(nsamp*nch,dim(cctfit$BUGSoutput$sims.array)[3]))  #this cctfit2 vs. cctfit difference is intentional
    
    cctfit2$BUGSoutput$sims.list[["th"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="th")]]], c(nsamp*nch,cctfit$n))
    cctfit2$BUGSoutput$sims.list[["g"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="g")]]], c(nsamp*nch,cctfit$n))
    if(cctfit$itemdiff==TRUE){cctfit2$BUGSoutput$sims.list[["lam"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lam")]]], c(nsamp*nch,cctfit$m))}
    
    cctfit2$BUGSoutput$sims.list[["Z"]] <- array(NA, c(nsamp*nch,cctfit$m,cctfit$V))
    if(cctfit$itemdiff==TRUE){
      cctfit2$BUGSoutput$sims.list[["lam"]] <- array(NA, c(nsamp*nch,cctfit$m,cctfit$V))
      cctfit2$BUGSoutput$sims.list[["lammu"]] <- array(NA, c(nsamp*nch,cctfit$V))
      cctfit2$BUGSoutput$sims.list[["lamtau"]] <- array(NA, c(nsamp*nch,cctfit$V))
    }
    cctfit2$BUGSoutput$sims.list[["thmu"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["thtau"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["gmu"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["gtau"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["p"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["pi"]] <- array(NA, c(nsamp*nch,cctfit$V))
    
    for(v in 1:V){
      cctfit2$BUGSoutput$sims.list[["Z"]][,,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Z")]][1:cctfit$m +((v-1)*cctfit$m)]], c(nsamp*nch,cctfit$m))
      if(cctfit$itemdiff==TRUE){
        cctfit2$BUGSoutput$sims.list[["lam"]][,,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lam")]][1:cctfit$m +((v-1)*cctfit$m)]], c(nsamp*nch,cctfit$m))
        cctfit2$BUGSoutput$sims.list[["lammu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lammu")]][v]], c(nsamp*nch))
        cctfit2$BUGSoutput$sims.list[["lamtau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lamtau")]][v]], c(nsamp*nch))
      }
      cctfit2$BUGSoutput$sims.list[["thmu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="thmu")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["thtau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="thtau")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["gmu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="gmu")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["gtau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="gtau")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["p"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="p")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["pi"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="pi")]][v]], c(nsamp*nch))
      
      cctfit2$BUGSoutput$mean$Z[,v] <- apply(cctfit2$BUGSoutput$sims.list[["Z"]][,,v],2,mean)
      if(cctfit$itemdiff==TRUE){
        cctfit2$BUGSoutput$mean$lam[,v] <- apply(cctfit2$BUGSoutput$sims.list[["lam"]][,,v],2,mean)
        cctfit2$BUGSoutput$mean$lammu[v] <- mean(cctfit2$BUGSoutput$sims.list[["lammu"]][,v])
        cctfit2$BUGSoutput$mean$lamtau[v] <- mean(cctfit2$BUGSoutput$sims.list[["lamtau"]][,v])
      }
      cctfit2$BUGSoutput$mean$thmu[v] <- mean(cctfit2$BUGSoutput$sims.list[["thmu"]][,v])
      cctfit2$BUGSoutput$mean$thtau[v] <- mean(cctfit2$BUGSoutput$sims.list[["thtau"]][,v])
      cctfit2$BUGSoutput$mean$gmu[v] <- mean(cctfit2$BUGSoutput$sims.list[["gmu"]][,v])
      cctfit2$BUGSoutput$mean$gtau[v] <- mean(cctfit2$BUGSoutput$sims.list[["gtau"]][,v])
      cctfit2$BUGSoutput$mean$p[v] <- mean(cctfit2$BUGSoutput$sims.list[["p"]][,v])
      cctfit2$BUGSoutput$mean$pi[v] <- mean(cctfit2$BUGSoutput$sims.list[["pi"]][,v])
    }
    
    if(nch!=1){
      cctfit2$BUGSoutput$summary[,1] <- apply(apply(cctfit2$BUGSoutput$sims.array,c(2,3),mean),2,mean)
      cctfit2$BUGSoutput$summary[,2] <- apply(apply(cctfit2$BUGSoutput$sims.array,c(2,3),sd),2,sd)
      cctfit2$BUGSoutput$summary[,3:7] <- t(apply(cctfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
      cctfit2$BUGSoutput$summary[,8] <- Rhat(cctfit2$BUGSoutput$sims.array)
      cctfit2$BUGSoutput$summary[,8][is.nan(cctfit2$BUGSoutput$summary[,8])] <- 1.000000
      
      dimnames(cctfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(cctfit$BUGSoutput$summary))
      dimnames(cctfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(cctfit$BUGSoutput$summary))
    }else{
      
      cctfit2$BUGSoutput$summary[,1] <- apply(cctfit2$BUGSoutput$sims.array,2,mean)
      cctfit2$BUGSoutput$summary[,2] <- apply(cctfit2$BUGSoutput$sims.array,2,sd)
      cctfit2$BUGSoutput$summary[,3:7] <- t(apply(cctfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
      cctfit2$BUGSoutput$summary <- cctfit2$BUGSoutput$summary[,-c(8,9)]
      dimnames(cctfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(cctfit$BUGSoutput$summary))
      dimnames(cctfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(cctfit$BUGSoutput$summary))
    }
    
    cctfit <- cctfit2; rm(cctfit2)
    
    return(cctfit)
  }
  
  ######################
  #Algorithm that corrects for label switching for the LTRM
  ######################
  labelswitchalgltrm <- function(cctfit,chnind=0){
    
    cctfit2 <- cctfit
    
    nch <- cctfit2$BUGSoutput$n.chains
    nsamp <- cctfit2$BUGSoutput$n.keep
    
    if(nch!=1){
      ntruths <- cctfit2$V
      truths <- array(NA,c(cctfit2$m,nch,ntruths))
      inds <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="T")]]
      inds <- matrix(inds,cctfit2$m,cctfit2$V)
      for(v in 1:cctfit$V){truths[,,v] <- t(apply(cctfit$BUGSoutput$sims.array[ ,,inds[,v]],c(2,3),mean))}
      
      V <- cctfit2$V
      
      chstart <- 1
      if(length(chnind)==1){ 
        chstart <- 2
        chnind <- array(NA,c(V,nch))
        chnind[1:V,1] <- 1:V
        
        for(v in 1:V){
          for(ch in chstart:nch){
            Tind <- c(1:V)[-chnind[,ch][!is.na(chnind[,ch])]]
            if(length(Tind)==0){Tind <- c(1:V)}
            
            chnind[v,ch] <- which(max(cor(truths[1:cctfit2$m, 1, v],truths[1:cctfit2$m, ch, Tind]))==cor(truths[1:cctfit2$m, 1, v],truths[1:cctfit2$m, ch, ]))
          }}
      }
      
      nsamp <- cctfit$BUGSoutput$n.keep
      
      if(cctfit$itemdiff==TRUE){
        inds2 <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lam")]]
        inds2 <- matrix(inds2,cctfit2$m,cctfit2$V)
        
        inds <- rbind(inds,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Tmu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Ttau")]],
                      inds2,
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lammu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lamtau")]],
                      matrix(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="gam")]],cctfit2$C-1,cctfit2$V),
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Emu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Etau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="amu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="atau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="bmu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="btau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="pi")]])
        
        rm(inds2)
      }else{
        inds <- rbind(inds,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Tmu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Ttau")]],
                      matrix(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="gam")]],cctfit2$C-1,cctfit2$V),
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Emu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Etau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="amu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="atau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="bmu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="btau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="pi")]])
      }
      
      einds <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]]
      tmpeinds <- cctfit$BUGSoutput$sims.array[ ,,einds]
      
      for(v in 1:cctfit$V){
        for(ch in chstart:nch){
          cctfit2$BUGSoutput$sims.array[ ,ch,inds[,v]] <- cctfit$BUGSoutput$sims.array[ ,ch,inds[,chnind[v,ch]]]  #this cctfit2 vs. cctfit difference is intentional
          if(chstart==2){cctfit2$BUGSoutput$sims.array[ ,ch,einds][cctfit$BUGSoutput$sims.array[ ,ch,einds]==chnind[v,ch]] <- chnind[v,1]}
          if(chstart==1){cctfit2$BUGSoutput$sims.array[ ,ch,einds][cctfit$BUGSoutput$sims.array[ ,ch,einds]==v] <- chnind[v,1]}
        }}
      
    }
    cctfit2$BUGSoutput$sims.list[["Om"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]]], c(nsamp*nch,cctfit$n))
    cctfit2$BUGSoutput$mean$Om <- apply(cctfit2$BUGSoutput$sims.list[["Om"]],2,mean)
    
    cctfit2$BUGSoutput$sims.matrix <- array(cctfit2$BUGSoutput$sims.array,c(nsamp*nch,dim(cctfit$BUGSoutput$sims.array)[3]))  #this cctfit2 vs. cctfit difference is intentional
    
    cctfit2$BUGSoutput$sims.list[["E"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="E")]]], c(nsamp*nch,cctfit$n))
    cctfit2$BUGSoutput$sims.list[["a"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="a")]]], c(nsamp*nch,cctfit$n))
    cctfit2$BUGSoutput$sims.list[["b"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="b")]]], c(nsamp*nch,cctfit$n))
    if(cctfit$itemdiff==TRUE){cctfit2$BUGSoutput$sims.list[["lam"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lam")]]], c(nsamp*nch,cctfit$m))}
    
    cctfit2$BUGSoutput$sims.list[["T"]] <- array(NA, c(nsamp*nch,cctfit$m,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["Tmu"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["Ttau"]] <- array(NA, c(nsamp*nch,cctfit$V))
    if(cctfit$itemdiff==TRUE){
      cctfit2$BUGSoutput$sims.list[["lam"]] <- array(NA, c(nsamp*nch,cctfit$m,cctfit$V))
      cctfit2$BUGSoutput$sims.list[["lammu"]] <- array(NA, c(nsamp*nch,cctfit$V))
      cctfit2$BUGSoutput$sims.list[["lamtau"]] <- array(NA, c(nsamp*nch,cctfit$V))
    }
    cctfit2$BUGSoutput$sims.list[["gam"]] <- array(NA, c(nsamp*nch,cctfit$C-1,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["Emu"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["Etau"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["amu"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["atau"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["bmu"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["btau"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["pi"]] <- array(NA, c(nsamp*nch,cctfit$V))
    
    if(cctfit$C==2 && length(dim(cctfit2$BUGSoutput$sims.list[["gam"]])) < 3 ){
      cctfit2$BUGSoutput$sims.list[["gam"]] <- array(cctfit2$BUGSoutput$sims.list[["gam"]], c(nsamp*nch,cctfit$C-1,cctfit$V))
    }
    
    for(v in 1:cctfit2$V){
      cctfit2$BUGSoutput$sims.list[["T"]][,,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="T")]][1:cctfit$m +((v-1)*cctfit$m)]], c(nsamp*nch,cctfit$m))
      cctfit2$BUGSoutput$sims.list[["Tmu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Tmu")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["Ttau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Ttau")]][v]], c(nsamp*nch))
      if(cctfit$itemdiff==TRUE){
        cctfit2$BUGSoutput$sims.list[["lam"]][,,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lam")]][1:cctfit$m +((v-1)*cctfit$m)]], c(nsamp*nch,cctfit$m))
        cctfit2$BUGSoutput$sims.list[["lammu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lammu")]][v]], c(nsamp*nch))
        cctfit2$BUGSoutput$sims.list[["lamtau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lamtau")]][v]], c(nsamp*nch))
      }
      cctfit2$BUGSoutput$sims.list[["gam"]][,,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="gam")]][1:(cctfit$C-1) +((v-1)*(cctfit$C-1))]], c(nsamp*nch,cctfit$C-1))
      cctfit2$BUGSoutput$sims.list[["Emu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Emu")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["Etau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Etau")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["amu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="amu")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["atau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="atau")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["bmu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="bmu")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["btau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="btau")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["pi"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="pi")]][v]], c(nsamp*nch))
      
      cctfit2$BUGSoutput$mean$T[,v] <- apply(cctfit2$BUGSoutput$sims.list[["T"]][,,v],2,mean)
      cctfit2$BUGSoutput$mean$Tmu[v] <- mean(cctfit2$BUGSoutput$sims.list[["Tmu"]][,v])
      cctfit2$BUGSoutput$mean$Ttau[v] <- mean(cctfit2$BUGSoutput$sims.list[["Ttau"]][,v])
      if(cctfit$itemdiff==TRUE){
        cctfit2$BUGSoutput$mean$lam[,v] <- apply(cctfit2$BUGSoutput$sims.list[["lam"]][,,v],2,mean)
        cctfit2$BUGSoutput$mean$lammu[v] <- mean(cctfit2$BUGSoutput$sims.list[["lammu"]][,v])
        cctfit2$BUGSoutput$mean$lamtau[v] <- mean(cctfit2$BUGSoutput$sims.list[["lamtau"]][,v])
      }
      if(cctfit$C==2){
        cctfit2$BUGSoutput$mean$gam[v] <- mean(cctfit2$BUGSoutput$sims.list[["gam"]][,,v])
      }else{cctfit2$BUGSoutput$mean$gam[,v] <- apply(cctfit2$BUGSoutput$sims.list[["gam"]][,,v],2,mean)}
      cctfit2$BUGSoutput$mean$Emu[v] <- mean(cctfit2$BUGSoutput$sims.list[["Emu"]][,v])
      cctfit2$BUGSoutput$mean$Etau[v] <- mean(cctfit2$BUGSoutput$sims.list[["Etau"]][,v])
      cctfit2$BUGSoutput$mean$amu[v] <- mean(cctfit2$BUGSoutput$sims.list[["amu"]][,v])
      cctfit2$BUGSoutput$mean$atau[v] <- mean(cctfit2$BUGSoutput$sims.list[["atau"]][,v])
      cctfit2$BUGSoutput$mean$bmu[v] <- mean(cctfit2$BUGSoutput$sims.list[["bmu"]][,v])
      cctfit2$BUGSoutput$mean$btau[v] <- mean(cctfit2$BUGSoutput$sims.list[["btau"]][,v])
      cctfit2$BUGSoutput$mean$pi[v] <- mean(cctfit2$BUGSoutput$sims.list[["pi"]][,v])
    }
    
    if(nch!=1){
      cctfit2$BUGSoutput$summary[,1] <- apply(apply(cctfit2$BUGSoutput$sims.array,c(2,3),mean),2,mean)
      cctfit2$BUGSoutput$summary[,2] <- apply(apply(cctfit2$BUGSoutput$sims.array,c(2,3),sd),2,sd)
      cctfit2$BUGSoutput$summary[,3:7] <- t(apply(cctfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
      cctfit2$BUGSoutput$summary[,8] <- Rhat(cctfit2$BUGSoutput$sims.array)
      cctfit2$BUGSoutput$summary[,8][is.nan(cctfit2$BUGSoutput$summary[,8])] <- 1.000000
      dimnames(cctfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(cctfit$BUGSoutput$summary))
      dimnames(cctfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(cctfit$BUGSoutput$summary))
    }else{
      
      cctfit2$BUGSoutput$summary[,1] <- apply(cctfit2$BUGSoutput$sims.array,2,mean)
      cctfit2$BUGSoutput$summary[,2] <- apply(cctfit2$BUGSoutput$sims.array,2,sd)
      cctfit2$BUGSoutput$summary[,3:7] <- t(apply(cctfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
      cctfit2$BUGSoutput$summary <- cctfit2$BUGSoutput$summary[,-c(8,9)]
      dimnames(cctfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(cctfit$BUGSoutput$summary))
      dimnames(cctfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(cctfit$BUGSoutput$summary))
    }
    
    cctfit <- cctfit2; rm(cctfit2)
    
    return(cctfit)
  }
  
  ######################
  #Algorithm that corrects for label switching for the CRM
  ######################
  labelswitchalgcrm <- function(cctfit,chnind=0){
    
    cctfit2 <- cctfit
    
    nch <- cctfit2$BUGSoutput$n.chains
    nsamp <- cctfit2$BUGSoutput$n.keep
    
    if(nch!=1){
      ntruths <- cctfit2$V
      truths <- array(NA,c(cctfit2$m,nch,ntruths))
      inds <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="T")]]
      inds <- matrix(inds,cctfit2$m,cctfit2$V)
      for(v in 1:cctfit$V){truths[,,v] <- t(apply(cctfit$BUGSoutput$sims.array[ ,,inds[,v]],c(2,3),mean))}
      
      V <- cctfit2$V
      
      chstart <- 1
      if(length(chnind)==1){ 
        chstart <- 2
        chnind <- array(NA,c(V,nch))
        chnind[1:cctfit2$V,1] <- 1:cctfit2$V
        
        for(v in 1:cctfit2$V){
          for(ch in chstart:nch){
            Tind <- c(1:cctfit2$V)[-chnind[,ch][!is.na(chnind[,ch])]]
            if(length(Tind)==0){Tind <- c(1:cctfit2$V)}
            
            chnind[v,ch] <- which(max(cor(truths[1:cctfit2$m, 1, v],truths[1:cctfit2$m, ch, Tind]))==cor(truths[1:cctfit2$m, 1, v],truths[1:cctfit2$m, ch, ]))
          }}
      }
      
      nsamp <- cctfit$BUGSoutput$n.keep
      
      if(cctfit$itemdiff==TRUE){
        inds2 <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lam")]]
        inds2 <- matrix(inds2,cctfit2$m,cctfit2$V)
        
        inds <- rbind(inds,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Tmu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Ttau")]],
                      inds2,
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lammu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lamtau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Emu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Etau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="amu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="atau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="bmu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="btau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="pi")]])
        
        rm(inds2)
      }else{
        inds <- rbind(inds,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Tmu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Ttau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Emu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Etau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="amu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="atau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="bmu")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="btau")]],
                      cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="pi")]])
      }
      einds <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]]
      tmpeinds <- cctfit$BUGSoutput$sims.array[ ,,einds]
      
      for(v in 1:cctfit$V){
        for(ch in chstart:nch){
          cctfit2$BUGSoutput$sims.array[ ,ch,inds[,v]] <- cctfit$BUGSoutput$sims.array[ ,ch,inds[,chnind[v,ch]]]  #this cctfit2 vs. cctfit difference is intentional
          
          if(chstart==2){cctfit2$BUGSoutput$sims.array[ ,ch,einds][cctfit$BUGSoutput$sims.array[ ,ch,einds]==chnind[v,ch]] <- chnind[v,1]}
          if(chstart==1){cctfit2$BUGSoutput$sims.array[ ,ch,einds][cctfit$BUGSoutput$sims.array[ ,ch,einds]==v] <- chnind[v,1]}
        }}
      
    }
    
    cctfit2$BUGSoutput$sims.list[["Om"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]]], c(nsamp*nch,cctfit$n))
    cctfit2$BUGSoutput$mean$Om <- apply(cctfit2$BUGSoutput$sims.list[["Om"]],2,mean)
    
    cctfit2$BUGSoutput$sims.matrix <- array(cctfit2$BUGSoutput$sims.array,c(nsamp*nch,dim(cctfit$BUGSoutput$sims.array)[3]))  #this cctfit2 vs. cctfit difference is intentional
    
    cctfit2$BUGSoutput$sims.list[["E"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="E")]]], c(nsamp*nch,cctfit$n))
    cctfit2$BUGSoutput$sims.list[["a"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="a")]]], c(nsamp*nch,cctfit$n))
    cctfit2$BUGSoutput$sims.list[["b"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="b")]]], c(nsamp*nch,cctfit$n))
    if(cctfit$itemdiff==TRUE){cctfit2$BUGSoutput$sims.list[["lam"]] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lam")]]], c(nsamp*nch,cctfit$m))}
    
    cctfit2$BUGSoutput$sims.list[["T"]] <- array(NA, c(nsamp*nch,cctfit$m,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["Tmu"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["Ttau"]] <- array(NA, c(nsamp*nch,cctfit$V))
    if(cctfit$itemdiff==TRUE){
      cctfit2$BUGSoutput$sims.list[["lam"]] <- array(NA, c(nsamp*nch,cctfit$m,cctfit$V))
      cctfit2$BUGSoutput$sims.list[["lammu"]] <- array(NA, c(nsamp*nch,cctfit$V))
      cctfit2$BUGSoutput$sims.list[["lamtau"]] <- array(NA, c(nsamp*nch,cctfit$V))
    }
    cctfit2$BUGSoutput$sims.list[["Emu"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["Etau"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["amu"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["atau"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["bmu"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["btau"]] <- array(NA, c(nsamp*nch,cctfit$V))
    cctfit2$BUGSoutput$sims.list[["pi"]] <- array(NA, c(nsamp*nch,cctfit$V))
    
    for(v in 1:cctfit2$V){
      cctfit2$BUGSoutput$sims.list[["T"]][,,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="T")]][1:cctfit$m +((v-1)*cctfit$m)]], c(nsamp*nch,cctfit$m))
      cctfit2$BUGSoutput$sims.list[["Tmu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Tmu")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["Ttau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Ttau")]][v]], c(nsamp*nch))
      if(cctfit$itemdiff==TRUE){
        cctfit2$BUGSoutput$sims.list[["lam"]][,,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lam")]][1:cctfit$m +((v-1)*cctfit$m)]], c(nsamp*nch,cctfit$m))
        cctfit2$BUGSoutput$sims.list[["lammu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lammu")]][v]], c(nsamp*nch))
        cctfit2$BUGSoutput$sims.list[["lamtau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="lamtau")]][v]], c(nsamp*nch))
      }
      cctfit2$BUGSoutput$sims.list[["Emu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Emu")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["Etau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Etau")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["amu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="amu")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["atau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="atau")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["bmu"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="bmu")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["btau"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="btau")]][v]], c(nsamp*nch))
      cctfit2$BUGSoutput$sims.list[["pi"]][,v] <- array(cctfit2$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="pi")]][v]], c(nsamp*nch))
      
      cctfit2$BUGSoutput$mean$T[,v] <- apply(cctfit2$BUGSoutput$sims.list[["T"]][,,v],2,mean)
      cctfit2$BUGSoutput$mean$Tmu[v] <- mean(cctfit2$BUGSoutput$sims.list[["Tmu"]][,v])
      cctfit2$BUGSoutput$mean$Ttau[v] <- mean(cctfit2$BUGSoutput$sims.list[["Ttau"]][,v])
      if(cctfit$itemdiff==TRUE){
        cctfit2$BUGSoutput$mean$lam[,v] <- apply(cctfit2$BUGSoutput$sims.list[["lam"]][,,v],2,mean)
        cctfit2$BUGSoutput$mean$lammu[v] <- mean(cctfit2$BUGSoutput$sims.list[["lammu"]][,v])
        cctfit2$BUGSoutput$mean$lamtau[v] <- mean(cctfit2$BUGSoutput$sims.list[["lamtau"]][,v])
      }
      cctfit2$BUGSoutput$mean$Emu[v] <- mean(cctfit2$BUGSoutput$sims.list[["Emu"]][,v])
      cctfit2$BUGSoutput$mean$Etau[v] <- mean(cctfit2$BUGSoutput$sims.list[["Etau"]][,v])
      cctfit2$BUGSoutput$mean$amu[v] <- mean(cctfit2$BUGSoutput$sims.list[["amu"]][,v])
      cctfit2$BUGSoutput$mean$atau[v] <- mean(cctfit2$BUGSoutput$sims.list[["atau"]][,v])
      cctfit2$BUGSoutput$mean$bmu[v] <- mean(cctfit2$BUGSoutput$sims.list[["bmu"]][,v])
      cctfit2$BUGSoutput$mean$btau[v] <- mean(cctfit2$BUGSoutput$sims.list[["btau"]][,v])
      cctfit2$BUGSoutput$mean$pi[v] <- mean(cctfit2$BUGSoutput$sims.list[["pi"]][,v])
    }
    
    if(nch!=1){
      cctfit2$BUGSoutput$summary[,1] <- apply(apply(cctfit2$BUGSoutput$sims.array,c(2,3),mean),2,mean)
      cctfit2$BUGSoutput$summary[,2] <- apply(apply(cctfit2$BUGSoutput$sims.array,c(2,3),sd),2,sd)
      cctfit2$BUGSoutput$summary[,3:7] <- t(apply(cctfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
      cctfit2$BUGSoutput$summary[,8] <- Rhat(cctfit2$BUGSoutput$sims.array)
      cctfit2$BUGSoutput$summary[,8][is.nan(cctfit2$BUGSoutput$summary[,8])] <- 1.000000
      dimnames(cctfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(cctfit$BUGSoutput$summary))
      dimnames(cctfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(cctfit$BUGSoutput$summary))
    }else{
      
      cctfit2$BUGSoutput$summary[,1] <- apply(cctfit2$BUGSoutput$sims.array,2,mean)
      cctfit2$BUGSoutput$summary[,2] <- apply(cctfit2$BUGSoutput$sims.array,2,sd)
      cctfit2$BUGSoutput$summary[,3:7] <- t(apply(cctfit2$BUGSoutput$sims.matrix,2,function(x) quantile(x,probs=c(.025,.25,.50,.75,.975))))
      cctfit2$BUGSoutput$summary <- cctfit2$BUGSoutput$summary[,-c(8,9)]
      dimnames(cctfit2$BUGSoutput$sims.array) <- list(NULL,NULL,rownames(cctfit$BUGSoutput$summary))
      dimnames(cctfit2$BUGSoutput$sims.matrix) <- list(NULL,rownames(cctfit$BUGSoutput$summary))
    }
    
    cctfit <- cctfit2; rm(cctfit2)
    
    return(cctfit)
  }
  
  
  ######################
  #This is run after the jags inference, we are still in the 'Apply CCT Model' function
  #- Reports the number of Rhats above 1.05 and 1.10   (before applying the algorithm that corrects for label-switching)
  #- Detects if a mixture model was applied
  #- If so, the appropriate label-correcting algorithm is applied
  #- Recalculates the Rhats, DIC, and statistics for the parameters
  #- Reports the number of Rhats above 1.05 and 1.10 
  #- Outputs the DIC that is calculated after the label-correcting algorithm (if applicable)
  ######################
#   message("\n ...Inference complete, data is saved as '",guidat$varname,"'")
#   if(clusters > 1){
#     message("\n    'cctfit$respmem' provides the respondent clustering")
#   }
  
  message("\n ...Performing final calculations")
  if(cctfit$BUGSoutput$n.chains > 1){
    cctfit$BUGSoutput$summary[,8] <- Rhat(cctfit$BUGSoutput$sims.array)
    cctfit$BUGSoutput$summary[,8][is.nan(cctfit$BUGSoutput$summary[,8])] <- 1.000000
    if(cctfit$whmodel== "GCM"){
       cctfit$Rhat$ncp <- length(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Z")]]),8])
       cctfit$Rhat$above110 <- sum(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Z")]]),8]>1.10)
       cctfit$Rhat$above105 <- sum(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Z")]]),8]>1.050)
      }else{
        cctfit$Rhat$ncp <- length(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]]),8])
        cctfit$Rhat$above110 <- sum(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]]),8]>1.10)
        cctfit$Rhat$above105 <- sum(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]]),8]>1.050)
    }
    message(paste("\nFor Continuous Parameters"))
    message(paste("Number of Rhats above 1.10 : ",cctfit$Rhat$above110,"/",cctfit$Rhat$ncp,"\nNumber of Rhats above 1.05 : ",cctfit$Rhat$above105,"/",cctfit$Rhat$ncp))
  }
  
  if(cctfit$V>1 && cctfit$BUGSoutput$n.chains==1){
    
    tmp <- unique(apply(cctfit$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]]],c(2),Mode))
    if(length(tmp)==1){
      message(paste("\n ...This chain has ",length(tmp)," culture rather than the ",cctfit$V," cultures requested", sep=""))
      message(paste("\n ...Try running the inference again",sep="" ))
    }
    if(length(tmp)!=1 && length(tmp) < cctfit$V){
      message(paste("\n ...This chain has ",tmp," cultures rather than the ",cctfit$V," cultures requested", sep=""))
      message(paste("\n ...Try running the inference again",sep="" ))
    }
  }
  
  if(cctfit$V>1 && cctfit$BUGSoutput$n.chains > 1){
    message("\n ...More than 1 culture applied with more than 1 chain")
    
    einds <- cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]]
    
    tmp <- apply(apply(cctfit$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]]],c(2,3),Mode),1,unique)
    
    chntoremove <- NULL
    if(is.list(tmp)){
      for(i in 1:cctfit$BUGSoutput$n.chains){
        if(length(tmp[[i]])!=cctfit$V){chntoremove <- c(chntoremove,i)}
      }
    }
    
    if(length(dim(tmp))==2){
      for(i in 1:cctfit$BUGSoutput$n.chains){
        if(length(tmp[,i])!=cctfit$V){chntoremove <- c(chntoremove,i)}
      }
    }
    
    if(is.null(dim(tmp)) && !is.list(tmp)){chntoremove <- c(1:cctfit$BUGSoutput$n.chains)}
    
    
    if(length(chntoremove)>0){
      if(length(chntoremove) < cctfit$BUGSoutput$n.chains){
        
        if(length(chntoremove)==1){
          message(paste("\n ...",length(chntoremove)," chain out of ",cctfit$BUGSoutput$n.chains," had fewer than ",cctfit$V," cultures requested",sep=""))
          message(paste("\n ...", "removing the ", length(chntoremove)," chain", sep=""))
        }else{
          message(paste("\n ...",length(chntoremove)," chains out of ",cctfit$BUGSoutput$n.chains," had fewer than ",cctfit$V," cultures requested",sep=""))
          message(paste("\n ...", "removing these ", length(chntoremove)," chains", sep=""))
        }
        
        cctfit$BUGSoutput$n.chains <- cctfit$BUGSoutput$n.chains-length(chntoremove)
        cctfit$BUGSoutput$n.chain <- cctfit$BUGSoutput$n.chains
        cctfit$BUGSoutput$n.sims <- cctfit$BUGSoutput$n.chains*cctfit$BUGSoutput$n.keep
        if(cctfit$BUGSoutput$n.chain==1){
          cctfit$BUGSoutput$sims.array <- array(cctfit$BUGSoutput$sims.array[,-chntoremove,], c(dim(cctfit$BUGSoutput$sims.array)[1],1,dim(cctfit$BUGSoutput$sims.array)[3]))
        }else{cctfit$BUGSoutput$sims.array <- cctfit$BUGSoutput$sims.array[,-chntoremove,]}
        
      }else{
        message(paste("\n ...All chains out of ",cctfit$BUGSoutput$n.chains," had fewer than ",cctfit$V," cultures requested", sep=""))
        message(paste("\n ...Try running the inference again",sep="" ))
      }
    }
    
    message("\n ...Computing the most-consistent labeling across chains")
    
    if(cctfit$whmodel=="GCM"){
      
      cctfit <- labelswitchalggcm(cctfit)
      
      cctfit$respmem <- apply(cctfit$BUGSoutput$sims.list$Om[,],2,Mode)
      tmeans <- cctfit$BUGSoutput$mean$Z
      tmeans[tmeans<.5] <- tmeans[tmeans<.5]+1
      ind <- rank(apply(abs(1-tmeans),2,mean))
    } 
    
    
    if(cctfit$whmodel=="LTRM"){
      
      cctfit <- labelswitchalgltrm(cctfit)
      
      cctfit$respmem <- apply(cctfit$BUGSoutput$sims.list$Om[,],2,Mode)
      tmeans <- cctfit$BUGSoutput$mean$T
      ind <- rank(-apply(tmeans,2,sd))
      
    }
    
    if(cctfit$whmodel=="CRM"){
      
      cctfit <- labelswitchalgcrm(cctfit)
      
      cctfit$respmem <- apply(cctfit$BUGSoutput$sims.list$Om[,],2,Mode)
      tmeans <- cctfit$BUGSoutput$mean$T
      ind <- rank(-apply(tmeans,2,sd))
    }
    
    if(cctfit$BUGSoutput$n.chains > 1){
      message(paste("\nFor Continuous Parameters:"))
      if(cctfit$whmodel== "GCM"){
        message(paste("Number of Rhats above 1.10 : ",sum(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Z")]]),8]>1.10),"/",length(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Z")]]),8]),"\nNumber of Rhats above 1.05 : ",sum(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Z")]]),8]>1.050),"/",length(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Z")]]),8]),sep=""))
      }else{
        message(paste("Number of Rhats above 1.10 : ",sum(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]]),8]>1.10),"/",length(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]]),8]),"\nNumber of Rhats above 1.05 : ",sum(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]]),8]>1.050),"/",length(cctfit$BUGSoutput$summary[-c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]]),8]),sep=""))
      }
    }
    if(cctfit$BUGSoutput$n.chains > 1){
      if(gui==TRUE){message(paste("\nFor Discrete Parameters:"))
                 message(paste("Use button 'Traceplot Discrete' to see their trace plots",sep=""))
                 }else{
                   message(paste("\nFor Discrete Parameters:"))
                   message(paste("Use function 'dtraceplot()' to see their trace plots",sep=""))  
                   message(paste("\nUse function 'cctmemb()' to see the respondent cluster assignments",sep=""))
                   
                 }
    }
    
  }else{
    cctfit$respmem <- rep(1,cctfit$n)
    if(cctfit$BUGSoutput$n.chains > 1){
      if(gui==TRUE){message(paste("\nFor Discrete Parameters:"))
                    message(paste("Use button 'Traceplot Discrete' to see their trace plots",sep=""))
      }else{
        message(paste("\nFor Discrete Parameters:"))
        message(paste("Use function 'dtraceplot()' to see their trace plots",sep=""))  
      }
    }
    
  }
  
  ######################
  #DIC is calculated for the model that was applied
  #The likelihood of the model is evaluated at each node of each sample and saved to cctfit$Lik
  ######################
  message("\n ...Calculating DIC")
  cctfit$BUGSoutput$pD <- var(cctfit$BUGSoutput$sims.list$deviance[,1])/2
  nsimstouse <- min(cctfit$BUGSoutput$n.sims,1000)
  ind <- sample(1:cctfit$BUGSoutput$n.sims,nsimstouse)
  
  storelik <- 0
  
  if(cctfit$whmodel=="GCM"){
    
    if(cctfit$V==1){
      cctfit$V <- 1;
      cctfit$BUGSoutput$sims.list$Om <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
      if(length(dim(cctfit$BUGSoutput$sims.list$Z)) < 3){
        cctfit$BUGSoutput$sims.list$Z <- array(cctfit$BUGSoutput$sims.list[["Z"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
      if(cctfit$itemdiff==TRUE){
        if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
          cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
      }
    }
    if(cctfit$itemdiff==FALSE){
      if(cctfit$V==1){
        cctfit$BUGSoutput$sims.list$lam <- array(.5,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
      }
      if(cctfit$V > 1){
        cctfit$BUGSoutput$sims.list$lam <- array(.5,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
      }
    }
    
    cctfit$BUGSoutput$sims.list$D <- array(-1,c(cctfit$BUGSoutput$n.sims,cctfit$n,cctfit$m))
    
    storefulllik <- 0
    if(storefulllik==1){
      cctfit$Lik  <- array(NA, c(cctfit$n,cctfit$m,cctfit$BUGSoutput$n.sims));
      for(samp in 1:cctfit$BUGSoutput$n.sims){
        cctfit$BUGSoutput$sims.list[["D"]][samp,,] <- (cctfit$BUGSoutput$sims.list[["th"]][samp,]*t(1-cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])) / 
          ( (cctfit$BUGSoutput$sims.list[["th"]][samp,]*t(1-cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])) + 
              ((1-cctfit$BUGSoutput$sims.list[["th"]][samp,])*t(cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])) ) 
        
        cctfit$Lik[,,samp] <- t( ( (t(cctfit$BUGSoutput$sims.list[["D"]][samp,,])+(t(1-cctfit$BUGSoutput$sims.list[["D"]][samp,,])%*%diag(cctfit$BUGSoutput$sims.list[["g"]][samp,])))^(t(cctfit$data[,])*cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]]))*(t(1-cctfit$BUGSoutput$sims.list[["D"]][samp,,])%*%diag(cctfit$BUGSoutput$sims.list[["g"]][samp,]))^(t(cctfit$data[,])*(1-cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]]))*(t(1-cctfit$BUGSoutput$sims.list[["D"]][samp,,])%*%diag(1-cctfit$BUGSoutput$sims.list[["g"]][samp,]))^((1-t(cctfit$data[,]))*cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])*(t(cctfit$BUGSoutput$sims.list[["D"]][samp,,])+t(1-cctfit$BUGSoutput$sims.list[["D"]][samp,,])%*%diag(1-cctfit$BUGSoutput$sims.list[["g"]][samp,]))^((1-t(cctfit$data[,]))*(1-cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]]))  )
        if(cctfit$mval==TRUE){cctfit$Lik[cbind(datob$thena,samp)] <- 1}
      }
      
      cctfit$BUGSoutput$sims.list$deviance <- array(-2*apply(log(cctfit$Lik),3,sum),c(cctfit$BUGSoutput$n.sims,1))
    }else{
      cctfit$LogLik  <- array(NA, c(cctfit$BUGSoutput$n.sims));
      for(samp in 1:cctfit$BUGSoutput$n.sims){
        Dtmp <- (cctfit$BUGSoutput$sims.list[["th"]][samp,]*t(1-cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])) / 
          ( (cctfit$BUGSoutput$sims.list[["th"]][samp,]*t(1-cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])) + 
              ((1-cctfit$BUGSoutput$sims.list[["th"]][samp,])*t(cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])) ) 
        
        Liktmp <- t( ( (t(Dtmp)+(t(1-Dtmp)%*%diag(cctfit$BUGSoutput$sims.list[["g"]][samp,])))^(t(cctfit$data[,])*cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]]))*(t(1-Dtmp)%*%diag(cctfit$BUGSoutput$sims.list[["g"]][samp,]))^(t(cctfit$data[,])*(1-cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]]))*(t(1-Dtmp)%*%diag(1-cctfit$BUGSoutput$sims.list[["g"]][samp,]))^((1-t(cctfit$data[,]))*cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])*(t(Dtmp)+t(1-Dtmp)%*%diag(1-cctfit$BUGSoutput$sims.list[["g"]][samp,]))^((1-t(cctfit$data[,]))*(1-cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]]))  )
        if(cctfit$mval==TRUE){Liktmp[datob$thena] <- 1}
        
        cctfit$LogLik[samp] <- sum(log(Liktmp))
      }
      
      cctfit$BUGSoutput$sims.list$deviance <- array(-2*(cctfit$LogLik),c(cctfit$BUGSoutput$n.sims,1))
    }
    
    cctfit$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="deviance")]]] <- array(cctfit$BUGSoutput$sims.list$deviance,c(cctfit$BUGSoutput$n.keep,cctfit$BUGSoutput$n.chains))
    cctfit$BUGSoutput$sims.matrix[,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="deviance")]]] <- cctfit$BUGSoutput$sims.list$deviance[,1]
    
    if(sum(cctfit$BUGSoutput$sims.list$deviance==Inf) >0){ 
      cctfit$BUGSoutput$mean$deviance <- mean(cctfit$BUGSoutput$sims.list$deviance[cctfit$BUGSoutput$sims.list$deviance!=Inf,1],na.rm=TRUE) #Dbar, also known as deviance
      cctfit$BUGSoutput$pD <- var(cctfit$BUGSoutput$sims.list$deviance[cctfit$BUGSoutput$sims.list$deviance!=Inf,1],na.rm=TRUE)/2  #pD, variance of the deviance, divided by 2
      cctfit$BUGSoutput$DIC <- cctfit$BUGSoutput$mean$deviance + cctfit$BUGSoutput$pD
    }else{
      cctfit$BUGSoutput$mean$deviance <- mean(cctfit$BUGSoutput$sims.list$deviance) #Dbar, also known as deviance
      cctfit$BUGSoutput$pD <- var(cctfit$BUGSoutput$sims.list$deviance[,1])/2  #pD, variance of the deviance, divided by 2
      cctfit$BUGSoutput$DIC <- cctfit$BUGSoutput$mean$deviance + cctfit$BUGSoutput$pD
    }
    
  }
  
  if(cctfit$whmodel=="LTRM"){
    if(cctfit$V==1){
      cctfit$V <- 1;
      cctfit$BUGSoutput$sims.list$Om <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
      if(length(dim(cctfit$BUGSoutput$sims.list$T)) < 3){
        cctfit$BUGSoutput$sims.list$T <- array(cctfit$BUGSoutput$sims.list[["T"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
      if(length(dim(cctfit$BUGSoutput$sims.list$gam))<3){cctfit$BUGSoutput$sims.list$gam <- array(cctfit$BUGSoutput$sims.list[["gam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$C-1,1))}
      if(cctfit$itemdiff==TRUE){
        if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
          cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
      }
    }
    if(cctfit$itemdiff==FALSE){
      if(cctfit$V==1){
        cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
      }
      if(cctfit$V > 1){
        cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
      }
    }
    
    
    storefulllik <- 0
    if(storefulllik==1){
      
      cctfit$BUGSoutput$sims.list$tau <- array(-1,c(cctfit$BUGSoutput$n.sims,cctfit$n,cctfit$m))
      cctfit$ppdelta <- array(NA, c(cctfit$n,cctfit$C-1,cctfit$BUGSoutput$n.sims))
      cctfit$ppdeltafull <- array(NA, c(cctfit$n,cctfit$C+1,cctfit$BUGSoutput$n.sims))
      cctfit$ppdeltafull[,1,] <- -1000000; cctfit$ppdeltafull[,cctfit$C+1,] <- 1000000
      cctfit$Lik  <- array(NA, c(cctfit$n,cctfit$m,cctfit$BUGSoutput$n.sims));
      
      for(samp in 1:cctfit$BUGSoutput$n.sims){
        cctfit$BUGSoutput$sims.list[["tau"]][samp,,] <- cctfit$BUGSoutput$sims.list[["E"]][samp,]*t(1/cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])
        cctfit$ppdeltafull[,2:(cctfit$C),samp] <- cctfit$BUGSoutput$sims.list[["a"]][samp,]*t(cctfit$BUGSoutput$sims.list[["gam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])+cctfit$BUGSoutput$sims.list[["b"]][samp,]
        
        for(i in 1:cctfit$n){
          cctfit$Lik[i,,samp] <-  pnorm(cctfit$ppdeltafull[i,cctfit$data[i,]+1,samp] ,cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,i]],cctfit$BUGSoutput$sims.list[["tau"]][samp,i,]^-.5)-pnorm(cctfit$ppdeltafull[i,cctfit$data[i,],samp] ,cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,i]],cctfit$BUGSoutput$sims.list[["tau"]][samp,i,]^-.5)
          if(cctfit$mval==TRUE){cctfit$Lik[cbind(datob$thena,samp)] <- 1}
        }
      }
      
      if(cctfit$C==2){cctfit$Lik[cctfit$Lik==0] <- 0.001}
      
      cctfit$BUGSoutput$sims.list$deviance <- array(-2*apply(log(cctfit$Lik),3,sum),c(cctfit$BUGSoutput$n.sims,1))
    }else{
      cctfit$LogLik  <- array(NA, c(cctfit$BUGSoutput$n.sims));
      Liktmp <- array(NA, c(cctfit$n,cctfit$m))
      deltatmp <- array(NA, c(cctfit$n,cctfit$C+1))
      deltatmp[,1] <- -100000; deltatmp[,cctfit$C+1] <- 1000000
      
      for(samp in 1:cctfit$BUGSoutput$n.sims){
        tautmp <- (cctfit$BUGSoutput$sims.list[["E"]][samp,]*t(exp(cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])))^-2
        deltatmp[,2:(cctfit$C)] <- cctfit$BUGSoutput$sims.list[["a"]][samp,]*t(cctfit$BUGSoutput$sims.list[["gam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])+cctfit$BUGSoutput$sims.list[["b"]][samp,]
        
        for(i in 1:cctfit$n){
          Liktmp[i,] <-  pnorm(deltatmp[i,cctfit$data[i,]+1],cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,i]],tautmp[i,]^-.5)-pnorm(deltatmp[i,cctfit$data[i,]],cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,i]],tautmp[i,]^-.5)
          if(cctfit$mval==TRUE){cctfit$Lik[datob$thena] <- 1}
        }
        if(cctfit$C==2){Liktmp[Liktmp==0] <- 0.001}
        
        cctfit$LogLik[samp] <- sum(log(Liktmp))
      }
      
      cctfit$BUGSoutput$sims.list$deviance <- array(-2*(cctfit$LogLik),c(cctfit$BUGSoutput$n.sims,1))
    }
    
    cctfit$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="deviance")]]] <- array(cctfit$BUGSoutput$sims.list$deviance,c(cctfit$BUGSoutput$n.keep,cctfit$BUGSoutput$n.chains))
    cctfit$BUGSoutput$sims.matrix[,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="deviance")]]] <- cctfit$BUGSoutput$sims.list$deviance[,1]
    
    
    if(sum(cctfit$BUGSoutput$sims.list$deviance==Inf) >0){ 
      cctfit$BUGSoutput$mean$deviance <- mean(cctfit$BUGSoutput$sims.list$deviance[cctfit$BUGSoutput$sims.list$deviance!=Inf,1],na.rm=TRUE) #Dbar, also known as deviance
      cctfit$BUGSoutput$pD <- var(cctfit$BUGSoutput$sims.list$deviance[cctfit$BUGSoutput$sims.list$deviance!=Inf,1],na.rm=TRUE)/2  #pD, variance of the deviance, divided by 2
      cctfit$BUGSoutput$DIC <- cctfit$BUGSoutput$mean$deviance + cctfit$BUGSoutput$pD
    }else{
      cctfit$BUGSoutput$mean$deviance <- mean(cctfit$BUGSoutput$sims.list$deviance) #Dbar, also known as deviance
      cctfit$BUGSoutput$pD <- var(cctfit$BUGSoutput$sims.list$deviance[,1])/2  #pD, variance of the deviance, divided by 2
      cctfit$BUGSoutput$DIC <- cctfit$BUGSoutput$mean$deviance + cctfit$BUGSoutput$pD
    }
    
  }
  
  if(cctfit$whmodel=="CRM"){
    if(cctfit$V==1){
      cctfit$V <- 1;
      cctfit$BUGSoutput$sims.list$Om <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
      if(length(dim(cctfit$BUGSoutput$sims.list$T)) < 3){
        cctfit$BUGSoutput$sims.list$T <- array(cctfit$BUGSoutput$sims.list[["T"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
      if(cctfit$itemdiff==TRUE){
        if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
          cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
      }
    }
    if(cctfit$itemdiff==FALSE){
      if(cctfit$V==1){
        cctfit$BUGSoutput$sims.list$lam <- array(0,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
      }
      if(cctfit$V > 1){
        cctfit$BUGSoutput$sims.list$lam <- array(0,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
      }
    }
    
    storefulllik <- 0
    if(storefulllik==1){
      cctfit$BUGSoutput$sims.list$tau <- array(-1,c(cctfit$BUGSoutput$n.sims,cctfit$n,cctfit$m))
      cctfit$Lik  <- array(NA, c(cctfit$n,cctfit$m,cctfit$BUGSoutput$n.sims));
      
      for(samp in 1:cctfit$BUGSoutput$n.sims){
        cctfit$BUGSoutput$sims.list[["tau"]][samp,,] <-  (cctfit$BUGSoutput$sims.list[["E"]][samp,]*t(exp(cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])))^-2
        cctfit$Lik[,,samp] <-  dnorm(cctfit$data,mean=(cctfit$BUGSoutput$sims.list[["a"]][samp,]*t(cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]]))+array(cctfit$BUGSoutput$sims.list[["b"]][samp,],c(cctfit$n,cctfit$m)),sd=cctfit$BUGSoutput$sims.list[["a"]][samp,]*(cctfit$BUGSoutput$sims.list[["tau"]][samp,,]^-.5))
        if(cctfit$mval==TRUE){cctfit$Lik[cbind(datob$thena,samp)] <- 1}
      }
      cctfit$BUGSoutput$sims.list$deviance <- array(-2*apply(log(cctfit$Lik),3,sum),c(cctfit$BUGSoutput$n.sims,1))
    }else{
      cctfit$LogLik  <- array(NA, c(cctfit$BUGSoutput$n.sims));
      
      for(samp in 1:cctfit$BUGSoutput$n.sims){
        tautmp <- (cctfit$BUGSoutput$sims.list[["E"]][samp,]*t(exp(cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])))^-2
        Liktmp <- dnorm(cctfit$data,mean=(cctfit$BUGSoutput$sims.list[["a"]][samp,]*t(cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]]))+array(cctfit$BUGSoutput$sims.list[["b"]][samp,],c(cctfit$n,cctfit$m)),sd=cctfit$BUGSoutput$sims.list[["a"]][samp,]*(tautmp^-.5))
        if(cctfit$mval==TRUE){Liktmp[datob$thena] <- 1}
        cctfit$LogLik[samp] <- sum(log(Liktmp))
      }
      cctfit$BUGSoutput$sims.list$deviance <- array(-2*(cctfit$LogLik),c(cctfit$BUGSoutput$n.sims,1))
    }
    
    cctfit$BUGSoutput$sims.array[,,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="deviance")]]] <- array(cctfit$BUGSoutput$sims.list$deviance,c(cctfit$BUGSoutput$n.keep,cctfit$BUGSoutput$n.chains))
    cctfit$BUGSoutput$sims.matrix[,cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="deviance")]]] <- cctfit$BUGSoutput$sims.list$deviance[,1]
    
    if(sum(cctfit$BUGSoutput$sims.list$deviance==Inf) >0){ 
      cctfit$BUGSoutput$mean$deviance <- mean(cctfit$BUGSoutput$sims.list$deviance[cctfit$BUGSoutput$sims.list$deviance!=Inf,1],na.rm=TRUE) #Dbar, also known as deviance
      cctfit$BUGSoutput$pD <- var(cctfit$BUGSoutput$sims.list$deviance[cctfit$BUGSoutput$sims.list$deviance!=Inf,1],na.rm=TRUE)/2  #pD, variance of the deviance, divided by 2
      cctfit$BUGSoutput$DIC <- cctfit$BUGSoutput$mean$deviance + cctfit$BUGSoutput$pD
    }else{
      cctfit$BUGSoutput$mean$deviance <- mean(cctfit$BUGSoutput$sims.list$deviance) #Dbar, also known as deviance
      cctfit$BUGSoutput$pD <- var(cctfit$BUGSoutput$sims.list$deviance[,1])/2  #pD, variance of the deviance, divided by 2
      cctfit$BUGSoutput$DIC <- cctfit$BUGSoutput$mean$deviance + cctfit$BUGSoutput$pD
    }
    
  }
  
  message(paste("DIC : ",round(cctfit$BUGSoutput$DIC,2),"   pD : ",round(cctfit$BUGSoutput$pD,2),sep=""))
  
  if(gui==TRUE){
    guidat <- get("guidat", pkg_globals)
    tkconfigure(guidat$plotresults.but, state="normal") 
    tkconfigure(guidat$doppc.but, state="normal") 
    tkconfigure(guidat$exportresults.but, state="normal")
    tkconfigure(guidat$printfit.but, state="normal")
    tkconfigure(guidat$summary.but, state="normal")
    tkconfigure(guidat$sendconsole.but, state="normal")
    tkconfigure(guidat$memb.but, state="normal")
    tkconfigure(guidat$mvest.but, state="normal")
    #if(guidat$mval == TRUE){tkconfigure(guidat$mvest.but, state="normal")}
    tkconfigure(guidat$traceplotdiscrete.but, state="normal")
    tkconfigure(guidat$traceplotall.but, state="normal")
    cctfit$guidat <- guidat
    assign("guidat", guidat, pkg_globals)
    class(cctfit) <- c("cct",class(cctfit))
    return(cctfit)
  }else{
    class(cctfit) <- c("cct",class(cctfit))
    return(cctfit)
  }
  
  }

######################
#Function for the 'Plot Results' Button
#- Creates a sophisticated plot of the posterior results, which includes
#    the posterior means of each parameters and their highest density intervals (HDIs, see Kruschke 2011)
#- Denotes a different symbol for the parameters pertaining to each culture
#- This function is also called during the file export, and .eps and .jpeg's of the plot are saved
######################
plotresultsfuncbutton <- function() {
  guidat <- get("guidat", pkg_globals);
  plotresultsfunc(cctfit=get(guidat$varname, pkg_globals),gui=TRUE)
}

plotresultsfunc <- function(cctfit,saveplots=0,savedir="",gui=FALSE) {
  
  hdi <- function(sampleVec, credMass=0.95){
    sortedPts = sort(sampleVec)
    ciIdxInc = floor(credMass * length( sortedPts) )
    nCIs = length( sortedPts) - ciIdxInc
    ciwidth = rep(0, nCIs)
    for(i in 1:nCIs){
      ciwidth[i] = sortedPts[i + ciIdxInc] - sortedPts[i]}
    HDImin = sortedPts[which.min(ciwidth)]
    HDImax = sortedPts[which.min(ciwidth)+ciIdxInc]
    HDIlim = c(HDImin,HDImax)
    return(HDIlim) 
  }
  
  Mode <- function(x) {ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
  cctfit$respmem <- apply(cctfit$BUGSoutput$sims.list$Om[,],2,Mode)
  
  if(saveplots==1){jpeg(file.path(gsub(".Rdata","results.jpg",savedir)),width = 6, height = 6, units = "in", pointsize = 12,quality=100,res=400)}
  if(saveplots==2){postscript(file=file.path(gsub(".Rdata","results.eps",savedir)), onefile=FALSE, horizontal=FALSE, width = 6, height = 6, paper="special", family="Times")}
  
  if(cctfit$whmodel=="GCM"){
    par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(2,2))
    
    sym <- c(21,22, 23, 24, 25, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14) # possible symbols: 
    bgcol <- c("black",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
    
    #if(cctfit$V==1){bgcol <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)}
    
    if(cctfit$V==1){
      if(length(dim(cctfit$BUGSoutput$sims.list[["Z"]])) < 3){
        cctfit$BUGSoutput$sims.list[["Z"]] <- array(cctfit$BUGSoutput$sims.list[["Z"]],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
      plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Truth (",Z[vk],")")),xlab="Item",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=21,bg="white")
      for(i in 1:dim(cctfit$BUGSoutput$sims.list[["Z"]])[3]){
        points(apply(cctfit$BUGSoutput$sims.list[["Z"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",Z[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg=bgcol[i])
      }}else{
        plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Truth (",Z[vk],") Per Culture")),xlab="Item",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=21,bg="white")
        for(i in 1:dim(cctfit$BUGSoutput$sims.list[["Z"]])[3]){
          points(apply(cctfit$BUGSoutput$sims.list[["Z"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",Z[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg=bgcol[i])
        }
      }
    
    if(cctfit$itemdiff==TRUE){
      if(cctfit$V==1){
        if(length(dim(cctfit$BUGSoutput$sims.list[["lam"]])) < 3){
          cctfit$BUGSoutput$sims.list[["lam"]] <- array(cctfit$BUGSoutput$sims.list[["lam"]],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
        hditmp <- apply(cctfit$BUGSoutput$sims.list[["lam"]][,,1],2,hdi)
        plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],")")),xlab="Item",ylim=c(0,1),ylab="Posterior Mean Value",las=1,pch=21,bg="white")
        for(i in 1:dim(cctfit$BUGSoutput$sims.list[["lam"]])[3]){
          points(apply(cctfit$BUGSoutput$sims.list[["lam"]][,,i],c(2),mean),xlab=expression(paste("Item Difficulty (",lambda[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg=bgcol[i])
        }
        segments(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,])
        arrows(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,],code=3,angle=90,length=.025)
      }else{plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],") Per Culture")),xlab="Item",ylim=c(0,1),ylab="Posterior Mean Value",las=1,pch=21)
            for(i in 1:dim(cctfit$BUGSoutput$sims.list[["lam"]])[3]){
              points(apply(cctfit$BUGSoutput$sims.list[["lam"]][,,i],c(2),mean),xlab=expression(paste("Item Difficulty (",lambda[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg=bgcol[i])
            }
      }
    }else{
      plot(c(1:cctfit$m),rep(.5,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],")")),xlab="Item",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[1],bg="white",col="grey")
      points(c(par("usr")[1],par("usr")[2]),c(par("usr")[3],par("usr")[4]),type="l",col="grey")
      points(c(par("usr")[1],par("usr")[2]),c(par("usr")[4],par("usr")[3]),type="l",col="grey")
    }
    
    plot(cctfit$BUGSoutput$mean$th,main=expression(paste("Respondent Competency (",theta[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[cctfit$respmem],bg=bgcol[cctfit$respmem])
    hditmp <- apply(cctfit$BUGSoutput$sims.list$th,2,hdi)
    segments(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,])
    arrows(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,],code=3,angle=90,length=.025)
    
    plot(cctfit$BUGSoutput$mean$g,main=expression(paste("Respondent Guessing Bias (",g[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[cctfit$respmem],bg=bgcol[cctfit$respmem])
    hditmp <- apply(cctfit$BUGSoutput$sims.list$g,2,hdi)
    segments(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,])
    arrows(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,],code=3,angle=90,length=.025)
  }
  
  if(cctfit$whmodel=="LTRM"){
    invlogit <- function(x){x <- 1 / (1 + exp(-x)); return(x)}
    
    par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(2,2))
    
    sym <- c(21,22, 23, 24, 25, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14) # possible symbols: 
    bgcol <- c("black",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
    
    if(cctfit$C==2){useinvlogit <- TRUE}else{useinvlogit <- FALSE}
    
    if(useinvlogit==TRUE){
      #tmp1 <- cctfit$BUGSoutput$sims.list[["T"]]; cctfit$BUGSoutput$sims.list[["T"]] <- invlogit(cctfit$BUGSoutput$sims.list[["T"]])
      #tmp2 <- cctfit$BUGSoutput$sims.list[["gam"]]; cctfit$BUGSoutput$sims.list[["gam"]] <- invlogit(cctfit$BUGSoutput$sims.list[["gam"]])
      #tmp3 <- cctfit$BUGSoutput$sims.list[["b"]]; cctfit$BUGSoutput$sims.list[["b"]] <- invlogit(cctfit$BUGSoutput$sims.list[["b"]])
      tmp1 <- cctfit$BUGSoutput$sims.list[["T"]]; cctfit$BUGSoutput$sims.list[["T"]] <- invlogit(cctfit$BUGSoutput$sims.list[["T"]])
      tmp2 <- cctfit$BUGSoutput$sims.list[["gam"]]; cctfit$BUGSoutput$sims.list[["gam"]] <- invlogit(cctfit$BUGSoutput$sims.list[["gam"]])
      tmp3 <- cctfit$BUGSoutput$sims.list[["b"]]; cctfit$BUGSoutput$sims.list[["b"]] <- invlogit(cctfit$BUGSoutput$sims.list[["b"]])
    }
    
    if(cctfit$V==1){
      newm <- ceiling(cctfit$m*1.092)
      if(length(dim(cctfit$BUGSoutput$sims.list[["T"]])) <3){
        cctfit$BUGSoutput$sims.list[["T"]] <- array(cctfit$BUGSoutput$sims.list[["T"]],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
      }
      if(length(dim(cctfit$BUGSoutput$sims.list[["gam"]])) <3){
        cctfit$BUGSoutput$sims.list[["gam"]] <- array(cctfit$BUGSoutput$sims.list[["gam"]],c(cctfit$BUGSoutput$n.sims,cctfit$C-1,1))
      }
      hditmp <- apply(cctfit$BUGSoutput$sims.list[["T"]][,,1],2,hdi)
      if(useinvlogit==TRUE){
        plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[vk],") and Thresholds (",gamma[vc],")")),xlab="Item",ylim=c(0,1),ylab="Posterior Mean Invlogit Value",las=1,pch=21,bg="white",axes=FALSE)
      }else{plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[vk],") and Thresholds (",gamma[vc],")")),xlab="Item",ylim=c(min(hditmp,min(apply(cctfit$BUGSoutput$sims.list[["gam"]],c(2,3),mean))),max(hditmp,max(apply(cctfit$BUGSoutput$sims.list[["gam"]],c(2,3),mean)))),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
      }
      for(i in 1:dim(cctfit$BUGSoutput$sims.list[["T"]])[3]){
        points(apply(cctfit$BUGSoutput$sims.list[["T"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",T[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg=bgcol[i])
        
        if(cctfit$C==2){
          text(newm,mean(cctfit$BUGSoutput$sims.list[["gam"]][,,i]),labels=sapply(c(1:(cctfit$C-1)),function(x) as.expression(substitute(list(gamma[x]),list(x=x)))))
          segments(1,mean(cctfit$BUGSoutput$sims.list[["gam"]][,,i]),max(cctfit$m+1,newm-2),mean(cctfit$BUGSoutput$sims.list[["gam"]][,,i]),lty=2);
        }
        else{
          text(newm,apply(cctfit$BUGSoutput$sims.list[["gam"]][,,i],c(2),mean),labels=sapply(c(1:(cctfit$C-1)),function(x) as.expression(substitute(list(gamma[x]),list(x=x)))))
          segments(1,apply(cctfit$BUGSoutput$sims.list[["gam"]][,,i],c(2),mean),max(cctfit$m+1,newm-2),apply(cctfit$BUGSoutput$sims.list[["gam"]][,,i],c(2),mean),lty=2);
        }
        segments(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,])
        arrows(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,],code=3,angle=90,length=.025)
        box()
        axis(2, labels = TRUE,las=1)
        axis(side = 1,labels=TRUE)
        axis(side = 1, at = newm, labels = expression(gamma[c]) )
      }}else{
        
        newm <- ceiling(cctfit$m*1.14)
        
        if(cctfit$C==2){
          
          if(useinvlogit==TRUE){
            plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[vk],") and Thresholds (",gamma[vc],") Per Culture")),xlab="Item",ylim=c(0,1),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
          }else{plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[vk],") and Thresholds (",gamma[vc],") Per Culture")),xlab="Item",ylim=c(min(min(apply(cctfit$BUGSoutput$sims.list[["T"]],c(2,3),mean))
                                                                                                                                                            ,min(apply(cctfit$BUGSoutput$sims.list[["gam"]],c(3),mean))),max(max(apply(cctfit$BUGSoutput$sims.list[["T"]],c(2,3),mean))
                                                                                                                                                                                                                             ,max(apply(cctfit$BUGSoutput$sims.list[["gam"]],c(3),mean)))),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
          }
          
          for(i in 1:dim(cctfit$BUGSoutput$sims.list[["T"]])[3]){
            points(apply(cctfit$BUGSoutput$sims.list[["T"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",T[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg=bgcol[i])
            text(cctfit$m+(((newm-cctfit$m)/dim(cctfit$BUGSoutput$sims.list[["T"]])[3])*i),mean(cctfit$BUGSoutput$sims.list[["gam"]][,,i]),labels=sapply(c(1:(cctfit$C-1)),function(x) as.expression(substitute(list(gamma[x]),list(x=x)))))
          }
        }else{
          if(useinvlogit==TRUE){
            plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[vk],") and Thresholds (",gamma[vc],") Per Culture")),xlab="Item",ylim=c(0,1),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
          }else{
            plot(1:newm,rep(NA,newm),main=expression(paste("Item Truth (",T[vk],") and Thresholds (",gamma[vc],") Per Culture")),xlab="Item",ylim=c(min(min(apply(cctfit$BUGSoutput$sims.list[["T"]],c(2,3),mean))
                                                                                                                                                        ,min(apply(cctfit$BUGSoutput$sims.list[["gam"]][,,],c(2,3),mean))),max(max(apply(cctfit$BUGSoutput$sims.list[["T"]],c(2,3),mean))
                                                                                                                                                                                                                               ,max(apply(cctfit$BUGSoutput$sims.list[["gam"]][,,],c(2,3),mean)))),ylab="Posterior Mean Value",las=1,pch=21,bg="white",axes=FALSE)
          }
          
          for(i in 1:dim(cctfit$BUGSoutput$sims.list[["T"]])[3]){
            points(apply(cctfit$BUGSoutput$sims.list[["T"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",T[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg=bgcol[i])
            text(cctfit$m+(((newm-cctfit$m)/dim(cctfit$BUGSoutput$sims.list[["T"]])[3])*i),apply(cctfit$BUGSoutput$sims.list[["gam"]][,,i],c(2),mean),labels=sapply(c(1:(cctfit$C-1)),function(x) as.expression(substitute(list(gamma[x]),list(x=x)))))
          }
        }
        box()
        axis(2, labels = TRUE,las=1)
        axis(side=1,at=axTicks(side=1)[axTicks(side=1)<= cctfit$m])
        axis(side = 1, at = cctfit$m+(((newm-cctfit$m)/dim(cctfit$BUGSoutput$sims.list[["T"]])[3])*(1:2)), 
             labels = sapply(1:2,function(x) as.expression(substitute(list(gamma[x*c]),list(x=x)))) )
      }
    
    if(cctfit$itemdiff==TRUE){
      if(cctfit$V==1){
        if(length(dim(cctfit$BUGSoutput$sims.list[["lam"]])) < 3){
          cctfit$BUGSoutput$sims.list[["lam"]] <- array(cctfit$BUGSoutput$sims.list[["lam"]],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
        hditmp <- apply(cctfit$BUGSoutput$sims.list[["lam"]][,,1],2,hdi)
        plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],")")),xlab="Item",ylim=c(min(hditmp),max(hditmp)),ylab="Posterior Mean Log Value",las=1,pch=21,bg="white")
        for(i in 1:dim(cctfit$BUGSoutput$sims.list[["lam"]])[3]){
          points(apply(cctfit$BUGSoutput$sims.list[["lam"]][,,i],c(2),mean),xlab=expression(paste("Item Difficulty (",lambda[vk],") Per Cluster")),ylab="Posterior Mean Log Value",las=1,ylim=c(0,1),pch=sym[i],bg=bgcol[i])
        }
        segments(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,])
        arrows(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,],code=3,angle=90,length=.025)
      }else{plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],") Per Culture")),xlab="Item",ylim=c(min(apply(cctfit$BUGSoutput$sims.list[["lam"]],c(2,3),mean)),max(apply(cctfit$BUGSoutput$sims.list[["lam"]],c(2,3),mean))),ylab="Posterior Mean Value",las=1,pch=21,bg="white")
            for(i in 1:dim(cctfit$BUGSoutput$sims.list[["lam"]])[3]){
              points(apply(cctfit$BUGSoutput$sims.list[["lam"]][,,i],c(2),mean),xlab=expression(paste("Item Difficulty (",lambda[vk],") Per Cluster")),ylab="Posterior Mean Log Value",las=1,ylim=c(0,1),pch=sym[i],bg=bgcol[i])
            }
      }
    }else{
      plot(c(1:cctfit$m),rep(.5,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],")")),xlab="Item",ylab="Posterior Mean Log Value",las=1,ylim=c(0,1),pch=sym[1],bg="white",col="grey")
      points(c(par("usr")[1],par("usr")[2]),c(par("usr")[3],par("usr")[4]),type="l",col="grey")
      points(c(par("usr")[1],par("usr")[2]),c(par("usr")[4],par("usr")[3]),type="l",col="grey")
    }
    
    
    hditmp <- apply(cctfit$BUGSoutput$sims.list$E,2,hdi)
    plot(cctfit$BUGSoutput$mean$E,main=expression(paste("Respondent Standard Error (",E[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,max(hditmp)),pch=sym[cctfit$respmem],bg=bgcol[cctfit$respmem])
    segments(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,])
    arrows(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,],code=3,angle=90,length=.025)
    
    if(useinvlogit==TRUE){
      plot(invlogit(cctfit$BUGSoutput$mean$b)-.5,main=expression(paste("Respondent Shift Bias (",b[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(-.5,.5),pch=sym[cctfit$respmem],bg=bgcol[cctfit$respmem])
      hditmp <- apply(cctfit$BUGSoutput$sims.list$b-.5,2,hdi) # because this is invprobit transformed still
      segments(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,])
      arrows(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,],code=3,angle=90,length=.025)
    }else{
      loga <- apply(log(cctfit$BUGSoutput$sims.list[["a"]]),2,mean)
      plot(-5000, -5000, main=expression(paste("Category Usage Bias (",a[i]," and ",b[i],")")),xlim=c(min(-1.6,-abs(cctfit$BUGSoutput$mean$b)),max(1.6,cctfit$BUGSoutput$mean$b)), ylim=c(min(-1,loga),max(1,loga)),xlab=expression(paste("Respondent Shift Bias (",b[i],")")),
           ylab=expression(paste("Log Respondent Scale Bias (",a[i],")")),las=1,pch=21)
      points(cbind(cctfit$BUGSoutput$mean$b,loga),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,max(hditmp)),pch=sym[cctfit$respmem],bg=bgcol[cctfit$respmem])
      segments(.625*par("usr")[1],0,.625*par("usr")[2],0)
      segments(0,.725*par("usr")[3],0,.725*par("usr")[4])
      text(0,.875*par("usr")[3],"Middle Categories",cex=.8)
      text(0,.875*par("usr")[4],"Outer Categories",cex=.8)
      text(.775*par("usr")[1],0,"Left \n Categories",cex=.8)
      text(.775*par("usr")[2],0,"Right \n Categories",cex=.8)
    }
    
    
    if(useinvlogit==TRUE){
      cctfit$BUGSoutput$sims.list[["T"]] <- tmp1 
      cctfit$BUGSoutput$sims.list[["gam"]] <- tmp2
      cctfit$BUGSoutput$sims.list[["b"]] <- tmp3
    }
    
  }
  
  if(cctfit$whmodel=="CRM"){
    
    par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(2,2))
    
    sym <- c(21,22, 23, 24, 25, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14) # possible symbols: 
    bgcol <- c("black",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
    
    if(cctfit$V==1){
      if(length(dim(cctfit$BUGSoutput$sims.list[["T"]])) < 3){
        cctfit$BUGSoutput$sims.list[["T"]] <- array(cctfit$BUGSoutput$sims.list[["T"]],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
      hditmp <- apply(cctfit$BUGSoutput$sims.list[["T"]][,,1],2,hdi)
      plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Truth (",T[vk],")")),xlab="Item",ylim=c(min(hditmp),max(hditmp)),ylab="Posterior Mean Value",las=1,pch=21,bg="white")
    }else{plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Truth (",T[vk],") Per Culture")),xlab="Item",ylim=c(min(apply(cctfit$BUGSoutput$sims.list[["T"]],c(2,3),mean)),max(apply(cctfit$BUGSoutput$sims.list[["T"]],c(2,3),mean))),ylab="Posterior Mean Value",las=1,pch=21,bg="white")}
    for(i in 1:dim(cctfit$BUGSoutput$sims.list[["T"]])[3]){
      points(apply(cctfit$BUGSoutput$sims.list[["T"]][,,i],c(2),mean),xlab=expression(paste("Item Truth (",T[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg=bgcol[i])
    }
    if(cctfit$V==1){
      segments(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,])
      arrows(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,],code=3,angle=90,length=.025)
    }
    
    if(cctfit$itemdiff==TRUE){
      if(cctfit$V==1){
        if(length(dim(cctfit$BUGSoutput$sims.list[["lam"]])) < 3){
          cctfit$BUGSoutput$sims.list[["lam"]] <- array(cctfit$BUGSoutput$sims.list[["lam"]],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
        hditmp <- apply(cctfit$BUGSoutput$sims.list[["lam"]][,,1],2,hdi)
        plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],")")),xlab="Item",ylim=c(min(hditmp),max(hditmp)),ylab="Posterior Mean Value",las=1,pch=21,bg="white")
        for(i in 1:dim(cctfit$BUGSoutput$sims.list[["lam"]])[3]){
          points(apply(cctfit$BUGSoutput$sims.list[["lam"]][,,i],c(2),mean),xlab=expression(paste("Item Difficulty (",lambda[vk],") Per Cluster")),ylab="Posterior Mean Value",las=1,ylim=c(0,1),pch=sym[i],bg=bgcol[i])
        }
        segments(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,])
        arrows(1:cctfit$m,hditmp[1,], 1:cctfit$m, hditmp[2,],code=3,angle=90,length=.025)
      }else{plot(1:cctfit$m,rep(NA,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],") Per Culture")),xlab="Item",ylim=c(min(apply(cctfit$BUGSoutput$sims.list[["lam"]],c(2,3),mean)),max(apply(cctfit$BUGSoutput$sims.list[["lam"]],c(2,3),mean))),ylab="Posterior Mean Log Value",las=1,pch=21,bg="white")
            for(i in 1:dim(cctfit$BUGSoutput$sims.list[["lam"]])[3]){
              points(apply(cctfit$BUGSoutput$sims.list[["lam"]][,,i],c(2),mean),xlab=expression(paste("Item Difficulty (",lambda[vk],") Per Cluster")),ylab="Posterior Mean Log Value",las=1,ylim=c(0,1),pch=sym[i],bg=bgcol[i])
            }
      }
    }else{
      plot(c(1:cctfit$m),rep(.5,cctfit$m),main=expression(paste("Item Difficulty (",lambda[vk],")")),xlab="Item",ylab="Posterior Mean Log Value",las=1,ylim=c(0,1),pch=sym[1],bg="white",col="grey")
      points(c(par("usr")[1],par("usr")[2]),c(par("usr")[3],par("usr")[4]),type="l",col="grey")
      points(c(par("usr")[1],par("usr")[2]),c(par("usr")[4],par("usr")[3]),type="l",col="grey")
    }
    
    hditmp <- apply(cctfit$BUGSoutput$sims.list$E,2,hdi)
    plot(cctfit$BUGSoutput$mean$E,main=expression(paste("Respondent Standard Error (",E[i],")")),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(min(hditmp),max(hditmp)),pch=sym[cctfit$respmem],bg=bgcol[cctfit$respmem])
    segments(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,])
    arrows(1:cctfit$n,hditmp[1,], 1:cctfit$n, hditmp[2,],code=3,angle=90,length=.025)
    #text(cctfit$n/4.25,1.3*par("usr")[3],"More Knowledgeable",cex=.8)
    #text(cctfit$n/4.25,.9*par("usr")[4],"Less Knowledgeable",cex=.8)
    
    
    loga <- apply(log(cctfit$BUGSoutput$sims.list[["a"]]),2,mean)
    plot(-5000, -5000, main=expression(paste("Category Usage Bias (",a[i]," and ",b[i],")")),xlim=c(min(-1.6,-abs(cctfit$BUGSoutput$mean$b)),max(1.6,cctfit$BUGSoutput$mean$b)), ylim=c(min(-1,loga),max(1,loga)),xlab=expression(paste("Respondent Shift Bias (",b[i],")")),
         ylab=expression(paste("Log Respondent Scale Bias (",a[i],")")),las=1,pch=21)
    points(cbind(cctfit$BUGSoutput$mean$b,loga),xlab="Respondent",ylab="Posterior Mean Value",las=1,ylim=c(0,max(hditmp)),pch=sym[cctfit$respmem],bg=bgcol[cctfit$respmem])
    segments(.625*par("usr")[1],0,.625*par("usr")[2],0)
    segments(0,.725*par("usr")[3],0,.725*par("usr")[4])
    text(0,.875*par("usr")[3],"Middle Categories",cex=.8)
    text(0,.875*par("usr")[4],"Outer Categories",cex=.8)
    text(.775*par("usr")[1],0,"Left \n Categories",cex=.8)
    text(.775*par("usr")[2],0,"Right \n Categories",cex=.8)
  }
  
  
  if(saveplots==1 || saveplots==2){dev.off()}
  
  saveplots <- 0
}

######################
#Function for the 'Run Checks' Button
#- Calculates the posterior predictive data for the model that was applied
#- Calculates 2 important posterior predictive checks and plots them
#- The VDI check percentiles are reported in the R console, and saved to cctfit$VDIperc
#- Performs factor analysis using fa() and polychoric() (if selected)
#- This function is also called in the file export and saves .eps and .jpegs of the plot
######################
ppcfuncbutton <- function(){
  guidat <- get("guidat", pkg_globals)
  cctfit <- get(guidat$varname, pkg_globals)
  
   recalc <- FALSE
   if(cctfit$whmodel=="LTRM" && cctfit$checksrun==TRUE){ 
      if(cctfit$polycor!=as.logical(as.numeric(tclvalue(guidat$polyvar)))){
        recalc <- TRUE
      }
    }
  
  assign(guidat$varname,
         ppcfunc(cctfit=cctfit,gui=TRUE,rerunchecks=recalc,polych=as.logical(as.numeric(tclvalue(guidat$polyvar)))),
         pkg_globals
  )
}

ppcfunc <- function(cctfit,saveplots=0,savedir="",gui=FALSE,polych=FALSE,rerunchecks=FALSE) {
  if(gui==TRUE){guidat <- get("guidat", pkg_globals)}
  
  if(cctfit$checksrun==FALSE || rerunchecks==TRUE){
    message("\n ...One moment, calculating posterior predictive checks")
    
    usesubset <- 1
    
    if(usesubset==1){
      subsetsize <- min(500,cctfit$BUGSoutput$n.sims)
      indices <- sample(cctfit$BUGSoutput$n.sims,min(subsetsize,cctfit$BUGSoutput$n.sims))
      
      if(cctfit$whmodel=="GCM"){
        if(cctfit$V==1){
          cctfit$V <- 1;
          cctfit$BUGSoutput$sims.list$Om <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
          if(length(dim(cctfit$BUGSoutput$sims.list$Z)) < 3){
            cctfit$BUGSoutput$sims.list$Z <- array(cctfit$BUGSoutput$sims.list[["Z"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
          if(cctfit$itemdiff==TRUE){
            if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
              cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
          }
        }
        if(cctfit$itemdiff==FALSE){
          if(cctfit$V==1){
            cctfit$BUGSoutput$sims.list$lam <- array(.5,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
          }
          if(cctfit$V > 1){
            cctfit$BUGSoutput$sims.list$lam <- array(.5,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
          } }
        
        cctfit$ppY <- array(NA, c(cctfit$n,cctfit$m,subsetsize));
        
        for(samp in indices){
          tautmp <- (cctfit$BUGSoutput$sims.list[["th"]][samp,]*t(1-cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])) / 
            ( (cctfit$BUGSoutput$sims.list[["th"]][samp,]*t(1-cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])) + 
                ((1-cctfit$BUGSoutput$sims.list[["th"]][samp,])*t(cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])) ) 
          
          cctfit$ppY[,,which(indices==samp)] <- matrix(rbinom((cctfit$n*cctfit$m),1,
                                                                (t(cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])*tautmp ) +
                                                                  t(t(1-tautmp)%*%diag(cctfit$BUGSoutput$sims.list[["g"]][samp,])) ), cctfit$n,cctfit$m)
        }
        
        if(cctfit$mval==TRUE){
          for(i in 1:dim(cctfit$datob$thena)[1]){cctfit$data[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2]] <- Mode(cctfit$ppY[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2],])}
          cctfit$MVest <- cbind(cctfit$datob$thena,cctfit$data[cctfit$datob$thena])
          colnames(cctfit$MVest) <- c("Pers","Item","Resp")
        }
        
        cctfit$mean$ppY <- rowMeans(cctfit$ppY[,,],dims=2)
        
        leigv <- 12
        options(warn=-3)
        ind <- indices 
        eigv <- matrix(-1,length(ind),leigv)
        tmp <- apply(cctfit$ppY,c(1,3),function(x) t(x))
        tmp2 <- apply(tmp[,,],3,function(x) cor(x))
        tmp3 <- array(tmp2,c(cctfit$n,cctfit$n,subsetsize))
        for(i in 1:length(ind)){
          suppressMessages(try(eigv[i,] <- fa(tmp3[,,i])$values[1:leigv],silent=TRUE))}
        wch <- -which(eigv[,1]==-1); 
        if(length(wch)==0){cctfit$ppeig <-  eigv}else{cctfit$ppeig <-  eigv[-which(eigv[,1]==-1),]}
        
        if(sum(apply(cctfit$data,1,function(x) sd(x,na.rm=TRUE))==0) > 0){
          tmp <- cctfit$data
          if(sum(apply(cctfit$datob$dat,1,function(x) sd(x,na.rm=TRUE))==0)==1){
            tmp[which(apply(tmp,1,function(x) sd(x,na.rm=TRUE))==0),1] <-  min(tmp[which(apply(tmp,1,function(x) sd(x,na.rm=TRUE))==0),1]+.01,.99)
          }else{
            tmp[which(apply(tmp,1,function(x) sd(x,na.rm=TRUE))==0),1] <- sapply(apply(tmp[which(apply(tmp,1,function(x) sd(x,na.rm=TRUE))==0),],1,mean),function(x) min(x+.01,.99))
          }
          cctfit$dateig <- suppressMessages(fa(cor(t(tmp)))$values[1:leigv])
        }else{
          cctfit$dateig <-  suppressMessages(fa(cor(t(cctfit$data)))$values[1:leigv])
        }
        
        options(warn=0)
        
        if(cctfit$V==1){
          varvec <- matrix(-1,length(ind),dim(cctfit$ppY)[2])
          vdi <- matrix(-1,dim(cctfit$ppY)[3],1)
          for(i in 1:length(ind)){
            varvec[i,] <- apply(cctfit$ppY[,,i],2,var)
          }
          vdi <- apply(varvec,1,var)
          cctfit$ppVDI <- vdi
          cctfit$datVDI <- var(apply(cctfit$data,2,var))
        }
        
        options(warn=-3)
        if(cctfit$V>1){
          varvec <- array(NA,c(length(ind),dim(cctfit$ppY)[2],cctfit$V))
          vdi <- matrix(NA,dim(cctfit$ppY)[3],cctfit$V)
          cctfit$datVDI  <- array(-1,c(1,cctfit$V))
          for(v in 1:cctfit$V){
            for(i in 1:length(ind)){
              if(sum(cctfit$BUGSoutput$sims.list$Om[ind[i],]==v)>1){
                suppressMessages(try(varvec[i,,v] <- apply(cctfit$ppY[cctfit$BUGSoutput$sims.list$Om[ind[i],]==v,,i],2,var),silent=TRUE))
              }else{
                suppressMessages(try(varvec[i,,v] <- var(cctfit$ppY[cctfit$BUGSoutput$sims.list$Om[ind[i],]==v,,i]),silent=TRUE))
              }
            }
            if(sum(cctfit$respmem==v)>1){
              cctfit$datVDI[v] <- var(apply(cctfit$data[cctfit$respmem==v,],2,var))
            }else{
              cctfit$datVDI[v] <- var(cctfit$data[cctfit$respmem==v,])
            }
          }
          vdi <- apply(varvec,c(1,3),var)
          cctfit$ppVDI <- vdi
          
          options(warn=0)
          if(sum(apply(cctfit$ppVDI,2,function(x) all(x==0,na.rm=TRUE)))>=1){
            for(i in which(apply(cctfit$ppVDI,2,function(x) all(x==0)))){
              cctfit$ppVDI[,i] <- NA
            }
          }
          for(v in 1:dim(cctfit$ppVDI)[2]){
            if(identical(unique(cctfit$ppVDI[,v]),c(NA,0)) || identical(unique(cctfit$ppVDI[,v]),c(0,NA))){cctfit$ppVDI[,v] <- NA}
          }
        }
        rm(vdi, varvec)
        rm(eigv,leigv,tmp,tmp2,tmp3,ind);
        
        cctfit$checksrun <- TRUE
      }
      
      if(cctfit$whmodel=="LTRM"){
        if(cctfit$V==1){
          cctfit$V <- 1;
          cctfit$BUGSoutput$sims.list$Om <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
          if(length(dim(cctfit$BUGSoutput$sims.list[["T"]])) < 3){
            cctfit$BUGSoutput$sims.list$T <- array(cctfit$BUGSoutput$sims.list[["T"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
          if(length(dim(cctfit$BUGSoutput$sims.list[["gam"]])) < 3){
            cctfit$BUGSoutput$sims.list$gam <- array(cctfit$BUGSoutput$sims.list[["gam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$C-1,1))}
          if(cctfit$itemdiff==TRUE){
            if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
              cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
          }
        }
        if(cctfit$itemdiff==FALSE){
          if(cctfit$V==1){
            cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
          }
          if(cctfit$V > 1){
            cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
          }
        }
        cctfit$ppY <- array(NA, c(cctfit$n,cctfit$m,subsetsize));
        cctfit$ppX <- array(NA, c(cctfit$n,cctfit$m,subsetsize));
        
        for(samp in indices){
          tautmp <- (cctfit$BUGSoutput$sims.list[["E"]][samp,]*t(exp(cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])))^-2
          
          cctfit$ppX[,,which(indices==samp)] <- matrix(
            rnorm((cctfit$n*cctfit$m),t(cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]]),(tautmp^(-.5))), 
            cctfit$n,cctfit$m)
          
          deltatmp <- cctfit$BUGSoutput$sims.list[["a"]][samp,]*t(cctfit$BUGSoutput$sims.list[["gam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])+cctfit$BUGSoutput$sims.list[["b"]][samp,]
          
          rc <- which(cctfit$ppX[,,which(indices==samp)] < deltatmp[,1],arr.ind=TRUE)
          cctfit$ppY[cbind(rc,array(which(indices==samp),dim(rc)[1]))] <- 1 
          for(c in 1:(cctfit$C-1)){
            rc <- which(cctfit$ppX[,,which(indices==samp)] > deltatmp[,c],arr.ind=TRUE)
            cctfit$ppY[cbind(rc,array(which(indices==samp),dim(rc)[1]))] <- (c+1) }
          
        }
        cctfit$mean$ppY <- rowMeans(cctfit$ppY[,,],dims=2)
        
        if(cctfit$mval==TRUE){
          for(i in 1:dim(cctfit$datob$thena)[1]){cctfit$data[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2]] <- Mode(cctfit$ppY[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2],])}
          cctfit$MVest <- cbind(cctfit$datob$thena,cctfit$data[cctfit$datob$thena])
          colnames(cctfit$MVest) <- c("Pers","Item","Resp")
        }
        
        leigv <- 12
        eigv <- matrix(-1,dim(cctfit$ppY)[3],leigv)
        options(warn=-3)
        ind <- sample(indices,min(cctfit$BUGSoutput$n.sims,250,subsetsize)) #250 is the number of samples
        tmp <- apply(cctfit$ppY,c(1,3),function(x) t(x))
        if(polych==TRUE){
          message("    note: utilizing polychoric correlations (time intensive)");
          tmp2 <- apply(tmp[,,],3,function(x) suppressMessages(polychoric(x,global=FALSE))$rho);
          tmp2 <- sapply(tmp2,function(x) array(x,c(cctfit$n,cctfit$n))) }else{
            tmp2 <- apply(tmp[,,],3,function(x) cor(x))}
        tmp3 <- array(tmp2,c(cctfit$n,cctfit$n,length(ind)))
        for(i in 1:length(ind)){
          suppressMessages(try(eigv[i,] <- fa(tmp3[,,i])$values[1:leigv],silent=TRUE))}
        wch <- -which(eigv[,1]==-1); 
        if(length(wch)==0){cctfit$ppeig <-  eigv}else{cctfit$ppeig <-  eigv[-which(eigv[,1]==-1),]}
        if(polych==TRUE){
          cctfit$dateig <-  suppressMessages(fa(polychoric(t(cctfit$data),global=FALSE)$rho)$values[1:leigv])}else{
            cctfit$dateig <-  suppressMessages(fa(cor(t(cctfit$data)))$values[1:leigv])
          }
        
        options(warn=0)
        
        if(cctfit$V==1){
          varvec <- matrix(NA,length(ind),dim(cctfit$ppY)[2])
          vdi <- matrix(NA,dim(cctfit$ppY)[3],1)
          for(i in 1:length(ind)){
            varvec[i,] <- apply(cctfit$ppY[,,i],2,var)
          }
          vdi <- apply(varvec,1,var)
          cctfit$ppVDI <- vdi
          cctfit$datVDI <- var(apply(cctfit$data,2,var))
        }
        options(warn=-3)
        if(cctfit$V>1){
          varvec <- array(NA,c(length(ind),dim(cctfit$ppY)[2],cctfit$V))
          vdi <- matrix(NA,dim(cctfit$ppY)[3],cctfit$V)
          cctfit$datVDI  <- array(-1,c(1,cctfit$V))
          for(v in 1:cctfit$V){
            for(i in 1:length(ind)){
              if(sum(cctfit$BUGSoutput$sims.list$Om[ind[i],]==v)>1){
                suppressMessages(try(varvec[i,,v] <- apply(cctfit$ppY[cctfit$BUGSoutput$sims.list$Om[ind[i],]==v,,i],2,var),silent=TRUE))
              }else{
                suppressMessages(try(varvec[i,,v] <- var(cctfit$ppY[cctfit$BUGSoutput$sims.list$Om[ind[i],]==v,,i]),silent=TRUE))
              }
            }
            if(sum(cctfit$respmem==v)>1){
              cctfit$datVDI[v] <- var(apply(cctfit$data[cctfit$respmem==v,],2,var))
            }else{
              cctfit$datVDI[v] <- var(cctfit$data[cctfit$respmem==v,])
            }
          }
          vdi <- apply(varvec,c(1,3),var)
          cctfit$ppVDI <- vdi
          
          options(warn=0)
          if(sum(apply(cctfit$ppVDI,2,function(x) all(x==0,na.rm=TRUE)))>=1){
            for(i in which(apply(cctfit$ppVDI,2,function(x) all(x==0)))){
              cctfit$ppVDI[,i] <- NA
            }
          }
          for(v in 1:dim(cctfit$ppVDI)[2]){
            if(identical(unique(cctfit$ppVDI[,v]),c(NA,0)) || identical(unique(cctfit$ppVDI[,v]),c(0,NA))){cctfit$ppVDI[,v] <- NA}
          }
        }
        cctfit$polycor <- polych
        
        rm(vdi, varvec)
        rm(eigv,leigv,tmp,tmp2,tmp3,ind);
        
        cctfit$checksrun <- TRUE
      }
      
      if(cctfit$whmodel=="CRM"){
        invlogit <- function(x){x <- 1 / (1 + exp(-x)); return(x)}
        logit <- function(x){x <- log(x/(1-x)); return(x)}
        
        usetransform <- 0
        
        version2 <- 1
        
        if(cctfit$V==1){
          cctfit$V <- 1;
          cctfit$BUGSoutput$sims.list$Om <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
          if(length(dim(cctfit$BUGSoutput$sims.list[["T"]])) < 3){
            cctfit$BUGSoutput$sims.list$T <- array(cctfit$BUGSoutput$sims.list[["T"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
          if(cctfit$itemdiff==TRUE){
            if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
              cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
          }
        }
        if(cctfit$itemdiff==FALSE){
          if(cctfit$V==1){
            cctfit$BUGSoutput$sims.list$lam <- array(0,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
          }
          if(cctfit$V > 1){
            cctfit$BUGSoutput$sims.list$lam <- array(0,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
          }
        }
        
        cctfit$ppY <- array(NA, c(cctfit$n,cctfit$m,subsetsize));
        cctfit$ppX <- array(NA, c(cctfit$n,cctfit$m,subsetsize));
        
        for(samp in indices){
          tautmp <- (cctfit$BUGSoutput$sims.list[["E"]][samp,]*t(exp(cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])))^-2
          
          cctfit$ppY[,,which(indices==samp)] <- matrix(
            rnorm((cctfit$n*cctfit$m),(cctfit$BUGSoutput$sims.list[["a"]][samp,]*t(cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]]))+matrix(rep(cctfit$BUGSoutput$sims.list[["b"]][samp,],cctfit$m),cctfit$n,cctfit$m),
                  (cctfit$BUGSoutput$sims.list[["a"]][samp,]*(tautmp^(-.5)))), 
            cctfit$n,cctfit$m)
          
          cctfit$ppX[,,which(indices==samp)] <- matrix(
            rnorm((cctfit$n*cctfit$m),t(cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]]),tautmp^(-.5)),c(cctfit$n,cctfit$m))
        }
        
        if(usetransform==1){
          if(sum(cctfit$ppX < logit(.001)) > 0){cctfit$ppX[cctfit$ppX < logit(.001)] <- logit(.001)}  
          if(sum(cctfit$ppY < logit(.001)) > 0){cctfit$ppY[cctfit$ppY < logit(.001)] <- logit(.001)}  
          
          if(sum(cctfit$ppX > logit(.999)) > 0){cctfit$ppX[cctfit$ppX > logit(.999)] <- logit(.999)}  
          if(sum(cctfit$ppY > logit(.999)) > 0){cctfit$ppY[cctfit$ppY > logit(.999)] <- logit(.999)}  
        }
        
        if(usetransform==2){
          if(sum(cctfit$ppX < logit(.0001)) > 0){cctfit$ppX[cctfit$ppX < logit(.0001)] <- logit(.0001)}  
          if(sum(cctfit$ppY < logit(.0001)) > 0){cctfit$ppY[cctfit$ppY < logit(.0001)] <- logit(.0001)}  
          
          if(sum(cctfit$ppX > logit(.9999)) > 0){cctfit$ppX[cctfit$ppX > logit(.9999)] <- logit(.9999)}  
          if(sum(cctfit$ppY > logit(.9999)) > 0){cctfit$ppY[cctfit$ppY > logit(.9999)] <- logit(.9999)}  
        }
        
        cctfit$mean$ppY <- rowMeans(cctfit$ppY[,,],dims=2)
        
        if(cctfit$mval==TRUE){
          for(i in 1:dim(cctfit$datob$thena)[1]){cctfit$data[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2]] <- mean(cctfit$ppY[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2],])}
          cctfit$MVest <- cbind(cctfit$datob$thena,cctfit$data[cctfit$datob$thena])
          colnames(cctfit$MVest) <- c("Pers","Item","Resp")
        }
        
        leigv <- 12
        options(warn=-3)
        ind <- indices #500 is the number of samples
        eigv <- matrix(-1,length(ind),leigv)
        tmp <- apply(cctfit$ppY,c(1,3),function(x) t(x))
        tmp2 <- apply(tmp[,,],3,function(x) cor(x))
        tmp3 <- array(tmp2,c(cctfit$n,cctfit$n,subsetsize))
        for(i in 1:length(ind)){
          suppressMessages(try(eigv[i,] <- fa(tmp3[,,i])$values[1:leigv],silent=TRUE))}
        wch <- -which(eigv[,1]==-1); 
        if(length(wch)==0){cctfit$ppeig <-  eigv}else{cctfit$ppeig <-  eigv[-which(eigv[,1]==-1),]}
        cctfit$dateig <-  suppressMessages(fa(cor(t(cctfit$data)))$values[1:leigv])
        options(warn=0)
        
        #VDI Check
        
        usetransform <- 0
        
        if(usetransform==0){
          invlogit <- function(x){return(x)}
        }
        
        
        if(cctfit$V==1){
          vdi <- matrix(-1,dim(cctfit$ppY)[3],1)
          varvec <- array(-1,c(length(ind),dim(cctfit$ppY)[2]))
          
          for(i in 1:length(ind)){
            varvec[i,] <- apply(invlogit(cctfit$ppY[,,i]),2,var)
            if(version2==1){
              varvec[i,] <- apply(invlogit(cctfit$ppX[,,i]),2,var)
            } 
          }
          vdi <- apply(varvec,1,var)
          cctfit$ppVDI <- apply(varvec,1,var) 
          
          cctfit$datVDI <- var(apply(invlogit(cctfit$data),2,var))  
          
          if(version2==1){
            cctfit$datVDI <- var(apply(invlogit((array(rep(1/cctfit$BUGSoutput$mean$a,times=cctfit$m),c(cctfit$n,cctfit$m))*cctfit$data)-matrix(rep(cctfit$BUGSoutput$mean$b,cctfit$m),cctfit$n,cctfit$m)),2,var))  
          } 
          
        }
        if(cctfit$V>1){
          varvec <- array(NA,c(length(ind),dim(cctfit$ppY)[2],cctfit$V))
          cctfit$datVDI  <- array(NA,c(1,cctfit$V))
          options(warn=-3)
          datest <- t(sapply(1:cctfit$n,function(i) ((1/mean(cctfit$BUGSoutput$sims.list$a[cctfit$BUGSoutput$sims.list$Om[,i]==cctfit$respmem[i],i]))*
                                                       cctfit$data[i,])-(mean(cctfit$BUGSoutput$sims.list$b[cctfit$BUGSoutput$sims.list$Om[,i]==cctfit$respmem[i],i]))
          ))
          
          for(v in 1:cctfit$V){
            for(i in 1:length(ind)){
              if(sum(cctfit$BUGSoutput$sims.list$Om[ind[i],]==v)>1){
                suppressMessages(try(varvec[i,,v] <- apply(invlogit(cctfit$ppY[cctfit$BUGSoutput$sims.list$Om[ind[i],]==v,,i]),2,var),silent=TRUE))
                if(version2==1){
                  suppressMessages(try(varvec[i,,v] <- apply(invlogit(cctfit$ppX[cctfit$BUGSoutput$sims.list$Om[ind[i],]==v,,i]),2,var),silent=TRUE))
                } 
              }else{
                suppressMessages(try(varvec[i,,v] <- var(invlogit(cctfit$ppY[cctfit$BUGSoutput$sims.list$Om[ind[i],]==v,,i])),silent=TRUE)) 
                if(version2==1){
                  suppressMessages(try(varvec[i,,v] <- var(invlogit(cctfit$ppX[cctfit$BUGSoutput$sims.list$Om[ind[i],]==v,,i])),silent=TRUE)) 
                } 
                
              }
            }
            if(sum(cctfit$respmem==v)>1){
              cctfit$datVDI[v] <- var(apply(datest[cctfit$respmem==v ,],2,var))
              
            }else{
              cctfit$datVDI[v] <- var(datest[cctfit$respmem==v ,])
              
            }
          }
          options(warn=0)
          vdi <- apply(varvec,c(1,3),function(x) var(x,na.rm=TRUE))
          cctfit$ppVDI <- vdi
          
          if(sum(apply(cctfit$ppVDI,2,function(x) all(x==0,na.rm=TRUE)))>=1){
            for(i in which(apply(cctfit$ppVDI,2,function(x) all(x==0)))){
              cctfit$ppVDI[,i] <- NA
            }
          }
          for(v in 1:dim(cctfit$ppVDI)[2]){
            if(identical(unique(cctfit$ppVDI[,v]),c(NA,0)) || identical(unique(cctfit$ppVDI[,v]),c(0,NA))){cctfit$ppVDI[,v] <- NA}
          }
        }
        if(usetransform==0){
          invlogit <-  function(x){x <- 1 / (1 + exp(-x)); return(x)}
        } 
        
        rm(vdi, varvec)
        rm(eigv,leigv,tmp,tmp2,tmp3,ind);
        
        cctfit$checksrun <- TRUE
      }
      
      
    }else{
      #########################################
      ##### Calculate Posterior Predictive Data From All Samples (memory intensive)
      #########################################
      if(cctfit$whmodel=="GCM"){
        if(cctfit$V==1){
          cctfit$V <- 1;
          cctfit$BUGSoutput$sims.list$Om <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
          if(length(dim(cctfit$BUGSoutput$sims.list$Z)) < 3){
            cctfit$BUGSoutput$sims.list$Z <- array(cctfit$BUGSoutput$sims.list[["Z"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
          if(cctfit$itemdiff==TRUE){
            if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
              cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
          }
        }
        if(cctfit$itemdiff==FALSE){
          if(cctfit$V==1){
            cctfit$BUGSoutput$sims.list$lam <- array(.5,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
          }
          if(cctfit$V > 1){
            cctfit$BUGSoutput$sims.list$lam <- array(.5,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
          } }
        
        cctfit$ppY <- array(NA, c(cctfit$n,cctfit$m,cctfit$BUGSoutput$n.sims));
        cctfit$BUGSoutput$sims.list$tau <- array(-1,c(cctfit$BUGSoutput$n.sims,cctfit$n,cctfit$m))
        
        for(samp in 1:cctfit$BUGSoutput$n.sims){
          cctfit$BUGSoutput$sims.list[["tau"]][samp,,] <- (cctfit$BUGSoutput$sims.list[["th"]][samp,]*t(1-cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])) / 
            ( (cctfit$BUGSoutput$sims.list[["th"]][samp,]*t(1-cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])) + 
                ((1-cctfit$BUGSoutput$sims.list[["th"]][samp,])*t(cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])) ) 
          
          cctfit$ppY[,,samp] <- matrix(rbinom((cctfit$n*cctfit$m),1,
                                              (t(cctfit$BUGSoutput$sims.list[["Z"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])*cctfit$BUGSoutput$sims.list[["tau"]][samp,,] ) +
                                                t(t(1-cctfit$BUGSoutput$sims.list[["tau"]][samp,,])%*%diag(cctfit$BUGSoutput$sims.list[["g"]][samp,])) ), cctfit$n,cctfit$m)
        }
        
        if(cctfit$mval==TRUE){
          for(i in 1:dim(cctfit$datob$thena)[1]){cctfit$data[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2]] <- Mode(cctfit$ppY[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2],])}
          cctfit$MVest <- cbind(cctfit$datob$thena,cctfit$data[cctfit$datob$thena])
          colnames(cctfit$MVest) <- c("Pers","Item","Resp")
        }
        
        cctfit$mean$ppY <- rowMeans(cctfit$ppY[,,],dims=2)
        
        leigv <- 12
        options(warn=-3)
        ind <- sample(cctfit$BUGSoutput$n.sims,min(500,cctfit$BUGSoutput$n.sims)) #500 is the number of samples
        eigv <- matrix(-1,length(ind),leigv)
        tmp <- apply(cctfit$ppY[,,ind],c(1,3),function(x) t(x))
        tmp2 <- apply(tmp[,,],3,function(x) cor(x))
        tmp3 <- array(tmp2,c(cctfit$n,cctfit$n,cctfit$BUGSoutput$n.sims))
        for(i in 1:length(ind)){
          suppressMessages(try(eigv[i,] <- fa(tmp3[,,i])$values[1:leigv],silent=TRUE))}
        wch <- -which(eigv[,1]==-1); 
        if(length(wch)==0){cctfit$ppeig <-  eigv}else{cctfit$ppeig <-  eigv[-which(eigv[,1]==-1),]}
        cctfit$dateig <-  suppressMessages(fa(cor(t(cctfit$data)))$values[1:leigv])
        options(warn=0)
        
        if(cctfit$V==1){
          varvec <- matrix(-1,length(ind),dim(cctfit$ppY)[2])
          vdi <- matrix(-1,dim(cctfit$ppY)[3],1)
          for(i in 1:length(ind)){
            varvec[i,] <- apply(cctfit$ppY[,,i],2,var)
          }
          vdi <- apply(varvec,1,var)
          cctfit$ppVDI <- vdi
          cctfit$datVDI <- var(apply(cctfit$data,2,var))
        }
        
        options(warn=-3)
        if(cctfit$V>1){
          varvec <- array(NA,c(length(ind),dim(cctfit$ppY)[2],cctfit$V))
          vdi <- matrix(NA,dim(cctfit$ppY)[3],cctfit$V)
          cctfit$datVDI  <- array(-1,c(1,cctfit$V))
          for(v in 1:cctfit$V){
            for(i in 1:length(ind)){
              if(sum(cctfit$BUGSoutput$sims.list$Om[ind[i],]==v)>1){
                suppressMessages(try(varvec[i,,v] <- apply(cctfit$ppY[cctfit$BUGSoutput$sims.list$Om[ind[i],]==v,,i],2,var),silent=TRUE))
              }else{
                suppressMessages(try(varvec[i,,v] <- var(cctfit$ppY[cctfit$BUGSoutput$sims.list$Om[ind[i],]==v,,i]),silent=TRUE))
              }
            }
            if(sum(cctfit$respmem==v)>1){
              cctfit$datVDI[v] <- var(apply(cctfit$data[cctfit$respmem==v,],2,var))
            }else{
              cctfit$datVDI[v] <- var(cctfit$data[cctfit$respmem==v,])
            }
          }
          vdi <- apply(varvec,c(1,3),var)
          cctfit$ppVDI <- vdi
          
          options(warn=0)
          if(sum(apply(cctfit$ppVDI,2,function(x) all(x==0,na.rm=TRUE)))>=1){
            for(i in which(apply(cctfit$ppVDI,2,function(x) all(x==0)))){
              cctfit$ppVDI[,i] <- NA
            }
          }
          for(v in 1:dim(cctfit$ppVDI)[2]){
            if(identical(unique(cctfit$ppVDI[,v]),c(NA,0)) || identical(unique(cctfit$ppVDI[,v]),c(0,NA))){cctfit$ppVDI[,v] <- NA}
          }
        }
        rm(vdi, varvec)
        rm(eigv,leigv,tmp,tmp2,tmp3,ind);
        
        cctfit$checksrun <- TRUE
      }
      
      if(cctfit$whmodel=="LTRM"){
        if(cctfit$V==1){
          cctfit$V <- 1;
          cctfit$BUGSoutput$sims.list$Om <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
          if(length(dim(cctfit$BUGSoutput$sims.list[["T"]])) < 3){
            cctfit$BUGSoutput$sims.list$T <- array(cctfit$BUGSoutput$sims.list[["T"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
          if(length(dim(cctfit$BUGSoutput$sims.list[["gam"]])) < 3){
            cctfit$BUGSoutput$sims.list$gam <- array(cctfit$BUGSoutput$sims.list[["gam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$C-1,1))}
          if(cctfit$itemdiff==TRUE){
            if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
              cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
          }
        }
        if(cctfit$itemdiff==FALSE){
          if(cctfit$V==1){
            cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
          }
          if(cctfit$V > 1){
            cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
          }
        }
        cctfit$ppY <- array(NA, c(cctfit$n,cctfit$m,cctfit$BUGSoutput$n.sims));
        cctfit$ppX <- array(NA, c(cctfit$n,cctfit$m,cctfit$BUGSoutput$n.sims));
        cctfit$ppdelta <- array(NA, c(cctfit$n,cctfit$C-1,cctfit$BUGSoutput$n.sims))
        
        cctfit$BUGSoutput$sims.list$tau <- array(-1,c(cctfit$BUGSoutput$n.sims,cctfit$n,cctfit$m))
        
        for(samp in 1:cctfit$BUGSoutput$n.sims){
          cctfit$BUGSoutput$sims.list[["tau"]][samp,,] <- cctfit$BUGSoutput$sims.list[["E"]][samp,]*t(1/cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])
          
          cctfit$ppX[,,samp] <- matrix(
            rnorm((cctfit$n*cctfit$m),t(cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]]),(cctfit$BUGSoutput$sims.list[["tau"]][samp,,]^(-.5))), 
            cctfit$n,cctfit$m)
          
          cctfit$ppdelta[,,samp] <- cctfit$BUGSoutput$sims.list[["a"]][samp,]*t(cctfit$BUGSoutput$sims.list[["gam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])+cctfit$BUGSoutput$sims.list[["b"]][samp,]
          
          rc <- which(cctfit$ppX[,,samp] < cctfit$ppdelta[,1,samp],arr.ind=TRUE)
          cctfit$ppY[cbind(rc,array(samp,dim(rc)[1]))] <- 1 
          for(c in 1:(cctfit$C-1)){
            rc <- which(cctfit$ppX[,,samp] > cctfit$ppdelta[,c,samp],arr.ind=TRUE)
            cctfit$ppY[cbind(rc,array(samp,dim(rc)[1]))] <- (c+1) }
          
        }
        cctfit$mean$ppY <- rowMeans(cctfit$ppY[,,],dims=2)
        
        if(cctfit$mval==TRUE){
          for(i in 1:dim(cctfit$datob$thena)[1]){cctfit$data[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2]] <- Mode(cctfit$ppY[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2],])}
          cctfit$MVest <- cbind(cctfit$datob$thena,cctfit$data[cctfit$datob$thena])
          colnames(cctfit$MVest) <- c("Pers","Item","Resp")
        }
        
        leigv <- 12
        eigv <- matrix(-1,dim(cctfit$ppY)[3],leigv)
        options(warn=-3)
        ind <- sample(cctfit$BUGSoutput$n.sims,min(250,cctfit$BUGSoutput$n.sims)) #250 is the number of samples
        tmp <- apply(cctfit$ppY[,,ind],c(1,3),function(x) t(x))
        if(polych==TRUE){tmp2 <- apply(tmp[,,],3,function(x) suppressMessages(polychoric(x,global=FALSE))$rho)}else{
          tmp2 <- apply(tmp[,,],3,function(x) cor(x))
        }
        tmp3 <- array(tmp2,c(cctfit$n,cctfit$n,length(ind))) 
        for(i in 1:length(ind)){
          suppressMessages(try(eigv[i,] <- fa(tmp3[,,i])$values[1:leigv],silent=TRUE))}
        wch <- -which(eigv[,1]==-1); 
        if(length(wch)==0){cctfit$ppeig <-  eigv}else{cctfit$ppeig <-  eigv[-which(eigv[,1]==-1),]}
        if(polych==TRUE){
          cctfit$dateig <-  suppressMessages(fa(polychoric(t(cctfit$data),global=FALSE)$rho)$values[1:leigv])}else{
            cctfit$dateig <-  suppressMessages(fa(cor(t(cctfit$data)))$values[1:leigv])
          }
        options(warn=0)
        
        if(cctfit$V==1){
          varvec <- matrix(NA,length(ind),dim(cctfit$ppY)[2])
          vdi <- matrix(NA,dim(cctfit$ppY)[3],1)
          for(i in 1:length(ind)){
            varvec[i,] <- apply(cctfit$ppY[,,i],2,var)
          }
          vdi <- apply(varvec,1,var)
          cctfit$ppVDI <- vdi
          cctfit$datVDI <- var(apply(cctfit$data,2,var))
        }
        options(warn=-3)
        if(cctfit$V>1){
          varvec <- array(NA,c(length(ind),dim(cctfit$ppY)[2],cctfit$V))
          vdi <- matrix(NA,dim(cctfit$ppY)[3],cctfit$V)
          cctfit$datVDI  <- array(-1,c(1,cctfit$V))
          for(v in 1:cctfit$V){
            for(i in 1:length(ind)){
              if(sum(cctfit$BUGSoutput$sims.list$Om[ind[i],]==v)>1){
                suppressMessages(try(varvec[i,,v] <- apply(cctfit$ppY[cctfit$BUGSoutput$sims.list$Om[ind[i],]==v,,i],2,var),silent=TRUE))
              }else{
                suppressMessages(try(varvec[i,,v] <- var(cctfit$ppY[cctfit$BUGSoutput$sims.list$Om[ind[i],]==v,,i]),silent=TRUE))
              }
            }
            if(sum(cctfit$respmem==v)>1){
              cctfit$datVDI[v] <- var(apply(cctfit$data[cctfit$respmem==v,],2,var))
            }else{
              cctfit$datVDI[v] <- var(cctfit$data[cctfit$respmem==v,])
            }
          }
          vdi <- apply(varvec,c(1,3),var)
          cctfit$ppVDI <- vdi
          
          options(warn=0)
          if(sum(apply(cctfit$ppVDI,2,function(x) all(x==0,na.rm=TRUE)))>=1){
            for(i in which(apply(cctfit$ppVDI,2,function(x) all(x==0)))){
              cctfit$ppVDI[,i] <- NA
            }
          }
          for(v in 1:dim(cctfit$ppVDI)[2]){
            if(identical(unique(cctfit$ppVDI[,v]),c(NA,0)) || identical(unique(cctfit$ppVDI[,v]),c(0,NA))){cctfit$ppVDI[,v] <- NA}
          }
        }
        cctfit$polycor <- polych
        
        rm(vdi, varvec)
        rm(eigv,leigv,tmp,tmp2,tmp3,ind);
        
        cctfit$checksrun <- TRUE
      }
      
      if(cctfit$whmodel=="CRM"){
        invlogit <- function(x){x <- 1 / (1 + exp(-x)); return(x)}
        logit <- function(x){x <- log(x/(1-x)); return(x)}
        
        if(cctfit$V==1){
          cctfit$V <- 1;
          cctfit$BUGSoutput$sims.list$Om <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$n));
          if(length(dim(cctfit$BUGSoutput$sims.list[["T"]])) < 3){
            cctfit$BUGSoutput$sims.list$T <- array(cctfit$BUGSoutput$sims.list[["T"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
          if(cctfit$itemdiff==TRUE){
            if(length(dim(cctfit$BUGSoutput$sims.list$lam)) < 3){
              cctfit$BUGSoutput$sims.list$lam <- array(cctfit$BUGSoutput$sims.list[["lam"]][,],c(cctfit$BUGSoutput$n.sims,cctfit$m,1))}
          }
        }
        if(cctfit$itemdiff==FALSE){
          if(cctfit$V==1){
            cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,1))
          }
          if(cctfit$V > 1){
            cctfit$BUGSoutput$sims.list$lam <- array(1,c(cctfit$BUGSoutput$n.sims,cctfit$m,cctfit$V))
          }
        }
        
        cctfit$ppY <- array(NA, c(cctfit$n,cctfit$m,cctfit$BUGSoutput$n.sims));
        cctfit$ppX <- array(NA, c(cctfit$n,cctfit$m,cctfit$BUGSoutput$n.sims));
        
        cctfit$BUGSoutput$sims.list$tau <- array(-1,c(cctfit$BUGSoutput$n.sims,cctfit$n,cctfit$m))
        
        for(samp in 1:cctfit$BUGSoutput$n.sims){
          cctfit$BUGSoutput$sims.list[["tau"]][samp,,] <- cctfit$BUGSoutput$sims.list[["E"]][samp,]*t(1/cctfit$BUGSoutput$sims.list[["lam"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]])
          
          cctfit$ppY[,,samp] <- matrix(
            rnorm((cctfit$n*cctfit$m),(cctfit$BUGSoutput$sims.list[["a"]][samp,]*t(cctfit$BUGSoutput$sims.list[["T"]][samp,,cctfit$BUGSoutput$sims.list[["Om"]][samp,]]))+matrix(rep(cctfit$BUGSoutput$sims.list[["b"]][samp,],cctfit$m),cctfit$n,cctfit$m),
                  (cctfit$BUGSoutput$sims.list[["a"]][samp,]*(cctfit$BUGSoutput$sims.list[["tau"]][samp,,]^(-.5)))), 
            cctfit$n,cctfit$m)
          
          cctfit$ppX[,,samp] <- (array(rep(1/cctfit$BUGSoutput$sims.list[["a"]][samp,],times=cctfit$m),c(cctfit$n,cctfit$m))*(cctfit$ppY[,,samp]))-matrix(rep(cctfit$BUGSoutput$sims.list[["b"]][samp,],times=cctfit$m),c(cctfit$n,cctfit$m))
        }
        cctfit$mean$ppY <- rowMeans(cctfit$ppY[,,],dims=2)
        
        if(cctfit$mval==TRUE){
          for(i in 1:dim(cctfit$datob$thena)[1]){cctfit$data[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2]] <- mean(cctfit$ppY[cctfit$datob$thena[i,1],cctfit$datob$thena[i,2],])}
          cctfit$MVest <- cbind(cctfit$datob$thena,cctfit$data[cctfit$datob$thena])
          colnames(cctfit$MVest) <- c("Pers","Item","Resp")
        }
        
        leigv <- 12
        options(warn=-3)
        ind <- sample(cctfit$BUGSoutput$n.sims,min(500,cctfit$BUGSoutput$n.sims)) #500 is the number of samples
        eigv <- matrix(-1,length(ind),leigv)
        tmp <- apply(cctfit$ppY[,,ind],c(1,3),function(x) t(x))
        tmp2 <- apply(tmp[,,],3,function(x) cor(x))
        tmp3 <- array(tmp2,c(cctfit$n,cctfit$n,cctfit$BUGSoutput$n.sims))
        for(i in 1:length(ind)){
          suppressMessages(try(eigv[i,] <- fa(tmp3[,,i])$values[1:leigv],silent=TRUE))}
        wch <- -which(eigv[,1]==-1); 
        if(length(wch)==0){cctfit$ppeig <-  eigv}else{cctfit$ppeig <-  eigv[-which(eigv[,1]==-1),]}
        cctfit$dateig <-  suppressMessages(fa(cor(t(cctfit$data)))$values[1:leigv])
        options(warn=0)
        
        if(cctfit$V==1){
          vdi <- matrix(-1,dim(cctfit$ppY)[3],1)
          varvec <- array(-1,c(length(ind),dim(cctfit$ppY)[2]))
          
          for(i in 1:length(ind)){
            #varvec[i,] <- apply(invlogit(cctfit$ppY[,,i]),2,var) 
            varvec[i,] <- apply(invlogit(cctfit$ppY[,,i]),2,var)
          }
          vdi <- apply(varvec,1,var)
          cctfit$ppVDI <- apply(varvec,1,var) 
          
          #cctfit$datVDI <- var(apply(invlogit(cctfit$data),2,var))  
          cctfit$datVDI <- var(apply(invlogit(cctfit$data),2,var))  
          
        }
        if(cctfit$V>1){
          varvec <- array(NA,c(length(ind),dim(cctfit$ppY)[2],cctfit$V))
          cctfit$datVDI  <- array(NA,c(1,cctfit$V))
          options(warn=-3)
          for(v in 1:cctfit$V){
            for(i in 1:length(ind)){
              if(sum(cctfit$BUGSoutput$sims.list$Om[ind[i],]==v)>1){
                #suppressMessages(try(varvec[i,,v] <- apply(invlogit(cctfit$ppY[cctfit$BUGSoutput$sims.list$Om[ind[i],]==v,,i]),2,var),silent=TRUE))
                suppressMessages(try(varvec[i,,v] <- apply(invlogit(cctfit$ppY[cctfit$BUGSoutput$sims.list$Om[ind[i],]==v,,i]),2,var),silent=TRUE))
              }else{
                #suppressMessages(try(varvec[i,,v] <- var(invlogit(cctfit$ppY[cctfit$BUGSoutput$sims.list$Om[ind[i],]==v,,i])),silent=TRUE)) 
                suppressMessages(try(varvec[i,,v] <- var(invlogit(cctfit$ppY[cctfit$BUGSoutput$sims.list$Om[ind[i],]==v,,i])),silent=TRUE)) 
              }
            }
            if(sum(cctfit$respmem==v)>1){
              #cctfit$datVDI[v] <- var(apply(invlogit(cctfit$data[cctfit$respmem==v,]),2,var)) 
              cctfit$datVDI[v] <- var(apply(invlogit(cctfit$data[cctfit$respmem==v,]),2,var)) 
            }else{
              #cctfit$datVDI[v] <- var(invlogit(cctfit$data[cctfit$respmem==v,])) 
              cctfit$datVDI[v] <- var(invlogit(cctfit$data[cctfit$respmem==v,])) 
            }
          }
          options(warn=0)
          vdi <- apply(varvec,c(1,3),function(x) var(x,na.rm=TRUE))
          cctfit$ppVDI <- vdi
          
          if(sum(apply(cctfit$ppVDI,2,function(x) all(x==0,na.rm=TRUE)))>=1){
            for(i in which(apply(cctfit$ppVDI,2,function(x) all(x==0)))){
              cctfit$ppVDI[,i] <- NA
            }
          }
          for(v in 1:dim(cctfit$ppVDI)[2]){
            if(identical(unique(cctfit$ppVDI[,v]),c(NA,0)) || identical(unique(cctfit$ppVDI[,v]),c(0,NA))){cctfit$ppVDI[,v] <- NA}
          }
        }
        
        rm(vdi, varvec)
        rm(eigv,leigv,tmp,tmp2,tmp3,ind);
        
        cctfit$checksrun <- TRUE
      }
      
    }
    
    if(cctfit$mval==TRUE && gui==FALSE){
      message(paste("\n ...Use function 'cctmvest()' to view model estimates for missing data",sep=""))
     # message(paste("\n ...Use  '",guidat$varname,"$MVest' to view posterior predictive estimates for missing data\n ...     '",guidat$varname,"$data' to view the full data matrix \n ...     '",guidat$varname,"$datamiss' to view the matrix with missing values",sep=""))
    }
    message("\n ...Posterior predictive checks complete")
  }
  
  if(saveplots==1){jpeg(file.path(gsub(".Rdata","ppc.jpg",savedir)),width = 6, height = 3, units = "in", pointsize = 12,quality=100,res=400)}
  if(saveplots==2){postscript(file=file.path(gsub(".Rdata","ppc.eps",savedir)), onefile=FALSE, horizontal=FALSE, width = 6, height = 6, paper="special", family="Times")}
  
  par(oma=c(0,0,0,0),mar=c(4,4,3,1),mgp=c(2.25,.75,0),mfrow=c(1,2))
  
  plot(NA, main="Culture Number Check",xlim=c(1,min(cctfit$n,8)),ylim=c(min(cctfit$dateig[min(cctfit$n,8)],cctfit$ppeig[,min(cctfit$n,8)]),ceiling(max(cctfit$ppeig[,1],cctfit$dateig[1]))),xlab="Eigenvalue",ylab="Value",las=1,pch=21,type="l",col="black")
  for(i in 1:dim(cctfit$ppeig)[1]){
    points(cctfit$ppeig[i,],col="grey",type="l") }
  points(cctfit$dateig,col="black",type="l")
  
  if(cctfit$V>1){
    cctfit$VDIperc <- array(NA,c(cctfit$V,1))
    color <- "black"; color2 <- "black"
    color <- c("black","dark grey","light grey",rainbow(cctfit$V))
    minmax <- array(NA,c(cctfit$V,2))
    minmax2 <- array(NA,c(cctfit$V,1))
    
    for(v in which(sort(which(!apply(cctfit$ppVDI,2,function(x) all(is.na(x))))) %in% sort(unique(cctfit$respmem)))){
      tmpdist <-ecdf(cctfit$ppVDI[,v])
      minmax[v,] <- c(quantile(tmpdist, .005),quantile(tmpdist, .995))
      minmax2[v] <- max(density(cctfit$ppVDI[,v],na.rm=TRUE)$y)
    }
    plot(NA,main="Item Difficulty Check", xlab= "VDI", xlim=c(min(cctfit$datVDI[!is.na(cctfit$datVDI)],minmax,na.rm=TRUE),max(cctfit$datVDI[!is.na(cctfit$datVDI)],minmax,na.rm=TRUE)), ylim=c(0,max(minmax2,na.rm=TRUE)), ylab="Value",las=1,col=color[1],type="l")
    for(v in which(sort(which(!apply(cctfit$ppVDI,2,function(x) all(is.na(x))))) %in% sort(unique(cctfit$respmem)))){
      tmpdist <-ecdf(cctfit$ppVDI[,v])
      points(density(cctfit$ppVDI[,v],na.rm=TRUE),lwd=1.5,main="Item Difficulty Check", xlab= "VDI", xlim=c(min(cctfit$datVDI[!is.na(cctfit$datVDI)],quantile(tmpdist, .005),na.rm=TRUE),max(cctfit$datVDI[!is.na(cctfit$datVDI)],quantile(tmpdist, .995),na.rm=TRUE)), ylim=c(0,max(density(cctfit$ppVDI,na.rm=TRUE)$y,na.rm=TRUE)), ylab="Value",las=1,col=color[v],type="l")
      if(saveplots!=1 && saveplots!=2){print(paste("VDI Culture ",v," : ",round(100*tmpdist(cctfit$datVDI[,v]),digits=2)," percentile"),col=color2) 
                                           cctfit$VDIperc[v,1] <- round(100*tmpdist(cctfit$datVDI[,v]),digits=2)
      }
    }
    rm(tmpdist)
    segments(cctfit$datVDI,0,cctfit$datVDI,par("usr")[4],col=color,lwd=1.5)
  }else{
    color <- "black"; color2 <- "black"
    tmpdist <-ecdf(cctfit$ppVDI)
    plot(density(cctfit$ppVDI,na.rm=TRUE),main="Item Difficulty Check", xlab= "VDI", xlim=c(min(cctfit$datVDI,quantile(tmpdist, .005)),max(cctfit$datVDI,quantile(tmpdist, .995))), ylim=c(0,max(density(cctfit$ppVDI,na.rm=TRUE)$y,na.rm=TRUE)), ylab="Value",las=1,col=color2)
    segments(cctfit$datVDI,0,cctfit$datVDI,par("usr")[4],col=color2)
    text(.7*par("usr")[2],.7*par("usr")[4],labels=paste(round(100*tmpdist(cctfit$datVDI),digits=2)," percentile"),col=color2) 
    if(saveplots!=1 && saveplots!=2){cctfit$VDIperc <- round(100*tmpdist(cctfit$datVDI),digits=2)
                                         print(paste("VDI Culture 1 : ",round(100*tmpdist(cctfit$datVDI),digits=2)," percentile"),col=color2)
    }; rm(tmpdist)
  }
  
  rm(color,color2)
  if(saveplots==1 || saveplots==2){dev.off()}
  saveplots <- 0
  
    return(cctfit)  
}

#######################
#Accessor Function for Missing Value Estimates
#######################
cctmvest <- function(cctfit){
  
  if(cctfit$mval==FALSE){
    message("\n ...Data has no missing values to estimate",sep="")
    return()
  }
  if(cctfit$checksrun==FALSE){
    message("\n ...Please use 'cctppc()' first",sep="")
  }else{
    return(cctfit$MVest)
  }
  
}

#######################
#Accessor Function for Cluster Memberships
#######################
cctmemb <- function(cctfit){
  return(cctfit$respmem)
}


#######################
#cctfit Summary Function
#######################
cctsum <- function(cctfit){
  tmp <- c(" ","X")
  if(cctfit$checksrun==FALSE){
  cctfit$sumtab <- t(array(
    c(paste("Data:"),cctfit$datob$datatype,paste("Model:"),cctfit$whmodel,
      paste("Respondents:"),cctfit$n,paste("Items:"),cctfit$m,
      paste("Observations:   "),paste(100*round(cctfit$datob$nobs / (cctfit$datob$nobs+cctfit$datob$nmiss),2),"%",sep=""),paste("Missing Values:   "),paste(100*round(cctfit$datob$nmiss / (cctfit$datob$nobs+cctfit$datob$nmiss),2),"%",sep=""),
      paste("Cultures:"),cctfit$V,paste("Item Difficulty"),paste("[",tmp[cctfit$itemdiff+1],"]",sep=""),
      paste("DIC:"),round(cctfit$BUGSoutput$DIC,2),paste("pD:"),round(cctfit$BUGSoutput$pD,2),
      paste("Rhats > 1.10:"),paste(cctfit$Rhat$above110,"/",cctfit$Rhat$ncp),paste("Rhats > 1.05:"),paste(cctfit$Rhat$above105,"/",cctfit$Rhat$ncp),
      paste("Samples:"),cctfit$BUGSoutput$n.iter,paste("Chains:"),cctfit$BUGSoutput$n.chains,
      paste("Burn-in:"),cctfit$BUGSoutput$n.burnin,paste("Thinning:"),cctfit$BUGSoutput$n.thin
      #,paste("PPC Calculated"),paste("[",tmp[cctfit$checksrun+1],"]",sep=""),paste(""),paste("")
    ),c(4,8)))
  }else{
      cctfit$sumtab <- t(array(
        c(paste("Data:"),cctfit$datob$datatype,paste("Model:"),cctfit$whmodel,
          paste("Respondents:"),cctfit$n,paste("Items:"),cctfit$m,
          paste("Observations:   "),paste(100*round(cctfit$datob$nobs / (cctfit$datob$nobs+cctfit$datob$nmiss),2),"%",sep=""),paste("Missing Values:   "),paste(100*round(cctfit$datob$nmiss / (cctfit$datob$nobs+cctfit$datob$nmiss),2),"%",sep=""),
          paste("Cultures:"),cctfit$V,paste("Item Difficulty"),paste("[",tmp[cctfit$itemdiff+1],"]",sep=""),
          paste("DIC:"),round(cctfit$BUGSoutput$DIC,2),paste("pD:"),round(cctfit$BUGSoutput$pD,2),
          paste("Rhats > 1.10:"),paste(cctfit$Rhat$above110,"/",cctfit$Rhat$ncp),paste("Rhats > 1.05:"),paste(cctfit$Rhat$above105,"/",cctfit$Rhat$ncp),
          paste("PPC Calculated"),paste("[",tmp[cctfit$checksrun+1],"]",sep=""),paste("VDI"),paste(round(c(cctfit$VDI)),collapse=', '),
          paste("Samples:"),cctfit$BUGSoutput$n.iter,paste("Chains:"),cctfit$BUGSoutput$n.chains,
          paste("Burn-in:"),cctfit$BUGSoutput$n.burnin,paste("Thinning:"),cctfit$BUGSoutput$n.thin
          ),c(4,9)))  
    }
  
  rownames(cctfit$sumtab) <- rep("", nrow(cctfit$sumtab))
  colnames(cctfit$sumtab) <- rep("", ncol(cctfit$sumtab))
  print(cctfit$sumtab,quote=F)
}  
######################
#Function for the 'Export Results' Button
#- Prompts user where to save, and the filename to use
#- Saves the inference results (cctfit) to an .Rdata file of that filename
#- Saves .eps and .jpeg plots of the scree plot, plot results button, and run checks button
#    if the checks were calculated (the checks button was pressed)
######################
exportfuncbutton <- function() {
  guidat <- get("guidat", pkg_globals)
  exportfunc(cctfit=get(guidat$varname, pkg_globals),gui=TRUE,guidat=guidat)
}

exportfunc <- function(cctfit,filename=paste(getwd(),"/CCTpackdata.Rdata",sep=""),gui=FALSE,guidat=NULL) {
  
  if(gui==TRUE){
    savedir <- tclvalue(tkgetSaveFile(initialfile = "CCTpackdata.Rdata",filetypes = "{{Rdata Files} {.Rdata}}"))
    if (!nchar(savedir)) {
      message("\n ...Export cancelled\n")
      return()}
  }else{
    savedir <- filename
  }
  
  if (savedir=="" && !file.exists("CCTpack")){dir.create(file.path(getwd(), "CCTpack"));
                                               savedir <- file.path(getwd(),"CCTpack","CCTpackdata.Rdata")
  } 
  
  if(substr(savedir, nchar(savedir)-6+1, nchar(savedir))!=".Rdata"){
    savedir <- paste(savedir,".Rdata",sep="")
  }
  
  message("\n ...Exporting results")
  write.csv(cctfit$BUGSoutput$summary,file.path(gsub(".Rdata","posterior.csv",savedir)))
  save(cctfit,file=file.path(savedir))
  
  if(cctfit$checksrun==TRUE){
    ppcfunc(cctfit=cctfit,saveplots=1,savedir); ppcfunc(cctfit=cctfit,saveplots=2,savedir)}
  plotresultsfunc(cctfit=cctfit,saveplots=1,savedir); plotresultsfunc(cctfit=cctfit,saveplots=2,savedir);
  screeplotfunc(datob=cctfit$datob,saveplots=1,savedir); screeplotfunc(datob=cctfit$datob,saveplots=2,savedir)
  
  message("\n ...Export complete\n")
}

######################
#Function for the 'Fit Summary' Button
#- Gives a basic summary of the fit
######################
summaryfuncbutton <- function() {
  guidat <- get("guidat", pkg_globals)

  summary(get(guidat$varname, pkg_globals))
}

######################
#Function for the 'Print Fit' Button
#- Prints the cctfit object in the console
######################
printfitfuncbutton <- function() {
  guidat <- get("guidat", pkg_globals)
  
  print(get(guidat$varname, pkg_globals))
}

######################
#Function for the 'Send to Console' Button
#- Gives user the code to load the cctfit object into the console
######################
sendconsolefuncbutton <- function() {
  guidat <- get("guidat", pkg_globals)

  message("\n ...Paste the following to send your fit object into the R console
           \n    ",guidat$varname," <- CCTpack:::pkg_globals$",guidat$varname,
          "\n",sep="")
}

######################
#Function for the 'Send to Console' Button
#- Gives user the code to load the cctfit object into the console
######################
sendconsolefuncbutton <- function() {
  guidat <- get("guidat", pkg_globals)
  message("",sep="")
  message("\n ...Paste the following to send your fit object into the R console
          \n    ",guidat$varname," <- CCTpack:::pkg_globals$",guidat$varname,
          "\n")
}

######################
#Function for the 'Resp Memberships' Button
#- Outputs the respondent memberships from the model fit
######################
membfuncbutton <- function(){
  guidat <- get("guidat", pkg_globals)
  message("\n ...Model-based cluster memberships for each respondent \n",sep="")
  print(get(guidat$varname, pkg_globals)$respmem)
  #print(get(guidat$varname, pkg_globals)$MVest)
}

######################
#Function for the 'NA Value Est' Button
#- Plot traceplots for all discrete parameters 3x3 plot size
######################
mvestfuncbutton <- function(){
  guidat <- get("guidat", pkg_globals)
  cctfit <- get(guidat$varname, pkg_globals)
  if(cctfit$mval==FALSE){
    message("\n ...Data has no missing values to estimate\n",sep="")
    return()
  }
  if(cctfit$checksrun==FALSE){
  message("\n ...Please use 'Run Checks' first",sep="")
  }else{
    message("\n ...Model estimates for missing values \n",sep="")
    print(cctfit$MVest)  
  }
}

######################
#Function for the 'Traceplot Discrete' Button
#- Plot traceplots for all discrete parameters 3x3 plot size
######################
traceplotdiscretefuncbutton <- function() {
  guidat <- get("guidat", pkg_globals)
  
  dtraceplot(get(guidat$varname, pkg_globals))
}

######################
#Function for the 'Traceplot Continuous' Button
#- Plot traceplots for all parameters 3x3 plot size
######################
traceplotallfuncbutton <- function() {
  guidat <- get("guidat", pkg_globals)
  
  traceplot(get(guidat$varname, pkg_globals),mfrow=c(3,3), ask = FALSE)
}


dtraceplot <- function(cctfit,ask=FALSE){
  if(cctfit$V==1 && cctfit$whmodel!="GCM"){
    message("\n ...There are no discrete nodes in this inference")
    return() 
  }

  if(cctfit$whmodel!="GCM"){
    inds <- c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]])
  }else{
    if(cctfit$V==1){
      inds <- c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Z")]])
    }
    if(cctfit$V>1){
      inds <- c(cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Z")]],cctfit$BUGSoutput$long.short[[which(cctfit$BUGSoutput$root.short=="Om")]])   
    }  
  }
   
  if(cctfit$BUGSoutput$n.chains > 1){
    cctfit$BUGSoutput$sims.array <- cctfit$BUGSoutput$sims.array[,,inds]  
  }else{
    dmnames <- dimnames(cctfit$BUGSoutput$sims.array); dmnames[[3]] <- dmnames[[3]][inds] 
    tmp <-  array(cctfit$BUGSoutput$sims.array[,,inds],c(cctfit$BUGSoutput$n.keep,length(inds),1))
    cctfit$BUGSoutput$sims.array <- aperm(tmp,c(1,3,2)); dimnames(cctfit$BUGSoutput$sims.array) <- dmnames;
  }
    
  traceplot(cctfit,ask=ask,mfrow=c(3,3))  
}

#######################
#End of Functions for the GUI and/or Standalone Functions
#######################

#######################
#Package Methods
#######################

setClass("cct",contains=c("rjags")); 
setGeneric("plot"); setGeneric("screeplot"); setGeneric("summary")
setMethod("plot",signature=c(x="cct"),definition=function(x) cctresults(x))
setMethod("summary",signature=c(object="cct"),definition=function(object) cctsum(object))
setMethod("screeplot",signature=c(x="data.frame"),definition=function(x,polych=F) cctscree(x,polych=polych))
setMethod("screeplot",signature=c(x="matrix"),definition=function(x,polych=F) cctscree(x,polych=polych))
setMethod("screeplot",signature=c(x="cct"),definition=function(x,polych=F) cctscree(x$datob$dat,polych=polych))

######################
#END CCTpack
######################