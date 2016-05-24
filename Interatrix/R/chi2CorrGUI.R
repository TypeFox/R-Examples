chi2CorrGUI <- function() {
 
	font1 <- tkfont.create(family = "courrier", size = 13, weight = "bold", slant = "italic")
	font2 <- tkfont.create(family = "courrier", size = 11)
 
	## Gets the path of the data file:
	getFilename <- function () {
    fildial <- paste("tk_getOpenFile", "-filetypes { \
	      {\"R binary files\" {.rda .rdata .Rdata .Rda}} \
    		{\"text files\" {.txt}} }")
	 	filename <- tclvalue(.Tcl(fildial))
    if(filename == "") {
      tkmessageBox(title = "Error", message = "No file was selected!")
      stop("Process killed because there is no data file dowloaded. Please restart.", call. = FALSE)
    } else
    	return(filename)
	}
  
	tkmessageBox(title = "Select source file", message = "Select the file contening data.")
	filedescpath <- getFilename()
  ext <- file_ext(filedescpath)
  if(ext == "txt")
		data.obs <- read.table(filedescpath, header = TRUE, sep = "\t")
  else
    data.obs <- get(load(filedescpath))
	data.obs <- na.omit(data.obs)
  lab.fac <- as.character(colnames(data.obs))
  
	## Precision of variables role:
	choose.variables <- tktoplevel()
	tkwm.title(choose.variables, "Settings of the analysis")
	ChooseOverall <- tkframe(choose.variables)
 
  ## Box to select qualitative factors for the model
	ChooseAtt <- tkframe(ChooseOverall)
	choose.att <- tklistbox(ChooseAtt, width = 25, height = ncol(data.obs), exportselection = FALSE, selectmode = "multiple", background = "white")
	tkgrid(tklabel(ChooseAtt, width = 30, text = "    Identify qualitative factor(s)   ", font = font2), choose.att)
  sapply(lab.fac, function(x) {tkinsert(choose.att, "end", x)})
	tkselection.set(choose.att, -1)
	tkpack(ChooseAtt)
  
  ## Box to select the two columns with the serological stutus 
	ChoosePar <- tkframe(ChooseOverall)
	choose.par <- tklistbox(ChoosePar, width = 25, height = ncol(data.obs), exportselection = FALSE, selectmode = "multiple", background = "white")
	tkgrid(tklabel(ChoosePar, width = 30, text = "       Select the two columns        
of serological statuses", font = font2), choose.par)
  sapply(lab.fac, function(x) {tkinsert(choose.par, "end", x)})
	tkselection.set(choose.par, -1)
	tkpack(ChoosePar)
  
  ## Field for the model expression
	Formulaframe <- tkframe(ChooseOverall)
	formula <- tclVar("")
	formul <- tkentry(Formulaframe, width = 25, textvariable = formula)
	tkgrid(tklabel(Formulaframe, width = 30, text = "Formula", font = font2), formul)
	tkpack(Formulaframe)
 
  ## Field for the number of repeated simulations
	NbSimuFrame <- tkframe(ChooseOverall)
	nsimu <- tclVar(500)
	sim <- tkentry(NbSimuFrame, width = 25, textvariable = nsimu)
	tkgrid(tklabel(NbSimuFrame, width = 30, text = "Simulation number", font = font2), sim)
	tkpack(NbSimuFrame) 

	tkpack(ChooseAtt, ChoosePar, Formulaframe, NbSimuFrame)
  
	OnOKChoice <- function() {
    simustart <- FALSE
		
    ## check variables
		nAtt <- as.numeric(tkcurselection(choose.att)) + 1
    nPar <- as.numeric(tkcurselection(choose.par)) + 1
    ParChoice <- lab.fac[nPar]
    formul <- as.character(tclvalue(formula))
    sim <- as.numeric(tclvalue(nsimu))
    
    namepara1 <- ParChoice[1]
    namepara2 <- ParChoice[2]
    
    if(length(ParChoice) != 2)
      tkmessageBox(message = "Please select the two columns of serologic statuses.")
    else if(class(formul) != "character")
      tkmessageBox(message = "Please enter a model formula.")
    else if(sim <= 0 | is.na(sim) | !is.numeric(sim))
      tkmessageBox(message = "Please enter a positive number of simulations.")
    else
      simustart <- TRUE

    
    ## start process   
    if(isTRUE(simustart)) {
			tkdestroy(choose.variables)
    
    	## data transformation for factors
      numCol <- c(nAtt, nPar)
      data.obs[, numCol] <- as.data.frame(lapply(data.obs[, numCol], as.factor))
   
      ## directory where to save simulation results
      SavePath <- tk_choose.dir(default = getwd(), caption = "Select directory to save results")
      if(is.na(SavePath)) {
        print(paste("The files are recorded in ", getwd(), sep = ""))
        SavePath <- getwd()
      }
      
      ## name and type of graphics
      getGraphName <- function () {
        fildial <- paste("tk_getSaveFile", "-filetypes { \
	                    {\".eps\" {.eps}} \
	                    {\".pdf\" {.pdf}} \
	                    {\".png\" {.png}} \
	                    {\".jpeg\" {.jpg}} }")
        filename <- tclvalue(.Tcl(fildial))
        return(filename)
      }
      tkmessageBox(title = "Graphics name and type", message = "Choose the name and the type of the graphics saved.")
      GraphPath <- getGraphName()
      
      waitingframe1 <- tktoplevel()
      tkwm.title(waitingframe1, "In progress")
      lab1 <- tklabel(waitingframe1, text = "Please wait...", width = 30, height = 10)
      tcl("update")
      tkgrid(lab1)
      tcl("update")
      
      graphics.off()
      ext <- substr(GraphPath, nchar(GraphPath) - 2, nchar(GraphPath))
      if(ext == "eps")
        postscript(file = GraphPath)
      else if(ext == "pdf")
        pdf(file = GraphPath)
      else if(ext == "png")
        png(filename = GraphPath)
      else if(ext == "jpg")
        jpeg(filename = GraphPath)

      ## Calculations
			res <- chi2Corr(formula = formul, data.obs = data.obs, namepara1 = namepara1, namepara2 = namepara2, nsimu = sim)
      dev.off()
      
      ## Close window
      tkdestroy(waitingframe1)
      
  		## Print in a file
      OutFileSaved <- paste(SavePath, "/Sim_", as.character(Sys.Date()), ".txt", sep = "")
      list2ascii(res, file = OutFileSaved)
      
  		## Print in the R consol
      res$chi2.corr.sim <- NULL   ## not print the 'sim' simulated chi2
      assign("simu_chi2corr", res, envir = parent.frame(n = 2))
      print(res)
      cat(substr(options("prompt")$prompt, 1, 2))
      
    } ## end simustart = TRUE
	}
  
  tkgrid(ChooseOverall)
	Choice.but <- tkbutton(choose.variables, text = "   OK   ", command = OnOKChoice, default = "active") 
	tkgrid(Choice.but)
	tkfocus(choose.variables)
}

