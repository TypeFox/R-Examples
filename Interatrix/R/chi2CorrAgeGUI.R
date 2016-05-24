chi2CorrAgeGUI <- function() {
  
	font1 <- tkfont.create(family = "courrier", size = 13, weight = "bold", slant = "italic")
	font2 <- tkfont.create(family = "courrier", size = 11)
  
	## Gets the path of the data file
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
	label.att <- tklabel(ChooseAtt, width = 40, text = "    Identify qualitative factor(s)   ", font = font2)
  subChooseAtt <- tkframe(ChooseAtt)
  yscr.att <- tkscrollbar(subChooseAtt, command = function(...) tkyview(choose.att, ...))
  choose.att <- tklistbox(subChooseAtt, width = 25, height = 8, exportselection = FALSE, selectmode = "multiple", background = "white", yscrollcommand = function(...) tkset(yscr.att, ...))
  sapply(lab.fac, function(x) {tkinsert(choose.att, "end", x)})
	tkselection.set(choose.att, -1)
  tkpack(choose.att, side = "left", fill = "both")
	tkpack(yscr.att, side = "right", fill = "y")
	tkpack(subChooseAtt, side = "right")
  tkpack(label.att, side = "left", fill = "x")
	tkpack(ChooseAtt)
  
  ## Box to select the two columns with the serological stutus 
	ChoosePar <- tkframe(ChooseOverall)
  label.par <- tklabel(ChoosePar, width = 40, text = "       Select the two columns       
of serological statuses", font = font2)
  subChoosePar <- tkframe(ChoosePar)
  yscr.par <- tkscrollbar(subChoosePar, command = function(...) tkyview(choose.par, ...))
  choose.par <- tklistbox(subChoosePar, width = 25, height = 8, exportselection = FALSE, selectmode = "multiple", background = "white", yscrollcommand = function(...) tkset(yscr.par, ...))
  sapply(lab.fac, function(x) {tkinsert(choose.par, "end", x)})
	tkselection.set(choose.par, -1)
  tkpack(choose.par, side = "left", fill = "both")
	tkpack(yscr.par, side = "right", fill = "y")
	tkpack(subChoosePar, side = "right")
  tkpack(label.par, side = "left", fill = "x")
	tkpack(ChoosePar)
  
  ## Field for w1
	W1frame <- tkframe(ChooseOverall)
	w1 <- tclVar("")
	W1 <- tkentry(W1frame, width = 27, textvariable = w1)
  W1label <- tklabel(W1frame, width = 40, text = "Antibodies' disappearance rate", font = font2) 
	tkgrid(W1label, W1)
	tkpack(W1frame)
  
  ## Field for w2
	W2frame <- tkframe(ChooseOverall)
	w2 <- tclVar("")
	W2 <- tkentry(W2frame, width = 27, textvariable = w2)
  W2label <- tklabel(W2frame, width = 40, text = "Antibodies' disappearance rate", font = font2) 
	tkgrid(W2label, W2)
	tkpack(W2frame)
  
  tkbind(choose.par, "<ButtonRelease-1>", function() {
  	      nPar <- as.numeric(tkcurselection(choose.par)) + 1
          ParChoice <- lab.fac[nPar]
          tkconfigure(W1label, text = paste(" Antibodies' disappearance rate of ", na.omit(ParChoice[1]), sep = ""))
          tkconfigure(W2label, text = paste(" Antibodies' disappearance rate of ", na.omit(ParChoice[2]), sep = ""))
	    })
  
  ## Box to select the column with the age class
  ChooseAge <- tkframe(ChooseOverall)
  label.age <- tklabel(ChooseAge, width = 40, text = "     Select the age factor     ", font = font2)
  subChooseAge <- tkframe(ChooseAge)
  yscr.age <- tkscrollbar(subChooseAge, command = function(...) tkyview(choose.age, ...))
  choose.age <- tklistbox(subChooseAge, width = 25, height = 8, exportselection = FALSE, selectmode = "single", background = "white", yscrollcommand = function(...) tkset(yscr.age, ...))
  sapply(lab.fac, function(x) {tkinsert(choose.age, "end", x)})
	tkselection.set(choose.age, -1)
  tkpack(choose.age, side = "left", fill = "both")
	tkpack(yscr.age, side = "right", fill = "y")
	tkpack(subChooseAge, side = "right")
  tkpack(label.age, side = "left", fill = "x")
	tkpack(ChooseAge)
  
  ## Field for mortality rate for each class of age
  Mortframe <- tkframe(ChooseOverall)
  mort <- tclVar("")
  Mort <- tkentry(Mortframe, width = 27, textvariable = mort)
  Mortlabel <- tklabel(Mortframe, width = 40, text = "Mortality rate by age class (separated by ';')", font = font2) 
  tkgrid(Mortlabel, Mort)
  tkpack(Mortframe)
  
  ## Field for mortality rate for each class of age
  BoundClassframe <- tkframe(ChooseOverall)
  boundclass <- tclVar("")
  BoundClass <- tkentry(BoundClassframe, width = 27, textvariable = boundclass)
  BoundClasslabel <- tklabel(BoundClassframe, width = 40, text = "Bounds of the age classes (separated by ';')", font = font2) 
  tkgrid(BoundClasslabel, BoundClass)
  tkpack(BoundClassframe)
  
  ## Field for the model expression
	Formulaframe <- tkframe(ChooseOverall)
	formula <- tclVar("")
	formul <- tkentry(Formulaframe, width = 27, textvariable = formula)
	tkgrid(tklabel(Formulaframe, width = 40, text = "Formula", font = font2), formul)
	tkpack(Formulaframe)
  
  ## Field for the number of repeated simulations
	NbSimuFrame <- tkframe(ChooseOverall)
	nsimu <- tclVar(500)
	sim <- tkentry(NbSimuFrame, width = 27, textvariable = nsimu)
	tkgrid(tklabel(NbSimuFrame, width = 40, text = "Simulation number", font = font2), sim)
	tkpack(NbSimuFrame) 
  
	tkpack(ChooseAtt, ChoosePar, W1frame, W2frame, ChooseAge, Mortframe, BoundClassframe, Formulaframe, NbSimuFrame)
  
	OnOKChoice <- function() {
    simustart <- FALSE
		
    ## collect variable values
		nAtt <- as.numeric(tkcurselection(choose.att)) + 1
    nPar <- as.numeric(tkcurselection(choose.par)) + 1
    ParChoice <- lab.fac[nPar]
    nAge <- as.numeric(tkcurselection(choose.age)) + 1
    if(length(nAge) != 1)
      tkmessageBox(message = "Please select one column for the class of age.")
    else
    	nClassAge <- nlevels(as.factor(data.obs[, nAge]))
    mortality <- as.character(tclvalue(mort))
    mortality <- as.numeric(unlist(strsplit(mortality, ";")))
    boundclassage <- as.character(tclvalue(boundclass))
    boundclassage <- as.numeric(unlist(strsplit(boundclassage, ";")))
    formul <- as.character(tclvalue(formula))
    sim <- as.numeric(tclvalue(nsimu))
    
    w1 <- as.numeric(tclvalue(w1))
    w2 <- as.numeric(tclvalue(w2))
    
    namepara1 <- ParChoice[1]
    namepara2 <- ParChoice[2]
    nameage <- lab.fac[nAge]
    
    ## check variables
    if(length(ParChoice) != 2)
      tkmessageBox(message = "Please select the two columns of serologic statuses.")
    else if((w1 < 0) | (w1 > 1) | (w2 < 0) | (w2 > 1))
      tkmessageBox(message = "The antibodies' disappearance rates must be comprised between 0 and 1.")
    else if(length(mortality) != nClassAge | !all(mortality >= 0) | !all(mortality <= 1))
      tkmessageBox(message = paste("There must be ", nClassAge, " mortality rates comprised between 0 and 1 and separated by ';'.", sep = ""))
    else if((length(boundclassage) != nClassAge + 1) | !all(boundclassage == sort(boundclassage)) | !all(boundclassage >= 0))
      tkmessageBox(message = paste("There must be ", nClassAge + 1, " bound(s) separated by ';' in increasing order.", sep = ""))
    else if(formul == "")
      tkmessageBox(message = "Please enter a model formula.")
    else if(sim <= 0 | is.na(sim))
      tkmessageBox(message = "Please enter a positive number of simulations.")
    else
      simustart <- TRUE
    
    ## start process   
    if(isTRUE(simustart)) {
			tkdestroy(choose.variables)
      
    	## data transformation for factors
      numCol <- c(nAtt, nPar, nAge)
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

      dialog <- tktoplevel()
			tkwm.title(dialog, "")
			frmbut <- tkframe(dialog)
      
      
      ## button to choose parallel computing 
			yesbut <- tkbutton(frmbut, text = "Yes", command = function() {
			  if(requireNamespace("foreach", quietly = TRUE) && requireNamespace("doParallel", quietly = TRUE)) {
          waitingframe1 <- tktoplevel()
          tkwm.title(waitingframe1, "In progress")
          lab1 <- tklabel(waitingframe1, text = "Please wait...", width = 30, height = 10)
          tcl("update")
          tkgrid(lab1)
          tcl("update")
          
          tkdestroy(dialog)
          
          ## Calculations
          ncores <- max(1, parallel::detectCores() - 1) ## cores available on the computer
          res <- chi2CorrAge(formula = formul, data.obs = data.obs, namepara1 = namepara1, namepara2 = namepara2, nameage = nameage, w1 = w1, w2 = w2, mort = mortality, a = boundclassage, nsimu = sim, nbcore = ncores)
          dev.off()
          
          ## Close window
          tkdestroy(waitingframe1)
          
		      ## Print in a file
          OutFileSaved <- paste(SavePath, "/Sim_", as.character(Sys.Date()), ".txt", sep = "")
          list2ascii(res, file = OutFileSaved)
          
		      ## Print in the R consol
          res$chi2.corr.sim <- NULL   ## not print the 'sim' simulated chi2 
          assign("simu_chi2corrage", res, envir = parent.frame(n = 2))
          print(res)
          cat(substr(options("prompt")$prompt, 1, 2))
        }
      })
    
      
    	## button to calculate without parallel computing 
			nobut <- tkbutton(frmbut, text = "No", command = function() {
        waitingframe2 <- tktoplevel()
        tkwm.title(waitingframe2, "In progress")
        lab1 <- tklabel(waitingframe2, text = "Please wait...", width = 30, height = 10)
        tcl("update")
        tkgrid(lab1)
        tcl("update")
        
        tkdestroy(dialog)
        
        ## Calculations
        res <- chi2CorrAge(formula = formul, data.obs = data.obs, namepara1 = namepara1, namepara2 = namepara2, nameage = nameage, w1 = w1, w2 = w2, mort = mortality, a = boundclassage, nsimu = sim, nbcore = 1)
        dev.off()
        
        ## Close window
        tkdestroy(waitingframe2)
        
		    ## Print in a file
        OutFileSaved <- paste(SavePath, "/Sim_", as.character(Sys.Date()), ".txt", sep = "")
        list2ascii(res, file = OutFileSaved)
        
		    ## Print in the R consol
        res$chi2.corr.sim <- NULL   ## not print the 'sim' simulated chi2 
        assign("simu_chi2corrage", res, envir = parent.frame(n = 2))
        print(res)
        cat(substr(options("prompt")$prompt, 1, 2))
      })
      
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
      
  		tkpack(tklabel(frmbut, text = ""), yesbut, side = "left")
  		tkpack(nobut, tklabel(frmbut, text = ""), side = "right")
  		tkpack(tklabel(dialog, text = "Simulations must run for a long time.\
            Do you want to use a parallel design (recommanded) ?"))
  		tkpack(tklabel(dialog, text = ""))
  		tkpack(frmbut, fill = "x")
  		tkpack(tklabel(dialog, text = ""))
      
    } ## end simustart=TRUE
	}
  
  tkgrid(ChooseOverall)
	Choice.but <- tkbutton(choose.variables, text = "   OK   ", command = OnOKChoice, default = "active") 
	tkgrid(Choice.but)
	tkfocus(choose.variables)
}

