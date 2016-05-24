.makeNewFactorAnalysisDialog<-function(){

	#make dialog
	dialog <- new(SimpleRDialog)
	dialog$setSize(500L,400L)
	dialog$setTitle("Factor Analysis")
	
	#add variable selector
	variableSelector <- new(VariableSelectorWidget)
	variableSelector$setTitle("data")
	addComponent(dialog,variableSelector,10,400,850,10)
	
	#add a list for the variables
	variableList<- new(VariableListWidget,variableSelector)
	variableList$setTitle("variables")
	addComponent(dialog, variableList,100,900,450, 420)
	
	#Add an 'Options' button
	JButton <- J("javax.swing.JButton")
	button <- new(JButton,"Options")
	addComponent(dialog,button,600,825,700,625)

	#make Options Dialog
	subDialog <- new(SimpleRSubDialog,dialog,"Factor Analysis: Options")
	setSize(subDialog,350,450)
	subDialog$setResizable(FALSE)
	
	radio <- new(ButtonGroupWidget,"Matrix",c("Correlations","Variance-covariance"))
	radio$setTitle("Matrix")
	addComponent(subDialog,radio,30,900,220,100)

	combo <- new(J("org.rosuda.deducer.widgets.ComboBoxWidget"),
					c("principal components analysis", "maximum likelihood",
					"principal axis factoring","minimum residual (OLS)",
					"weighted least squares (WLS)","generalized weighted least squares (GLS)"))
	combo$setTitle("Extraction",TRUE)
	addComponent(subDialog, combo,250,900,390, 100)	
	combo$setDefaultModel("maximum likelihood")


	combo1 <-  new(J("org.rosuda.deducer.widgets.ComboBoxWidget"),"Rotation",
					c("none", "varimax", "quartimax", "bentlerT", "geominT", 
					"promax", "oblimin", "simplimax", "bentlerQ",  "geominQ", "cluster"))
	addComponent(subDialog, combo1,420,900,560, 100)
	combo1$setDefaultModel(c("promax"))
	
	# n of factors 
	textArea <- new(TextFieldWidget,"# of factors")
	addComponent(dialog, textArea,725,825,850, 625)
	textArea$setDefaultModel(c("2"))
	textArea$setInteger(TRUE)
	textArea$setLowerBound(1)

	#options for sorting factor loadings
	transBoxesp <- new(CheckBoxesWidget,c("Sorted by size","Cut loadings less than:","Save scores"))
	transBoxesp$setTitle("Loadings")
	addComponent(subDialog, transBoxesp,590,625,820,100)

	# value of loadings suppressed
	textArea0 <- new(TextFieldWidget)
	addComponent(subDialog, textArea0,670,770,740,630)
	textArea0$setDefaultModel(c("0.2"))
	textArea0$setNumeric(TRUE)
	textArea0$setLowerBound(0)
	
	#Listen for the button to be pressed
	ActionListener <- J("org.rosuda.deducer.widgets.event.RActionListener")
	actionFunction <- function(cmd,ActionEvent){
		subDialog$setLocationRelativeTo(button)
		subDialog$run()
	}
	listener <- new(ActionListener)
	listener$setFunction(toJava(actionFunction))
	button$addActionListener(listener)

	#Add an 'Plots' button
	JButton1 <- J("javax.swing.JButton")
	button1 <- new(JButton1,"Plots")
	addComponent(dialog,button1,475,825,575,625)

	#make Plots Dialog
	subDialog1 <- new(SimpleRSubDialog,dialog,"Factor Analysis: Plots")
	setSize(subDialog1,280,200)
	subDialog1$setResizable(FALSE)

	#options for plotting the variables
	transBoxes <- new(CheckBoxesWidget,"Plots",c("path diagram","parallel analysis","biplot"))
	addComponent(subDialog1, transBoxes,50,800,700, 200)
	

	#Listen for the button1 to be pressed
	ActionListener1 <- J("org.rosuda.deducer.widgets.event.RActionListener")
	actionFunction1 <- function(cmd,ActionEvent){
		subDialog1$setLocationRelativeTo(button1)
		subDialog1$run()
	}
	listener1 <- new(ActionListener1)
	listener1$setFunction(toJava(actionFunction1))
	button1$addActionListener(listener1)

	# Help button:
	dialog$addHelpButton("pmwiki.php?n=Main.FactorAnalysis")

	#required library ###############################################
	status<-.deducer$requirePackage("GPArotation")
	if(status=="installed"){
		execute("library('GPArotation')")
		}else if(status=="not installed"){
			stop("package GPArotation required")}
					
	status1<-.deducer$requirePackage("psych")					
	if(status1=="installed"){
		execute("library('psych')")
		}else if(status1=="not installed"){
			stop("package psych required")}
			
	status2<-.deducer$requirePackage("mvnormtest")					
	if(status2=="installed"){
		execute("library('mvnormtest')")
		}else if(status2=="not installed"){
			stop("package mvnormtest required")}
	###################################################################
	
	
	
	
	#make sure at least two variables are selected###################			
	.factorAnalysisCheckFunction <- function(state){
		if(length(state$variables)<2)
			return("Please select at least two variables")
		return("")
	}
	###################################################################
	####################required package#######
	status1<-.deducer$requirePackage("psych")					
	if(status1=="installed"){
		execute("library('psych')")
	}else if(status1=="not installed"){
		stop("package psych required")}


	.factorAnalysisRunFunction<- function(state) {
		 #print(state) #a print statement is useful for debugging
		
		
		
		
		# data frame ###################	
		ddf<-state$variables
		ddga<-state$data
		ddfa<-paste(ddga, "$", sep="")	 
		ddf<-paste(ddfa, ddf, sep="")
		form <- paste(ddf," ")
		form<-paste(form)
		forma<-paste(form[1])
		for( var in form[-1])
			forma<-paste(forma,", ", var)
		cmd<- paste("pr.model<-data.frame(", forma) 
		cmd<-paste(cmd, ")")
		cmd <- paste (cmd)
		
		
		# state rotation ###################
		if("none" %in%state$Rotation)
				rota <- "\"none\""
		if("varimax" %in%state$Rotation)
				rota <- "\"varimax\""
		if("quartimax" %in%state$Rotation)
				rota <- "\"quartimax\""
		if("bentlerT" %in%state$Rotation)
				rota <- "\"bentlerT\""
		if("geominT" %in%state$Rotation)
				rota <- "\"geominT\""
		if("promax" %in%state$Rotation)
				rota <- "\"promax\""
		if("oblimin" %in%state$Rotation)
				rota <- "\"oblimin\""
		if("simplimax" %in%state$Rotation)
				rota <- "\"simplimax\""
		if("bentlerQ" %in%state$Rotation)
				rota <- "\"bentlerQ\""
		if("geominQ" %in%state$Rotation)
				rota <- "\"geominQ\""
		if("cluster" %in%state$Rotation)
				rota <- "\"cluster\""
		
		# state extraction ###################
		if("maximum likelihood" %in%state$Extraction)
				extra <- "\"ml\""
		if("principal axis factoring" %in%state$Extraction)
				extra <- "\"pa\""
		if("minimum residual (OLS)" %in%state$Extraction)
				extra <- "\"minres\""
		if("weighted least squares (WLS)" %in%state$Extraction)
				extra <- "\"wls\""
		if("generalized weighted least squares (GLS)" %in%state$Extraction)
				extra <- "\"gls\""
		if("principal components analysis" %in%state$Extraction)
				extra <- "\"gls\""
		
		# state matrix ###################
		if("Correlations" %in% state$Matrix)
				cova <- "FALSE"
		if("Variance-covariance" %in% state$Matrix)
				cova <- "TRUE"
		
		# PCA ###################
		cmd7<- {cmd7<- paste("pr.model1<-principal(pr.model,") 
		cmd7<-paste(cmd7, "nfactors=", textArea$getModel(), ",")
		cmd7<-paste(cmd7, "rotate=", rota, ",")
		cmd7<-paste(cmd7, "covar=", cova,",")
		cmd7<-paste(cmd7, "scores=TRUE")
		cmd7<-paste(cmd7, ")")}
		
		# Factor analysis ###################
		cmd8<- {cmd8<- paste("pr.model1<-fa(pr.model,") 
		cmd8<-paste(cmd8, "nfactors=", textArea$getModel(), ",")
		cmd8<-paste(cmd8, "rotate=", rota, ",")
		cmd8<-paste(cmd8, "fm=", extra, ",")
		cmd8<-paste(cmd8, "covar=", cova)
		cmd8<-paste(cmd8, ")")}
		
		
		# print function ###################
		ifelse("principal components analysis" %in%state$Extraction, cmd11<- paste (cmd7),cmd11 <- paste (cmd8))
		
		# Loadings options ###################
		cmd12 <- {cmd12<- paste("print.psych(pr.model1")
		if ("Sorted by size" %in% state$Loadings) cmd12<-paste(cmd12, ", sort=T")
		if ("Cut loadings less than:" %in% state$Loadings) cmd12<-paste(cmd12, ", cut =",textArea0$getModel())
		cmd12<-paste(cmd12,")")}
		cm12.b <- if ("Save scores" %in% state$Loadings) paste(ddga,".1 <- cbind(pr.model, pr.model1$scores)", sep="") else character()
		
		# Plot options ###################
		nplots <- 0L
		cmd1 <- character()
		if ("path diagram" %in% state$Plots)
			{nplots <- nplots+1L
			cmd1<-paste("pr.model2<- fa.diagram(pr.model1)")}
		if ("parallel analysis" %in% state$Plots)
			{nplots <- nplots+1L
			cmd1<-paste(cmd1,"pr.model3<- fa.parallel(pr.model)",sep="\n")}
		if (nplots > 0L)
			{plot1<-paste("JavaGD(width=800, height=600, ps=12)")
			if (nplots == 2L) plot1<-paste(plot1,"par(mfrow=c(2,1))",sep="\n")}
		if ("biplot" %in% state$Plots)
			{plot2<-paste("JavaGD(width=500, height=500, ps=12)")
			cmd2<-paste("biplot(pr.model1)")}
		
		
		cmd000<-paste("print(mshapiro.test(t(na.omit(as.matrix(pr.model)))))")
		
		
		# remove data frame ###################
		cmd13<-paste("rm(pr.model)")
		
		# execute sequence ###################
		if (nplots > 0L){
			if ("biplot" %in% state$Plots){
				comm<-paste(cmd,"\n", cmd11,"\n",cmd12,"\n",cm12.b,"\n",plot1,"\n",cmd1,"\n",plot2,"\n",cmd2,"\n",cmd000,"\n",cmd13)
			}else{
				comm<-paste(cmd,"\n", cmd11,"\n",cmd12,"\n",cm12.b,"\n",plot1,"\n",cmd1,"\n",cmd000,"\n",cmd13)}
		}else{
			if ("biplot" %in% state$Plot){
				comm<-paste(cmd,"\n", cmd11,"\n",cmd12,"\n",cm12.b,"\n",plot2,"\n",cmd2,"\n",cmd000,"\n",cmd13)
			}else{
				comm<-paste(cmd,"\n", cmd11,"\n",cmd12,"\n",cm12.b,"\n",cmd000,"\n",cmd13)}}
		
		execute(comm)
	}

	dialog$setCheckFunction(toJava(.factorAnalysisCheckFunction))
	dialog$setRunFunction(toJava(.factorAnalysisRunFunction))
	dialog
}


###################################### 

#            Reliabiltity 	         #

######################################	


.makeNewReliabilityDialog<-function(){

	#make dialog
	dialog <- new(SimpleRDialog)
	dialog$setSize(500L,400L)
	dialog$setTitle("Cronbach Alpha")
	
	#add variable selector
	variableSelector <- new(VariableSelectorWidget)
	variableSelector$setTitle("data")
	addComponent(dialog,variableSelector,10,400,850,10)
	
	#add a list for the variables
	variableList<- new(VariableListWidget,variableSelector)
	variableList$setTitle("variables")
	addComponent(dialog, variableList,100,900,450, 420)
 
	# Help button:
	dialog$addHelpButton("pmwiki.php?n=Main.ReliabilityAnalysis")

	#make sure at least two variables are selected###################			
	.ReliabilityCheckFunction <- function(state){
		if(length(state$variables)<2)
			return("Please select at least two variables")
		return("")
	}



	.ReliabilityRunFunction <- function(state){
		
		### data.frame ########
		ddf<-state$variables
		ddga<-state$data
		ddfa<-paste(ddga, "$", sep="")		 
		ddf<-paste(ddfa, ddf, sep="")
		form <- paste(ddf," ")
		form<-paste(form)
		forma<-paste(form[1])
		for( var in form[-1])
		forma<-paste(forma,", ", var)
		
		cmd<- paste("pr.model<-data.frame(", forma) 
		
		cmd<-paste(cmd, ")")
		
		cmd <- paste (cmd)
		
		
		
		## model ########
		cmd1<- paste("pr.model1<-psych::alpha(pr.model") 
		
		cmd1<-paste(cmd1, ")")
		
		cmd1 <- paste (cmd1,"\n","print(pr.model1)")
		
		# remove data frame ###################
		cmd13<-paste("rm(pr.model)")
		
		comm<-paste(cmd,"\n",cmd1,"\n",cmd13)
		
		execute(comm)
	}

	dialog$setCheckFunction(toJava(.ReliabilityCheckFunction))
	dialog$setRunFunction(toJava(.ReliabilityRunFunction ))
	dialog
}




########################################################################
#
#				Intraclass correlation coefficient
#
########################################################################
.makeICCDialog <- function(){
	
	
	RFunction <- J("org.rosuda.deducer.widgets.param.RFunction")
	RFunctionDialog <- J("org.rosuda.deducer.widgets.param.RFunctionDialog")
	ParamCharacter <- J("org.rosuda.deducer.widgets.param.ParamCharacter")
	ParamNumeric <- J("org.rosuda.deducer.widgets.param.ParamNumeric")
	ParamMultipleVariables <- J("org.rosuda.deducer.widgets.param.ParamMultipleVariables")
	
	iccFunc <- new(RFunction,"icc")
	
	ratingsParam <- new(ParamMultipleVariables,"ratings")
	iccFunc$add(ratingsParam)
	
	modelParam <- new(ParamCharacter,"model","oneway")
	modelParam$setOptions(c("oneway", "twoway"))
	modelParam$setViewType(modelParam$VIEW_COMBO)
	iccFunc$add(modelParam)
	
	typeParam <- new(ParamCharacter,"type","consistency")
	typeParam$setOptions(c("consistency", "agreement"))
	typeParam$setViewType(typeParam$VIEW_COMBO)
	iccFunc$add(typeParam)
	
	unitParam <- new(ParamCharacter,"unit","single")
	unitParam$setOptions(c("single", "average"))
	unitParam$setViewType(unitParam$VIEW_COMBO)
	iccFunc$add(unitParam)
	
	nullParam <- new(ParamNumeric,"r0")
	nullParam$setTitle("H0")
	nullParam$setValue(0)
	nullParam$setLowerBound(0)
	iccFunc$add(nullParam)
	
	conf <- new(ParamNumeric,"conf.level",.95)
	conf$setLowerBound(0)
	conf$setUpperBound(1)
	iccFunc$add(conf)
	
	rfd <- new(RFunctionDialog, iccFunc)
	rfd$setSize(425L,500L)
	rfd$setLocationRelativeTo(.jnull())
	rfd
}

##########################################################

#            Linear Discriminant Analysis 	         #

##########################################################

.makeNewLDADialog <- function(){
	dialog <- new(SimpleRDialog)
	dialog$setTitle("Linear Discriminant Analysis")
	dialog$setSize(500L,450L)
	# Variable selector
	variableSelector <- new(VariableSelectorWidget)
	variableSelector$setTitle("data")
	addComponent(dialog,variableSelector,20,300,790,20)
	# Outcome (groups)
	groups <- new(SingleVariableWidget,variableSelector)
	groups$setTitle("Groups",show=TRUE)
	addComponent(dialog,groups,40,750,220,320)
	# Predictors list
	predictorsList <- new(VariableListWidget,variableSelector)
	predictorsList$setTitle("Predictors",show=TRUE)
	addComponent(dialog,predictorsList,240,750,790,320)
	# Assumption icons
	AssumptionIcon <- J("org.rosuda.deducer.toolkit.AssumptionIcon")
	largeAssump <- new(AssumptionIcon,"/icons/N_or_norm_assump.png","Large Sample or Normal",NULL,"Large Sample or Normal")
	addComponent(dialog,largeAssump,805,140,862,40)
	outliersAssump <- new(AssumptionIcon,"/icons/outlier_assump.png","No outliers",NULL,"No outliers")
	addComponent(dialog,outliersAssump,805,195,862,145)
	eqvarAssump <- new(AssumptionIcon,"/icons/eqvar_assump.png","Equal variances",NULL,"Equal variances")
	addComponent(dialog,eqvarAssump,805,250,862,200)
	colinearAssump <- new(AssumptionIcon,"/icons/func_assump.png","No collinearity",NULL,"No collinearity")
	addComponent(dialog,colinearAssump,805,305,862,255)
	# Buttons: options, variable selection, show, export
	JButton <- J("javax.swing.JButton")
	optionsButton <- new(JButton,"Options")
	addComponent(dialog,optionsButton,40,960,120,770)
	selectionButton <- new(JButton,"Selection")
	addComponent(dialog,selectionButton,140,960,220,770)
	showButton <- new(JButton,"Show")
	addComponent(dialog,showButton,240,960,320,770)
	exportButton <- new(JButton,"Export")
	addComponent(dialog,exportButton,340,960,420,770)
	# Help button:
	dialog$addHelpButton("pmwiki.php?n=Main.LinearDiscriminantAnalysis")
	# Options subdialog
	optionsSubDialog <- new(SimpleRSubDialog,dialog,"LDA Options")
	optionsSubDialog$setSize(250L,275L)
	priorButtons <- new(ButtonGroupWidget,"Prior probabilites",c("Observed","Equal"))
	priorButtons$setTitle("prior")
	addComponent(optionsSubDialog,priorButtons,20,960,320,20)
	methodComboBox <- new(ComboBoxWidget,"Method",c("Moment","MLE","MVE","Robust (t)"))
	methodComboBox$setTitle("method")
	methodComboBox$setDefaultModel("Moment")
	addComponent(optionsSubDialog,methodComboBox,340,680,550,20)
	nuLabel <- new(JLabel,"Df")
	addComponent(optionsSubDialog,nuLabel,340,960,420,705)
	nuTextField <- new(TextFieldWidget)
	nuTextField$setTitle("nu")
	nuTextField$setDefaultModel("5")
	nuTextField$getTextField()$setEnabled(FALSE)
	addComponent(optionsSubDialog,nuTextField,420,960,520,700)
	tolLabel <- new(JLabel,"Tolerance to singularity")
	addComponent(optionsSubDialog,tolLabel,570,800,640,50)
	tolTextField <- new(TextFieldWidget)
	tolTextField$setTitle("tol")
	tolTextField$setDefaultModel("0.0001")
	addComponent(optionsSubDialog,tolTextField,640,500,750,50)
	# Selection subdialog
	selectionSubDialog <- new(SimpleRSubDialog,dialog,"LDA: Selection of variables")
	selectionSubDialog$setSize(250L,150L)
	selectionButtons <- new(ButtonGroupWidget,c("All","Stepwise"))
	selectionButtons$setTitle("selection")
	addComponent(selectionSubDialog,selectionButtons,20,450,600,60)
	niveauLabel <- new(JLabel,"Max. p-value")
	addComponent(selectionSubDialog,niveauLabel,180,920,300,500)
	niveauTextField <- new(TextFieldWidget)
	niveauTextField$setTitle("niveau")
	niveauTextField$setDefaultModel("0.2")
	niveauTextField$getTextField()$setEnabled(FALSE)
	addComponent(selectionSubDialog,niveauTextField,350,920,550,500)
	# Show subdialog
	showSubDialog <- new(SimpleRSubDialog,dialog,"LDA: Show additional results")
	showSubDialog$setSize(400L,250L)
	ctablesCheckBoxes <- new(CheckBoxesWidget,"Classification tables",c("Original data","Cross-validation"),2L)
	ctablesCheckBoxes$setTitle("ctables")
	addComponent(showSubDialog,ctablesCheckBoxes,20,960,280,20)
	statsCheckBoxes <- new(CheckBoxesWidget,"Discriminant functions",c("Centroids","Statistics","Structure matrix","Plot"),2L)
	statsCheckBoxes$setTitle("dfs")
	addComponent(showSubDialog,statsCheckBoxes,320,960,740,20)
	# Export subdialog
	exportSubDialog <- new(SimpleRSubDialog,dialog,"LDA: Export values")
	exportSubDialog$setSize(250L,250L)
	exportCheckBoxes <- new(CheckBoxesWidget,c("Classification","A posteriori probabilities","Discriminant functions"))
	exportCheckBoxes$setTitle("export")
	addComponent(exportSubDialog,exportCheckBoxes,20,960,400,60)
	# Look for data frames in the workspace
	newdataLabel <- new(JLabel,"Apply to data frame:")
	addComponent(exportSubDialog,newdataLabel,420,920,550,60)
	newdataSelector <- new(VariableSelectorWidget)
	newdataSelector$setTitle("newdata")
	newdataComboBox <- newdataSelector$getJComboBox()
	addComponent(exportSubDialog,newdataComboBox,550,900,680,60)
	addComponent(exportSubDialog,newdataSelector,1,2,2,1)
	newdataSelector$setVisible(FALSE)
	# Listeners to show subdialogs
	ActionListener <- J("org.rosuda.deducer.widgets.event.RActionListener")
	linkSubDialog <- function(button,subdialog){
		listenerFunction <- function(cmd,ActionEvent){
			subdialog$setLocationRelativeTo(button)
			subdialog$run()
		}
		listener <- new(ActionListener)
		listener$setFunction(toJava(listenerFunction))
		button$addActionListener(listener)
	}
	linkSubDialog(optionsButton,optionsSubDialog)
	linkSubDialog(selectionButton,selectionSubDialog)
	linkSubDialog(showButton,showSubDialog)
	linkSubDialog(exportButton,exportSubDialog)
	# Listeners of methodComboBox and selectionButtons to show nuTextArea and niveauTextArea, respectively
	methodListenerFunction <- function(cmd,ActionEvent){
		method <- methodComboBox$getModel()
		nuTextField$getTextField()$setEnabled(method == "Robust (t)")
		ctablesCheckBoxes$getBoxes()$get(1L)$setEnabled(method %in% c("Moment", "MLE"))
	}
	methodListener <- new(ActionListener)
	methodListener$setFunction(toJava(methodListenerFunction))
	methodComboBox$addListener(methodListener)
	selectionListenerFunction <- function(cmd,ActionEvent){
		if (selectionButtons$getModel()=="Stepwise") niveauTextField$getTextField()$setEnabled(TRUE) else niveauTextField$getTextField()$setEnabled(FALSE)
	}
	selectionListener <- new(ActionListener)
	selectionListener$setFunction(toJava(selectionListenerFunction))
	selectionButtons$addListener(selectionListener)
	
	# Check function
	.ldaCheck <- function(state){
		if (length(state$Groups)==0) return("A grouping factor and at least one predictor must be selected")
		dataframe <-  get(state$data,.GlobalEnv)
		if (!is.factor(dataframe[[state$Groups]])) return("The grouping variable must be a factor")
		if (!all(sapply(dataframe[state$Predictors], "is.numeric"))) return("The predictors must be numeric variables")
		if (nlevels(dataframe[[state$Groups]]) < 2) return("The grouping factor must have at least two levels")
		return("")
	}

	# LDA Run function
	.ldaRun <- function(state){
		cmd <- character()
		# Define formula
		ldaFormula <- paste(state$Groups, "~", paste(state$Predictors,collapse="+"), sep="")
		# Options
		optionsList <- character()
		if (state$prior == "Equal"){
			ng <- nlevels(get(state$data,.GlobalEnv)[[state$Groups]])
			optionsList <- c(optionsList, paste("prior=rep(1,", ng, ")/", ng, sep=""))
		}
		optionsList <- c(optionsList, paste("method=\"", switch(state$method,
			Moment="moment", MLE="mle", MVE="mve", "Robust (t)"="t"), "\"", sep=""))
		if (state$method == "Robust (t)") optionsList <- c(optionsList, paste("nu=", state$nu, sep=""))
		optionsList <- c(optionsList, paste("tol=", state$tol, sep=""))
		optionsList <- paste(optionsList, collapse=", ")
		# Selection method
		if (state$selection == "Stepwise"){
			require(klaR)
			gwcmd <- paste("varSelection <- greedy.wilks(",	ldaFormula, ", data=", state$data,
				", niveau=", state$niveau,")\n", sep="")
			gwcmd <- paste(gwcmd, "varSelection\n")
			cmd <- paste(cmd, gwcmd)
			ldacmd <- paste("ldaModel <- lda(formula(varSelection), data=", state$data,
			if (length(optionsList)>0) ", ", optionsList, ")\n", sep="")
		} else {
			ldacmd <- paste("ldaModel <- lda(", ldaFormula, ", data=", state$data,
			if (length(optionsList)>0) ", ", optionsList, ")\n", sep="")
		}
		# Lda command
		cmd <- paste(cmd, ldacmd)
		cmd <- paste(cmd, "ldaModel\n")
		# Self-classification needed
		if (any(c( "Original data" %in% state$ctables,
			c("Centroids", "Statistics", "Structre matrix") %in% state$dfs,
			c("Classification","A posteriori probabilities","Discriminant functions") %in% state$export ))){
			cmd <- paste(cmd, paste("ldaPredict <- predict(ldaModel, ", state$data, ")\n"))
			cmd <- paste(cmd, "ldaFunctions <- ldaPredict$x\n")
			cmd <- paste(cmd, "ldaGroups <- model.frame(ldaModel)[,1]\n")
			ldaClassif <- TRUE
		}else ldaClassif <- FALSE
		if (any(c("Original data", "Cross-validation") %in% state$ctables)){
			cmd <- paste(cmd, "ldaClassification <- list()\n")
			if ("Original data" %in% state$ctables){
				cmd <- paste(cmd, "ldaClassification$\"From original data\" <- table(ldaGroups, ldaPredict$class)\n")
				cmd <- paste(cmd, "ldaClassification$\"From original data - proportions\" <- round(addmargins(prop.table(ldaClassification$\"From original data\")), 2)\n")
				cmd <- paste(cmd, "ldaClassification$\"From original data\" <- addmargins(ldaClassification$\"From original data\")\n")
			}
			if ("Cross-validation" %in% state$ctables && methodComboBox$getModel() %in% c("Moment","MLE")){
				cmd <- paste(cmd, "ldaCV <- update(ldaModel, CV=TRUE)\n")
				cmd <- paste(cmd, "ldaClassification$\"Cross-validation\" <- table(ldaGroups, ldaCV$class)\n")
				cmd <- paste(cmd, "ldaClassification$\"Cross-validation - proportions\" <- round(addmargins(prop.table(ldaClassification$\"Cross-validation\")), 2)\n")
				cmd <- paste(cmd, "ldaClassification$\"Cross-validation\" <- addmargins(ldaClassification$\"Cross-validation\")\n")
				cmd <- paste(cmd, "rm(ldaCV)\n")
			}
			cmd <- paste(cmd, "print(ldaClassification)\n")
			cmd <- paste(cmd, "rm(ldaClassification)\n")
		}
		if ("Centroids" %in% state$dfs){
			cmd <- paste(cmd, "ldaCentroids <- by(ldaFunctions, ldaGroups, colMeans)\n")
			cmd <- paste(cmd, "ldaCentroids\n")
			cmd <- paste(cmd, "rm(ldaCentroids)\n")
		}
		if ("Statistics" %in% state$dfs){
			cmd <- paste(cmd, "summary.aov(lm(ldaFunctions ~ ldaGroups))\n")
			cmd <- paste(cmd, "ldaMSS <- ldaModel$svd^2 * (length(ldaModel$lev)-1) / (ldaModel$N-length(ldaModel$lev))\n")
			cmd <- paste(cmd, paste("ldaStatTable <- matrix(, nrow=2, ncol=ncol(ldaFunctions),\n",
				"dimnames=list(c(\"Wilks' lambda\", \"Canonical correlation\"), colnames(ldaFunctions)))\n"))
			cmd <- paste(cmd, "ldaStatTable[1,] <- 1/(1+ldaMSS)\n")
			cmd <- paste(cmd, "ldaStatTable[2,] <- sqrt(1-ldaStatTable[1,])\n")
			cmd <- paste(cmd, "ldaStatTable\nrm(ldaMSS,ldaStatTable)\n")
		}
		if ("Structure matrix" %in% state$dfs){
			cmd <- paste(cmd, "ldaStructureMatrix <- cor(model.frame(ldaModel)[,-1],ldaFunctions)\n")
			cmd <- paste(cmd, "ldaStructureMatrix\n")
			cmd <- paste(cmd, "rm(ldaStructureMatrix)\n")
		}
		if ("Plot" %in% state$dfs){
			if (state$selection == "Stepwise") eval(parse(text=gwcmd))
			ldaModel <- eval(parse(text=ldacmd))
			n <- length(ldaModel$svd)
			cmd <- paste(cmd, switch(n,
				"1" = "ldahist(ldaFunctions, ldaGroups, type=\"both\")\n",
				"2" = "scatterplot(ldaFunctions[,2] ~ ldaFunctions[,1] | ldaGroups, smooth=FALSE, reg.line=FALSE, xlab=colnames(ldaFunctions)[1], ylab=colnames(ldaFunctions)[2])\n",
				"pairs(ldaModel,type=\"trellis\")\n"))
		}
		if (ldaClassif) cmd <- paste(cmd, "rm(ldaPredict, ldaFunctions, ldaGroups)\n")
		if (any(c("Classification","A posteriori probabilities","Discriminant functions") %in% state$export)){
			cmd <- paste(cmd, paste("ldaExport <- predict(ldaModel, ", state$newdata, ")\n", sep=""))
			if ("Classification" %in% state$export){
				cmdbind <- paste(state$newdata, " <- cbind(", state$newdata, ", Classification = ldaExport$class)\n", sep="")
				cmd <- paste(cmd, cmdbind)
			}
			if ("A posteriori probabilities" %in% state$export){
				cmdbind <- paste(state$newdata, " <- cbind(", state$newdata, ", ldaExport$posterior)\n", sep="")
				cmd <- paste(cmd, cmdbind)
			}
			if ("Discriminant functions" %in% state$export){
				cmdbind <- paste(state$newdata, " <- cbind(", state$newdata, ", ldaExport$x)\n", sep="")
				cmd <- paste(cmd, cmdbind)
			}
			cmd <- paste(cmd, "rm(ldaExport)\n")
		}
		cmd <- paste(cmd, "rm(ldaModel)\n")
		execute(cmd)
	}

	# Link check and run functions
	dialog$setCheckFunction(toJava(.ldaCheck))
	dialog$setRunFunction(toJava(.ldaRun))
	return(dialog)
}

.onLoad <- function(libname, pkgname) { 

	#if deducer gui is not running, do minimal load
	deducerLoaded <- try(.deducer == .jnull(),silent=TRUE)
	if(inherits(deducerLoaded,"try-error") || deducerLoaded)
		return(NULL)

	.registerDialog("Factor analysis",.makeNewFactorAnalysisDialog)
	.registerDialog("Reliability",.makeNewReliabilityDialog)
	.registerDialog("Intraclass correlation",.makeICCDialog)
	.registerDialog("Linear discriminant analysis",.makeNewLDADialog)
	
	deducer.addMenu("Psych")
	deducer.addMenuItem("Factor analysis",,".getDialog('Factor analysis')$run()","Psych")
	deducer.addMenuItem("Reliability",,".getDialog('Reliability')$run()","Psych")	
	deducer.addMenuItem("Intraclass correlation",,".getDialog('Intraclass correlation')$run()","Psych")	
	deducer.addMenuItem("Linear discriminant analysis",,".getDialog('Linear discriminant analysis')$run()","Psych")	
	if(.windowsGUI){
		winMenuAdd("Psych")
		winMenuAddItem("Psych", "Factor analysis", "deducer('Factor analysis')")
		winMenuAddItem("Psych", "Reliability", "deducer('Reliability')")
		winMenuAddItem("Psych", "Intraclass correlation", "deducer('Intraclass correlation')")
		winMenuAddItem("Psych", "Linear discriminant analysis", "deducer('Linear discriminant analysis')")
	}else if(.jgr){
		if(exists("jgr.getMenuNames") && exists("jgr.insertMenu")){
			menus <- jgr.getMenuNames()
			index <- which(menus=="Packages & Data")
			if(length(index)==0) 
				index <- 1
			jgr.insertMenu("Psych",index)
		}else
			jgr.addMenu("Psych")
		
		jgr.addMenuItem("Psych", "Factor analysis", "deducer('Factor analysis')")
		jgr.addMenuItem("Psych", "Reliability", "deducer('Reliability')")
		jgr.addMenuItem("Psych", "Intraclass correlation", "deducer('Intraclass correlation')")
		jgr.addMenuItem("Psych", "Linear discriminant analysis", "deducer('Linear discriminant analysis')")
	}

}

