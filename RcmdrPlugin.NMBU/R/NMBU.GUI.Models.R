# $Id: NMBU.GUI.Models.R 35 2014-01-10 21:17:26Z khliland $

##
## GUI functions for the Model menu
##



## GUI tips
#
# Usual code structure:
#    1. Intitialise window and prepare graphical elements
#    2. onOK function contianing actions to perform
#       2.1 Collect values from GUI
#       2.2 Test if combination of values is usable
#       2.3 Perform main calculations, print, update models/data etc.
#    3. Set up GUI.
#       - tkgrid() adds elements. Explicit placement and width/heigth by colum, row, columnspan and rowspan
#       - Frames with graphical elements are safer than direct placement of elements due to version conflicts.
#       - dialogSuffix() defines the final size of the grid used for elements.

#########################################################
# Prediction for LDA/QDA, PCR/PLSR and linear regression
predictRegressionModel <- function(){
  # Boks som tar imot verdiene til alle preditorene, kommaseparert (navnerekkef?lge m? framkomme, f.eks. med formula).
  # Konfidensgrad for CI og PI.
  .activeModel <- ActiveModel()
  model.class <- justDoIt(paste("class(", .activeModel, ")", sep=""))[1]
  if(model.class == "mvr"){
    initializeDialog(title=gettextRcmdr("Prediction in multivariate regression"))
    ff <- justDoIt(paste("fparse(formula(", .activeModel, "))", sep=""))
    ncomp <- justDoIt(paste(.activeModel, "$ncomp", sep=""))
  } else {
    if(model.class == "lda" || model.class == "qda"){
      initializeDialog(title=gettextRcmdr("Prediction in discriminant analysis"))
      ff <- justDoIt(paste("colnames(", .activeModel, "$means)", sep=""))
    } else {
      if(model.class == "lm" || model.class == "lmm"){
        initializeDialog(title=gettextRcmdr("Prediction in linear regression"))
        ff <- justDoIt(paste("fparse(formula(", .activeModel, "))", sep=""))
      } else {
        initializeDialog(title=gettextRcmdr("Prediction in generalized linear model"))
        ff <- justDoIt(paste("fparse(formula(", .activeModel, "))", sep=""))
      }
    }}
  
  n.eff <- length(ff)
  effNames  <- list()
  effFrames <- list()
  effs      <- list()
  for(i in 1:n.eff){
    effNames[[i]] <- tclVar('0')
    effFrames[[i]] <- tkframe(top)
    effs[[i]] <- ttkentry(effFrames[[i]], width="7", textvariable=effNames[[i]])
  }
  if(model.class == "lm" || model.class == "lmm"){
    levelName   <- tclVar("0.95")
    levelFrame  <- tkframe(top)
    level       <- ttkentry(levelFrame, width="5", textvariable=levelName)
    extrapFrame <- tkframe(top)
  }
  if(model.class == "mvr"){
    compName   <- tclVar(ncomp)
    compFrame  <- tkframe(top)
    comp       <- ttkentry(compFrame, width="5", textvariable=compName)
    typeFrame  <- tkframe(top)
  }
  
  onOK <- function(){ # Actions to perform
    closeDialog()
    the.effs <- list()
    for(i in 1:n.eff){
      the.effs[[i]] <- tclvalue(effNames[[i]])
	  if(grepl(',',the.effs[[i]],fixed=TRUE)){
		errorCondition(recall=predictRegressionModel, message=gettextRcmdr("Use period (.) instead of comma (,) as decimal mark."))
	  }
    }
    if(model.class == "lm" || model.class == "lmm"){
      levs <- tclvalue(levelName)
      if(trim.blanks(levs) == gettextRcmdr("")){
        errorCondition(recall=predictRegressionModel, message=gettextRcmdr("Please specify confidence level."))
      }
    }
    if(model.class == "mvr"){
      the.comps <- tclvalue(compName)
      if(trim.blanks(the.comps) == gettextRcmdr("")){
        errorCondition(recall=predictRegressionModel, message=gettextRcmdr("Please specify the number of components to use."))
      }
      the.type <- as.character(tclvalue(typeVariable))
    }
    
    doItAndPrint(paste("attach(", ActiveDataSet(), ")", sep=""))
    command <- paste("predict(", .activeModel, ", data.frame(", sep="")
    for(i in 1:n.eff){
      command <- paste(command, "'", ff[[i]], "'", "=", the.effs[[i]], sep="")
      if(i<n.eff)
        command <- paste(command, ", ", sep="")
    }
    if(model.class == "mvr"){
      command <- paste(command, "), ncomp=", the.comps, ", type='", the.type, "')", sep="")
    } else {
      if(model.class == "lm" || model.class == "lmm"){
        command <- gsub("predict(",  "predict_CI_PI(", command, fixed=TRUE)
        extrap <- tclvalue(extrapVariable)
        if(extrap == gettextRcmdr("1")){
          command <- paste(command, "), level=", levs, ", xXXx=TRUE)", sep="")
        } else {
          command <- paste(command, "), level=", levs, ")", sep="")
        }
      } else {
        if(model.class == "glm"){
          command <- gsub("predict(",  "predict_link_response(", command, fixed=TRUE)
          command <- paste(command, "))", sep="")				
        } else {
          command <- paste(command, "))[1:2]", sep="")
        }
      }			
    }
    doItAndPrint(command)
    if((model.class == "lm" || model.class == "lmm") && extrap == gettextRcmdr("1")){
      doItAndPrint(paste("max(hatvalues(", .activeModel, ")) # h_max", sep=""))}
    # if(model.class == "lm"){
    # extrap <- tclvalue(extrapVariable)
    # if(extrap == gettextRcmdr("1")){
    # command <- paste("X.tmp <- model.matrix(",.activeModel, ");   x.tmp <- data.frame(", sep="")
    # for(i in 1:n.eff){
    # command <- paste(command, "'", ff[[i]], "'", "=", the.effs[[i]], sep="")
    # if(i<n.eff)
    # command <- paste(command, ", ", sep="")
    # }
    # justDoIt(paste(command, ")", sep=""))
    # justDoIt(paste("x.tmp <- cbind('(Intercept)'=1, as.matrix(x.tmp))", sep=""))
    # doItAndPrint(paste("xXXx(X.tmp,x.tmp) # h_00", sep=""))
    # doItAndPrint(paste("max(hatvalues(", .activeModel, ")) # h_max", sep=""))
    # justDoIt("rm(c('X.tmp','x.tmp'))")
    # }
    # }
    doItAndPrint(paste("detach(", ActiveDataSet(), ")", sep=""))
    tkfocus(CommanderWindow())
  }
  
  # Set up GUI
  OKCancelHelp(helpSubject="predict", model=TRUE)
  tkgrid(labelRcmdr(top,text='Specify value of explanatory variable(s)', fg="blue"),sticky="w", row=1, column=1, columnspan=max(n.eff,2))
  for(i in 1:n.eff){
    tkgrid(labelRcmdr(effFrames[[i]],text=gettextRcmdr(ff[i])),sticky="w")
    tkgrid(effs[[i]], sticky="w")
    tkgrid(effFrames[[i]], row=2, column=i)
  }
  tkgrid(labelRcmdr(top,text=''),sticky="w", row=3, column=1, columnspan=max(n.eff,2))
  
  n.row <- 4
  if(model.class == "lm" || model.class == "lmm"){
    tkgrid(labelRcmdr(levelFrame, text=gettextRcmdr("Confidence level:")), level, sticky="w")
    tkgrid(levelFrame, sticky="n", row=4, column=1, columnspan=max(floor(n.eff/2),1))
    checkBoxes(frame="extrapFrame", boxes=c("extrap"), initialValues=c("0"),
               labels=gettextRcmdr(c("Detect hidden extrapolation")))
    tkgrid(extrapFrame, sticky="n", row=4, column=max(2, ceiling(n.eff/2)), columnspan=n.eff-floor(n.eff/2))
    n.row <- 5
  }
  if(model.class == "mvr"){
    tkgrid(labelRcmdr(compFrame, text=gettextRcmdr("Number of components:")), comp, sticky="w")
    tkgrid(compFrame, sticky="n", row=4, column=1, columnspan=max(n.eff,2)) 
    radioButtons(name="type", buttons=c("response", "scores"), values=c("response", "scores"), initialValue="response",
                 labels=gettextRcmdr(c("Response", "Scores")), title=gettextRcmdr("Type of prediction"))
    tkgrid(typeFrame, sticky="n", row=5, column=1, columnspan=max(n.eff,2)) 
    n.row <- 6
  }
  tkgrid(buttonsFrame, stick="s", row=n.row, column=1, columnspan=max(n.eff,2))
  tkgrid.configure(helpButton, sticky="se")
  dialogSuffix(rows=n.row, columns=max(n.eff,2))
}

###################################
## Stepwise and subset methods
forwardAdd <- function(){
  .activeModel <- ActiveModel()
  availableTerms <- justDoIt(paste("attr(",.activeModel,"$terms,'term.labels')",sep=""))
  initializeDialog(title=gettextRcmdr("Stepwise forward selection"))
  
  # Boks for signifikans ved ny variabel
  alphaName <- tclVar("0.2")
  alphaFrame <- tkframe(top)
  fullFrame <- tkframe(top)
  alpha <- ttkentry(alphaFrame, width="10", textvariable=alphaName)
  xBox <- variableListBox(top, availableTerms, selectmode="multiple", title=gettextRcmdr("Compulsory variables (pick zero or more)"))
  onOK <- function(){ # Actions to perform
    x <- getSelection(xBox)
    closeDialog()
    the.alpha <- tclvalue(alphaName)
    if(trim.blanks(the.alpha) == gettextRcmdr("")){
      warning('Please specify the alpha level for entering')
    }
    fullmodel <- tclvalue(fullmodelVariable)
    force.in <- 'NULL'
    if(length(x)>0){
      force.in <- paste("'",paste(x, collapse="+"),"'", sep="")}
    command <- paste("forward(", .activeModel, ", alpha=", the.alpha, ", full=",ifelse(fullmodel == gettextRcmdr("1"), 'TRUE', 'FALSE'),", force.in=",force.in,")", sep="")
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  
  # Set up GUI
  OKCancelHelp(helpSubject="add", model=TRUE)
  tkgrid(labelRcmdr(alphaFrame, text=gettextRcmdr("Alpha to enter:")), alpha, sticky="w")
  tkgrid(alphaFrame, sticky="n")
  checkBoxes(frame="fullFrame", boxes=c("fullmodel"), initialValues=c("0"),
             labels=gettextRcmdr(c("Extended output")))
  tkgrid(fullFrame, sticky="n")
  tkgrid(getFrame(xBox), sticky = "n")
  tkgrid(buttonsFrame, stick="s")
  tkgrid.configure(helpButton, sticky="se")
  dialogSuffix(rows=3, columns=2)
}

backwardDrop <- function(){
  .activeModel <- ActiveModel()
  availableTerms <- justDoIt(paste("attr(",.activeModel,"$terms,'term.labels')",sep=""))
  initializeDialog(title=gettextRcmdr("Stepwise backward eliminiation"))
  
  # Boks for signifikans ved ny variabel
  alphaName <- tclVar("0.2")
  alphaFrame <- tkframe(top)
  fullFrame <- tkframe(top)
  alpha <- ttkentry(alphaFrame, width="10", textvariable=alphaName)
  xBox <- variableListBox(top, availableTerms, selectmode="multiple", title=gettextRcmdr("Compulsory variables (pick zero or more)"))
  onOK <- function(){ # Actions to perform
    x <- getSelection(xBox)
    closeDialog()
    the.alpha <- tclvalue(alphaName)
    if(trim.blanks(the.alpha) == gettextRcmdr("")){
      warning('Please specify the alpha level for elimination')
    }
    fullmodel <- tclvalue(fullmodelVariable)
    force.in <- 'NULL'
    if(length(x)>0){
      force.in <- paste("c('",paste(x, collapse="','"),"')", sep="")}
    command <- paste("backward(", .activeModel, ", alpha=", the.alpha, ", full=",ifelse(fullmodel == gettextRcmdr("1"), 'TRUE', 'FALSE'),", force.in=",force.in,")", sep="")
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  
  # Set up GUI
  OKCancelHelp(helpSubject="drop", model=TRUE)
  tkgrid(labelRcmdr(alphaFrame, text=gettextRcmdr("Alpha to remove:")), alpha, sticky="w")
  tkgrid(alphaFrame, sticky="n")
  checkBoxes(frame="fullFrame", boxes=c("fullmodel"), initialValues=c("0"),
             labels=gettextRcmdr(c("Extended output")))
  tkgrid(fullFrame, sticky="n")
  tkgrid(getFrame(xBox), sticky = "n")
  tkgrid(buttonsFrame, stick="s")
  tkgrid.configure(helpButton, sticky="se")
  dialogSuffix(rows=3, columns=2)
}

forwardBackward <- function(){
  .activeModel <- ActiveModel()
  initializeDialog(title=gettextRcmdr("Stepwise forward/backward selection"))
  
  # Boks for signifikans ved ny variabel
  alphaName <- tclVar("0.15")
  alphaFrame <- tkframe(top)
  alpha1Name <- tclVar("0.15")
  alpha1Frame <- tkframe(top)
  fullFrame <- tkframe(top)
  alpha <- ttkentry(alphaFrame, width="10", textvariable=alphaName)
  alpha1 <- ttkentry(alpha1Frame, width="10", textvariable=alpha1Name)
  onOK <- function(){ # Actions to perform
    closeDialog()
    the.alpha <- tclvalue(alphaName)
    if(trim.blanks(the.alpha) == gettextRcmdr("")){
      warning('Please specify the alpha level for entering')
    }
    the.alpha1 <- tclvalue(alpha1Name)
    if(trim.blanks(the.alpha1) == gettextRcmdr("")){
      warning('Please specify the alpha level for removing')
    }
    fullmodel <- tclvalue(fullmodelVariable)
    if(fullmodel == gettextRcmdr("1")){
      command <- paste("stepWise(", .activeModel, ", alpha.enter=", the.alpha, ", alpha.remove=", the.alpha1, ", full=TRUE)", sep="")
      doItAndPrint(command)
	} else{
      command <- paste("stepWise(", .activeModel, ", alpha.enter=", the.alpha, ", alpha.remove=", the.alpha1, ", full=FALSE)", sep="")
      doItAndPrint(command)}
    
    tkfocus(CommanderWindow())
  }
  
  # Set up GUI
  OKCancelHelp(helpSubject="add", model=TRUE)
  tkgrid(labelRcmdr(alphaFrame, text=gettextRcmdr("Alpha to enter:")), alpha, sticky="w")
  tkgrid(alphaFrame, sticky="n")
  tkgrid(labelRcmdr(alpha1Frame, text=gettextRcmdr("Alpha to remove:")), alpha1, sticky="w")
  tkgrid(alpha1Frame, sticky="n")
  checkBoxes(frame="fullFrame", boxes=c("fullmodel"), initialValues=c("0"),
             labels=gettextRcmdr(c("Extended output")))
  tkgrid(fullFrame, sticky="n")
  tkgrid(buttonsFrame, stick="s")
  tkgrid.configure(helpButton, sticky="se")
  dialogSuffix(rows=4, columns=2)
}

backwardForward <- function(){
  .activeModel <- ActiveModel()
  initializeDialog(title=gettextRcmdr("Stepwise backward/forward selection"))
  
  # Boks for signifikans ved ny variabel
  alphaName <- tclVar("0.15")
  alphaFrame <- tkframe(top)
  alpha1Name <- tclVar("0.15")
  alpha1Frame <- tkframe(top)
  fullFrame <- tkframe(top)
  alpha <- ttkentry(alphaFrame, width="10", textvariable=alphaName)
  alpha1 <- ttkentry(alpha1Frame, width="10", textvariable=alpha1Name)
  onOK <- function(){ # Actions to perform
    closeDialog()
    the.alpha <- tclvalue(alphaName)
    if(trim.blanks(the.alpha) == gettextRcmdr("")){
      warning('Please specify the alpha level for removing')
    }
    the.alpha1 <- tclvalue(alpha1Name)
    if(trim.blanks(the.alpha1) == gettextRcmdr("")){
      warning('Please specify the alpha level for entering')
    }
    fullmodel <- tclvalue(fullmodelVariable)
    if(fullmodel == gettextRcmdr("1")){
      command <- paste("stepWiseBack(", .activeModel, ", alpha.remove=", the.alpha, ", alpha.enter=", the.alpha1, ", full=TRUE)", sep="")
      doItAndPrint(command)}
    else{
      command <- paste("stepWiseBack(", .activeModel, ", alpha.remove=", the.alpha, ", alpha.enter=", the.alpha1, ", full=FALSE)", sep="")
      doItAndPrint(command)}
    
    tkfocus(CommanderWindow())
  }
  
  # Set up GUI
  OKCancelHelp(helpSubject="drop", model=TRUE)
  tkgrid(labelRcmdr(alphaFrame, text=gettextRcmdr("Alpha to remove:")), alpha, sticky="w")
  tkgrid(alphaFrame, sticky="n")
  tkgrid(labelRcmdr(alpha1Frame, text=gettextRcmdr("Alpha to enter:")), alpha1, sticky="w")
  tkgrid(alpha1Frame, sticky="n")
  checkBoxes(frame="fullFrame", boxes=c("fullmodel"), initialValues=c("0"),
             labels=gettextRcmdr(c("Extended output")))
  tkgrid(fullFrame, sticky="n")
  tkgrid(buttonsFrame, stick="s")
  tkgrid.configure(helpButton, sticky="se")
  dialogSuffix(rows=4, columns=2)
}

bestSubsets <- function(){
  .activeModel <- ActiveModel()
  availableTerms <- justDoIt(paste("attr(",.activeModel,"$terms,'term.labels')",sep=""))
  initializeDialog(title=gettextRcmdr("Best subset regression"))
  
  # Bokser for antall
  nbestName <- tclVar("3")
  nbestFrame <- tkframe(top)
  nvmaxName <- tclVar("3")
  nvmaxFrame <- tkframe(top)
  nbest <- ttkentry(nbestFrame, width="10", textvariable=nbestName)
  nvmax <- ttkentry(nvmaxFrame, width="10", textvariable=nvmaxName)
  # Valg av faste variabler
  xBox <- variableListBox(top, availableTerms, selectmode="multiple", title=gettextRcmdr("Compulsory variables (pick zero or more)"))
  onOK <- function(){ # Actions to perform
    x <- getSelection(xBox)
    closeDialog()
    the.nbest <- tclvalue(nbestName)
    if(trim.blanks(the.nbest) == gettextRcmdr("")){
      warning('Please specify the number of models per model size')
    }
    the.nvmax <- tclvalue(nvmaxName)
    if(trim.blanks(the.nvmax) == gettextRcmdr("")){
      warning('Please specify the maximum model size')
    }
    force.in <- 'NULL'
    if(length(x)>0){
      force.in <- paste("'c(",paste(is.element(availableTerms,x), collapse=","),")'", sep="")}
    command <- paste("best.subsets(", .activeModel, ", nbest=", the.nbest, ", nvmax=", the.nvmax, ", force.in=", force.in,")", sep="")
    doItAndPrint(command)
    
    tkfocus(CommanderWindow())
  }
  
  # Set up GUI
  OKCancelHelp(helpSubject="regsubsets", model=TRUE)
  tkgrid(labelRcmdr(nbestFrame, text=gettextRcmdr("Number of models per model size:")), nbest, sticky="w")
  tkgrid(nbestFrame, sticky="n")
  tkgrid(labelRcmdr(nvmaxFrame, text=gettextRcmdr("Maximum model size:")), nvmax, sticky="w")
  tkgrid(nvmaxFrame, sticky="n")
  tkgrid(getFrame(xBox), sticky = "n")
  tkgrid(buttonsFrame, stick="s")
  tkgrid.configure(helpButton, sticky="se")
  dialogSuffix(rows=4, columns=2)
}

################################
## 2D discrinimant plot
discriminantPlot <- function(){
  initializeDialog(title=gettextRcmdr("Two-dimensional discriminant plot"))
  boxFrame <- tkframe(top)
  sliderFrame <- tkframe(top)
  sliderValue <- tclVar("100")
  componentsSlider <- tkscale(sliderFrame, from=50, to=500, showvalue=FALSE, variable=sliderValue,
                              resolution=50, orient="horizontal")
  componentsShow <- labelRcmdr(sliderFrame, textvariable=sliderValue, width=3, justify="right")
  
  onOK <- function(){ # Actions to perform
    the.slider <- as.numeric(tclvalue(sliderValue))
    closeDialog()
    the.regions <- tclvalue(regionsVariable)
    the.contours <- tclvalue(contoursVariable)
    if(the.regions == gettextRcmdr("1") && the.contours == gettextRcmdr("1")){
      command <- paste("plotDA(regions=TRUE, contours=TRUE, resolution=", the.slider, ")", sep="")
    }
    if(the.regions == gettextRcmdr("0") && the.contours == gettextRcmdr("1")){
      command <- paste("plotDA(regions=FALSE, contours=TRUE, resolution=", the.slider, ")", sep="")
    }
    if(the.regions == gettextRcmdr("1") && the.contours == gettextRcmdr("0")){
      command <- paste("plotDA(regions=TRUE, contours=FALSE, resolution=", the.slider, ")", sep="")
    }
    if(the.regions == gettextRcmdr("0") && the.contours == gettextRcmdr("0")){
      command <- paste("plotDA(regions=FALSE, contours=FALSE, resolution=", the.slider, ")", sep="")
    }
    logger(command)
    justDoIt(command)
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="plot", model=TRUE)
  tkgrid(labelRcmdr(top, text=gettextRcmdr("Plot resolution:"), fg="blue"), sticky="w")    
  tkgrid(componentsSlider, componentsShow, sticky="nw")
  tkgrid(sliderFrame, sticky="w")
  checkBoxes(frame="boxFrame", boxes=c("regions","contours"), initialValues=c("1","0"),
             labels=gettextRcmdr(c("Show discriminant regions", "Show binormal contours")))
  tkgrid(boxFrame, sticky="w")
  tkgrid(buttonsFrame, stick="s")
  tkgrid.configure(helpButton, sticky="se")
  dialogSuffix(rows=5, columns=2)
}

#####################################
# ANOVA contrasts
contrastGUI <- function(){
  .activeModel <- ActiveModel()
  levs <- justDoIt(paste(.activeModel, "$xlevels", sep=""))
  if(length(levs)==0 || length(levs)>2){
    return()
  }
  if(length(levs)>1){
    initializeDialog(title=gettextRcmdr("Test contrasts in ANOVA (choose factor)"))
    xFrame <- tkframe(top)
    xBox <- variableListBox(xFrame, names(levs),
                            title=gettextRcmdr("Factors (pick one)"))
    onOK <- function(){ # Actions to perform
      x <- getSelection(xBox)
      closeDialog()
      if (length(x)==0) {
        errorCondition(recall=contrastGUI, 
                       message=gettextRcmdr("No variables selected."))
        return()
      }
      contrastGUI2(which(names(levs)==x),levs)
    }
    OKCancelHelp(helpSubject="fit.contrast", model=TRUE)
    tkgrid(getFrame(xBox), sticky="nw")
    tkgrid(xFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    tkgrid.configure(helpButton, sticky="se")
    dialogSuffix(rows=2, columns=1)
  } else{
    contrastGUI2(1,levs)
  }
}
contrastGUI2 <- function(varNumb,levs){
  # Boks som tar imot verdiene til alle prediktorene, kommaseparert (navnerekkef?lge m? framkomme, f.eks. med formula).
  # Konfidensgrad for CI og PI.
  initializeDialog(title=gettextRcmdr("Test contrasts in ANOVA"))
  
  .activeModel <- ActiveModel()
  ff <- justDoIt(paste(.activeModel, "$xlevels[[", varNumb, "]]", sep=""))
  fn <- justDoIt(paste("names(", .activeModel, "$xlevels)[[", varNumb, "]]", sep=""))
  
  n.eff <- length(ff)
  effNames  <- list()
  effFrames <- list()
  effs      <- list()
  if(n.eff>2){
    for(i in 1:n.eff){
      effNames[[i]] <- tclVar('0')
      effFrames[[i]] <- tkframe(top)
      effs[[i]] <- ttkentry(effFrames[[i]], width=5, textvariable=effNames[[i]])
    }
  } else{
    effNames[[1]] <- tclVar('1')
    effNames[[2]] <- tclVar('-1')
    effFrames[[1]] <- tkframe(top)
    effFrames[[2]] <- tkframe(top)
    effs[[1]] <- ttkentry(effFrames[[1]], width=5, textvariable=effNames[[1]])
    effs[[2]] <- ttkentry(effFrames[[2]], width=5, textvariable=effNames[[2]])
  }
  
  levelName <- tclVar("0.95")
  levelFrame <- tkframe(top)
  level <- ttkentry(levelFrame, width="10", textvariable=levelName)
  
  onOK <- function(){ # Actions to perform
    closeDialog()
    effValues <- numeric(n.eff)
    for(i in 1:n.eff){
      effValues[i] <- try(eval(parse(text=(tclvalue(effNames[[i]])))))
    }
    if(abs(sum(effValues)) >= 0.1){
      errorCondition(recall=contrastGUI2, message=gettextRcmdr("Contrasts must sum to 0 (total deviance allowed=0.1)."))
      return()
    }
    levelValue <- as.numeric(tclvalue(levelName))
    if(levelValue <= 0 || levelValue >= 1){
      errorCondition(recall=contrastGUI2, message=gettextRcmdr("Confidence level must be >0 and <1."))
      return()
    }
    # level.name <- justDoIt(paste("names(",.activeModel, "$xlevels)", sep=""))
	Library("gmodels")
    command <- paste("print(fit.contrast(model=", ActiveModel(), ", varname='",fn, "', df=TRUE, coeff=c(", effValues[1], sep="")
    for(i in 2:n.eff){
      command <- paste(command, ", ", effValues[i], sep="")
    }
    command <- paste(command, "), conf.int=", levelValue,"))", sep="")
    
    # If interaction, test per level, not across levels
    has.interaction <- justDoIt(paste("max(apply(attr(",.activeModel, "$terms,'factors'),1,sum))>1", sep=""))
    if(length(levs)>1 && has.interaction){
      the.levels <- justDoIt(paste("levels(",ActiveDataSet(),"$",names(levs[3-varNumb]),")",sep=""))
      for(i in 1:length(the.levels)){
        doItAndPrint(paste(paste("cat('Contrast for ",names(levs[3-varNumb]), " = ", the.levels[1], "\n')", sep=""),command, sep="; "))
        #				doItAndPrint(command)
        the.levels <- c(the.levels[-1],the.levels[1])
        doItAndPrint(paste(ActiveDataSet(),"$",names(levs[3-varNumb]), " <- factor(", ActiveDataSet(),"$",names(levs[3-varNumb]), ", levels=c('",paste(the.levels,sep="",collapse="','"),"'))",sep=""))
      }
    } else {
      doItAndPrint(command)
    }
    #logger(command)
    tkfocus(CommanderWindow())
  }
  
  # Set up GUI
  OKCancelHelp(helpSubject="fit.contrast", model=TRUE)
  tkgrid(labelRcmdr(top,text=paste('Specify contrast for ', fn, sep=""), fg="blue"),sticky="w", row=1, column=1, columnspan=n.eff)
  for(i in 1:n.eff){
    tkgrid(labelRcmdr(effFrames[[i]],text=gettextRcmdr(ff[i])),sticky="w")
    tkgrid(effs[[i]], sticky="w")
    tkgrid(effFrames[[i]], row=2, column=i)
  }
  tkgrid(labelRcmdr(top,text=''),sticky="w", row=3, column=1, columnspan=n.eff)
  tkgrid(labelRcmdr(levelFrame, text=gettextRcmdr("Confidence level:")), level, sticky="w")
  tkgrid(levelFrame, sticky="w", row=4, column=1, columnspan=n.eff)
  tkgrid(buttonsFrame, stick="s", row=5, column=1, columnspan=n.eff)
  tkgrid.configure(helpButton, sticky="se")
  dialogSuffix(rows=5, columns=n.eff)
}

#####################################
# Post hoc tests
postHocGUI <- function(){
  initializeDialog(title=gettextRcmdr("Post hoc pair-wise tests"))
  .activeModel <- ActiveModel()
  effects <- eval(parse(text=paste("attr(terms(formula(",.activeModel,")),'term.labels')", sep="")))
#  if(glmP()){
#	effects <- effects[!grepl(":",effects)]
#  }
  xFrame <- tkframe(top)
#  if(glmP()){
    xBox <- variableListBox(xFrame, effects,
                            title=gettextRcmdr("Effects (pick one)"))
 # } else {
  #  xBox <- variableListBox(xFrame, effects, selectmode="multiple",
   #                         title=gettextRcmdr("Effects (pick one or more)"))
  #}
  comboFrame <- tkframe(top)
  levs <- justDoIt(paste(.activeModel, "$xlevels", sep=""))
  levNames <- names(levs)
  values <- NULL
  valuesLookUp <- NULL
  comboVar <- tclVar()
  for(i in 1:length(levNames)){
    values <- c(values, paste(levNames[i], levs[[i]], sep=", "))
    valuesLookUp <- rbind(valuesLookUp, cbind(rep(levNames[i],length(levs[[i]])),levs[[i]]))
  }
  onOK <- function(){
	x <- getSelection(xBox)
	selectedDunnett <- which(values==tclvalue(comboVar))
	if(length(x)==0 && length(selectedDunnett)==0){
      errorCondition(recall=postHocGUI, message=gettextRcmdr("Pick at least one effect for pair-wise comparisons other than Dunnett."))
      return()
    }
	if(length(x)==1){
		selected <- paste("'",x,"'",sep="")
	} else {
		selected <- paste("c('",paste(x,sep="",collapse="','"),"')", sep="")
	}

    the.tukey <- tclvalue(tukeyName)
    if(tclvalue(tukeyTestsVariable)== gettextRcmdr("1")){
		command <- paste("print(simple.glht(", ActiveModel(), ",",selected,", level=", the.tukey,"))",sep="")
#	  if(glmP()){
#		command <- paste("summary(glht(", ActiveModel(), ", linfct=mcp(", x, "='Tukey')))", sep="")
#	  } else {
#        command <- paste("TukeyHSD(", ActiveModel(), ",",selected,", conf.level=",the.tukey,")",sep="")
#	  }
      doItAndPrint(command)
    } 
    if(tclvalue(tukeyGroupsVariable)== gettextRcmdr("1")){
		command <- paste("cld(simple.glht(", ActiveModel(), ",",selected,", level=", the.tukey,"))",sep="")
#	  if(glmP()){
#		command <- paste("multcomp:::print.cld(multcomp::cld(glht(", ActiveModel(), ", linfct=mcp(", x, "='Tukey'))))", sep="")
#	  } else {
#		command <- paste("cld(TukeyHSD(", ActiveModel(), ", ", selected, "), ", 1-as.numeric(the.tukey),")", sep="")
#	  }
      doItAndPrint(command)
    }

#	if(!glmP()){
		the.bonferroni <- tclvalue(bonferroniName)
		if(tclvalue(bonferroniTestsVariable)== gettextRcmdr("1")){
		  command <- paste("print(simple.glht(", ActiveModel(), ",",selected,", corr='Bonferroni', level=", the.tukey,"))",sep="")
#		  command <- paste("Bonferroni(", ActiveModel(), ",",selected,", conf.level=",the.bonferroni,")",sep="")
		  doItAndPrint(command)
		}
		if(tclvalue(bonferroniGroupsVariable)== gettextRcmdr("1")){
		  command <- paste("cld(simple.glht(", ActiveModel(), ",",selected,", corr='Bonferroni', level=", the.tukey,"))",sep="")
		  doItAndPrint(command)
#		  doItAndPrint(paste("cld(Bonferroni(", ActiveModel(), ", ", selected, "), ", 1-as.numeric(the.bonferroni), ")", sep=""))
		}

		the.fisher <- tclvalue(fisherName)
		if(tclvalue(fisherTestsVariable)== gettextRcmdr("1")){
		  command <- paste("print(simple.glht(", ActiveModel(), ",",selected,", corr='Fisher', level=", the.tukey,"))",sep="")
#		  command <- paste("Fisher(", ActiveModel(), ",",selected,", conf.level=",the.fisher,")",sep="")
		  doItAndPrint(command)
		}
		if(tclvalue(fisherGroupsVariable)== gettextRcmdr("1")){
		  command <- paste("cld(simple.glht(", ActiveModel(), ",",selected,", corr='Fisher', level=", the.tukey,"))",sep="")
		  doItAndPrint(command)
#		  doItAndPrint(paste("cld(Fisher(", ActiveModel(), ", ", selected, "), ", 1-as.numeric(the.fisher), ")", sep=""))
		}
		
		if(length(selectedDunnett)>0){
		  effect <- valuesLookUp[selectedDunnett,1]
		  level  <- valuesLookUp[selectedDunnett,2]
		  command <- paste("levels(", ActiveModel(), "$model[,'", effect, "'])", sep="")
		  effLevs.orig <- effLevs <- justDoIt(command)
		  if(effLevs[1]!=level){
			effLevs <- c(level, effLevs[!is.element(effLevs,level)])
			command <- paste(ActiveDataSet(), "$", effect, " <- factor(", ActiveDataSet(), "$", effect, ", levels=c('", paste(effLevs,sep="", collapse="', '"), "'))", sep="")
			doItAndPrint(command)
			doItAndPrint(paste(ActiveModel(), " <- update(",ActiveModel(),")",sep=""))
		  }
		  command <- paste("summary(glht(", ActiveModel(), ", linfct = mcp(", effect, " = 'Dunnett')))", sep="")
		  doItAndPrint(command)
		  if(effLevs.orig[1]!=level){
			command <- paste(ActiveDataSet(), "$", effect, " <- factor(", ActiveDataSet(), "$", effect, ", levels=c('", paste(effLevs.orig,sep="", collapse="', '"), "'))", sep="")
			doItAndPrint(command)
			doItAndPrint(paste(ActiveModel(), " <- update(",ActiveModel(),")",sep=""))
		  }
		}
#	}
    closeDialog()
    return()
  }
  OKCancelHelp(helpSubject="simple.glht")
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(xFrame, row=1, column=1, columnspan=2, sticky="w")
  
  topFrame <- tkframe(top)
  tkgrid(labelRcmdr(topFrame, text=gettextRcmdr("Tukey HSD"), fg="blue"),sticky="w")
  tkgrid(topFrame, row=2, column=1, columnspan=1, sticky="w")
  checkBoxes(frame="tukeyFrame", boxes=c("tukeyTests","tukeyGroups"), initialValues=c("0","0"),
             labels=gettextRcmdr(c("Tests","Groups")))
  tkgrid(tukeyFrame, row=3, column=1, columnspan=1, rowspan=1, sticky="w")
  tukeyName  <- tclVar("0.95")
  tukeyFrameGr <- tkframe(top)
  tukey <- ttkentry(tukeyFrameGr, width="6", textvariable=tukeyName)
  tkgrid(labelRcmdr(tukeyFrameGr, text=gettextRcmdr("Confidence")), tukey, sticky="w")
  tkgrid(tukeyFrameGr, sticky="w", row=3, column=2)
  
#  if(!glmP()){
	  topFrame2 <- tkframe(top)
	  tkgrid(labelRcmdr(topFrame2, text=gettextRcmdr("Bonferroni"), fg="blue"),sticky="w")
	  tkgrid(topFrame2, row=4, column=1, columnspan=1, sticky="w")
	  checkBoxes(frame="bonferroniFrame", boxes=c("bonferroniTests","bonferroniGroups"), initialValues=c("0","0"),
				 labels=gettextRcmdr(c("Tests","Groups")))
	  tkgrid(bonferroniFrame, row=5, column=1, columnspan=1, rowspan=1, sticky="w")
	  bonferroniName  <- tclVar("0.95")
	  bonferroniFrameGr <- tkframe(top)
	  bonferroni <- ttkentry(bonferroniFrameGr, width="6", textvariable=bonferroniName)
	  tkgrid(labelRcmdr(bonferroniFrameGr, text=gettextRcmdr("Confidence")), bonferroni, sticky="w")
	  tkgrid(bonferroniFrameGr, sticky="w", row=5, column=2)
	  
	  topFrame3 <- tkframe(top)
	  tkgrid(labelRcmdr(topFrame3, text=gettextRcmdr("Fisher"), fg="blue"),sticky="w")
	  tkgrid(topFrame3, row=6, column=1, columnspan=1, sticky="w")
	  checkBoxes(frame="fisherFrame", boxes=c("fisherTests","fisherGroups"), initialValues=c("0","0"),
				 labels=gettextRcmdr(c("Tests","Groups")))
	  tkgrid(fisherFrame, row=7, column=1, columnspan=1, rowspan=1, sticky="w")
	  fisherName  <- tclVar("0.95")
	  fisherFrameGr <- tkframe(top)
	  fisher <- ttkentry(fisherFrameGr, width="6", textvariable=fisherName)
	  tkgrid(labelRcmdr(fisherFrameGr, text=gettextRcmdr("Confidence")), fisher, sticky="w")
	  tkgrid(fisherFrameGr, sticky="w", row=7, column=2)
	  
	  topFrame1 <- tkframe(top)
	  tkgrid(labelRcmdr(topFrame1, text=gettextRcmdr("Dunnet"), fg="blue"),sticky="w")
	  tkgrid(topFrame1, row=8, column=1, columnspan=1, sticky="w")
	  combo <- ttkcombobox(comboFrame, values=values, textvariable=comboVar)
	  tkgrid(labelRcmdr(comboFrame, text=gettextRcmdr("Choose comparison:")), combo, sticky="w")
	  tkgrid(comboFrame, sticky="w", column=1, row=9, columnspan=2)
#  }
  tkgrid(buttonsFrame, sticky="w", column=1, row=10, columnspan=2)
  dialogSuffix(rows=10, columns=2)
}


#####################################
# Plot mixture experiment
mixtureGUI <- function(){
  Library("lattice")
  initializeDialog(title=gettextRcmdr("Plot response surface for mixture design"))
  .numeric <- Numeric()
  comboVar <- tclVar()
  comboFrame <- tkframe(top)
  comboVar2 <- tclVar()
  comboFrame2 <- tkframe(top)
  formatFrame <- tkframe(top)
  zoomFrame <- tkframe(top)
  onOK <- function(){ # Actions to perform
    mix.format <- as.character(tclvalue(label.formatVariable))
    zoomed <- tclvalue(zoomedVariable)
    if(zoomed == gettextRcmdr("1")){
      zoomed <- "TRUE"
    } else {
      zoomed <- "FALSE"
    }
    show.points <- tclvalue(samplesVariable)
    if(show.points == gettextRcmdr("1")){
      show.points <- "TRUE"
    } else {
      show.points <- "FALSE"
    }
    logen <- ifelse(tclvalue(logVariable)== gettextRcmdr("1"), TRUE, FALSE)
    n.ticks <- which(as.character(2:10)==tclvalue(comboVar))+1
    if(length(n.ticks)==0) n.ticks <- "6"
    n.grade <- which(as.character(seq(5,25,5))==tclvalue(comboVar2))*5
    if(length(n.grade)==0) n.grade <- "15"
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    .activeModel <- ActiveModel()
    formula1 <- justDoIt(paste("formula(",.activeModel, ")", sep=""))
    if(logen){
      doItAndPrint(paste("mixture.contour(", .activeDataSet, ", ", paste(formula1[2],formula1[1],formula1[3],sep=" "), ", n.tick=", n.ticks, ", n.grade=", n.grade,", FUN='exp', FUN2='log', mix.format='",mix.format,"', show.points=",show.points,", zoomed=",zoomed,ifelse(show.points,", pch=21, cex=1.25, col.points='black', fill.points='white'",""),")", sep=""))
    } else {
      doItAndPrint(paste("mixture.contour(", .activeDataSet, ", ", paste(formula1[2],formula1[1],formula1[3],sep=" "), ", n.tick=", n.ticks, ", n.grade=", n.grade,", mix.format='",mix.format,"', show.points=",show.points,", zoomed=",zoomed,ifelse(show.points,", pch=21, cex=1.25, col.points='black', fill.points='white'",""),")", sep=""))
    }
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="lm")
  combo <- ttkcombobox(comboFrame, values=as.character(2:10), textvariable=comboVar, width=3)
  combo2 <- ttkcombobox(comboFrame2, values=as.character(seq(5,25,5)), textvariable=comboVar2, width=3)
  tkgrid(labelRcmdr(comboFrame, text=gettextRcmdr("Plot ticks (default=6):")), combo, sticky="w")
  tkgrid(labelRcmdr(comboFrame2, text=gettextRcmdr("Plot gradings (approximate, default=15):")), combo2, sticky="w")
  tkgrid(comboFrame, sticky="w", column=1, row=1, columnspan=1)
  tkgrid(comboFrame2, sticky="w", column=1, row=2, columnspan=1)
  checkBoxes(frame="logFrame", boxes=c("log"), initialValues=c("0"),
             labels=gettextRcmdr(c("log transformed response")))
  tkgrid(logFrame, row=3, column=1, columnspan=1, rowspan=1, sticky="w")
  radioButtonsNMBU(formatFrame,name="label.format", buttons=c("frac", "dec"), values=c("frac", "dec"), initialValue = "frac",
                  labels=gettextRcmdr(c("Fraction", "Decimal")))
  tkgrid(formatFrame, row=4, column=1, columnspan=1, rowspan=1, sticky="w")
  checkBoxes(frame="zoomFrame", boxes=c("samples","zoomed"), initialValues=c("0"), labels=gettextRcmdr(c("Show samples", "Zoom on samples")))
  tkgrid(zoomFrame, row=5, column=1, columnspan=1, rowspan=1, sticky="w")
  tkgrid(buttonsFrame, sticky="w", row=6, column=1, columnspan=1)
  dialogSuffix(rows=6, columns=1)
}

#####################################
# ANOVA for regression
anova_reg_GUI <- function(){
  doItAndPrint(paste("anova_reg(", ActiveModel(), ")", sep=""))
}

#####################################
# PRESS
PRESS.GUI <- function(){
  doItAndPrint(paste("PRESS(", ActiveModel(), ")", sep=""))
  doItAndPrint(paste("R2pred(", ActiveModel(), ")", sep=""))
  doItAndPrint(paste("rmsep(", ActiveModel(), ")", sep=""))
}

######################################
## Customized observation statistics
addObservationStatisticsNMBU <- function(){
  .activeDataSet <- ActiveDataSet()
  .activeModel <- ActiveModel()
  if (is.null(.activeModel)) return()
  addVariable <- function(name, type=NULL){
    variable <- paste(name, ".", ifelse(is.null(type),"",paste(type,".",sep="")), .activeModel, sep="")
    if (is.element(variable, .variables)) {
      ans <- checkReplace(variable)
      if (tclvalue(ans) == "no") return()
    }
    command <- paste(name, "(", .activeModel, ifelse(is.null(type),"",paste(", type='",type,"'",sep="")), ")", sep="")
    if(is.matrix(justDoIt(command))){
      justDoIt(paste(.activeDataSet, " <- cbind(", .activeDataSet ,", extend.colnames(", command,",'", variable, "'))", sep=""))
      logger(paste(.activeDataSet, " <- cbind(", .activeDataSet ,", extend.colnames(", command,",'", variable, "'))", sep=""))
    } else {
      justDoIt(paste(.activeDataSet, "$", variable, " <- ", command, sep=""))
      logger(paste(.activeDataSet, "$", variable, " <- ", command, sep=""))
    }
  }
  if (getRcmdr("modelWithSubset")){
    Message(message=
              gettextRcmdr("Observation statistics not available\nfor a model fit to a subset of the data."),
            type="error")
    tkfocus(.commander)
    return()
  }
  initializeDialog(title=gettextRcmdr("Add Observation Statistics to Data"))
  .activeModel <- ActiveModel()
  .activeDataSet <- ActiveDataSet()
  .variables <- Variables()
  obsNumberExists <- is.element("obsNumber", .variables)
  activate1 <- c( checkMethod("fitted", .activeModel, default=TRUE, reportError=FALSE),
                  checkMethod("residuals", .activeModel, default=TRUE, reportError=FALSE))
  activate2 <- c( checkMethod("hatvalues", .activeModel, reportError=FALSE),
                  TRUE, # dffits
                  checkMethod("cooks.distance", .activeModel, reportError=FALSE))
  activate3 <- c( TRUE,TRUE) # stdres,  studres
  activate4 <- c( TRUE,TRUE) # pearson, spearson
  activate5 <- c( TRUE, TRUE) # PRESS.res, PRESS.pred
  checkBoxes(frame="selectFrame1", boxes=c(c("fitted", "residuals")[activate1]),
             initialValues=c(rep(0, sum(activate1))),
             labels=c(gettextRcmdr(c("Fitted values", "Residuals"))[activate1]))
  checkBoxes(frame="selectFrame2", boxes=c(c("hatvalues", "dffits", "cookd")[activate2]),
             initialValues=c(rep(0, sum(activate2))),
             labels=c(gettextRcmdr(c("Hat-values", "DFFITS", "Cook's distances"))[activate2]))
  checkBoxes(frame="selectFrame3", boxes=c(c("stdres", "studres")[activate3]),
             initialValues=c(rep(0, sum(activate3)), "0"),
             labels=c(gettextRcmdr(c("Standardized residuals", "Studentized residuals"))[activate3]))
  checkBoxes(frame="selectFrame4", boxes=c(c("pearson", "spearson")[activate4]),
             initialValues=c(rep(0, sum(activate4))),
             labels=c(gettextRcmdr(c("Pearson residuals", "Standardized Pearson residuals"))[activate4]))
  checkBoxes(frame="selectFrame5", boxes=c(c("pressres","presspred")[activate5]),
             initialValues=c(rep(0, sum(activate5))),
             labels=c(gettextRcmdr(c("CV residual", "CV predictions"))[activate5]))
  checkBoxes(frame="selectFrame6", boxes=c("obsNumbers"),
             initialValues=c("0"),
             labels=c( gettextRcmdr("Observation indices")))
  onOK <- function(){
    closeDialog()
    if (activate1[1] && tclvalue(fittedVariable) == 1) addVariable("fitted")
    if (activate1[2] && tclvalue(residualsVariable) == 1) addVariable("residuals")
    if (activate3[1] && tclvalue(stdresVariable) == 1) addVariable("stdres")
    if (activate3[2] && tclvalue(studresVariable) == 1) addVariable("studres")
    if (activate4[1] && tclvalue(pearsonVariable) == 1) addVariable("residuals","pearson")
    if (activate4[2] && tclvalue(spearsonVariable) == 1) addVariable("spearson")
    if (activate2[1] && tclvalue(hatvaluesVariable) == 1) addVariable("hatvalues")
    if (activate2[2] && tclvalue(dffitsVariable) == 1) addVariable("dffits")
    if (activate2[3] && tclvalue(cookdVariable) == 1) addVariable("cooks.distance")
    if (activate5[1] && tclvalue(pressresVariable) == 1) addVariable("PRESS.res")
    if (activate5[2] && tclvalue(presspredVariable) == 1) addVariable("PRESS.pred")
    if (tclvalue(obsNumbersVariable) == 1){
      proceed <- if (obsNumberExists) tclvalue(checkReplace("obsNumber")) else "yes"
      if (proceed == "yes") {
        command <- paste(.activeDataSet, "$obsNumber <- 1:nrow(", .activeDataSet, ")", sep="")
        justDoIt(command)
        logger(command)
      }
    }
    activeDataSet(.activeDataSet, flushModel=FALSE)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="influence.measures")
  
  tkgrid(selectFrame1, sticky="w", column=1, row=1, rowspan=2)
  tkgrid(selectFrame2, sticky="w", column=1, row=4, rowspan=3)
  tkgrid(selectFrame3, sticky="w", column=2, row=1, rowspan=2)
  tkgrid(selectFrame4, sticky="w", column=2, row=4, rowspan=2)
  tkgrid(selectFrame5, sticky="w", column=2, row=7, rowspan=2)
  tkgrid(selectFrame6, sticky="w", column=1, row=8, rowspan=1)

  spaceFrame1 <- tkframe(top)
  tkgrid(labelRcmdr(spaceFrame1, text=gettextRcmdr(" ")), sticky="w")
  tkgrid(spaceFrame1, sticky="w", row=3, column=1, columnspan=1)
  spaceFrame2 <- tkframe(top)
  tkgrid(labelRcmdr(spaceFrame2, text=gettextRcmdr(" ")), sticky="w")
  tkgrid(spaceFrame2, sticky="w", row=3, column=2, columnspan=1)
  spaceFrame3 <- tkframe(top)
  tkgrid(labelRcmdr(spaceFrame3, text=gettextRcmdr(" ")), sticky="w")
  tkgrid(spaceFrame3, sticky="w", row=6, column=2, columnspan=1)
  spaceFrame4 <- tkframe(top)
  tkgrid(labelRcmdr(spaceFrame4, text=gettextRcmdr(" ")), sticky="w")
  tkgrid(spaceFrame4, sticky="w", row=7, column=1, columnspan=1)
  
  tkgrid(buttonsFrame, sticky="w", row=9, column=1, columnspan=2)
  dialogSuffix(rows=9, columns=1)
}


#####################################
# Goodness of fit tests
deviance_tests_GUI <- function(){
  doItAndPrint("deviance_tests()")
}


#####################################
# (X'X)^-1 = solve(crossprod(X))
invXtXGUI <- function(){
  initializeDialog(title=gettextRcmdr("Store and/or display (X'X)^-1 (inverse regressor covariance)"))
  mainFrame <- tkframe(top)
  onOK <- function(){ # Actions to perform
    command <- paste("solve(crossprod(model.matrix(",ActiveModel(),")))",sep="")
    stored <- tclvalue(storeVariable)
    if(stored == gettextRcmdr("1")){
      doItAndPrint(paste("invXtX <- ", command, sep=""))
    }
    displayed <- tclvalue(displayVariable)
    if(displayed == gettextRcmdr("1")){
      doItAndPrint(command)
    }
    closeDialog()
    tkfocus(CommanderWindow())
  }
  # Set up GUI
  OKCancelHelp(helpSubject="lm")
  checkBoxes(frame="mainFrame", boxes=c("store","display"), initialValues=c("1","1"),
             labels=gettextRcmdr(c("Store matrix as invXtX","Display matrix")))
  tkgrid(mainFrame, row=1, column=1, columnspan=1, rowspan=1, sticky="w")
  tkgrid(buttonsFrame, sticky="w", row=2, column=1, columnspan=1)
  dialogSuffix(rows=2, columns=1)
}

#####################################
# Display regression coefficients
coefNMBU <- function(){
  .activeModel <- ActiveModel()
  model.class <- justDoIt(paste("class(", .activeModel, ")", sep=""))[1]
  if(model.class=="mvr"){
	ncomp <- justDoIt(paste(.activeModel, "$ncomp", sep=""))
	doItAndPrint(paste("coef(", ActiveModel(), ", ncomp=", ncomp, ")", sep=""))
  } else {
    doItAndPrint(paste("coef(", ActiveModel(), ")", sep=""))
  }
}

#####################################
# CI for the model grand
CIgrandMeanNMBU <- function(){
	doItAndPrint(paste("CIgrandMean(", ActiveModel(), ")", sep=""))
}