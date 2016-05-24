#' A graphical user interface for the package rankFD
#' 
#' This function provides a graphical user interface for calculating rank-based 
#' statistical tests in general factorial designs.
#' 
#' The function produces a GUI for the calculation of the test statistics and 
#' for plotting. Data can be loaded via the "load data" button. The formula 
#' and the significance level alpha (default: 0.05) need to be specified.
#' One can choose between two different null hypotheses to be tested as well as
#' weighted or unweighted effects.
#' If the plot option is chosen, an additional window opens containing
#' information on the plots.
#' 
#' 
#' @export

# GUI
calculateGUI <- function() {
  ## Run on "Load"
  requireNamespace("RGtk2", quietly=TRUE)
  if(!("package:RGtk2" %in% search())){attachNamespace("RGtk2")}
  getDirectory <- function(button, user.data){
    directory = file.choose()
    RGtk2::gtkEntrySetText(filename,directory)
  }  
  ## Run on "OK"
  performStatistics <- function(button, user.data) {
    res <- NULL
    d <- NULL
    error <- NULL
    warning <- NULL
    # Get the information about data and the file
    the.file <- filename$getText()
    the.formula <- formula(filename1$getText())
    the.alpha <- as.numeric(filename3$getText())
    the.plot <- toPlot$active
    the.effect <-c("weighted","unweighted")[comboboxeffect$active+1]
    the.hypothesis <-c("H0F","H0p")[comboboxhypothesis$active+1]
    the.sep <- sepEntry$getText()
    the.headers <- headersEntry$active
    the.dec <- decEntry$getText()
    the.factors = toFactors$active
    
    
    d <- read.table(the.file, sep = the.sep, header = the.headers,
                    dec = the.dec)
    
    dat.ModelTEST <- model.frame(the.formula,d)
    #--------Numbers of factors---------#
    
    for(i in length(dat.ModelTEST):2) dat.ModelTEST[,i] <-as.factor(dat.ModelTEST[,i])
    nfTEST = ncol(dat.ModelTEST)-1
    
    n.levelsTEST = c()
    for (i in 2:(nfTEST+1)){
      n.levelsTEST[i-1]= nlevels(dat.ModelTEST[,i])
    }
    if(!(nfTEST==1 && n.levelsTEST ==2)){
      res <- rangtests(the.formula, d,plot_CI = the.plot,alpha=the.alpha,effect=the.effect,
                       hypothesis=the.hypothesis, Factor.Information=the.factors)
      print(res)
    }
    
    if(nfTEST == 1 && n.levelsTEST == 2){
      ###########################################################
      
      #-------------GUI Function!!--------#
      calculateGUItwosamples <- function() {
        twosamples <- function(button, user.data) {
          
          error <- NULL
                    
          if (!is.null(error)) {
            hbox <- RGtk2::gtkHBoxNew()
            vbox$packStart(hbox,FALSE,FALSE,0)
            label <- RGtk2::gtkLabel(error)
            hbox$packStart(label,FALSE,FALSE,0)
          }
          
          
          the.alternative <-c("two.sided","less","greater")[comboboxalternative$active+1]
          the.nperm <-as.numeric(filename14$getText())
          the.alpha <- as.numeric(filename16$getText())
          the.shift <- toPlot$active
          the.output <-toOutput$active
          the.wilcoxon <-c("asymptotic","exact")[comboboxwilcoxon$active+1]
          the.method = c("normal","t.app","logit","probit","permu")[combobox$active+1]
          res<-rank.two.samples(the.formula, d, conf.level = 1-the.alpha, plot.simci=the.plot,
                                alternative = the.alternative, method =the.method,  info = the.output, 
                                shift.int=the.shift,nperm = the.nperm,wilcoxon=the.wilcoxon) 
          
          print(res)
          ## Run on "OK" 
        }
        
        ##################################################
        # Create window
        window <- RGtk2::gtkWindow()
        # Add title
        window["title"] <- "Two-sample Rank Tests"
        
        # Add a frame
        frame <- RGtk2::gtkFrameNew("Please specify which method should be computed")
        window$add(frame)
        
        # Create vertical container for file name entry
        vbox <- RGtk2::gtkVBoxNew(FALSE, 8)
        vbox$setBorderWidth(24)
        frame$add(vbox)
        # Add horizontal container for every widget line
        hbox <- RGtk2::gtkHBoxNew(FALSE, 8)
        vbox$packStart(hbox, FALSE, FALSE, 0)
        
        hbox <- RGtk2::gtkHBoxNew(FALSE,8)
        vbox$packStart(hbox, FALSE, FALSE, 0)
        labeltest <- RGtk2::gtkLabelNewWithMnemonic("__#----------Tests and Confidence Intervals for Relative Effects-------#")
        hbox$packStart(labeltest,FALSE,FALSE,0)    
        
        ############################################################
        
        # Add an horizontal container to specify parameters
        
        hbox <- RGtk2::gtkHBoxNew(FALSE,8)
        vbox$packStart(hbox, FALSE, FALSE, 0)
        labeltest <- RGtk2::gtkLabelNewWithMnemonic("Method")
        hbox$packStart(labeltest,FALSE,FALSE,0)
        model<-RGtk2::rGtkDataFrame(c("normal","t.app","logit","probit","permu"))
        combobox <- RGtk2::gtkComboBox(model)
        #combobox allowing to decide whether we want result as integer or double
        crt <- RGtk2::gtkCellRendererText()
        combobox$packStart(crt)
        combobox$addAttribute(crt, "text", 0)
        
        RGtk2::gtkComboBoxSetActive(combobox,0)
        hbox$packStart(combobox)
        
        labelAlternative <- RGtk2::gtkLabelNewWithMnemonic("Alternative")
        hbox$packStart(labelAlternative,FALSE,FALSE,0)
        alternative<-RGtk2::rGtkDataFrame(c("two.sided","less","greater"))
        comboboxalternative <- RGtk2::gtkComboBox(alternative)
        #combobox allowing to decide whether we want result as integer or double
        crtalternative <- RGtk2::gtkCellRendererText()
        comboboxalternative$packStart(crtalternative)
        comboboxalternative$addAttribute(crtalternative, "text", 0)
        RGtk2::gtkComboBoxSetActive(comboboxalternative,0)
        hbox$packStart(comboboxalternative)
        
        label14 <- RGtk2::gtkLabelNewWithMnemonic("nperm")
        hbox$packStart(label14,FALSE,FALSE,0)
        # Add entry in the second column; named "filename14"
        filename14 <- RGtk2::gtkEntryNew()
        filename14$setWidthChars(10)
        filename14$setText(10000)
        label14$setMnemonicWidget(filename14)
        hbox$packStart(filename14,FALSE,FALSE,0)
        
        hbox <- RGtk2::gtkHBoxNew(FALSE,8)
        vbox$packStart(hbox, FALSE, FALSE, 0)
        
        # Add separator
        vbox$packStart(RGtk2::gtkHSeparatorNew(), FALSE, FALSE, 0)
                
        labelWilcoxon <- RGtk2::gtkLabelNewWithMnemonic("Wilcoxon Test")
        hbox$packStart(labelWilcoxon,FALSE,FALSE,0)
        wilcoxon<-RGtk2::rGtkDataFrame(c("asymptotic","exact"))
        comboboxwilcoxon <- RGtk2::gtkComboBox(wilcoxon)
        crtwilcoxon <- RGtk2::gtkCellRendererText()
        comboboxwilcoxon$packStart(crtwilcoxon)
        comboboxwilcoxon$addAttribute(crtwilcoxon, "text", 0)
        
        RGtk2::gtkComboBoxSetActive(comboboxwilcoxon,0)
        hbox$packStart(comboboxwilcoxon)        
        
        label16 <- RGtk2::gtkLabelNewWithMnemonic("Alpha")
        hbox$packStart(label16,FALSE,FALSE,0)
        # Add entry in the second column; named "filename16"
        filename16 <- RGtk2::gtkEntryNew()
        filename16$setWidthChars(10)
        filename16$setText(0.05)
        label16$setMnemonicWidget(filename16)
        hbox$packStart(filename16,FALSE,FALSE,0)
        
        # Add Shift Effect Option
        hbox <- RGtk2::gtkHBoxNew(FALSE,8)
        vbox$packStart(hbox, FALSE, FALSE, 0)
        label <- RGtk2::gtkLabelNewWithMnemonic("Shift Effects?")
        hbox$packStart(label,FALSE,FALSE,0)
        toPlot <- RGtk2::gtkCheckButton()
        hbox$packStart(toPlot,FALSE,FALSE,0)
        
        # Add Output Information Option
        hbox <- RGtk2::gtkHBoxNew(FALSE,8)
        vbox$packStart(hbox, FALSE, FALSE, 0)
        label <- RGtk2::gtkLabelNewWithMnemonic("Output Information?")
        hbox$packStart(label,FALSE,FALSE,0)
        toOutput <- RGtk2::gtkCheckButton()
        hbox$packStart(toOutput,FALSE,FALSE,0)
        
        ############################################################
        
        # Add button
        the.buttons <- RGtk2::gtkHButtonBoxNew()
        the.buttons$setBorderWidth(5)
        vbox$add(the.buttons)
        the.buttons$setLayout("spread")
        the.buttons$setSpacing(40)
        buttonOK <- RGtk2::gtkButtonNewFromStock("gtk-ok")
        RGtk2::gSignalConnect(buttonOK, "clicked", twosamples)
        the.buttons$packStart(buttonOK,fill=F)
        
        ##################################################
        
        buttonCancel <- RGtk2::gtkButtonNewFromStock("gtk-close")
        RGtk2::gSignalConnect(buttonCancel, "clicked", window$destroy)
        the.buttons$packStart(buttonCancel,fill=F)
        
      }
      ###################################
      
      #end of two samples!!
      calculateGUItwosamples()
            
    }
  }
  
  # Create window
  window <- RGtk2::gtkWindow()
  # Add title
  window["title"] <- "Rank Methods for Factorial Designs"
  
  # Add a frame
  frame <- RGtk2::gtkFrameNew("Specify data location and formula...")
  window$add(frame)
  
  # Create vertical container for file name entry
  vbox <- RGtk2::gtkVBoxNew(FALSE, 8)
  vbox$setBorderWidth(24)
  frame$add(vbox)
  # Add horizontal container for every widget line
  hbox <- RGtk2::gtkHBoxNew(FALSE, 8)
  vbox$packStart(hbox, FALSE, FALSE, 0)
  
  label <- RGtk2::gtkLabelNewWithMnemonic("_File name")
  hbox$packStart(label,FALSE,FALSE,0)
  # Add entry in the second column; named "filename"
  filename <- RGtk2::gtkEntryNew()
  filename$setWidthChars(50)
  label$setMnemonicWidget(filename)
  hbox$packStart(filename,FALSE,FALSE,0)
  
  ############################################################
  # Add label in first column
  label1 <- RGtk2::gtkLabelNewWithMnemonic("_Formula")
  hbox$packStart(label1,FALSE,FALSE,0)
  # Add entry in the second column; named "filename1"
  filename1 <- RGtk2::gtkEntryNew()
  filename1$setWidthChars(50)
  label1$setMnemonicWidget(filename1)
  hbox$packStart(filename1,FALSE,FALSE,0)
  
  ############################################################
  
  label3 <- RGtk2::gtkLabelNewWithMnemonic("_alpha")
  hbox$packStart(label3,FALSE,FALSE,0)
  # Add entry in the second column; named "filename3"
  filename3 <- RGtk2::gtkEntryNew()
  filename3$setWidthChars(10)
  filename3$setText(0.05)
  label3$setMnemonicWidget(filename3)
  hbox$packStart(filename3,FALSE,FALSE,0)
    
  #####################################################
  # Add separator
  hbox <- RGtk2::gtkHBoxNew(FALSE,8)
  vbox$packStart(hbox, FALSE, FALSE, 0)
  vbox$packStart(RGtk2::gtkHSeparatorNew(), FALSE, FALSE, 0)
  
  labelhypothesis <- RGtk2::gtkLabelNewWithMnemonic("Hypothesis")
  hbox$packStart(labelhypothesis,FALSE,FALSE,0)
  hypothesis<-RGtk2::rGtkDataFrame(c("Distribution Functions H_0^F","Relative Effects H_0^p"))
  comboboxhypothesis <- RGtk2::gtkComboBox(hypothesis)
  crthypothesis <- RGtk2::gtkCellRendererText()
  comboboxhypothesis$packStart(crthypothesis)
  comboboxhypothesis$addAttribute(crthypothesis, "text", 0)
  RGtk2::gtkComboBoxSetActive(comboboxhypothesis,0)
  hbox$packStart(comboboxhypothesis)
  
  labeleffect <- RGtk2::gtkLabelNewWithMnemonic("Effects")
  hbox$packStart(labeleffect,FALSE,FALSE,0)
  effect<-RGtk2::rGtkDataFrame(c("weighted (n_i/N)","unweighted (1/d)"))
  comboboxeffect <- RGtk2::gtkComboBox(effect)
  crteffect <- RGtk2::gtkCellRendererText()
  comboboxeffect$packStart(crteffect)
  comboboxeffect$addAttribute(crteffect, "text", 0)
  
  RGtk2::gtkComboBoxSetActive(comboboxeffect,0)
  hbox$packStart(comboboxeffect)
  
  #####################################################
  
  ############################################################
  
  # Add an horizontal container to specify input file options
  # are headers included in the file?
  hbox <- RGtk2::gtkHBoxNew(FALSE,8)
  vbox$packStart(hbox, FALSE, FALSE, 0)
  label <- RGtk2::gtkLabelNewWithMnemonic("_Headers?")
  hbox$packStart(label,FALSE,FALSE,0)
  headersEntry <- RGtk2::gtkCheckButton()
  headersEntry$active <- TRUE
  hbox$packStart(headersEntry,FALSE,FALSE,0)
  label$setMnemonicWidget(headersEntry)
  
  # are headers included in the file?
  label <- RGtk2::gtkLabelNewWithMnemonic("Col. _Separator?")
  hbox$packStart(label,FALSE,FALSE,0)
  sepEntry <- RGtk2::gtkEntryNew()
  sepEntry$setWidthChars(1)
  sepEntry$setText("")
  hbox$packStart(sepEntry,FALSE,FALSE,0)
  label$setMnemonicWidget(sepEntry)
  
  # what's the character used for decimal points?
  label <- RGtk2::gtkLabelNewWithMnemonic("_Dec. character?")
  hbox$packStart(label,FALSE,FALSE,0)
  decEntry <- RGtk2::gtkEntryNew()
  decEntry$setWidthChars(1)
  decEntry$setText(".")
  hbox$packStart(decEntry,FALSE,FALSE,0)
  label$setMnemonicWidget(decEntry)
  
  # Add separator
  vbox$packStart(RGtk2::gtkHSeparatorNew(), FALSE, FALSE, 0)
  
  # Add plot-option
  hbox <- RGtk2::gtkHBoxNew(FALSE,8)
  vbox$packStart(hbox, FALSE, FALSE, 0)
  label <- RGtk2::gtkLabelNewWithMnemonic("Plot _Results?")
  hbox$packStart(label,FALSE,FALSE,0)
  toPlot <- RGtk2::gtkCheckButton()
  hbox$packStart(toPlot,FALSE,FALSE,0)
  
  # Add Factor information
  #hbox <- RGtk2::gtkHBoxNew(FALSE,8)
  #vbox$packStart(hbox, FALSE, FALSE, 0)
  label <- RGtk2::gtkLabelNewWithMnemonic("Individual Factor Results?")
  hbox$packStart(label,FALSE,FALSE,0)
  toFactors <- RGtk2::gtkCheckButton()
  hbox$packStart(toFactors,FALSE,FALSE,0)
  
  # Add button
  the.buttons <- RGtk2::gtkHButtonBoxNew()
  the.buttons$setBorderWidth(5)
  vbox$add(the.buttons)
  the.buttons$setLayout("spread")
  the.buttons$setSpacing(40)
  buttonOK <- RGtk2::gtkButtonNewFromStock("gtk-ok")
  buttonLoad <- RGtk2::gtkButtonNewFromStock("Load Data")
  RGtk2::gSignalConnect(buttonOK, "clicked", performStatistics)
  RGtk2::gSignalConnect(buttonLoad, "clicked", getDirectory)
  the.buttons$packStart(buttonOK,fill=F)
  the.buttons$packStart(buttonLoad,fill=F)
  
  ##################################################
  
  buttonCancel <- RGtk2::gtkButtonNewFromStock("gtk-close")
  RGtk2::gSignalConnect(buttonCancel, "clicked", window$destroy)
  the.buttons$packStart(buttonCancel,fill=F)
}



rangtests <- function(formula, data,alpha=0.05, CI.method=c("Logit","Normal"),
                      plot_CI=FALSE,effect,hypothesis=c("H0F","H0p"), 
                      Factor.Information=FALSE){
  
  requireNamespace("RGtk2", quietly=TRUE)
  
  #-------- Bestimme das gegebene Modell------#
  
  dat.Model0 <- model.frame(formula, data)
  
  #--------Numbers of factors---------#
  
  for(i in length(dat.Model0):2) dat.Model0[,i] <-as.factor(dat.Model0[,i])
  nf = ncol(dat.Model0)-1
  
  n.levels = c()
  names.levels=list()
  
  for (i in 2:(nf+1)){n.levels[i-1]= nlevels(dat.Model0[,i])
                      names.levels[[i-1]] = levels(dat.Model0[,i])
  }
  names(names.levels) <- names(dat.Model0[,2:(nf+1)])
  
  #########################################################################
  #########################################################################
  #########################################################################
  
  #-------------Hypothesenmatrizen--------#
  
  perm_names <- t(attr(terms(formula), "factors")[-1, ])
  nr_hypo <- attr(terms(formula), "factors")
  fac_names <- colnames(nr_hypo)
  Hypotheses0 = HC(n.levels,"Hyp",perm_names,fac_names)
  Hypotheses=Hypotheses0[[1]]
  CI.Matrices = HC(n.levels,"CI",perm_names,fac_names)[[1]]
  
  n.hypotheses = length(Hypotheses)
  n.levels.prod=prod(n.levels)
  #Output.names <- attr(terms(formula), "term.labels")
  Output.names <- Hypotheses0[[2]]
  
  #-------Sortiere die Daten nach den Faktoren------#
  for (i in length(dat.Model0):2) {
    dat.Model0 <- dat.Model0[order(dat.Model0[,i]),]}
  
  #-----------------Fuege Pseudo Faktor dem Datensatz zu-----#
  dat.Model0$Interaction = interaction(dat.Model0[,2:length(dat.Model0)],sep=":")
  dat.response <- dat.Model0[,1]
  
  #----------------Compute means and variances etc---------#
  n <- aggregate(formula,data=dat.Model0,length)
  for(i in (length(n)-1):1) {
    n <-n[order(n[,i]),]
  }
  colnames(n)[ncol(n)]<-"Size"
  
  #-------------Compute Inference Methods------#
  dat.Model0$INum <- as.factor(rep(1:n.levels.prod, n$Size))
  
  #--------------Compute the relative Effects----#
  H0pW <-Effects(dat.response, dat.Model0$INum,effect)
  n$pd <- c(H0pW$pd)
  n$Var <- c(diag(H0pW$VBF))
  
  CI <- Limits(c(H0pW$pd),H0pW$VBF,alpha,c(H0pW$N))
  n$L.Normal <- CI$Normal[,1]
  n$U.Normal <- CI$Normal[,2]
  n$L.Logit <- CI$Logit[,1]
  n$U.Logit <- CI$Logit[,2]
  
  WTS = matrix(0,n.hypotheses,3)
  ATS = matrix(0,n.hypotheses, 4)
  ATSp = matrix(0,n.hypotheses,4)
  Descriptive.Factors = list()
  Levels.Factors = list()
  
  
  for(i in 1:n.hypotheses){
        
    WTS[i,] = Wald(c(H0pW$pd),Hypotheses[[i]],H0pW$VH0F)
    ATS[i,] =ANOVATYP(c(H0pW$pd),Hypotheses[[i]],H0pW$VH0F,n$Size)
    ATSp[i,]= ANOVATYPH0P(c(H0pW$pd),Hypotheses[[i]],H0pW$VBF,n$Size,H0pW$dfATS)
    CILimits <-Limits(c(CI.Matrices[[i]]%*%c(H0pW$pd)),CI.Matrices[[i]]%*%H0pW$VBF%*%t(CI.Matrices[[i]]),alpha,H0pW$N)
    
    Descriptives <-data.frame(pd=CI.Matrices[[i]]%*%n$pd,
                              Var= c(diag(CI.Matrices[[i]]%*%H0pW$VBF%*%t(CI.Matrices[[i]]))),
                              L.Normal=CILimits$Normal[,1],
                              U.Normal = CILimits$Normal[,2],
                              L.Logit = CILimits$Logit[,1],
                              U.Logit = CILimits$Logit[,2])
    Output.namesi <-Output.names[i] 
    formula.act <- as.formula(paste(names(dat.Model0)[1], Output.namesi, sep=" ~ "))
    aha <- data.frame(aggregate(formula.act,data=dat.Model0,mean))
    
    for(ii in (length(aha)-1):1) {aha <-aha[order(aha[,ii]),]}
    Descriptive.Factors[[i]] <-data.frame(aha,Descriptives)
    Descriptive.Factors[[i]] <-Descriptive.Factors[[i]][,-length(aha)]
    if(length(grep(":", Output.names[i]))<1){
      pos <- which(names(dat.Model0)==Output.names[i])
      Levels.Factors[[i]] = data.frame(X=levels(dat.Model0[,pos]))
    }
    
    if(length(grep(":", Output.names[i]))>=1){
      facs.singles <- c(strsplit(Output.names[i], ":")[[1]])
      
      Levels.Factors[[i]] = data.frame(n[,facs.singles])
      
      
    }
  }
  
  names(Descriptive.Factors) <- Output.names
  
  rownames(WTS) <- Output.names
  rownames(ATS) <- Output.names
  rownames(ATSp)<- Output.names
  colnames(WTS) <- c("Statistic", "df", "p-Value")
  colnames(ATS) <- c("Statistic", "df1", "df2", "p-Value")
  colnames(ATSp) <- c("Statistic", "df1", "df2", "p-Value")

  
  if(plot_CI==TRUE){
    
    calculateGUIplot <- function() {
      plotting <- function(button, user.data) {
        
        ######################################################################
        # PLOT GUI for One way!!!
        if(nf ==1){
          Faktor = fac_names[1]}
        
        if(nf > 1){
          Faktor <- filename$getText()}
        
        Fak.split <- strsplit(Faktor,":")[[1]]
        l.Fak.split <- length(Fak.split)
        error <- NULL
        
        Title <- filename2$getText()
        CI.method = c("Logit","Normal")[comboboxcimethod$active+1]
        line_width <- as.numeric(filename3$getText())
        
        if (!(Faktor %in% fac_names)) {
          error <- "Please enter a valid factor name"
        }
        
        if (!is.null(error)) {
          hbox <- RGtk2::gtkHBoxNew()
          vbox$packStart(hbox,FALSE,FALSE,0)
          label <- RGtk2::gtkLabel(error)
          hbox$packStart(label,FALSE,FALSE,0)
        }
        
        for(i in 1:n.hypotheses){
          
          if(names(Descriptive.Factors)[i] == Faktor){
            posP <- which(names(Descriptive.Factors)[i] == Faktor)
            DatenPlot <- data.frame(Descriptive.Factors[[i]])
            
            
            switch(CI.method, Logit={
              upper = DatenPlot$U.Logit 
              lower = DatenPlot$L.Logit},
              Normal={
                upper =DatenPlot$U.Normal
                lower = DatenPlot$L.Normal})
            
            
            if (l.Fak.split==1){
              print(xyplot(pd ~ DatenPlot[,1], group=DatenPlot[,1],data = DatenPlot, 
                           type = 'p',ylim=c(0,1),col=1,pch=7,cex=1.3,xlab=paste(names(DatenPlot[1])),
                           ylab="",upper = upper,main=Title,lwd=line_width,
                           lower = lower,
                           panel = function(x, y, ...){
                             panel.superpose(x, y,
panel.groups = function(x, y, upper, lower,subscripts, ..., font, fontface) {
upper <- upper[subscripts]
lower <- lower[subscripts]
panel.arrows(x, lower, x, upper,code=4,lwd=4)   
panel.points(x, lower,pch="_",cex=3,lwd=4,col=1)
panel.points(x, upper,pch="_",cex=3,lwd=4,col=1)
} , ...)
                             panel.xyplot(x, y, ...)}))}
            
            if (l.Fak.split==2){
              print(xyplot(pd ~ DatenPlot[,1]|DatenPlot[,2], 
                           group=DatenPlot[,1],data = DatenPlot, type = 'p',ylim=c(0,1),
                           pch=7,cex=1.3,lwd=line_width,xlab=paste(names(DatenPlot[1])),col=1,main=Title,
                           upper = upper,
                           lower = lower,
                           panel = function(x, y, ...){
                             panel.superpose(x, y,
panel.groups = function(x, y, upper, lower,subscripts, ..., font, fontface) {
upper <- upper[subscripts]
lower <- lower[subscripts]
panel.arrows(x, lower, x, upper,code=4,lwd=4)   
panel.points(x, lower,pch="_",cex=3,lwd=4,col=1)
panel.points(x, upper,pch="_",cex=3,lwd=4,col=1)
} , ...)
                             panel.xyplot(x, y, ...)
                           }))
              
            }
            
            
            if (l.Fak.split==3){
              print(xyplot(pd ~ DatenPlot[,1]|DatenPlot[,2]*DatenPlot[,3], 
                           group=DatenPlot[,1],data = DatenPlot, type = 'p',ylim=c(0,1),col=1,
                           pch=7,cex=1.3,lwd=line_width,xlab=paste(names(DatenPlot[1])),
                           main=Title,
                           upper = upper,
                           lower = lower,
                           panel = function(x, y, ...){
                             panel.superpose(x, y,
panel.groups = function(x, y, upper, lower,subscripts, ..., font, fontface) {
upper <- upper[subscripts]
lower <- lower[subscripts]
panel.arrows(x, lower, x, upper,code=4,lwd=4)   
panel.points(x, lower,pch="_",cex=3,lwd=4,col=1)
panel.points(x, upper,pch="_",cex=3,lwd=4,col=1)
} , ...)
     panel.xyplot(x, y, ...)
                           }))
              
            }
            
            if (l.Fak.split>=4){
              stop("4 and higher way interactions cannot be plotted!")
            }
            
          }
          
        }
        
      }
      
      # Create window
      window <- RGtk2::gtkWindow()
      # Add title
      window["title"] <- "Plot"
      
      # Add a frame
      frame <- RGtk2::gtkFrameNew("Please choose the factor you wish to plot (for interaction type something like group1:group2).")
      window$add(frame)
      
      # Create vertical container for file name entry
      vbox <- RGtk2::gtkVBoxNew(FALSE, 8)
      vbox$setBorderWidth(24)
      frame$add(vbox)
      # Add horizontal container for every widget line
      hbox <- RGtk2::gtkHBoxNew(FALSE, 8)
      vbox$packStart(hbox, FALSE, FALSE, 0)
      
      # Add label in first column
      if(nf >1){
        label <- RGtk2::gtkLabelNewWithMnemonic("_Factor")
        hbox$packStart(label,FALSE,FALSE,0)
        # Add entry in the second column; named "filename"
        filename <- RGtk2::gtkEntryNew()
        filename$setWidthChars(50)
        label$setMnemonicWidget(filename)
        hbox$packStart(filename,FALSE,FALSE,0)
      } 
      ############################################################
      
      # Add an horizontal container to specify parameters
      hbox <- RGtk2::gtkHBoxNew(FALSE,8)
      vbox$packStart(hbox, FALSE, FALSE, 0)
      
      label2 <- RGtk2::gtkLabelNewWithMnemonic("_Title")
      hbox$packStart(label2,FALSE,FALSE,0)
      # Add entry in the second column; named "filename2"
      filename2 <- RGtk2::gtkEntryNew()
      filename2$setWidthChars(10)
      label2$setMnemonicWidget(filename2)
      hbox$packStart(filename2,FALSE,FALSE,0)
      
      label3 <- RGtk2::gtkLabelNewWithMnemonic("_lwd")
      hbox$packStart(label3,FALSE,FALSE,0)
      # Add entry in the second column; named "filename3"
      filename3 <- RGtk2::gtkEntryNew()
      filename3$setWidthChars(10)
      filename3$setText(2)
      label3$setMnemonicWidget(filename3)
      hbox$packStart(filename3,FALSE,FALSE,0)
      
      labelcimethod <- RGtk2::gtkLabelNewWithMnemonic("CI Method")
      hbox$packStart(labelcimethod,FALSE,FALSE,0)
      cimethod<-RGtk2::rGtkDataFrame(c("Logit","Normal"))
      comboboxcimethod <- RGtk2::gtkComboBox(cimethod)
      #combobox allowing to decide whether we want result as integer or double
      crtcimethod <- RGtk2::gtkCellRendererText()
      comboboxcimethod$packStart(crtcimethod)
      comboboxcimethod$addAttribute(crtcimethod, "text", 0)
      RGtk2::gtkComboBoxSetActive(comboboxcimethod,0)
      hbox$packStart(comboboxcimethod)
      
      ############################################################
      
      # Add button
      the.buttons <- RGtk2::gtkHButtonBoxNew()
      the.buttons$setBorderWidth(5)
      vbox$add(the.buttons)
      the.buttons$setLayout("spread")
      the.buttons$setSpacing(40)
      buttonOK <- RGtk2::gtkButtonNewFromStock("gtk-ok")
      RGtk2::gSignalConnect(buttonOK, "clicked", plotting)
      the.buttons$packStart(buttonOK,fill=F)
      
      ##################################################
      
      buttonCancel <- RGtk2::gtkButtonNewFromStock("gtk-close")
      RGtk2::gSignalConnect(buttonCancel, "clicked", window$destroy)
      the.buttons$packStart(buttonCancel,fill=F)
    }
    
    
    #end of plot_ci==TRUE
    calculateGUIplot()
    
  }
  
  if(hypothesis=="H0F"){
    result <- list(Descriptive=n, Wald.Type.Statistic = WTS, ANOVA.Type.Statistic=ATS)
  }
  
  if(hypothesis=="H0p"){
    result <- list(Descriptive=n,  ANOVA.Type.Statistic=ATSp)
  } 
  
  
  
  
  
  if (Factor.Information==TRUE){
    print(Descriptive.Factors)
  }
  return(result)
}


