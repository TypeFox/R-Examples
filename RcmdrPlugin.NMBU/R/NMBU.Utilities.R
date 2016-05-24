# $Id: NMBU.Utilities.R 35 2014-01-10 21:17:26Z khliland $

## Utility functions for RcmdrPlugin.NMBU

###################
# Nothing (to make menus work)
.dc <- function(){}

################################
# Model and data availability
lmmP <- function() activeModelP() && any(class(get(ActiveModel()))[1] == "lmm")
aovP <- function() lmP()
# {
#	.activeModel <- ActiveModel()
#	ff <- justDoIt(paste(.activeModel, "$xlevels", sep=""))
#	activeModelP() && (length(ff)==1 || length(ff)==2)
#	}
clustP <- function() length(listHclustSolutions())>0
pcaP <- function() activeModelP() && any(class(get(ActiveModel()))[1] == c('prcomp'))
plsP <- function() activeModelP() && any(class(get(ActiveModel()))[1] == c('mvr'))
daP <- function() activeModelP() && (any(class(get(ActiveModel()))[1] == c('lda')) || any(class(get(ActiveModel()))[1] == c('qda')))
contP <- function(){ # Check if model contains continuous variable
  .activeModel <- ActiveModel()
  effects <- fparse(formula(.activeModel))
  for(i in 1:length(effects)){
    nu <- eval(parse(text=paste("class(", ActiveDataSet(), '$', effects[[i]],")")))
    if(any(nu%in%c("numeric","integer"))){
      return(FALSE)
    }
  }
  TRUE
}
variablesP <- function(n=1) activeDataSetP() && length(listVariables()) >= n
mixP <- function(){
  activeModelP() && any(class(get(ActiveModel()))[1] == c('lm'))
}

###########################################################
# Colour function for numerical and catergorical variables
make.colours <- function(object) {
  if(is.factor(object))
    my.colours <- as.numeric(object)
  else {
    object <- object-min(object)
    object <- 4*object/max(object)
    my.colours <- rgb(0,0,rep(0,length(object)))
    my.colours[object<=1] <- rgb(0,object[object<=1],1)
    my.colours[object<=2 & object>1] <- rgb(0,1,2-object[object<=2 & object>1])
    my.colours[object<=3 & object>2] <- rgb(object[object<=3 & object>2]-2,1,0)
    my.colours[object<=4 & object>3] <- rgb(1,4-object[object<=4 & object>3],0)
  }
  my.colours
}

#################################
# Estimated beta CI for PLSR/PCR
confint.mvr <- function (object, parm, level=0.95, ...){
  print("Normal approximated jackknife estimated confidence intervals!")
  print(var.jack(object)[,1,1] %o% c(qnorm((1-level)/2),qnorm(1-(1-level)/2)) + object$coefficients[,1,object$ncomp])
}
Confint.mvr <- function(object, parm, level = 0.95, ...) confint (object, parm=parm, level=0.95, ...)

##############################
# Dummy coding of categories
dummy <- function(y){
  if(!is.factor(y)){
    y <- factor(y)
  }
  n <- length(y)
  groups  <- levels(y)
  g       <- length(groups)
  Y.dummy <- matrix(0,n,g)
  for(i in 1:g){
    Y.dummy[,i]  <- y==groups[i]
  }
  colnames(Y.dummy) <- groups
  Y.dummy
}


#################################################
# Create dummy variables and dummy data.frames
dummify <- function (y,n,name) {
  groups <- levels(y)
  g <- length(groups)-1
  Y.dummy <- matrix(0, n, g)
  for (i in 1:g) {
    Y.dummy[, i] <- y == groups[i]
  }
  Y.dummy
}
Dummify <- function(data, main.effects, response){
  n <- length(main.effects)
  N <- dim(data)[1]
  data.list <- list(data[,response])
  for(i in 1:n){
    data.list[[i+1]] <- dummify(data[,main.effects[i]],N,main.effects[i])
  }
  names(data.list) <- c(response,main.effects)
  out <- eval(parse(text =paste("model.frame(", response, " ~ ", paste(main.effects,sep="",collapse=" + "), ", data=data.list)", sep="")))
  out
}
#############################
# List of Hclust solutions
listHclustSolutions <- function(envir=.GlobalEnv, ...) {
  objects <- ls(envir=envir, ...)
  if (length(objects) == 0) NULL
  else objects[sapply(objects,
                      function(.x) "hclust" == class(get(.x, envir=envir))[1]) ]
}

############################
# Model formula adapted to Linear model (extended)
modelFormula2 <- defmacro(frame=top, hasLhs=TRUE, expr={
  checkAddOperator <- function(rhs){
    rhs.chars <- rev(strsplit(rhs, "")[[1]])
    if (length(rhs.chars) < 1) return(FALSE)
    check.char <- if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
      rhs.chars[1] else rhs.chars[2]
    !is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%"))
  }
  .variables <- Variables()
  word <- paste("\\[", gettextRcmdr("factor"), "\\]", sep="")
  variables <- paste(.variables,
                     ifelse(is.element(.variables, Factors()), paste("[", gettextRcmdr("factor"), "]", sep=""), ""))
  xBox <- variableListBox(frame, variables, selectmode="multiple", title=gettextRcmdr("Variables (double-click to formula)"))
  onDoubleClick <- if (!hasLhs){
    function(){
      var <- getSelection(xBox)
      tkselection.clear(xBox$listbox, "0", "end")            		
      if (length(grep(word, var)) == 1) var <- sub(word, "",  var)
      tkfocus(rhsEntry)
      rhs <- tclvalue(rhsVariable)
      rhs.chars <- rev(strsplit(rhs, "")[[1]])
      check.char <- if (length(rhs.chars) > 0){
        if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
          rhs.chars[1] else rhs.chars[2]
      }
      else ""
      tclvalue(rhsVariable) <- if (rhs == "" ||
                                   is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%")))
        paste(rhs, var, sep="")
      else paste(rhs, "+", var)
      tkicursor(rhsEntry, "end")
      tkxview.moveto(rhsEntry, "1")
    }
  }
  else{
    function(){
      var <- getSelection(xBox)
      which <- tkcurselection(xBox$listbox)
      tkselection.clear(xBox$listbox, "0", "end")
      if (length(grep(word, var)) == 1) var <- sub(word, "",  var)
      lhs <- tclvalue(lhsVariable)
      if (lhs == "" || tclvalue(tkselection.present(lhsEntry)) == "1"){
        tclvalue(lhsVariable) <- var
        tkselection.clear(lhsEntry)
        tkfocus(rhsEntry)
      }
      else {
        tkfocus(rhsEntry)
        rhs <- tclvalue(rhsVariable)
        rhs.chars <- rev(strsplit(rhs, "")[[1]])
        check.char <- if (length(rhs.chars) > 0){
          if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
            rhs.chars[1] else rhs.chars[2]
        }
        else ""
        tclvalue(rhsVariable) <- if (rhs == "" ||
                                     is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%")))
          paste(rhs, var, sep="")
        else paste(rhs, "+", var)
      }
      tkicursor(rhsEntry, "end")
      tkxview.moveto(rhsEntry, "1")
    }
  }
  tkbind(xBox$listbox, "<Double-ButtonPress-1>", onDoubleClick)
  onPlus <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")										
    if ((check <- !checkAddOperator(rhs)) && length(var) == 0) return()
    if (length(var) > 1){
      if (length(grep(word, var)) > 0) var <- sub(word, "",  var)
      if (length(var) > 1) var <- paste(var, collapse=" + ")
    }
    tclvalue(rhsVariable) <- paste(rhs, if (!check) " + ", var, sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onTimes <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")						
    if ((check <- !checkAddOperator(rhs)) && length(var) == 0) return()
    if (length(var) > 1){
      if (length(grep(word, var)) > 0) var <- sub(word, "",  var)
      var <- trim.blanks(var)
      if (length(var) > 1) var <- paste(var, collapse="*")
      tclvalue(rhsVariable) <- paste(rhs, if (!check) " + ", var, sep="")
    }
    else tclvalue(rhsVariable) <- paste(rhs, if (!check) "*", sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onColon <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")						
    if ((check <- !checkAddOperator(rhs)) && length(var) == 0) return()
    if (length(var) > 1){
      if (length(grep(word, var)) > 0) var <- sub(word, "",  var)
      var <- trim.blanks(var)
      if (length(var) > 1) var <- paste(var, collapse=":")
      tclvalue(rhsVariable) <- paste(rhs, if (!check) " + ", var, sep="")
    }
    else tclvalue(rhsVariable) <- paste(rhs, if (!check) ":", sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onSlash <- function(){
    rhs <- tclvalue(rhsVariable)
    if (!checkAddOperator(rhs)) return()
    tclvalue(rhsVariable) <- paste(rhs, "/",  sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onIn <- function(){
    rhs <- tclvalue(rhsVariable)
    if (!checkAddOperator(rhs)) return()
    tclvalue(rhsVariable) <- paste(rhs, "%in% ")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onMinus <- function(){
    rhs <- tclvalue(rhsVariable)
    if (!checkAddOperator(rhs)) return()
    tclvalue(rhsVariable) <- paste(rhs, "- ")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onPower <- function(){
    rhs <- tclvalue(rhsVariable)
    if (!checkAddOperator(rhs)) return()
    tclvalue(rhsVariable) <- paste(rhs, "^", sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onLeftParen <- function(){
    tkfocus(rhsEntry)
    rhs <- tclvalue(rhsVariable)
    tclvalue(rhsVariable) <- paste(rhs, "(", sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onRightParen <- function(){
    rhs <- tclvalue(rhsVariable)
    if (!checkAddOperator(rhs)) return()
    tclvalue(rhsVariable) <- paste(rhs, ")", sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  outerOperatorsFrame <- tkframe(frame)
  operatorsFrame <- tkframe(outerOperatorsFrame)
  splinePolyFrame <- tkframe(outerOperatorsFrame)
  plusButton <- buttonRcmdr(operatorsFrame, text="+", width="3", command=onPlus)
  timesButton <- buttonRcmdr(operatorsFrame, text="*", width="3", command=onTimes)
  colonButton <- buttonRcmdr(operatorsFrame, text=":", width="3", command=onColon)
  slashButton <- buttonRcmdr(operatorsFrame, text="/", width="3", command=onSlash)
  inButton <- buttonRcmdr(operatorsFrame, text="%in%", width="5", command=onIn)
  minusButton <- buttonRcmdr(operatorsFrame, text="-", width="3", command=onMinus)
  powerButton <- buttonRcmdr(operatorsFrame, text="^", width="3", command=onPower)
  leftParenButton <- buttonRcmdr(operatorsFrame, text="(", width="3", command=onLeftParen)
  rightParenButton <- buttonRcmdr(operatorsFrame, text=")", width="3", command=onRightParen)
  onBSpline <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")
    if (length(var) == 0) var <- " "
    if (grepl("\\[factor\\]", var)){
      Message("spline requires a numeric variable", type="error")
      return()
    }
    if (length(var) > 1){
      Message("cannot select more than one variable", type="error")
      return()
    }
    check <- !checkAddOperator(rhs)
    tclvalue(rhsVariable) <- paste(rhs, 
                                   if (!check) paste(" + bs(", var, ", df=", tclvalue(dfSplineVar), ")", sep="") 
                                   else paste(" bs(", var, ", df=", tclvalue(dfSplineVar), ")", sep=""),
                                   sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onNatSline <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")
    if (length(var) == 0) var <- " "
    if (grepl("\\[factor\\]", var)){
      Message("spline requires a numeric variable", type="error")
      return()
    }
    if (length(var) > 1){
      Message("cannot select more than one variable", type="error")
      return()
    }
    check <- !checkAddOperator(rhs)
    tclvalue(rhsVariable) <- paste(rhs, 
                                   if (!check) paste(" + ns(", var, ", df=", tclvalue(dfSplineVar), ")", sep="") 
                                   else paste(" ns(", var, ", df=", tclvalue(dfSplineVar), ")", sep=""),
                                   sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onPoly <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")
    if (length(var) == 0) var <- " "
    if (grepl("\\[factor\\]", var)){
      Message("polynomial requires a numeric variable", type="error")
      return()
    }
    if (length(var) > 1){
      Message("cannot select more than one variable", type="error")
      return()
    }
    check <- !checkAddOperator(rhs)
    tclvalue(rhsVariable) <- paste(rhs, 
                                   if (!check) paste(" + poly(", var, ", degree=", tclvalue(degPolyVar), ")", sep="") 
                                   else paste(" poly(", var, ", degree=", tclvalue(degPolyVar), ")", sep=""),
                                   sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onRawPoly <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")
    if (length(var) == 0) var <- " "
    if (grepl("\\[factor\\]", var)){
      Message("polynomial requires a numeric variable", type="error")
      return()
    }
    if (length(var) > 1){
      Message("cannot select more than one variable", type="error")
      return()
    }
    check <- !checkAddOperator(rhs)
    tclvalue(rhsVariable) <- paste(rhs, 
                                   if (!check) paste(" + poly(", var, ", degree=", tclvalue(degPolyVar), ", raw=TRUE)", sep="") 
                                   else paste(" poly(", var, ", degree=", tclvalue(degPolyVar), ", raw=TRUE)", sep=""),
                                   sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  bsplineButton <- buttonRcmdr(splinePolyFrame, text=gettextRcmdr("B-spline\n"), width="10", command=onBSpline)
  nsplineButton <- buttonRcmdr(splinePolyFrame, text=gettextRcmdr("natural\nspline"), width="10", command=onNatSline)
  polyButton <- buttonRcmdr(splinePolyFrame, text=gettextRcmdr("orthogonal\npolynomial"), width="10", command=onPoly)
  RawPolyButton <- buttonRcmdr(splinePolyFrame, text=gettextRcmdr("raw\npolynomial"), width="10", command=onRawPoly)
  dfSplineVar <- tclVar("5")
  degPolyVar <- tclVar("2")
  dfDegFrame <- tkframe(outerOperatorsFrame)
  dfSplineSpin <- tkspinbox(dfDegFrame, textvariable=dfSplineVar, state="readonly", from=2, to=10, width=2)
  degPolySpin <- tkspinbox(dfDegFrame, textvariable=degPolyVar, state="readonly", from=2, to=5, width=2)
  tkgrid(plusButton, timesButton, colonButton, slashButton, inButton, minusButton,
         powerButton, leftParenButton, rightParenButton, sticky="w")
  tkgrid(labelRcmdr(dfDegFrame, text=gettextRcmdr("df for splines: ")), dfSplineSpin,  sticky="se")
  tkgrid(labelRcmdr(dfDegFrame, text=gettextRcmdr("deg. for polynomials: ")), degPolySpin, sticky="se")
  formulaFrame <- tkframe(frame)
  formulaFrameMain <- tkframe(formulaFrame)
  onFormulaHelp <- function () print(help("formula"))
  formulaHelpButton <- buttonRcmdr(formulaFrame, text=gettextRcmdr("Model formula\nhelp"), command=onFormulaHelp,
                                   image="::image::helpIcon", compound="left")
  if (hasLhs){
    tkgrid(labelRcmdr(outerOperatorsFrame, text=gettextRcmdr("Model Formula"), 
                      fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w")
    tkgrid(labelRcmdr(outerOperatorsFrame, text="Operators (click to formula):  "), operatorsFrame, sticky="nw")
    tkgrid(bsplineButton, nsplineButton, polyButton, RawPolyButton, sticky="nw")
    tkgrid(labelRcmdr(outerOperatorsFrame, text=gettextRcmdr("Splines/Polynomials:\n(select variable and click)")), 
           splinePolyFrame, dfDegFrame, sticky="nw")
    lhsVariable <- if (currentModel) tclVar(currentFields$lhs) else tclVar("")
    rhsVariable <- if (currentModel) tclVar(currentFields$rhs) else tclVar("")
    rhsEntry <- ttkentry(formulaFrameMain, width="75", textvariable=rhsVariable)
    rhsXscroll <- ttkscrollbar(formulaFrameMain,
                               orient="horizontal", command=function(...) tkxview(rhsEntry, ...))
    tkconfigure(rhsEntry, xscrollcommand=function(...) tkset(rhsXscroll, ...))
    lhsEntry <- ttkentry(formulaFrameMain, width="10", textvariable=lhsVariable)
    lhsScroll <- ttkscrollbar(formulaFrameMain,
                              orient="horizontal", command=function(...) tkxview(lhsEntry, ...))
    tkconfigure(lhsEntry, xscrollcommand=function(...) tkset(lhsScroll, ...))
    tkgrid(lhsEntry, labelRcmdr(formulaFrameMain, text=" ~    "), rhsEntry, sticky="w")
    tkgrid(lhsScroll, labelRcmdr(formulaFrameMain, text=""), rhsXscroll, sticky="w")
    tkgrid.configure(lhsScroll, sticky="ew")
    tkgrid(formulaFrameMain, labelRcmdr(formulaFrame, text="  "), formulaHelpButton, sticky="nw")
  }
  else{
    rhsVariable <- if (currentModel) tclVar(currentFields$rhs) else tclVar("")
    rhsEntry <- ttkentry(formulaFrameMain, width="75", textvariable=rhsVariable)
    rhsXscroll <- ttkscrollbar(formulaFrameMain,
                               orient="horizontal", command=function(...) tkxview(rhsEntry, ...))
    tkconfigure(rhsEntry, xscrollcommand=function(...) tkset(rhsXscroll, ...))
    tkgrid(labelRcmdr(formulaFrameMain, text="   ~ "), rhsEntry, sticky="w")
    tkgrid(labelRcmdr(formulaFrameMain, text=""), rhsXscroll, sticky="w")
    tkgrid(formulaFrameMain, labelRcmdr(formulaFrame, text="  "), formulaHelpButton, sticky="nw")
  }
  tkgrid.configure(rhsXscroll, sticky="ew")
})

## ... and including random factor r(
modelFormula3 <- defmacro(frame=top, hasLhs=TRUE, expr={
  checkAddOperator <- function(rhs){
    rhs.chars <- rev(strsplit(rhs, "")[[1]])
    if (length(rhs.chars) < 1) return(FALSE)
    check.char <- if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
      rhs.chars[1] else rhs.chars[2]
    !is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%"))
  }
  .variables <- Variables()
  word <- paste("\\[", gettextRcmdr("factor"), "\\]", sep="")
  variables <- paste(.variables,
                     ifelse(is.element(.variables, Factors()), paste("[", gettextRcmdr("factor"), "]", sep=""), ""))
  xBox <- variableListBox(frame, variables, selectmode="multiple", title=gettextRcmdr("Variables (double-click to formula)"))
  onDoubleClick <- if (!hasLhs){
    function(){
      var <- getSelection(xBox)
      tkselection.clear(xBox$listbox, "0", "end")            		
      if (length(grep(word, var)) == 1) var <- sub(word, "",  var)
      tkfocus(rhsEntry)
      rhs <- tclvalue(rhsVariable)
      rhs.chars <- rev(strsplit(rhs, "")[[1]])
      check.char <- if (length(rhs.chars) > 0){
        if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
          rhs.chars[1] else rhs.chars[2]
      }
      else ""
      tclvalue(rhsVariable) <- if (rhs == "" ||
                                   is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%")))
        paste(rhs, var, sep="")
      else paste(rhs, "+", var)
      tkicursor(rhsEntry, "end")
      tkxview.moveto(rhsEntry, "1")
    }
  }
  else{
    function(){
      var <- getSelection(xBox)
      which <- tkcurselection(xBox$listbox)
      tkselection.clear(xBox$listbox, "0", "end")
      if (length(grep(word, var)) == 1) var <- sub(word, "",  var)
      lhs <- tclvalue(lhsVariable)
      if (lhs == "" || tclvalue(tkselection.present(lhsEntry)) == "1"){
        tclvalue(lhsVariable) <- var
        tkselection.clear(lhsEntry)
        tkfocus(rhsEntry)
      }
      else {
        tkfocus(rhsEntry)
        rhs <- tclvalue(rhsVariable)
        rhs.chars <- rev(strsplit(rhs, "")[[1]])
        check.char <- if (length(rhs.chars) > 0){
          if ((rhs.chars[1] != " ") || (length(rhs.chars) == 1))
            rhs.chars[1] else rhs.chars[2]
        }
        else ""
        tclvalue(rhsVariable) <- if (rhs == "" ||
                                     is.element(check.char, c("+", "*", ":", "/", "-", "^", "(", "%", "r(")))
          paste(rhs, var, sep="")
        else paste(rhs, "+", var)
      }
      tkicursor(rhsEntry, "end")
      tkxview.moveto(rhsEntry, "1")
    }
  }
  tkbind(xBox$listbox, "<Double-ButtonPress-1>", onDoubleClick)
  onPlus <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")										
    if ((check <- !checkAddOperator(rhs)) && length(var) == 0) return()
    if (length(var) > 1){
      if (length(grep(word, var)) > 0) var <- sub(word, "",  var)
      if (length(var) > 1) var <- paste(var, collapse=" + ")
    }
    tclvalue(rhsVariable) <- paste(rhs, if (!check) " + ", var, sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onTimes <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")						
    if ((check <- !checkAddOperator(rhs)) && length(var) == 0) return()
    if (length(var) > 1){
      if (length(grep(word, var)) > 0) var <- sub(word, "",  var)
      var <- trim.blanks(var)
      if (length(var) > 1) var <- paste(var, collapse="*")
      tclvalue(rhsVariable) <- paste(rhs, if (!check) " + ", var, sep="")
    }
    else tclvalue(rhsVariable) <- paste(rhs, if (!check) "*", sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onColon <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")						
    if ((check <- !checkAddOperator(rhs)) && length(var) == 0) return()
    if (length(var) > 1){
      if (length(grep(word, var)) > 0) var <- sub(word, "",  var)
      var <- trim.blanks(var)
      if (length(var) > 1) var <- paste(var, collapse=":")
      tclvalue(rhsVariable) <- paste(rhs, if (!check) " + ", var, sep="")
    }
    else tclvalue(rhsVariable) <- paste(rhs, if (!check) ":", sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onSlash <- function(){
    rhs <- tclvalue(rhsVariable)
    if (!checkAddOperator(rhs)) return()
    tclvalue(rhsVariable) <- paste(rhs, "/",  sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onIn <- function(){
    rhs <- tclvalue(rhsVariable)
    if (!checkAddOperator(rhs)) return()
    tclvalue(rhsVariable) <- paste(rhs, "%in% ")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onMinus <- function(){
    rhs <- tclvalue(rhsVariable)
    if (!checkAddOperator(rhs)) return()
    tclvalue(rhsVariable) <- paste(rhs, "- ")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onPower <- function(){
    rhs <- tclvalue(rhsVariable)
    if (!checkAddOperator(rhs)) return()
    tclvalue(rhsVariable) <- paste(rhs, "^", sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onLeftParen <- function(){
    tkfocus(rhsEntry)
    rhs <- tclvalue(rhsVariable)
    tclvalue(rhsVariable) <- paste(rhs, "(", sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onRightParen <- function(){
    rhs <- tclvalue(rhsVariable)
    if (!checkAddOperator(rhs)) return()
    tclvalue(rhsVariable) <- paste(rhs, ")", sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onRandom <- function(){
    tkfocus(rhsEntry)
    rhs <- tclvalue(rhsVariable)
    tclvalue(rhsVariable) <- paste(rhs, "r(", sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  outerOperatorsFrame <- tkframe(frame)
  operatorsFrame <- tkframe(outerOperatorsFrame)
  splinePolyFrame <- tkframe(outerOperatorsFrame)
  plusButton <- buttonRcmdr(operatorsFrame, text="+", width="3", command=onPlus)
  timesButton <- buttonRcmdr(operatorsFrame, text="*", width="3", command=onTimes)
  colonButton <- buttonRcmdr(operatorsFrame, text=":", width="3", command=onColon)
  slashButton <- buttonRcmdr(operatorsFrame, text="/", width="3", command=onSlash)
  inButton <- buttonRcmdr(operatorsFrame, text="%in%", width="5", command=onIn)
  minusButton <- buttonRcmdr(operatorsFrame, text="-", width="3", command=onMinus)
  powerButton <- buttonRcmdr(operatorsFrame, text="^", width="3", command=onPower)
  leftParenButton <- buttonRcmdr(operatorsFrame, text="(", width="3", command=onLeftParen)
  rightParenButton <- buttonRcmdr(operatorsFrame, text=")", width="3", command=onRightParen)
  randomButton <- buttonRcmdr(operatorsFrame, text="r(", width="3", command=onRandom)
  onBSpline <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")
    if (length(var) == 0) var <- " "
    if (grepl("\\[factor\\]", var)){
      Message("spline requires a numeric variable", type="error")
      return()
    }
    if (length(var) > 1){
      Message("cannot select more than one variable", type="error")
      return()
    }
    check <- !checkAddOperator(rhs)
    tclvalue(rhsVariable) <- paste(rhs, 
                                   if (!check) paste(" + bs(", var, ", df=", tclvalue(dfSplineVar), ")", sep="") 
                                   else paste(" bs(", var, ", df=", tclvalue(dfSplineVar), ")", sep=""),
                                   sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onNatSline <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")
    if (length(var) == 0) var <- " "
    if (grepl("\\[factor\\]", var)){
      Message("spline requires a numeric variable", type="error")
      return()
    }
    if (length(var) > 1){
      Message("cannot select more than one variable", type="error")
      return()
    }
    check <- !checkAddOperator(rhs)
    tclvalue(rhsVariable) <- paste(rhs, 
                                   if (!check) paste(" + ns(", var, ", df=", tclvalue(dfSplineVar), ")", sep="") 
                                   else paste(" ns(", var, ", df=", tclvalue(dfSplineVar), ")", sep=""),
                                   sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onPoly <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")
    if (length(var) == 0) var <- " "
    if (grepl("\\[factor\\]", var)){
      Message("polynomial requires a numeric variable", type="error")
      return()
    }
    if (length(var) > 1){
      Message("cannot select more than one variable", type="error")
      return()
    }
    check <- !checkAddOperator(rhs)
    tclvalue(rhsVariable) <- paste(rhs, 
                                   if (!check) paste(" + poly(", var, ", degree=", tclvalue(degPolyVar), ")", sep="") 
                                   else paste(" poly(", var, ", degree=", tclvalue(degPolyVar), ")", sep=""),
                                   sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  onRawPoly <- function(){
    rhs <- tclvalue(rhsVariable)
    var <- getSelection(xBox)
    tkselection.clear(xBox$listbox, "0", "end")
    if (length(var) == 0) var <- " "
    if (grepl("\\[factor\\]", var)){
      Message("polynomial requires a numeric variable", type="error")
      return()
    }
    if (length(var) > 1){
      Message("cannot select more than one variable", type="error")
      return()
    }
    check <- !checkAddOperator(rhs)
    tclvalue(rhsVariable) <- paste(rhs, 
                                   if (!check) paste(" + poly(", var, ", degree=", tclvalue(degPolyVar), ", raw=TRUE)", sep="") 
                                   else paste(" poly(", var, ", degree=", tclvalue(degPolyVar), ", raw=TRUE)", sep=""),
                                   sep="")
    tkicursor(rhsEntry, "end")
    tkxview.moveto(rhsEntry, "1")
  }
  bsplineButton <- buttonRcmdr(splinePolyFrame, text=gettextRcmdr("B-spline\n"), width="10", command=onBSpline)
  nsplineButton <- buttonRcmdr(splinePolyFrame, text=gettextRcmdr("natural\nspline"), width="10", command=onNatSline)
  polyButton <- buttonRcmdr(splinePolyFrame, text=gettextRcmdr("orthogonal\npolynomial"), width="10", command=onPoly)
  RawPolyButton <- buttonRcmdr(splinePolyFrame, text=gettextRcmdr("raw\npolynomial"), width="10", command=onRawPoly)
  dfSplineVar <- tclVar("5")
  degPolyVar <- tclVar("2")
  dfDegFrame <- tkframe(outerOperatorsFrame)
  dfSplineSpin <- tkspinbox(dfDegFrame, textvariable=dfSplineVar, state="readonly", from=2, to=10, width=2)
  degPolySpin <- tkspinbox(dfDegFrame, textvariable=degPolyVar, state="readonly", from=2, to=5, width=2)
  tkgrid(plusButton, timesButton, colonButton, slashButton, inButton, minusButton,
         powerButton, leftParenButton, rightParenButton, randomButton, sticky="w")
  tkgrid(labelRcmdr(dfDegFrame, text=gettextRcmdr("df for splines: ")), dfSplineSpin,  sticky="se")
  tkgrid(labelRcmdr(dfDegFrame, text=gettextRcmdr("deg. for polynomials: ")), degPolySpin, sticky="se")
  formulaFrame <- tkframe(frame)
  formulaFrameMain <- tkframe(formulaFrame)
  onFormulaHelp <- function () print(help("formula"))
  formulaHelpButton <- buttonRcmdr(formulaFrame, text=gettextRcmdr("Model formula\nhelp"), command=onFormulaHelp,
                                   image="::image::helpIcon", compound="left")
  if (hasLhs){
    tkgrid(labelRcmdr(outerOperatorsFrame, text=gettextRcmdr("Model Formula"), 
                      fg=getRcmdr("title.color"), font="RcmdrTitleFont"), sticky="w")
    tkgrid(labelRcmdr(outerOperatorsFrame, text="Operators (click to formula):  "), operatorsFrame, sticky="nw")
    tkgrid(bsplineButton, nsplineButton, polyButton, RawPolyButton, sticky="nw")
    tkgrid(labelRcmdr(outerOperatorsFrame, text=gettextRcmdr("Splines/Polynomials:\n(select variable and click)")), 
           splinePolyFrame, dfDegFrame, sticky="nw")
    lhsVariable <- if (currentModel) tclVar(currentFields$lhs) else tclVar("")
    rhsVariable <- if (currentModel) tclVar(currentFields$rhs) else tclVar("")
    rhsEntry <- ttkentry(formulaFrameMain, width="75", textvariable=rhsVariable)
    rhsXscroll <- ttkscrollbar(formulaFrameMain,
                               orient="horizontal", command=function(...) tkxview(rhsEntry, ...))
    tkconfigure(rhsEntry, xscrollcommand=function(...) tkset(rhsXscroll, ...))
    lhsEntry <- ttkentry(formulaFrameMain, width="10", textvariable=lhsVariable)
    lhsScroll <- ttkscrollbar(formulaFrameMain,
                              orient="horizontal", command=function(...) tkxview(lhsEntry, ...))
    tkconfigure(lhsEntry, xscrollcommand=function(...) tkset(lhsScroll, ...))
    tkgrid(lhsEntry, labelRcmdr(formulaFrameMain, text=" ~    "), rhsEntry, sticky="w")
    tkgrid(lhsScroll, labelRcmdr(formulaFrameMain, text=""), rhsXscroll, sticky="w")
    tkgrid.configure(lhsScroll, sticky="ew")
    tkgrid(formulaFrameMain, labelRcmdr(formulaFrame, text="  "), formulaHelpButton, sticky="nw")
  }
  else{
    rhsVariable <- if (currentModel) tclVar(currentFields$rhs) else tclVar("")
    rhsEntry <- ttkentry(formulaFrameMain, width="75", textvariable=rhsVariable)
    rhsXscroll <- ttkscrollbar(formulaFrameMain,
                               orient="horizontal", command=function(...) tkxview(rhsEntry, ...))
    tkconfigure(rhsEntry, xscrollcommand=function(...) tkset(rhsXscroll, ...))
    tkgrid(labelRcmdr(formulaFrameMain, text="   ~ "), rhsEntry, sticky="w")
    tkgrid(labelRcmdr(formulaFrameMain, text=""), rhsXscroll, sticky="w")
    tkgrid(formulaFrameMain, labelRcmdr(formulaFrame, text="  "), formulaHelpButton, sticky="nw")
  }
  tkgrid.configure(rhsXscroll, sticky="ew")
})


#################################################
# Radio buttons inside custom frame
radioButtonsNMBU <- defmacro(the.frame, name, buttons, values=NULL, initialValue=..values[1], labels, 
                             title="", title.color="blue", right.buttons=TRUE,
                             expr={
                               ..values <- if (is.null(values)) buttons else values
                               ..frame <- paste(name, "Frame", sep="")
                               assign(..frame, the.frame)
                               ..variable <- paste(name, "Variable", sep="")
                               assign(..variable, tclVar(initialValue))
                               if(title != ""){
                                 tkgrid(labelRcmdr(eval(parse(text=..frame)), text=title, foreground=title.color), columnspan=2, sticky="w")
                               }
                               for (i in 1:length(buttons)) {
                                 ..button <- paste(buttons[i], "Button", sep="")
                                 assign(..button,
                                        ttkradiobutton(eval(parse(text=..frame)), variable=eval(parse(text=..variable)), value=..values[i]))
                                 if (right.buttons) tkgrid(labelRcmdr(eval(parse(text=..frame)), text=labels[i], justify="left"), eval(parse(text=..button)), sticky="w")
                                 else  tkgrid(eval(parse(text=..button)), labelRcmdr(eval(parse(text=..frame)), text=labels[i], justify="left"), sticky="w")
                               }
                             }
)


#################################################
# Automatic import of data
auto.import <- function(){
  MAC <- Sys.info()[1] == "Darwin"
  Dataset <- list()
  if(!MAC){
    try(Dataset[[1]] <- read.table("clipboard", header=FALSE, sep="\t", na.strings="NA", 
                                   dec=".", strip.white=FALSE))
  } else {
    try(Dataset[[1]] <- read.table(pipe("pbpaste"), header=FALSE, sep="\t", na.strings="NA", 
                                   dec=".", strip.white=FALSE))
  }
  if(length(Dataset)==0 || dim(Dataset[[1]])[2]==1){
    sep <- ""
    strip.white <- TRUE
    Dataset <- list()
    if(!MAC){
      Dataset[[1]] <- read.table("clipboard", header=FALSE, sep="", na.strings="NA", 
                                 dec=".", strip.white=TRUE)
      Dataset[[2]] <- read.table("clipboard", header=FALSE, sep="", na.strings="NA", 
                                 dec=",", strip.white=TRUE)
      Dataset[[3]] <- read.table("clipboard", header=TRUE, sep="", na.strings="NA", 
                                 dec=".", strip.white=TRUE)
      Dataset[[4]] <- read.table("clipboard", header=TRUE, sep="", na.strings="NA", 
                                 dec=",", strip.white=TRUE)
    } else {
      Dataset[[1]] <- read.table(pipe("pbpaste"), header=FALSE, sep="", na.strings="NA", 
                                 dec=".", strip.white=TRUE)
      Dataset[[2]] <- read.table(pipe("pbpaste"), header=FALSE, sep="", na.strings="NA", 
                                 dec=",", strip.white=TRUE)
      Dataset[[3]] <- read.table(pipe("pbpaste"), header=TRUE, sep="", na.strings="NA", 
                                 dec=".", strip.white=TRUE)
      Dataset[[4]] <- read.table(pipe("pbpaste"), header=TRUE, sep="", na.strings="NA", 
                                 dec=",", strip.white=TRUE)
    }
  } else {
    sep <- "\\t"
    strip.white <- FALSE
    if(!MAC){
      Dataset[[2]] <- read.table("clipboard", header=FALSE, sep="\t", na.strings="NA", 
                                 dec=",", strip.white=FALSE)
      Dataset[[3]] <- read.table("clipboard", header=TRUE, sep="\t", na.strings="NA", 
                                 dec=".", strip.white=FALSE)
      Dataset[[4]] <- read.table("clipboard", header=TRUE, sep="\t", na.strings="NA", 
                                 dec=",", strip.white=FALSE)
    } else {
      Dataset[[2]] <- read.table(pipe("pbpaste"), header=FALSE, sep="\t", na.strings="NA", 
                                 dec=",", strip.white=FALSE)
      Dataset[[3]] <- read.table(pipe("pbpaste"), header=TRUE, sep="\t", na.strings="NA", 
                                 dec=".", strip.white=FALSE)
      Dataset[[4]] <- read.table(pipe("pbpaste"), header=TRUE, sep="\t", na.strings="NA", 
                                 dec=",", strip.white=FALSE)
    }
  }
  
  D <- dim(Dataset[[1]])
  numbers <- matrix(0,ncol=4,nrow=D[2])
  for(i in 1:D[2]){
    if(class(Dataset[[1]][,i])=="integer" || class(Dataset[[1]][,i])=="numeric")
      numbers[i,1] <- 1
    if(class(Dataset[[2]][,i])=="integer" || class(Dataset[[2]][,i])=="numeric")
      numbers[i,2] <- 1
    if(class(Dataset[[3]][,i])=="integer" || class(Dataset[[3]][,i])=="numeric")
      numbers[i,3] <- 1
    if(class(Dataset[[4]][,i])=="integer" || class(Dataset[[4]][,i])=="numeric")
      numbers[i,4] <- 1
  }
  best <- which.max(apply(numbers,2,sum))
  Dataset <- Dataset[[best]]
  header <- ifelse(best>=3,TRUE,FALSE)
  dec    <- ifelse(best%%2==1,".",",")
  list(Dataset, sep=sep, strip.white=strip.white, dec=dec, header=header)
}


#################################################
# Adapted formulaFields including mer class compatibility
formulaFields2 <- function(model, hasLhs=TRUE, glm=FALSE){
  if(class(model)[1]=="mer"){
    call <- model@call
  } else {
    call <- model$call
  }
  formula <- as.character(call$formula)
  if (hasLhs){
    lhs <- formula[2]
    rhs <- formula[3]
  } else {
    lhs <- NULL
    rhs <- formula[2]
  }
  data <- as.character(call$data)
  which.subset <- which("subset" == names(call))
  subset <- if (0 == length(which.subset)) ""
  else as.character(call)[[which.subset]]
  if (glm) {
    fam <- as.character(call$family)
    family <- fam[1]
    link <- fam[2]
  }
  else {
    family <- NULL
    link <- NULL
  }
  list(lhs=lhs, rhs=rhs, data=data, subset=subset, family=family, link=link)
}
