
##== R-Commander Plug-in for Meta-Analysis  ==##

#require(MAd)
#require(metafor)

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

if (getRversion() >= '2.15.1') globalVariables(c('CatCompcmd',
                                                 'RemoveRows',
                                                 'modelC1Variable',
                                                 'modelC2Variable',
                                                 'modelFRVariable',
                                                 'subsetVariable',
                                                 'linearRegressionModel',
                                                 'modelFRVariable',
                                                 'modelNVariable',                                                 
                                                 "CatModGraphcmd",
                                                 "CatModcmd",
                                                 "ComplDatacmd",
                                                 "CorAttencmd",
                                                 "ForestPlotcmd",
                                                 "FunnelPlotcmd",
                                                 "Kappacmd",
                                                 "MAreg2cmd",
                                                 "MAregGraphcmd",
                                                 "MeanDiffcmd",
                                                 "MeanDiffdcmd",
                                                 "MeanDiffgcmd",
                                                 "MetaCcmd",
                                                 "MetaGcmd",
                                                 "MultiModGraphcmd",
                                                 "OmnibusEScmd",
                                                 "PubBiascmd",
                                                 "Wifuncmd",
                                                 "ancova_to_d1cmd",
                                                 "ancova_to_d2cmd",
                                                 "d_to_gcmd",
                                                 "f.ancova_to_dcmd",
                                                 "f_to_dcmd",
                                                 "factscmd",
                                                 "influencecmd",
                                                 "lor_to_dcmd",
                                                 "mean_to_d2cmd",
                                                 "mean_to_dcmd",
                                                 "p.ancova_to_d1cmd",
                                                 "p.ancova_to_d2cmd",
                                                 "p_to_d1cmd",
                                                 "p_to_d2cmd",
                                                 "prop_to_dcmd",
                                                 "prop_to_orcmd",
                                                 "qqnormcmd",
                                                 "r_from_chicmd",
                                                 "r_from_d1cmd",
                                                 "r_from_dcmd",
                                                 "r_from_tcmd",
                                                 "r_to_dcmd",
                                                 "radialcmd",
                                                 "rankcmd",
                                                 "regcmd",
                                                 "robustSEcmd",
                                                 "rstandardcmd",
                                                 "rstudentcmd",
                                                 "t_to_dcmd",
                                                 "trimcmd",
                                                 "tt.ancova_to_dcmd"
))

### ABOUT ####
#aboutcmd <- citation('RcmdrPlugin.MA')

aboutcmd <- function(){
  doItAndPrint(paste("This meta-analysis package provides a menu-driven, graphical user interface environment (e.g., SPSS, CMA) for conducting a meta-analysis. It uses recommended procedures as described in The Handbook of Research Synthesis and Meta-Analysis (Cooper, Hedges, &	Valentine, 2009). For more detail, see: http://acdelre.weebly.com \n 
  To cite this packages (and supporting packages), see output in main R console window (or type commands below into the console). \n"))  
  doItAndPrint(paste("citation('RcmdrPlugin.MA')"))
  doItAndPrint(paste("citation('MAd')"))
  doItAndPrint(paste("citation('metafor')"))

}

### START ####

startcmd <- function(){  
  doItAndPrint(paste("vignette('tutorial', package='RcmdrPlugin.MA')"))
  
}


### ADD MAc R PACKAGE FUNCTION 'agg' [HUNTER & SCHMIDT (2004) WITHIN-STUDY AGGREGATION FOR CORRELATION ES]

aggrs <- function (r, cor = 0.5) 
{
  k <- length(r)
  rbar <- cor
  r.xY <- sum(r)/(1 * sqrt(k + k * (k - 1) * rbar))
  return(r.xY)
}

aggC <- function (id, r, n, cor = 0.5, mod = NULL, data) 
{
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("id", "r", "n", "mod", "cor", "data"), names(mf), 
                0)
  mf <- mf[c(1, args)]
  mf[[1]] <- as.name("model.frame")
  mf.id <- mf[[match("id", names(mf))]]
  id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  mf.r <- mf[[match("r", names(mf))]]
  r <- eval(mf.r, data, enclos = sys.frame(sys.parent()))
  mf.n <- mf[[match("n", names(mf))]]
  n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
  mf.mod <- mf[[match("mod", names(mf))]]
  mod <- eval(mf.mod, data, enclos = sys.frame(sys.parent()))
  if (is.null(mod)) {
    st <- unique(id)
    out <- data.frame(id = st)
    for (i in 1:length(st)) {
      out$id[i] <- st[i]
      out$r[i] <- aggrs(r = r[id == st[i]], cor)
      out$n[i] <- round(mean(n[id == st[i]]), 0)
    }
  }
  if (!is.null(mod)) {
    st <- as.factor(id)
    st <- unique(st)
    um <- unique(mod)
    out <- data.frame(id = rep(st, rep(length(um), length(st))))
    out$mod <- rep(um, length(st))
    for (i in 1:length(st)) {
      for (j in 1:length(um)) {
        ro <- (i - 1) * length(um) + j
        m1 <- match(id, st[i], nomatch = 0)
        m2 <- match(mod, um[j], nomatch = 0)
        num <- sum(m1 * m2)
        out$r[ro] <- ifelse(num == 0, NA, aggrs(r = r[id == 
                                                        st[i] & mod == um[j]], cor))
        out$n[ro] <- round(mean(n[id == st[i] & mod ==  um[j]]), 0)
      }
    }
    out <- out[is.na(out$r) == 0, ]
  }
  return(out)
}


#== New Functions (3.07.10) ==##

# facts 3.18.10

factscmd <- function(){
  dataSet <- activeDataSet()
  initializeDialog(title=gettextRcmdr("Convert to Categorical Variable"))
  variablesFrame <- tkframe(top)
  .variable <- Variables()
  xBox <- variableListBox(variablesFrame, .variable, selectmode="multiple",
                          title=gettextRcmdr("variable to convert to categorical (pick one)"))
  newDataSetName <- tclVar(gettextRcmdr("<same as active data set>"))
  dataSetNameFrame <- tkframe(top)
  dataSetNameEntry <- ttkentry(dataSetNameFrame, width="25", textvariable=newDataSetName)
  onOK <- function(){
    value <- trim.blanks(tclvalue(newDataSetName))
    if (value == gettextRcmdr("<same as active data set>")) value <- ActiveDataSet()
    if (!is.valid.name(value)){
      errorCondition(recall=factscmd,
                     message=paste('"', value, '" ', gettextRcmdr("is not a valid name."), sep=""))
      return()
    }
    if (is.element(value, listDataSets())) {
      if ("no" == tclvalue(checkReplace(value, type=gettextRcmdr("Data set")))){
        closeDialog()
        RemoveRows()
        return()
      }
    }
    x <- getSelection(xBox)
    closeDialog()
    meta <- dataSet
    #modelN <- as.character(tclvalue(modelNVariable)) 
    #command <- paste("ComplData(", meta, ", ", paste(x, collapse=","),
    #                 ", type= '",modelN,"')", sep="")
    #command <- paste(value, " <- ", ActiveDataSet(), "[", removeRows, ",]", sep="")
    #logger(command)
    #result <- justDoIt(command)
    #command <- paste("facts(", meta, ",", meta, "$", x, ")", sep="")
    command <- paste("facts(", meta, ",'", x, "')", sep="")
#     logger(paste(value, " <- ", command, sep=""))
#     assign(value, justDoIt(command), envir=.GlobalEnv)
#     doItAndPrint(value)
    command <- paste(value, " <- ", command, sep="")
    doItAndPrint(command)
    #command2 <- (paste(value))
    #doItAndPrint(command2)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="facts", model=TRUE)
  tkgrid(labelRcmdr(dataSetNameFrame, text=gettextRcmdr("Name for new data set")), sticky="w")
  tkgrid(dataSetNameEntry, sticky="w")
  tkgrid(dataSetNameFrame, sticky="w")
  tkgrid(labelRcmdr(variablesFrame, text="    "), getFrame(xBox), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  tkgrid.configure(helpButton, sticky="e")
  dialogSuffix(rows=4, columns=2)
}


# ancova to d1

ancova_to_d1cmd <- function(){
  initializeDialog(title=gettextRcmdr("ancova to mean diff (adj. SD)"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  tmtVar <- tclVar(gettextRcmdr(" "))
  tmtFrame <- tkframe(labelsFrame)
  tmtEntry <- ttkentry(tmtFrame, width="8", textvariable=tmtVar)
  tkgrid(labelRcmdr(tmtFrame, text=gettextRcmdr("adj. tmt grp mean (ANCOVA)"), fg="blue"), sticky="w")
  tkgrid(tmtEntry, sticky="w")
  tkgrid(tmtFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  conVar <- tclVar(gettextRcmdr(" "))
  conFrame <- tkframe(labelsFrame)
  conEntry <- ttkentry(conFrame, width="8", textvariable=conVar)
  tkgrid(labelRcmdr(conFrame, text=gettextRcmdr("adj comp grp mean (ANCOVA)"), fg="blue"), sticky="w")
  tkgrid(conEntry, sticky="w")
  tkgrid(conFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  sdVar <- tclVar(gettextRcmdr(" "))
  sdFrame <- tkframe(labelsFrame)
  sdEntry <- ttkentry(sdFrame, width="8", textvariable=sdVar)
  tkgrid(labelRcmdr(sdFrame, text=gettextRcmdr("adj SD"), fg="blue"), sticky="w")
  tkgrid(sdEntry, sticky="w")
  tkgrid(sdFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("n of tmt grp"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n2Var <- tclVar(gettextRcmdr(" "))
  n2Frame <- tkframe(labelsFrame)
  n2Entry <- ttkentry(n2Frame, width="8", textvariable=n2Var)
  tkgrid(labelRcmdr(n2Frame, text=gettextRcmdr("n of comp grp"), fg="blue"), sticky="w")
  tkgrid(n2Entry, sticky="w")
  tkgrid(n2Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  corVar <- tclVar(gettextRcmdr(" "))
  corFrame <- tkframe(labelsFrame)
  corEntry <- ttkentry(corFrame, width="8", textvariable=corVar)
  tkgrid(labelRcmdr(corFrame, text=gettextRcmdr("covar/multi cor"), fg="blue"), sticky="w")
  tkgrid(corEntry, sticky="w")
  tkgrid(corFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  covVar <- tclVar(gettextRcmdr(" "))
  covFrame <- tkframe(labelsFrame)
  covEntry <- ttkentry(covFrame, width="8", textvariable=covVar)
  tkgrid(labelRcmdr(covFrame, text=gettextRcmdr("n of covar."), fg="blue"), sticky="w")
  tkgrid(covEntry, sticky="w")
  tkgrid(covFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  onOK <- function(){
    closeDialog()
    tmt <- trim.blanks(tclvalue(tmtVar))
    tmt <- paste (' ', tmt, '', sep="")   # paste(' stat, ') 
    
    con <- trim.blanks(tclvalue(conVar))
    con <- paste(', ', con, '', sep="")
    
    sd <- trim.blanks(tclvalue(sdVar))
    sd <- paste(', ', sd, '', sep="")
    
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    n2 <- trim.blanks(tclvalue(n2Var))
    n2 <- paste(', ', n2, '', sep="")
    
    cor <- trim.blanks(tclvalue(corVar))
    cor <- paste(', ', cor, '', sep="")
    
    cov <- trim.blanks(tclvalue(covVar))
    cov <- paste(', ', cov, '', sep="")
    
    doItAndPrint(paste("ancova_to_d1(", tmt, con, sd, n1, n2, cor, cov,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="ancova_to_d1")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}



# ancova to d2

ancova_to_d2cmd <- function(){
  initializeDialog(title=gettextRcmdr("ancova to mean diff (pooled SD)"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  tmtVar <- tclVar(gettextRcmdr(" "))
  tmtFrame <- tkframe(labelsFrame)
  tmtEntry <- ttkentry(tmtFrame, width="8", textvariable=tmtVar)
  tkgrid(labelRcmdr(tmtFrame, text=gettextRcmdr("adj tmt grp mean (ANCOVA)"), fg="blue"), sticky="w")
  tkgrid(tmtEntry, sticky="w")
  tkgrid(tmtFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  conVar <- tclVar(gettextRcmdr(" "))
  conFrame <- tkframe(labelsFrame)
  conEntry <- ttkentry(conFrame, width="8", textvariable=conVar)
  tkgrid(labelRcmdr(conFrame, text=gettextRcmdr("adj comp grp mean (ANCOVA)"), fg="blue"), sticky="w")
  tkgrid(conEntry, sticky="w")
  tkgrid(conFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  sdVar <- tclVar(gettextRcmdr(" "))
  sdFrame <- tkframe(labelsFrame)
  sdEntry <- ttkentry(sdFrame, width="8", textvariable=sdVar)
  tkgrid(labelRcmdr(sdFrame, text=gettextRcmdr("pooled SD"), fg="blue"), sticky="w")
  tkgrid(sdEntry, sticky="w")
  tkgrid(sdFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("n of tmt grp"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n2Var <- tclVar(gettextRcmdr(" "))
  n2Frame <- tkframe(labelsFrame)
  n2Entry <- ttkentry(n2Frame, width="8", textvariable=n2Var)
  tkgrid(labelRcmdr(n2Frame, text=gettextRcmdr("n of comp grp"), fg="blue"), sticky="w")
  tkgrid(n2Entry, sticky="w")
  tkgrid(n2Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  corVar <- tclVar(gettextRcmdr(" "))
  corFrame <- tkframe(labelsFrame)
  corEntry <- ttkentry(corFrame, width="8", textvariable=corVar)
  tkgrid(labelRcmdr(corFrame, text=gettextRcmdr("covar/multi cor"), fg="blue"), sticky="w")
  tkgrid(corEntry, sticky="w")
  tkgrid(corFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  covVar <- tclVar(gettextRcmdr(" "))
  covFrame <- tkframe(labelsFrame)
  covEntry <- ttkentry(covFrame, width="8", textvariable=covVar)
  tkgrid(labelRcmdr(covFrame, text=gettextRcmdr("n of covar"), fg="blue"), sticky="w")
  tkgrid(covEntry, sticky="w")
  tkgrid(covFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  onOK <- function(){
    closeDialog()
    tmt <- trim.blanks(tclvalue(tmtVar))
    tmt <- paste (' ', tmt, '', sep="")   # paste(' stat, ') 
    
    con <- trim.blanks(tclvalue(conVar))
    con <- paste(', ', con, '', sep="")
    
    sd <- trim.blanks(tclvalue(sdVar))
    sd <- paste(', ', sd, '', sep="")
    
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    n2 <- trim.blanks(tclvalue(n2Var))
    n2 <- paste(', ', n2, '', sep="")
    
    cor <- trim.blanks(tclvalue(corVar))
    cor <- paste(', ', cor, '', sep="")
    
    cov <- trim.blanks(tclvalue(covVar))
    cov <- paste(', ', cov, '', sep="")
    
    doItAndPrint(paste("ancova_to_d2(", tmt, con, sd, n1, n2, cor, cov,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="ancova_to_d2")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}


# d to g

d_to_gcmd <- function(){
  initializeDialog(title=gettextRcmdr("d to unbiased g"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  tmtVar <- tclVar(gettextRcmdr(" "))
  tmtFrame <- tkframe(labelsFrame)
  tmtEntry <- ttkentry(tmtFrame, width="8", textvariable=tmtVar)
  tkgrid(labelRcmdr(tmtFrame, text=gettextRcmdr("d"), fg="blue"), sticky="w")
  tkgrid(tmtEntry, sticky="w")
  tkgrid(tmtFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  conVar <- tclVar(gettextRcmdr(" "))
  conFrame <- tkframe(labelsFrame)
  conEntry <- ttkentry(conFrame, width="8", textvariable=conVar)
  tkgrid(labelRcmdr(conFrame, text=gettextRcmdr("variance of d"), fg="blue"), sticky="w")
  tkgrid(conEntry, sticky="w")
  tkgrid(conFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("n of tmt grp"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n2Var <- tclVar(gettextRcmdr(" "))
  n2Frame <- tkframe(labelsFrame)
  n2Entry <- ttkentry(n2Frame, width="8", textvariable=n2Var)
  tkgrid(labelRcmdr(n2Frame, text=gettextRcmdr("n of comp grp"), fg="blue"), sticky="w")
  tkgrid(n2Entry, sticky="w")
  tkgrid(n2Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  
  onOK <- function(){
    closeDialog()
    tmt <- trim.blanks(tclvalue(tmtVar))
    tmt <- paste (' ', tmt, '', sep="")   
    
    con <- trim.blanks(tclvalue(conVar))
    con <- paste(', ', con, '', sep="")
    
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    n2 <- trim.blanks(tclvalue(n2Var))
    n2 <- paste(', ', n2, '', sep="")
    
    doItAndPrint(paste("d_to_g(", tmt, con, n1, n2,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="d_to_g")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}



# f to d

f_to_dcmd <- function(){
  initializeDialog(title=gettextRcmdr("f-value to d"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  tmtVar <- tclVar(gettextRcmdr(" "))
  tmtFrame <- tkframe(labelsFrame)
  tmtEntry <- ttkentry(tmtFrame, width="8", textvariable=tmtVar)
  tkgrid(labelRcmdr(tmtFrame, text=gettextRcmdr("f-value"), fg="blue"), sticky="w")
  tkgrid(tmtEntry, sticky="w")
  tkgrid(tmtFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("n of tmt grp"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n2Var <- tclVar(gettextRcmdr(" "))
  n2Frame <- tkframe(labelsFrame)
  n2Entry <- ttkentry(n2Frame, width="8", textvariable=n2Var)
  tkgrid(labelRcmdr(n2Frame, text=gettextRcmdr("n of comp grp"), fg="blue"), sticky="w")
  tkgrid(n2Entry, sticky="w")
  tkgrid(n2Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  
  onOK <- function(){
    closeDialog()
    tmt <- trim.blanks(tclvalue(tmtVar))
    tmt <- paste (' ', tmt, '', sep="")   
    
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    n2 <- trim.blanks(tclvalue(n2Var))
    n2 <- paste(', ', n2, '', sep="")
    
    doItAndPrint(paste("f_to_d(", tmt, n1, n2,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="f_to_d")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}


# f ancova to d

f.ancova_to_dcmd <- function(){
  initializeDialog(title=gettextRcmdr("f-value (ANCOVA) to d"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  tmtVar <- tclVar(gettextRcmdr(" "))
  tmtFrame <- tkframe(labelsFrame)
  tmtEntry <- ttkentry(tmtFrame, width="8", textvariable=tmtVar)
  tkgrid(labelRcmdr(tmtFrame, text=gettextRcmdr("adj tmt grp mean (ANCOVA)"), fg="blue"), sticky="w")
  tkgrid(tmtEntry, sticky="w")
  tkgrid(tmtFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("n of tmt grp"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n2Var <- tclVar(gettextRcmdr(" "))
  n2Frame <- tkframe(labelsFrame)
  n2Entry <- ttkentry(n2Frame, width="8", textvariable=n2Var)
  tkgrid(labelRcmdr(n2Frame, text=gettextRcmdr("n of comp grp"), fg="blue"), sticky="w")
  tkgrid(n2Entry, sticky="w")
  tkgrid(n2Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  corVar <- tclVar(gettextRcmdr(" "))
  corFrame <- tkframe(labelsFrame)
  corEntry <- ttkentry(corFrame, width="8", textvariable=corVar)
  tkgrid(labelRcmdr(corFrame, text=gettextRcmdr("covar/multi cor"), fg="blue"), sticky="w")
  tkgrid(corEntry, sticky="w")
  tkgrid(corFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  covVar <- tclVar(gettextRcmdr(" "))
  covFrame <- tkframe(labelsFrame)
  covEntry <- ttkentry(covFrame, width="8", textvariable=covVar)
  tkgrid(labelRcmdr(covFrame, text=gettextRcmdr("n of covar"), fg="blue"), sticky="w")
  tkgrid(covEntry, sticky="w")
  tkgrid(covFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  onOK <- function(){
    closeDialog()
    tmt <- trim.blanks(tclvalue(tmtVar))
    tmt <- paste (' ', tmt, '', sep="")   # paste(' stat, ') 
    
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    n2 <- trim.blanks(tclvalue(n2Var))
    n2 <- paste(', ', n2, '', sep="")
    
    cor <- trim.blanks(tclvalue(corVar))
    cor <- paste(', ', cor, '', sep="")
    
    cov <- trim.blanks(tclvalue(covVar))
    cov <- paste(', ', cov, '', sep="")
    
    doItAndPrint(paste("f.ancova_to_d(", tmt, n1, n2, cor, cov,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="f.ancova_to_d")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}



# lor to d

lor_to_dcmd <- function(){
  initializeDialog(title=gettextRcmdr("log odds ratio to d"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  tmtVar <- tclVar(gettextRcmdr(" "))
  tmtFrame <- tkframe(labelsFrame)
  tmtEntry <- ttkentry(tmtFrame, width="8", textvariable=tmtVar)
  tkgrid(labelRcmdr(tmtFrame, text=gettextRcmdr("log odds ratio (LOR)"), fg="blue"), sticky="w")
  tkgrid(tmtEntry, sticky="w")
  tkgrid(tmtFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("variance of LOR"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  onOK <- function(){
    closeDialog()
    tmt <- trim.blanks(tclvalue(tmtVar))
    tmt <- paste (' ', tmt, '', sep="")   
    
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    doItAndPrint(paste("lor_to_d(", tmt, n1,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="lor_to_d")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}


# mean to d

mean_to_dcmd <- function(){
  initializeDialog(title=gettextRcmdr("raw means and SDs to d"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  tmtVar <- tclVar(gettextRcmdr(" "))
  tmtFrame <- tkframe(labelsFrame)
  tmtEntry <- ttkentry(tmtFrame, width="8", textvariable=tmtVar)
  tkgrid(labelRcmdr(tmtFrame, text=gettextRcmdr("tmt grp mean"), fg="blue"), sticky="w")
  tkgrid(tmtEntry, sticky="w")
  tkgrid(tmtFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  conVar <- tclVar(gettextRcmdr(" "))
  conFrame <- tkframe(labelsFrame)
  conEntry <- ttkentry(conFrame, width="8", textvariable=conVar)
  tkgrid(labelRcmdr(conFrame, text=gettextRcmdr("comp grp mean"), fg="blue"), sticky="w")
  tkgrid(conEntry, sticky="w")
  tkgrid(conFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  sdVar <- tclVar(gettextRcmdr(" "))
  sdFrame <- tkframe(labelsFrame)
  sdEntry <- ttkentry(sdFrame, width="8", textvariable=sdVar)
  tkgrid(labelRcmdr(sdFrame, text=gettextRcmdr("SD of tmt grp"), fg="blue"), sticky="w")
  tkgrid(sdEntry, sticky="w")
  tkgrid(sdFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("n of tmt grp"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n2Var <- tclVar(gettextRcmdr(" "))
  n2Frame <- tkframe(labelsFrame)
  n2Entry <- ttkentry(n2Frame, width="8", textvariable=n2Var)
  tkgrid(labelRcmdr(n2Frame, text=gettextRcmdr("n of comp grp"), fg="blue"), sticky="w")
  tkgrid(n2Entry, sticky="w")
  tkgrid(n2Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  corVar <- tclVar(gettextRcmdr(" "))
  corFrame <- tkframe(labelsFrame)
  corEntry <- ttkentry(corFrame, width="8", textvariable=corVar)
  tkgrid(labelRcmdr(corFrame, text=gettextRcmdr("SD of comp grp"), fg="blue"), sticky="w")
  tkgrid(corEntry, sticky="w")
  tkgrid(corFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  onOK <- function(){
    closeDialog()
    tmt <- trim.blanks(tclvalue(tmtVar))
    tmt <- paste (' ', tmt, '', sep="")   # paste(' stat, ') 
    
    con <- trim.blanks(tclvalue(conVar))
    con <- paste(', ', con, '', sep="")
    
    sd <- trim.blanks(tclvalue(sdVar))
    sd <- paste(', ', sd, '', sep="")
    
    cor <- trim.blanks(tclvalue(corVar))
    cor <- paste(', ', cor, '', sep="")
    
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    n2 <- trim.blanks(tclvalue(n2Var))
    n2 <- paste(', ', n2, '', sep="")
    
    
    
    
    doItAndPrint(paste("mean_to_d(", tmt, con, sd, cor, n1, n2,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="mean_to_d")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}

# mean to d II

mean_to_d2cmd <- function(){
  initializeDialog(title=gettextRcmdr("means with pooled SD to d"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  tmtVar <- tclVar(gettextRcmdr(" "))
  tmtFrame <- tkframe(labelsFrame)
  tmtEntry <- ttkentry(tmtFrame, width="8", textvariable=tmtVar)
  tkgrid(labelRcmdr(tmtFrame, text=gettextRcmdr("tmt grp mean"), fg="blue"), sticky="w")
  tkgrid(tmtEntry, sticky="w")
  tkgrid(tmtFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  conVar <- tclVar(gettextRcmdr(" "))
  conFrame <- tkframe(labelsFrame)
  conEntry <- ttkentry(conFrame, width="8", textvariable=conVar)
  tkgrid(labelRcmdr(conFrame, text=gettextRcmdr("comp grp mean"), fg="blue"), sticky="w")
  tkgrid(conEntry, sticky="w")
  tkgrid(conFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  sdVar <- tclVar(gettextRcmdr(" "))
  sdFrame <- tkframe(labelsFrame)
  sdEntry <- ttkentry(sdFrame, width="8", textvariable=sdVar)
  tkgrid(labelRcmdr(sdFrame, text=gettextRcmdr("pooled SD"), fg="blue"), sticky="w")
  tkgrid(sdEntry, sticky="w")
  tkgrid(sdFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("n of tmt grp"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n2Var <- tclVar(gettextRcmdr(" "))
  n2Frame <- tkframe(labelsFrame)
  n2Entry <- ttkentry(n2Frame, width="8", textvariable=n2Var)
  tkgrid(labelRcmdr(n2Frame, text=gettextRcmdr("n of comp grp"), fg="blue"), sticky="w")
  tkgrid(n2Entry, sticky="w")
  tkgrid(n2Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  onOK <- function(){
    closeDialog()
    tmt <- trim.blanks(tclvalue(tmtVar))
    tmt <- paste (' ', tmt, '', sep="")   # paste(' stat, ') 
    
    con <- trim.blanks(tclvalue(conVar))
    con <- paste(', ', con, '', sep="")
    
    sd <- trim.blanks(tclvalue(sdVar))
    sd <- paste(', ', sd, '', sep="")
    
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    n2 <- trim.blanks(tclvalue(n2Var))
    n2 <- paste(', ', n2, '', sep="")
    doItAndPrint(paste("mean_to_d2(", tmt, con, sd, n1, n2,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="mean_to_d2")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}


# mean diff (d & g)

MeanDiffcmd <- function(){
  initializeDialog(title=gettextRcmdr("Automated mean differences (g)"))
  variablesFrame <- tkframe(top)
  UpdateModelNumber()
  modelName <- tclVar(paste("meandiff.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  subsetBox()
  onOK <- function(){ 
    value <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(value)){
      UpdateModelNumber(-1)
      errorCondition(recall=MeanDiffcmd, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), value))
      return()
    }
    closeDialog() 
    meta <- ActiveDataSet()
    modelN <- as.character(tclvalue(modelNVariable)) 
    command <- paste(paste("MeanDiff(", meta, ", denom= '",modelN,"')", sep=""))
#     logger(paste(value, " <- ", command, sep=""))
#     assign(value, justDoIt(command), envir=.GlobalEnv)
#     doItAndPrint(value)
    command <- paste(value, " <- ", command, sep="")
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="MeanDiff", model=TRUE)
  radioButtons(name="modelN", buttons=c("pooled.sd", "control.sd"), values=c("pooled.sd", "control.sd"),   
               labels=gettextRcmdr(c("pooled.sd", "control.sd")), title=gettextRcmdr("demoninator:"))
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for data:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
  tkgrid.configure(helpButton, sticky="e")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(modelNFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=4, columns=2)
}  


# mean diff (d)

MeanDiffdcmd <- function(){
  initializeDialog(title=gettextRcmdr("Automated mean differences (d)"))
  variablesFrame <- tkframe(top)
  UpdateModelNumber()
  modelName <- tclVar(paste("meandiff.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  subsetBox()
  onOK <- function(){ 
    value <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(value)){
      UpdateModelNumber(-1)
      errorCondition(recall=MeanDiffdcmd, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), value))
      return()
    }
    closeDialog() 
    meta <- ActiveDataSet()
    modelN <- as.character(tclvalue(modelNVariable)) 
    command <- paste(paste("MeanDiffd(", meta, ", denom= '",modelN,"')", sep=""))
#     logger(paste(value, " <- ", command, sep=""))
#     assign(value, justDoIt(command), envir=.GlobalEnv)
#     doItAndPrint(value)
    command <- paste(value, " <- ", command, sep="")
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="MeanDiffd", model=TRUE)
  radioButtons(name="modelN", buttons=c("pooled.sd", "control.sd"), values=c("pooled.sd", "control.sd"),   
               labels=gettextRcmdr(c("pooled.sd", "control.sd")), title=gettextRcmdr("demoninator:"))
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for data:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
  tkgrid.configure(helpButton, sticky="e")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(modelNFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=4, columns=2)
}  


# mean diff (g)

MeanDiffgcmd <- function(){
  initializeDialog(title=gettextRcmdr("Automated mean differences (g)"))
  variablesFrame <- tkframe(top)
  UpdateModelNumber()
  modelName <- tclVar(paste("meandiff.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  #subsetBox()
  onOK <- function(){ 
    value <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(value)){
      UpdateModelNumber(-1)
      errorCondition(recall=MeanDiffdcmd, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), value))
      return()
    }
    closeDialog() 
    meta <- ActiveDataSet()
    #modelN <- as.character(tclvalue(modelNVariable)) 
    command <- paste(paste("MeanDiffd(", meta,")", sep=""))
#     logger(paste(value, " <- ", command, sep=""))
#     assign(value, justDoIt(command), envir=.GlobalEnv)
#     doItAndPrint(value)
    command <- paste(value, " <- ", command, sep="")
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="MeanDiffg", model=TRUE)
  #radioButtons(name="modelN", buttons=c("pooled.sd", "control.sd"), values=c("pooled.sd", "control.sd"),   
  #             labels=gettextRcmdr(c("pooled.sd", "control.sd")), title=gettextRcmdr("demoninator:"))
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for data:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
  tkgrid.configure(helpButton, sticky="e")
  tkgrid(variablesFrame, sticky="w")
  #tkgrid(modelNFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=4, columns=2)
}  

MeanDiffgcmd <- function(){
  initializeDialog(title=gettextRcmdr("Automated mean differences (g)"))
  labelsFrame <- tkframe(top)
  variablesFrame <- tkframe(top)
  .variable <- Variables()
  .numeric <- Numeric()
  statFrame <- tkframe(labelsFrame)
  
  xBox <- variableListBox(variablesFrame, .variable, 
                          title=gettextRcmdr("Grp 1 N variable (pick one)"))
  yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Grp 1 post-mean value (pick one)"))
  zBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Grp 1 std. deviation of mean value (pick one)"))
  aBox <- variableListBox(variablesFrame, .variable, 
                          title=gettextRcmdr("Grp 2 N variable (pick one)"))
  bBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Grp 2 post-mean value (pick one)"))
  cBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Grp 2 std. deviation of mean value (pick one)"))
  
  UpdateModelNumber()
  #statVar <- tclVar(gettextRcmdr(".50 "))  
  #statFrame <- tkframe(labelsFrame)
  #statEntry <- ttkentry(statFrame, width="5", textvariable=statVar)
  #tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("estimated correlation btwn outcome measures"), fg="blue"), sticky="w")
  #tkgrid(statEntry, sticky="w")
  #tkgrid(statFrame, labelRcmdr(labelsFrame, text=" Default is .50  (Wampold, 1997)  "), sticky="w")
  
  modelName <- tclVar(paste("dat.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  
  onOK <- function(){
    x <- getSelection(xBox)
    y <- getSelection(yBox)
    z <- getSelection(zBox)
    a <- getSelection(aBox)
    b <- getSelection(bBox)
    c <- getSelection(cBox)
    closeDialog()
    if (0 == length(x)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=gettextRcmdr("No  sample (n.1) variable selected."))
      return()
    }
    if (0 == length(y)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=gettextRcmdr("No  post mean (m.1) variable selected."))
      return()
    }
    if (0 == length(z)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=gettextRcmdr("No sd (sd.1) variable selected."))
      return()
    }
    if (0 == length(a)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=gettextRcmdr("No  sample (n.2) variable selected."))
      return()
    }
    if (0 == length(b)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=gettextRcmdr("No  post mean (m.2) variable selected."))
      return()
    }
    if (0 == length(c)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=gettextRcmdr("No sd (sd.2) variable selected."))
      return()
    }
    value <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(value)){
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), value))
      return()
    }
    modelN <- as.character(tclvalue(modelNVariable)) 
    #modelFR <- as.character(tclvalue(modelFRVariable)) 
    #stat <- trim.blanks(tclvalue(statVar))
    #stat <- paste (' ', stat, '', sep="")
    command <- paste("compute_dgs(" , x, ", ", y,", ", z,", ", a,", ", b,", ", c,",
                     denom= '",modelN,"', data=", ActiveDataSet(), ")", sep="")
#     logger(paste(value, " <- ", command, sep=""))
#     assign(value, justDoIt(command), envir=.GlobalEnv)
    command <- paste(value, " <- ", command, sep="")
    justDoIt(command)
  }
  OKCancelHelp(helpSubject="compute_dgs", model=TRUE)
  radioButtons(name="modelN", buttons=c("pooled.sd", "control.sd"), values=c("pooled.sd", "control.sd"),   
               labels=gettextRcmdr(c("pooled.sd", "control.sd")), title=gettextRcmdr("demoninator:"))  
  #tkgrid(modelFRFrame, sticky="w")
  tkgrid(modelNFrame, sticky="w")
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for data.frame:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
  tkgrid(getFrame(xBox), labelRcmdr(variablesFrame, text="    "), getFrame(yBox),getFrame(zBox), 
         getFrame(aBox),getFrame(bBox),getFrame(cBox),sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  tkgrid.configure(helpButton, sticky="e")
  dialogSuffix(rows=7, columns=7)
}

# p ancova to d

p.ancova_to_d1cmd <- function(){
  initializeDialog(title=gettextRcmdr("one-tailed p-value (ANCOVA) to d"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  tmtVar <- tclVar(gettextRcmdr(" "))
  tmtFrame <- tkframe(labelsFrame)
  tmtEntry <- ttkentry(tmtFrame, width="8", textvariable=tmtVar)
  tkgrid(labelRcmdr(tmtFrame, text=gettextRcmdr("p-value from an ANCOVA"), fg="blue"), sticky="w")
  tkgrid(tmtEntry, sticky="w")
  tkgrid(tmtFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("n of tmt grp"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n2Var <- tclVar(gettextRcmdr(" "))
  n2Frame <- tkframe(labelsFrame)
  n2Entry <- ttkentry(n2Frame, width="8", textvariable=n2Var)
  tkgrid(labelRcmdr(n2Frame, text=gettextRcmdr("n of comp grp"), fg="blue"), sticky="w")
  tkgrid(n2Entry, sticky="w")
  tkgrid(n2Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  corVar <- tclVar(gettextRcmdr(" "))
  corFrame <- tkframe(labelsFrame)
  corEntry <- ttkentry(corFrame, width="8", textvariable=corVar)
  tkgrid(labelRcmdr(corFrame, text=gettextRcmdr("covar/multi cor"), fg="blue"), sticky="w")
  tkgrid(corEntry, sticky="w")
  tkgrid(corFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  covVar <- tclVar(gettextRcmdr(" "))
  covFrame <- tkframe(labelsFrame)
  covEntry <- ttkentry(covFrame, width="8", textvariable=covVar)
  tkgrid(labelRcmdr(covFrame, text=gettextRcmdr("n of covar"), fg="blue"), sticky="w")
  tkgrid(covEntry, sticky="w")
  tkgrid(covFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  onOK <- function(){
    closeDialog()
    tmt <- trim.blanks(tclvalue(tmtVar))
    tmt <- paste (' ', tmt, '', sep="")   # paste(' stat, ') 
    
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    n2 <- trim.blanks(tclvalue(n2Var))
    n2 <- paste(', ', n2, '', sep="")
    
    cor <- trim.blanks(tclvalue(corVar))
    cor <- paste(', ', cor, '', sep="")
    
    cov <- trim.blanks(tclvalue(covVar))
    cov <- paste(', ', cov, '', sep="")
    
    doItAndPrint(paste("p.ancova_to_d1(", tmt, n1, n2, cor, cov,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="p.ancova_to_d1")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}

# p ancova to d II

p.ancova_to_d2cmd <- function(){
  initializeDialog(title=gettextRcmdr("two-tailed p-value (ANCOVA) to d"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  tmtVar <- tclVar(gettextRcmdr(" "))
  tmtFrame <- tkframe(labelsFrame)
  tmtEntry <- ttkentry(tmtFrame, width="8", textvariable=tmtVar)
  tkgrid(labelRcmdr(tmtFrame, text=gettextRcmdr("two-tailed p-value from an ANCOVA"), fg="blue"), sticky="w")
  tkgrid(tmtEntry, sticky="w")
  tkgrid(tmtFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("n of tmt grp"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n2Var <- tclVar(gettextRcmdr(" "))
  n2Frame <- tkframe(labelsFrame)
  n2Entry <- ttkentry(n2Frame, width="8", textvariable=n2Var)
  tkgrid(labelRcmdr(n2Frame, text=gettextRcmdr("n of comp grp"), fg="blue"), sticky="w")
  tkgrid(n2Entry, sticky="w")
  tkgrid(n2Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  corVar <- tclVar(gettextRcmdr(" "))
  corFrame <- tkframe(labelsFrame)
  corEntry <- ttkentry(corFrame, width="8", textvariable=corVar)
  tkgrid(labelRcmdr(corFrame, text=gettextRcmdr("covar/multi cor"), fg="blue"), sticky="w")
  tkgrid(corEntry, sticky="w")
  tkgrid(corFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  covVar <- tclVar(gettextRcmdr(" "))
  covFrame <- tkframe(labelsFrame)
  covEntry <- ttkentry(covFrame, width="8", textvariable=covVar)
  tkgrid(labelRcmdr(covFrame, text=gettextRcmdr("n of covar"), fg="blue"), sticky="w")
  tkgrid(covEntry, sticky="w")
  tkgrid(covFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  onOK <- function(){
    closeDialog()
    tmt <- trim.blanks(tclvalue(tmtVar))
    tmt <- paste (' ', tmt, '', sep="")   # paste(' stat, ') 
    
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    n2 <- trim.blanks(tclvalue(n2Var))
    n2 <- paste(', ', n2, '', sep="")
    
    cor <- trim.blanks(tclvalue(corVar))
    cor <- paste(', ', cor, '', sep="")
    
    cov <- trim.blanks(tclvalue(covVar))
    cov <- paste(', ', cov, '', sep="")
    
    doItAndPrint(paste("p.ancova_to_d2(", tmt, n1, n2, cor, cov,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="p.ancova_to_d2")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}


# p to d1

p_to_d1cmd <- function(){
  initializeDialog(title=gettextRcmdr("one-tailed p-value to d"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  tmtVar <- tclVar(gettextRcmdr(" "))
  tmtFrame <- tkframe(labelsFrame)
  tmtEntry <- ttkentry(tmtFrame, width="8", textvariable=tmtVar)
  tkgrid(labelRcmdr(tmtFrame, text=gettextRcmdr("p-value"), fg="blue"), sticky="w")
  tkgrid(tmtEntry, sticky="w")
  tkgrid(tmtFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("n of tmt grp"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n2Var <- tclVar(gettextRcmdr(" "))
  n2Frame <- tkframe(labelsFrame)
  n2Entry <- ttkentry(n2Frame, width="8", textvariable=n2Var)
  tkgrid(labelRcmdr(n2Frame, text=gettextRcmdr("n of comp grp"), fg="blue"), sticky="w")
  tkgrid(n2Entry, sticky="w")
  tkgrid(n2Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  
  onOK <- function(){
    closeDialog()
    tmt <- trim.blanks(tclvalue(tmtVar))
    tmt <- paste (' ', tmt, '', sep="")   
    
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    n2 <- trim.blanks(tclvalue(n2Var))
    n2 <- paste(', ', n2, '', sep="")
    
    doItAndPrint(paste("p_to_d1(", tmt, n1, n2,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="p_to_d1")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}

# p to d2

p_to_d2cmd <- function(){
  initializeDialog(title=gettextRcmdr("two-tailed p-value to d"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  tmtVar <- tclVar(gettextRcmdr(" "))
  tmtFrame <- tkframe(labelsFrame)
  tmtEntry <- ttkentry(tmtFrame, width="8", textvariable=tmtVar)
  tkgrid(labelRcmdr(tmtFrame, text=gettextRcmdr("p-value"), fg="blue"), sticky="w")
  tkgrid(tmtEntry, sticky="w")
  tkgrid(tmtFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("n of tmt grp"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n2Var <- tclVar(gettextRcmdr(" "))
  n2Frame <- tkframe(labelsFrame)
  n2Entry <- ttkentry(n2Frame, width="8", textvariable=n2Var)
  tkgrid(labelRcmdr(n2Frame, text=gettextRcmdr("n of comp grp"), fg="blue"), sticky="w")
  tkgrid(n2Entry, sticky="w")
  tkgrid(n2Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  
  onOK <- function(){
    closeDialog()
    tmt <- trim.blanks(tclvalue(tmtVar))
    tmt <- paste (' ', tmt, '', sep="")   
    
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    n2 <- trim.blanks(tclvalue(n2Var))
    n2 <- paste(', ', n2, '', sep="")
    
    doItAndPrint(paste("p_to_d2(", tmt, n1, n2,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="p_to_d2")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}


# prop to d

prop_to_dcmd <- function(){
  initializeDialog(title=gettextRcmdr("proportions to d"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  tmtVar <- tclVar(gettextRcmdr(" "))
  tmtFrame <- tkframe(labelsFrame)
  tmtEntry <- ttkentry(tmtFrame, width="8", textvariable=tmtVar)
  tkgrid(labelRcmdr(tmtFrame, text=gettextRcmdr("proportion one"), fg="blue"), sticky="w")
  tkgrid(tmtEntry, sticky="w")
  tkgrid(tmtFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("n of tmt grp"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n2Var <- tclVar(gettextRcmdr(" "))
  n2Frame <- tkframe(labelsFrame)
  n2Entry <- ttkentry(n2Frame, width="8", textvariable=n2Var)
  tkgrid(labelRcmdr(n2Frame, text=gettextRcmdr("n of comp grp"), fg="blue"), sticky="w")
  tkgrid(n2Entry, sticky="w")
  tkgrid(n2Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  corVar <- tclVar(gettextRcmdr(" "))
  corFrame <- tkframe(labelsFrame)
  corEntry <- ttkentry(corFrame, width="8", textvariable=corVar)
  tkgrid(labelRcmdr(corFrame, text=gettextRcmdr("proportion two"), fg="blue"), sticky="w")
  tkgrid(corEntry, sticky="w")
  tkgrid(corFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  onOK <- function(){
    closeDialog()
    tmt <- trim.blanks(tclvalue(tmtVar))
    tmt <- paste (' ', tmt, '', sep="")   # proportion 1 
    
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    n2 <- trim.blanks(tclvalue(n2Var))
    n2 <- paste(', ', n2, '', sep="")
    
    cor <- trim.blanks(tclvalue(corVar)) #prop 2
    cor <- paste(', ', cor, '', sep="")
    
    
    doItAndPrint(paste("prop_to_d(", tmt, cor, n1, n2,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="prop_to_d")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}


# prop to or

prop_to_orcmd <- function(){
  initializeDialog(title=gettextRcmdr("proportions to odds ratio"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  tmtVar <- tclVar(gettextRcmdr(" "))
  tmtFrame <- tkframe(labelsFrame)
  tmtEntry <- ttkentry(tmtFrame, width="8", textvariable=tmtVar)
  tkgrid(labelRcmdr(tmtFrame, text=gettextRcmdr("proportion one"), fg="blue"), sticky="w")
  tkgrid(tmtEntry, sticky="w")
  tkgrid(tmtFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("total sample size for group A and B"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n2Var <- tclVar(gettextRcmdr(" "))
  n2Frame <- tkframe(labelsFrame)
  n2Entry <- ttkentry(n2Frame, width="8", textvariable=n2Var)
  tkgrid(labelRcmdr(n2Frame, text=gettextRcmdr("total sample size for group C and D"), fg="blue"), sticky="w")
  tkgrid(n2Entry, sticky="w")
  tkgrid(n2Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  corVar <- tclVar(gettextRcmdr(" "))
  corFrame <- tkframe(labelsFrame)
  corEntry <- ttkentry(corFrame, width="8", textvariable=corVar)
  tkgrid(labelRcmdr(corFrame, text=gettextRcmdr("proportion two"), fg="blue"), sticky="w")
  tkgrid(corEntry, sticky="w")
  tkgrid(corFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  onOK <- function(){
    closeDialog()
    tmt <- trim.blanks(tclvalue(tmtVar))
    tmt <- paste (' ', tmt, '', sep="")   # proportion 1 
    
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    n2 <- trim.blanks(tclvalue(n2Var))
    n2 <- paste(', ', n2, '', sep="")
    
    cor <- trim.blanks(tclvalue(corVar)) #prop 2
    cor <- paste(', ', cor, '', sep="")
    
    
    doItAndPrint(paste("prop_to_or(", tmt, cor, n1, n2,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="prop_to_or")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}


# pub bias

PubBiascmd <- function(){
  initializeDialog(title=gettextRcmdr("Publication bias"))
  variablesFrame <- tkframe(top)
  UpdateModelNumber()
  modelName <- tclVar(paste("meandiff.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  #subsetBox()
  onOK <- function(){ 
    value <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(value)){
      UpdateModelNumber(-1)
      errorCondition(recall=MeanDiffdcmd, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), value))
      return()
    }
    closeDialog() 
    meta <- ActiveDataSet()
    command <- paste(paste("PubBias(", meta,")", sep=""))
#     logger(paste(value, " <- ", command, sep=""))
#     assign(value, justDoIt(command), envir=.GlobalEnv)
#     doItAndPrint(value)
    command <- paste(value, " <- ", command, sep="")
    doItAndPrint(command)
    command2 <- (paste(value))
    doItAndPrint(command2)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="PubBias", model=TRUE)
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for data:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
  tkgrid.configure(helpButton, sticky="e")
  tkgrid(variablesFrame, sticky="w")
  #tkgrid(modelNFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=4, columns=2)
}  


# r to d

r_to_dcmd <- function(){
  initializeDialog(title=gettextRcmdr("correlation to d"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("sample size"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  corVar <- tclVar(gettextRcmdr(" "))
  corFrame <- tkframe(labelsFrame)
  corEntry <- ttkentry(corFrame, width="8", textvariable=corVar)
  tkgrid(labelRcmdr(corFrame, text=gettextRcmdr("correlation value"), fg="blue"), sticky="w")
  tkgrid(corEntry, sticky="w")
  tkgrid(corFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  onOK <- function(){
    closeDialog()
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    cor <- trim.blanks(tclvalue(corVar)) #prop 2
    cor <- paste( ' ', cor, '', sep="")   # removed from after paste( : ',
    
    
    doItAndPrint(paste("r_to_d(", cor, n1,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="r_to_d")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}

# t ancova to d

tt.ancova_to_dcmd <- function(){
  initializeDialog(title=gettextRcmdr("t-value (ANCOVA) to d"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  tmtVar <- tclVar(gettextRcmdr(" "))
  tmtFrame <- tkframe(labelsFrame)
  tmtEntry <- ttkentry(tmtFrame, width="8", textvariable=tmtVar)
  tkgrid(labelRcmdr(tmtFrame, text=gettextRcmdr("t-value (ANCOVA)"), fg="blue"), sticky="w")
  tkgrid(tmtEntry, sticky="w")
  tkgrid(tmtFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("n of tmt grp"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n2Var <- tclVar(gettextRcmdr(" "))
  n2Frame <- tkframe(labelsFrame)
  n2Entry <- ttkentry(n2Frame, width="8", textvariable=n2Var)
  tkgrid(labelRcmdr(n2Frame, text=gettextRcmdr("n of comp grp"), fg="blue"), sticky="w")
  tkgrid(n2Entry, sticky="w")
  tkgrid(n2Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  corVar <- tclVar(gettextRcmdr(" "))
  corFrame <- tkframe(labelsFrame)
  corEntry <- ttkentry(corFrame, width="8", textvariable=corVar)
  tkgrid(labelRcmdr(corFrame, text=gettextRcmdr("covar/multi cor"), fg="blue"), sticky="w")
  tkgrid(corEntry, sticky="w")
  tkgrid(corFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  covVar <- tclVar(gettextRcmdr(" "))
  covFrame <- tkframe(labelsFrame)
  covEntry <- ttkentry(covFrame, width="8", textvariable=covVar)
  tkgrid(labelRcmdr(covFrame, text=gettextRcmdr("n of covar"), fg="blue"), sticky="w")
  tkgrid(covEntry, sticky="w")
  tkgrid(covFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  onOK <- function(){
    closeDialog()
    tmt <- trim.blanks(tclvalue(tmtVar))
    tmt <- paste (' ', tmt, '', sep="")   # paste(' stat, ') 
    
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    n2 <- trim.blanks(tclvalue(n2Var))
    n2 <- paste(', ', n2, '', sep="")
    
    cor <- trim.blanks(tclvalue(corVar))
    cor <- paste(', ', cor, '', sep="")
    
    cov <- trim.blanks(tclvalue(covVar))
    cov <- paste(', ', cov, '', sep="")
    
    doItAndPrint(paste("tt.ancova_to_d(", tmt, n1, n2, cor, cov,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="tt.ancova_to_d")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}


# t to d

t_to_dcmd <- function(){
  initializeDialog(title=gettextRcmdr("t-value to d"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  
  tmtVar <- tclVar(gettextRcmdr(" "))
  tmtFrame <- tkframe(labelsFrame)
  tmtEntry <- ttkentry(tmtFrame, width="8", textvariable=tmtVar)
  tkgrid(labelRcmdr(tmtFrame, text=gettextRcmdr("t-value"), fg="blue"), sticky="w")
  tkgrid(tmtEntry, sticky="w")
  tkgrid(tmtFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n1Var <- tclVar(gettextRcmdr(" "))
  n1Frame <- tkframe(labelsFrame)
  n1Entry <- ttkentry(n1Frame, width="8", textvariable=n1Var)
  tkgrid(labelRcmdr(n1Frame, text=gettextRcmdr("n of tmt group"), fg="blue"), sticky="w")
  tkgrid(n1Entry, sticky="w")
  tkgrid(n1Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  n2Var <- tclVar(gettextRcmdr(" "))
  n2Frame <- tkframe(labelsFrame)
  n2Entry <- ttkentry(n2Frame, width="8", textvariable=n2Var)
  tkgrid(labelRcmdr(n2Frame, text=gettextRcmdr("n of comp group"), fg="blue"), sticky="w")
  tkgrid(n2Entry, sticky="w")
  tkgrid(n2Frame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  
  onOK <- function(){
    closeDialog()
    tmt <- trim.blanks(tclvalue(tmtVar))
    tmt <- paste (' ', tmt, '', sep="")   # proportion 1 
    
    n1 <- trim.blanks(tclvalue(n1Var))
    n1 <- paste(', ', n1, '', sep="")
    
    n2 <- trim.blanks(tclvalue(n2Var))
    n2 <- paste(', ', n2, '', sep="")
    
    doItAndPrint(paste("t_to_d(", tmt, n1, n2,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="t_to_d")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}


# weight function

Wifuncmd <- function(){
  initializeDialog(title=gettextRcmdr("Add weights to dataset"))
  variablesFrame <- tkframe(top)
  UpdateModelNumber()
  modelName <- tclVar(paste("df.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  #subsetBox()
  onOK <- function(){ 
    value <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(value)){
      UpdateModelNumber(-1)
      errorCondition(recall=MeanDiffdcmd, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), value))
      return()
    }
    closeDialog() 
    meta <- ActiveDataSet()
    command <- paste(paste("Wifun(", meta,")", sep=""))
#     logger(paste(value, " <- ", command, sep=""))
#     assign(value, justDoIt(command), envir=.GlobalEnv)
#     doItAndPrint(value)
    command <- paste(value, " <- ", command, sep="")
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="Wifun", model=TRUE)
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for data:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
  tkgrid.configure(helpButton, sticky="e")
  tkgrid(variablesFrame, sticky="w")
  #tkgrid(modelNFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=4, columns=2)
}  



##==== Calculate ES ====##


# r_from_chi

r_from_chicmd <- function(){
  initializeDialog(title=gettextRcmdr("r from chi-squared"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  statVar <- tclVar(gettextRcmdr(" "))
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="8", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("reported statistic"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  modnameVar <- tclVar(gettextRcmdr(" "))
  modnameFrame <- tkframe(labelsFrame)
  modnameEntry <- ttkentry(modnameFrame, width="8", textvariable=modnameVar)
  tkgrid(labelRcmdr(modnameFrame, text=gettextRcmdr("sample size"), fg="blue"), sticky="w")
  tkgrid(modnameEntry, sticky="w")
  tkgrid(modnameFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  onOK <- function(){
    closeDialog()
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="")   # paste(' stat, ') 
    n <- trim.blanks(tclvalue(modnameVar))
    n <- paste(', ', n, '', sep="")
    doItAndPrint(paste("r_from_chi(", stat, n,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="r_from_chi")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}



# r_from_d 

r_from_dcmd <- function(){
  initializeDialog(title=gettextRcmdr("r from mean difference"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  statVar <- tclVar(gettextRcmdr(" "))
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="8", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("reported d statistic"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  modnameVar <- tclVar(gettextRcmdr(" "))
  modnameFrame <- tkframe(labelsFrame)
  modnameEntry <- ttkentry(modnameFrame, width="8", textvariable=modnameVar)
  tkgrid(labelRcmdr(modnameFrame, text=gettextRcmdr("variance of d"), fg="blue"), sticky="w")
  tkgrid(modnameEntry, sticky="w")
  tkgrid(modnameFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  onOK <- function(){
    closeDialog()
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="")   # paste(' stat, ') 
    var.d <- trim.blanks(tclvalue(modnameVar))
    var.d <- paste(', ', var.d, '', sep="")
    doItAndPrint(paste("r_from_d(", stat, var.d,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="r_from_d")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}


# r_from_d1 

r_from_d1cmd <- function(){
  initializeDialog(title=gettextRcmdr("r from mean difference II"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  statVar <- tclVar(gettextRcmdr(" "))
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="8", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("reported d statistic"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  modnameVar <- tclVar(gettextRcmdr(" "))
  modnameFrame <- tkframe(labelsFrame)
  modnameEntry <- ttkentry(modnameFrame, width="8", textvariable=modnameVar)
  tkgrid(labelRcmdr(modnameFrame, text=gettextRcmdr("n of 1st group"), fg="blue"), sticky="w")
  tkgrid(modnameEntry, sticky="w")
  tkgrid(modnameFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  ylimVar <- tclVar(gettextRcmdr(" "))
  ylimFrame <- tkframe(labelsFrame)
  ylimEntry <- ttkentry(ylimFrame, width="8", textvariable=ylimVar)
  tkgrid(labelRcmdr(ylimFrame, text=gettextRcmdr("n of 2nd group"), fg="blue"), sticky="w")
  tkgrid(ylimEntry, sticky="w")
  tkgrid(ylimFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  varVar <- tclVar(gettextRcmdr(" "))
  varFrame <- tkframe(labelsFrame)
  varEntry <- ttkentry(varFrame, width="8", textvariable=varVar)
  tkgrid(labelRcmdr(varFrame, text=gettextRcmdr("variance of d"), fg="blue"), sticky="w")
  tkgrid(varEntry, sticky="w")
  tkgrid(varFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  onOK <- function(){
    closeDialog()
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="")   # paste(' stat, ') 
    n1 <- trim.blanks(tclvalue(modnameVar))
    n1 <- paste(', ', n1, '', sep="")
    n2 <- trim.blanks(tclvalue(ylimVar))
    n2 <- paste(', ', n2, '', sep="")
    var <- trim.blanks(tclvalue(varVar))
    var <- paste(', ', var, '', sep="")
    doItAndPrint(paste("r_from_d1(", stat, n1, n2, var,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="r_from_d1")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}

# r_from_t

r_from_tcmd <- function(){
  initializeDialog(title=gettextRcmdr("r from t-test"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  statVar <- tclVar(gettextRcmdr(" "))
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="8", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("reported statistic"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  modnameVar <- tclVar(gettextRcmdr(" "))
  modnameFrame <- tkframe(labelsFrame)
  modnameEntry <- ttkentry(modnameFrame, width="25", textvariable=modnameVar)
  tkgrid(labelRcmdr(modnameFrame, text=gettextRcmdr("sample size"), fg="blue"), sticky="w")
  tkgrid(modnameEntry, sticky="w")
  tkgrid(modnameFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  onOK <- function(){
    closeDialog()
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="")  
    n <- trim.blanks(tclvalue(modnameVar))
    n <- paste(', ', n, '', sep="")
    doItAndPrint(paste("r_from_t(", stat, n,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="r_from_t")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(variablesFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=8, columns=2)
}



## DEPENDENCY FUNCTIONS

# ROBUST SE

robustSEcmd <- function(){
  initializeDialog(title=gettextRcmdr("Robust Standard Errors"))
  labelsFrame <- tkframe(top)
  variablesFrame <- tkframe(top)
  .variable <- Variables()
  .numeric <- Numeric()
  statFrame <- tkframe(labelsFrame)
  
  xBox <- variableListBox(variablesFrame, .variable, 
                          title=gettextRcmdr("Cluster variable where the dependencies are present \\ (Be certain active dataset is the same one used when creating the model object)"))
#   yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Effect size (ES) (pick one)"))
#   zBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("ES Variance (pick one)"))
#   aBox <- variableListBox(variablesFrame, .variable, 
#                           title=gettextRcmdr("Grp 1 N variable (pick one)"))
#   bBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Grp 2 N variable (pick one)"))
#   
  UpdateModelNumber()
  statVar <- tclVar(gettextRcmdr("m1"))  
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="5", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("Regression model object name for robust SE:"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text="e.g., m1"), sticky="w")
  
  modelName <- tclVar(paste("robust", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  
  onOK <- function(){
    x <- getSelection(xBox)
#     y <- getSelection(yBox)
#     z <- getSelection(zBox)
#     a <- getSelection(aBox)
#     b <- getSelection(bBox)
    closeDialog()
    if (0 == length(x)) {
      UpdateModelNumber(-1)
      errorCondition(recall=robustSEcmd, message=gettextRcmdr("You must select a cluster variable."))
      return()
    }
#     if (0 == length(a)) {
#       UpdateModelNumber(-1)
#       errorCondition(recall=MetaGcmd, message=gettextRcmdr("No  sample (n.1) variable selected."))
#       return()
#     }
#     if (0 == length(b)) {
#       UpdateModelNumber(-1)
#       errorCondition(recall=MetaGcmd, message=gettextRcmdr("No sample (n.2) variable selected."))
#       return()
#     }
#     if (0 == length(y)) {
#       UpdateModelNumber(-1)
#       errorCondition(recall=MetaGcmd, message=gettextRcmdr("You must select a ES variable."))
#       return()
#     }
#     if (0 == length(z)) {
#       UpdateModelNumber(-1)
#       errorCondition(recall=MetaGcmd, message=gettextRcmdr("No  variance (var) variable selected."))
#       return()
#     }
    
    value <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(value)){
      UpdateModelNumber(-1)
      errorCondition(recall=robustSEcmd, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), value))
      return()
    }
    
    #modelFR <- as.character(tclvalue(modelFRVariable)) 
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="")
    command <- paste("robustSE(model = ",stat, ", cluster =" , ActiveDataSet(), "$", x, ")", sep="")
    #     logger(paste(value, " <- ", command, sep=""))
    #     assign(value, justDoIt(command), envir=.GlobalEnv)
    command <- paste(value, " <- ", command, sep="")
    doItAndPrint(command)
    command2 <- paste(value)
    doItAndPrint(command2) 
  }
  OKCancelHelp(helpSubject="robustSE", model=TRUE)
  
  #tkgrid(modelFRFrame, sticky="w")
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
#   tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(zBox),getFrame(xBox), 
#          getFrame(aBox),getFrame(bBox),sticky="nw")
  tkgrid(getFrame(xBox), labelRcmdr(variablesFrame, text=""), sticky="w")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  tkgrid.configure(helpButton, sticky="e")
  dialogSuffix(rows=7, columns=7)
}





##==== Within-Study Aggregation

# agg function (for Cohen's d or Hedges g)

MetaGcmd <- function(){
  initializeDialog(title=gettextRcmdr("Within-study aggregation for mean difference ESs"))
  labelsFrame <- tkframe(top)
  variablesFrame <- tkframe(top)
  .variable <- Variables()
  .numeric <- Numeric()
  statFrame <- tkframe(labelsFrame)
  
  xBox <- variableListBox(variablesFrame, .variable, 
                          title=gettextRcmdr("id variable (pick one)"))
  yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Effect size (ES) (pick one)"))
  zBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("ES Variance (pick one)"))
  aBox <- variableListBox(variablesFrame, .variable, 
                          title=gettextRcmdr("Grp 1 N variable (pick one)"))
  bBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Grp 2 N variable (pick one)"))
  
  UpdateModelNumber()
  statVar <- tclVar(gettextRcmdr(".50 "))  
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="5", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("estimated correlation btwn outcome measures"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text=" Default is .50  (Wampold, 1997)  "), sticky="w")
  
  modelName <- tclVar(paste("ag.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  
  onOK <- function(){
    x <- getSelection(xBox)
    y <- getSelection(yBox)
    z <- getSelection(zBox)
    a <- getSelection(aBox)
    b <- getSelection(bBox)
    closeDialog()
    if (0 == length(x)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=gettextRcmdr("You must select an id variable."))
      return()
    }
    if (0 == length(a)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=gettextRcmdr("No  sample (n.1) variable selected."))
      return()
    }
    if (0 == length(b)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=gettextRcmdr("No sample (n.2) variable selected."))
      return()
    }
    if (0 == length(y)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=gettextRcmdr("You must select a ES variable."))
      return()
    }
    if (0 == length(z)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=gettextRcmdr("No  variance (var) variable selected."))
      return()
    }
    
    value <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(value)){
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), value))
      return()
    }
    
    #modelFR <- as.character(tclvalue(modelFRVariable)) 
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="")
    command <- paste("agg(" , x, ", ", y,", ", z,", ", a,", ", b,", cor = ",stat," , data=", ActiveDataSet(), ")", sep="")
#     logger(paste(value, " <- ", command, sep=""))
#     assign(value, justDoIt(command), envir=.GlobalEnv)
    command <- paste(value, " <- ", command, sep="")
    doItAndPrint(command)
  }
  OKCancelHelp(helpSubject="agg", model=TRUE)
  
  #tkgrid(modelFRFrame, sticky="w")
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for data.frame:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
  tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(zBox),getFrame(xBox), 
         getFrame(aBox),getFrame(bBox),sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  tkgrid.configure(helpButton, sticky="e")
  dialogSuffix(rows=7, columns=7)
}


# agg function (for correlations)

MetaCcmd <- function(){
  initializeDialog(title=gettextRcmdr("Within-study aggregation for correlation ESs"))
  labelsFrame <- tkframe(top)
  variablesFrame <- tkframe(top)
  .variable <- Variables()
  .numeric <- Numeric()
  statFrame <- tkframe(labelsFrame)
  
  xBox <- variableListBox(variablesFrame, .variable, 
                          title=gettextRcmdr("id variable (pick one)"))
  yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Effect size (ES) (pick one)"))
  #zBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("ES Variance (pick one)"))
  aBox <- variableListBox(variablesFrame, .variable, 
                          title=gettextRcmdr("N variable (pick one)"))
  #bBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Grp 2 N variable (pick one)"))
  
  UpdateModelNumber()
  statVar <- tclVar(gettextRcmdr(".50 "))  
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="5", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("estimated correlation btwn outcome measures"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text=" Default is .50  (Wampold, 1997)  "), sticky="w")
  
  modelName <- tclVar(paste("ag.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  
  onOK <- function(){
    x <- getSelection(xBox)
    y <- getSelection(yBox)
   # z <- getSelection(zBox)
    a <- getSelection(aBox)
   # b <- getSelection(bBox)
    closeDialog()
    if (0 == length(x)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=gettextRcmdr("You must select an id variable."))
      return()
    }
    if (0 == length(a)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=gettextRcmdr("No sample (N) variable selected."))
      return()
    }
    #if (0 == length(b)) {
    #  UpdateModelNumber(-1)
    #  errorCondition(recall=MetaGcmd, message=gettextRcmdr("No sample (n.2) variable selected."))
    #  return()
    #}
    if (0 == length(y)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=gettextRcmdr("You must select a ES variable."))
      return()
    }
    #if (0 == length(z)) {
    #  UpdateModelNumber(-1)
    #  errorCondition(recall=MetaGcmd, message=gettextRcmdr("No  variance (var) variable selected."))
    #  return()
    #}
    
    value <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(value)){
      UpdateModelNumber(-1)
      errorCondition(recall=MetaGcmd, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), value))
      return()
    }
    
    #modelFR <- as.character(tclvalue(modelFRVariable)) 
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="")
    #command <- paste("agg(" , x, ", ", y,", ", z,", ", a,", ", b,", cor = ",stat," , data=", ActiveDataSet(), ")", sep="")
    command <- paste("aggC(" , x, ", ", y,", ", a,", cor = ",stat," , data=", ActiveDataSet(), ")", sep="")
    #     logger(paste(value, " <- ", command, sep=""))
    #     assign(value, justDoIt(command), envir=.GlobalEnv)
    command <- paste(value, " <- ", command, sep="")
    doItAndPrint(command)
  }
  OKCancelHelp(helpSubject="agg", model=TRUE)
  
  #tkgrid(modelFRFrame, sticky="w")
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for data.frame:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
  #tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(zBox),getFrame(xBox), 
        # getFrame(aBox),getFrame(bBox),sticky="nw")
  tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(xBox), 
         getFrame(aBox),sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  tkgrid.configure(helpButton, sticky="e")
  dialogSuffix(rows=7, columns=7)
}









ComplDatacmd <- function(){
  initializeDialog(title=gettextRcmdr("Complete Dataset"))
  variablesFrame <- tkframe(top)
  labelsFrame <- tkframe(top)
  .variable <- Variables()
  xBox <- variableListBox(variablesFrame, .variable, selectmode="multiple",
                          title=gettextRcmdr("variables for complete data (pick one)"))
  UpdateModelNumber()
  statVar <- tclVar(gettextRcmdr(".50 "))  
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="5", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("estimated correlation btwn outcome measures"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text=" Default is .50  (Wampold, 1997)  "), sticky="w")
  modelName <- tclVar(paste("data.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  subsetBox()
  onOK <- function(){
    x <- getSelection(xBox)
    closeDialog()
    if (0 == length(x)) {
      UpdateModelNumber(-1)
      errorCondition(recall=ComplDatacmd , message=gettextRcmdr("No variables selected."))
      return()
    }
    value <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(value)){
      UpdateModelNumber(-1)
      errorCondition(recall=ComplDatacmd, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), value))
      return()
    }
    meta <- ActiveDataSet()
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="") 
    modelN <- as.character(tclvalue(modelNVariable)) 
    command <- paste("ComplData(", meta, ", ", meta, "$", x, ", cor = ",stat," , type= '",modelN,"')", sep="")
#     logger(paste(value, " <- ", command, sep=""))
#     assign(value, justDoIt(command), envir=.GlobalEnv)
#     doItAndPrint(value)
    command <- paste(value, " <- ", command, sep="")
    doItAndPrint(command)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="ComplData", model=TRUE)
  radioButtons(name="modelN", buttons=c("independent", "dependent"), values=c("independent", "dependent"),   
               labels=gettextRcmdr(c("independent", "dependent")), title=gettextRcmdr("Type of aggregation:"))
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(modelFrame, sticky="w")
  tkgrid(labelRcmdr(variablesFrame, text="    "), getFrame(xBox), sticky="nw")
  tkgrid(modelNFrame, sticky="w")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  tkgrid.configure(helpButton, sticky="e")
  dialogSuffix(rows=4, columns=2)
}

##==== Omnibus Analysis ====##


OmnibusEScmd <- function(){
  initializeDialog(title=gettextRcmdr("Omnibus Analysis"))
  variablesFrame <- tkframe(top)
  .variable <- Variables()
  .numeric <- Numeric()
  #xBox <- variableListBox(variablesFrame, .variable, selectmode="multiple",
  #                        title=gettextRcmdr("Moderator (pick one or more)"))
  yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Effect size (ES) (pick one)"))
  zBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("ES Variance (pick one)"))
  
  UpdateModelNumber()
  modelName <- tclVar(paste("omn", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  subsetBox()
  onOK <- function(){
    # x <- getSelection(xBox)
    y <- getSelection(yBox)
    z <- getSelection(zBox)
    closeDialog()
    if (0 == length(y)) {
      UpdateModelNumber(-1)
      errorCondition(recall=OmnibusEScmd, message=gettextRcmdr("You must select a response variable."))
      return()
    }
    #  if (0 == length(x)) {
    #    UpdateModelNumber(-1)
    #    errorCondition(recall=OmnibusEScmd, message=gettextRcmdr("No explanatory variables selected."))
    #    return()
    #}
    if (0 == length(z)) {
      UpdateModelNumber(-1)
      errorCondition(recall=OmnibusEScmd, message=gettextRcmdr("No variance variables selected."))
      return()
    }
    if (is.element(y, z)) {
      UpdateModelNumber(-1)
      errorCondition(recall=OmnibusEScmd, message=gettextRcmdr("Response and variance variables must be         
                                                               different."))
      return()
    }
    subset <- tclvalue(subsetVariable)
    if (trim.blanks(subset) == gettextRcmdr("<all valid cases>") || trim.blanks(subset) == ""){
      subset <- ""
      putRcmdr("modelWithSubset", FALSE)
    }
    else{
      subset <- paste(", subset=", subset, sep="")
      putRcmdr("modelWithSubset", TRUE)
    }
    value <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(value)){
      UpdateModelNumber(-1)
      errorCondition(recall=OmnibusEScmd, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), value))
      return()
    }
    if (is.element(value, listLinearModels())) {
      if ("no" == tclvalue(checkReplace(value, type=gettextRcmdr("Model")))){
        UpdateModelNumber(-1)
        linearRegressionModel()
        return()
      }
    }
    modelFR <- as.character(tclvalue(modelFRVariable)) 
    modelFRb <- as.character(tclvalue(modelFRbVariable))
    meta <-ActiveDataSet()     
    command <- paste("mareg(", y, "~", " 1",
                     ", var=", z,", data=", ActiveDataSet(), subset, ",method='",modelFR,"',knha=",modelFRb,")", sep="")
#     logger(paste(value, " <- ", command, sep=""))
#     assign(value, justDoIt(command), envir=.GlobalEnv)
    command <- paste(value, " <- ", command, sep="")
    doItAndPrint(command)
    command2 <- (paste("summary(", value, ")"))
     doItAndPrint(command2) 
     
  }
  OKCancelHelp(helpSubject="omni", model=TRUE)
  radioButtons(name="modelFR", buttons=c("REML","FE", "HE", "DL", "SJ", "ML",  "EB"), 
               values=c("REML","FE", "HE", "DL", "SJ", "ML",  "EB"),
               labels=gettextRcmdr(c("Restricted maximum-likelihood",
                                     "Fixed effects", 
                                     "Hedges", 
                                     "DerSimonian-Laird", 
                                     "Sidik-Jonkman", 
                                     "Maximum-likelihood",  
                                     "Empirical Bayes")), 
               title=gettextRcmdr("Model"))  
  radioButtons(name="modelFRb", buttons=c("FALSE","TRUE"), 
               values=c("FALSE","TRUE"),
               labels=gettextRcmdr(c("No",
                                     "Yes")), 
               title=gettextRcmdr("Knapp and Hartung adjustment")) 
  tkgrid(modelFRFrame, sticky="w")
    tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
  tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(zBox), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(subsetFrame, sticky="w")
  tkgrid(modelFRbFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  tkgrid.configure(helpButton, sticky="e")
  dialogSuffix(rows=8, columns=3)
  }


##==== Moderator ====#

# macat function


CatModcmd <- function(){
  initializeDialog(title=gettextRcmdr("Categorical moderation (single predictor)"))
  variablesFrame <- tkframe(top)
  .variable <- Variables()
  .numeric <- Numeric()
  .factor <- Factors()
  xBox <- variableListBox(variablesFrame, .factor, 
                          title=gettextRcmdr("Moderator (pick one)"))
  yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Effect size (ES) (pick one)"))
  zBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("ES Variance (pick one)"))
  
  UpdateModelNumber()
  modelName <- tclVar(paste("catmod.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  subsetBox()
  onOK <- function(){
    x <- getSelection(xBox)
    y <- getSelection(yBox)
    z <- getSelection(zBox)
    closeDialog()
    if (0 == length(y)) {
      UpdateModelNumber(-1)
      errorCondition(recall=CatModcmd, message=gettextRcmdr("You must select an effect size variable."))
      return()
    }
    if (0 == length(x)) {
      UpdateModelNumber(-1)
      errorCondition(recall=CatModcmd, message=gettextRcmdr("No moderator variables selected."))
      return()
    }
    if (0 == length(z)) {
      UpdateModelNumber(-1)
      errorCondition(recall=CatModcmd, message=gettextRcmdr("No variance variables selected."))
      return()
    }
    if (is.element(y, x)) {
      UpdateModelNumber(-1)
      errorCondition(recall=CatModcmd, message=gettextRcmdr("Response and explanatory variables must be         
                                                            different."))
      return()
    }
    subset <- tclvalue(subsetVariable)
    if (trim.blanks(subset) == gettextRcmdr("<all valid cases>") || trim.blanks(subset) == ""){
      subset <- ""
      putRcmdr("modelWithSubset", FALSE)
    }
    else{
      subset <- paste(", subset=", subset, sep="")
      putRcmdr("modelWithSubset", TRUE)
    }
    value <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(value)){
      UpdateModelNumber(-1)
      errorCondition(recall=CatModcmd, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), value))
      return()
    }
    if (is.element(value, listLinearModels())) {
      if ("no" == tclvalue(checkReplace(value, type=gettextRcmdr("Model")))){
        UpdateModelNumber(-1)
        linearRegressionModel()
        return()
      }
    }
    
    modelFR <- as.character(tclvalue(modelFRVariable)) 
    meta <-ActiveDataSet()     
    command <- paste("macat(", y,",  var=", z,",  mod=", x,", data=", ActiveDataSet(), ",method='",modelFR,"')", sep="")
#     logger(paste(value, " <- ", command, sep=""))
#     assign(value, doItAndPrint(command), envir=.GlobalEnv)
    command <- paste(value, " <- ", command, sep="")
    doItAndPrint(command)
    command2 <- (paste(value))
    doItAndPrint(command2)
  }
  OKCancelHelp(helpSubject="macat", model=TRUE)
  radioButtons(name="modelFR", buttons=c("Fixed", "Random"), 
               values=c("fixed", "random"),
               labels=gettextRcmdr(c("fixed", "random")), 
               title=gettextRcmdr("model"))   
  tkgrid(modelFRFrame, sticky="w")
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
  tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(zBox), getFrame(xBox), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  tkgrid.configure(helpButton, sticky="e")
  dialogSuffix(rows=5, columns=1)
  }



#CatComp

CatCompcmd <- function(){
  initializeDialog(title=gettextRcmdr("Direct Categorical Moderator Comparison"))
  variablesFrame <- tkframe(top)
  .factor <- Factors()
  .numeric <- Numeric()
  xBox <- variableListBox(variablesFrame, .factor, title=gettextRcmdr("moderator variable (pick one)"))
  yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Effect size (ES) (pick one)"))
  zBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("ES Variance (pick one)"))
  
  UpdateModelNumber()
  modelName <- tclVar(paste("modcompdata.", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  subsetBox()
  onOK <- function(){
    x <- getSelection(xBox)
    y <- getSelection(yBox)
    z <- getSelection(zBox)  
    closeDialog()
    if (0 == length(y)) {
      UpdateModelNumber(-1)
      errorCondition(recall=CatCompcmd, message=gettextRcmdr("You must select a moderator variable."))
      return()
    }
    value <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(value)){
      UpdateModelNumber(-1)
      errorCondition(recall=CatCompcmd, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), value))
      return()
    }
    meta <- ActiveDataSet()
    modelN <- as.character(tclvalue(modelNVariable)) 
    modelC1 <- as.character(tclvalue(modelC1Variable)) 
    modelC2 <- as.character(tclvalue(modelC2Variable))
    command <- paste("macatC(x1=",modelC1,", x2=",modelC2,", g=", y, ", var=", z, ",mod=", paste(x, collapse=","),
                     ",  type= '",modelN,"', data=", ActiveDataSet(),")", sep="")
#     logger(paste(value, " <- ", command, sep=""))
#     assign(value, justDoIt(command), envir=.GlobalEnv)
    command <- paste(value, " <- ", command, sep="")
    doItAndPrint(command)
    command2 <- (paste(value))
    doItAndPrint(command2)
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="CatComp", model=TRUE)
  radioButtons(name="modelN", buttons=c("post.hoc", "planned"), values=c("post.hoc", "planned"),   
               labels=gettextRcmdr(c("post.hoc", "planned")), title=gettextRcmdr("method"))
  radioButtons(name="modelC1", buttons=c("one", "two", "three","four", "five","six"), 
               values=c("1", "2", "3","4", "5","6"),   
               labels=gettextRcmdr(c("one", "two", "three","four", "five","six")),
               title=gettextRcmdr("choose 1st levels of factor to compare"))    
  radioButtons(name="modelC2", buttons=c("one", "two", "three","four", "five","six"), 
               values=c("1", "2", "3","4", "5","6"),   
               labels=gettextRcmdr(c("one", "two", "three","four", "five","six")), 
               title=gettextRcmdr("choose 2nd levels of factor to compare"))    
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
  tkgrid(labelRcmdr(variablesFrame, text="    "), getFrame(xBox), getFrame(yBox),getFrame(zBox),sticky="nw")
  tkgrid(modelNFrame, sticky="w")
  tkgrid(modelC1Frame, sticky="w")
  tkgrid(modelC2Frame, sticky="w")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(subsetFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  #tkgrid(buttonsFrame, stick="w")
  tkgrid.configure(helpButton, sticky="e")
  dialogSuffix(rows=4, columns=2)
}


#mareg  

MAreg2cmd <- function(){
  initializeDialog(title=gettextRcmdr("Meta-Regression"))
  variablesFrame <- tkframe(top)
  .variable <- Variables()
  .numeric <- Numeric()
  xBox <- variableListBox(variablesFrame, .variable, selectmode="multiple",
                          title=gettextRcmdr("Moderator (pick one or more)"))
  yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("Effect size (ES) (pick one)"))
  zBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("ES Variance (pick one)"))
  
  UpdateModelNumber()
  modelName <- tclVar(paste("m", getRcmdr("modelNumber"), sep=""))
  modelFrame <- tkframe(top)
  model <- ttkentry(modelFrame, width="20", textvariable=modelName)
  subsetBox()
  onOK <- function(){
    x <- getSelection(xBox)
    y <- getSelection(yBox)
    z <- getSelection(zBox)
    closeDialog()
    if (0 == length(y)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MAreg2cmd, message=gettextRcmdr("You must select a response variable."))
      return()
    }
    if (0 == length(x)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MAreg2cmd, message=gettextRcmdr("No explanatory variables selected."))
      return()
    }
    if (0 == length(z)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MAreg2cmd, message=gettextRcmdr("No variance variables selected."))
      return()
    }
    if (is.element(y, x)) {
      UpdateModelNumber(-1)
      errorCondition(recall=MAreg2cmd, message=gettextRcmdr("Response and explanatory variables must be         
                                                            different."))
      return()
    }
    subset <- tclvalue(subsetVariable)
    if (trim.blanks(subset) == gettextRcmdr("<all valid cases>") || trim.blanks(subset) == ""){
      subset <- ""
      putRcmdr("modelWithSubset", FALSE)
    }
    else{
      subset <- paste(", subset=", subset, sep="")
      putRcmdr("modelWithSubset", TRUE)
    }
    value <- trim.blanks(tclvalue(modelName))
    if (!is.valid.name(value)){
      UpdateModelNumber(-1)
      errorCondition(recall=MAreg2cmd, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), value))
      return()
    }
    if (is.element(value, listLinearModels())) {
      if ("no" == tclvalue(checkReplace(value, type=gettextRcmdr("Model")))){
        UpdateModelNumber(-1)
        linearRegressionModel()
        return()
      }
    }
    modelFR <- as.character(tclvalue(modelFRVariable))
    modelFRb <- as.character(tclvalue(modelFRbVariable))
    meta <-ActiveDataSet()     
    command <- paste("mareg(", y, "~", paste(x, collapse="+"),
                     ", var=", z,", data=", ActiveDataSet(), subset, ",method='",modelFR,"',knha=",modelFRb,")", sep="")
#     logger(paste(value, " <- ", command, sep=""))
#     assign(value, justDoIt(command), envir=.GlobalEnv)
    command <- paste(value, " <- ", command, sep="")
    doItAndPrint(command)
    command2 <- (paste("summary(", value, ")"))
    doItAndPrint(command2)
    
  }
  OKCancelHelp(helpSubject="mareg", model=TRUE)
  radioButtons(name="modelFR", buttons=c("REML","FE", "HE", "DL", "SJ", "ML",  "EB"), 
               values=c("REML","FE", "HE", "DL", "SJ", "ML",  "EB"),
               labels=gettextRcmdr(c("Restricted maximum-likelihood",
                                     "Fixed effects", 
                                     "Hedges", 
                                     "DerSimonian-Laird", 
                                     "Sidik-Jonkman", 
                                     "Maximum-likelihood",  
                                     "Empirical Bayes")), 
               title=gettextRcmdr("model")) 
  radioButtons(name="modelFRb", buttons=c("FALSE","TRUE"), 
               values=c("FALSE","TRUE"),
               labels=gettextRcmdr(c("No",
                                     "Yes")), 
               title=gettextRcmdr("Knapp and Hartung adjustment"))
  tkgrid(modelFRFrame, sticky="w")
  tkgrid(labelRcmdr(modelFrame, text=gettextRcmdr("Enter name for model:")), model, sticky="w")
  tkgrid(modelFrame, sticky="w")
  tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(zBox),getFrame(xBox), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  tkgrid(subsetFrame, sticky="w")
  tkgrid(modelFRbFrame, sticky="w")
  tkgrid(buttonsFrame, stick="w")
  tkgrid.configure(helpButton, sticky="e")
  dialogSuffix(rows=5, columns=1)
  }


#### DIAGNOSTIC  ####
# radial(temp)
# qqnorm()
## assess for funnel plot asymmetry
# regression test (Egger et al., 1997)
# regtest(m0.del, model="lm")  # indicates NO pub bias
# regtest(temp, model="rma")
# ranktest (Begg & Mazumdar (1994)). another method to assess pub bias
# ranktest(m0.del) 

influencecmd <- function(){
  initializeDialog(title=gettextRcmdr("Detect influential studies in model"))
  labelsFrame <- tkframe(top)
  statVar <- tclVar(gettextRcmdr(" "))  
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="25", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("Type name of model to diagnose"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text=" (Should be name of existing model object, e.g, m1)  "), sticky="w")
  onOK <- function(){
    closeDialog()
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="") 
    doItAndPrint(paste("influence(", stat, " )", sep=""))
    doItAndPrint(paste("plot(influence(", stat, "), plotdfb=TRUE)", sep=""))
  }
  OKCancelHelp(helpSubject="influence", model=TRUE)
  
  tkgrid(labelsFrame, sticky="w") 
  tkgrid.configure(helpButton, sticky="e")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=5, columns=3)
}

rstudentcmd <- function(){
  initializeDialog(title=gettextRcmdr("Extract externally standardized residuals"))
  labelsFrame <- tkframe(top)
  statVar <- tclVar(gettextRcmdr(" "))  
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="25", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("Type name of model to diagnose"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text=" (Should be name of existing model object, e.g, m1)  "), sticky="w")
  onOK <- function(){
    closeDialog()
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="") 
    doItAndPrint(paste("rstudent(", stat, " )", sep=""))
  }
  OKCancelHelp(helpSubject="rstudent")
  
  tkgrid(labelsFrame, sticky="w") 
  tkgrid.configure(helpButton, sticky="e")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=5, columns=3)
}

rstandardcmd <- function(){
  initializeDialog(title=gettextRcmdr("Extract internally standardized residuals (less conservative thn external)"))
  labelsFrame <- tkframe(top)
  statVar <- tclVar(gettextRcmdr(" "))  
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="25", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("Type name of model to diagnose"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text=" (Should be name of existing model object, e.g, m1)  "), sticky="w")
  onOK <- function(){
    closeDialog()
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="") 
    doItAndPrint(paste("rstandard(", stat, " )", sep=""))
  }
  OKCancelHelp(helpSubject="rstandard")
  
  tkgrid(labelsFrame, sticky="w") 
  tkgrid.configure(helpButton, sticky="e")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=5, columns=3)
}

trimcmd <- function(){
  initializeDialog(title=gettextRcmdr("Trim and fill (assess publication bias)"))
  labelsFrame <- tkframe(top)
  statVar <- tclVar(gettextRcmdr(" "))  
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="25", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("Type name of model to diagnose"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text=" (Should be name of existing model object, e.g, omn1)  "), sticky="w")
  onOK <- function(){
    closeDialog()
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="") 
    doItAndPrint(paste("trimfill(", stat, " )", sep=""))
    doItAndPrint(paste("funnel(trimfill(", stat, " ))", sep=""))
  }
  OKCancelHelp(helpSubject="trimfill")
  
  tkgrid(labelsFrame, sticky="w") 
  tkgrid.configure(helpButton, sticky="e")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=5, columns=3)
}

regcmd <- function(){
  initializeDialog(title=gettextRcmdr("Regression test (funnel plot asymmetry)"))
  labelsFrame <- tkframe(top)
  statVar <- tclVar(gettextRcmdr(" "))  
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="25", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("Type name of model to diagnose"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text=" (Should be name of existing model object, e.g, omn1)  "), sticky="w")
  onOK <- function(){
    closeDialog()
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="") 
    doItAndPrint(paste("regtest(", stat, " )", sep=""))
  }
  OKCancelHelp(helpSubject="regtest")
  
  tkgrid(labelsFrame, sticky="w") 
  tkgrid.configure(helpButton, sticky="e")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=5, columns=3)
}

rankcmd <- function(){
  initializeDialog(title=gettextRcmdr("Rank correlation test (assess publication bias)"))
  labelsFrame <- tkframe(top)
  statVar <- tclVar(gettextRcmdr(" "))  
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="25", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("Type name of model to diagnose"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text=" (Should be name of existing model object, e.g, omn1)  "), sticky="w")
  onOK <- function(){
    closeDialog()
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="") 
    doItAndPrint(paste("ranktest(", stat, " )", sep=""))
  }
  OKCancelHelp(helpSubject="ranktest")
  
  tkgrid(labelsFrame, sticky="w") 
  tkgrid.configure(helpButton, sticky="e")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=5, columns=3)
}

qqnormcmd <- function(){
  initializeDialog(title=gettextRcmdr("Normal Q-Q plot"))
  labelsFrame <- tkframe(top)
  statVar <- tclVar(gettextRcmdr(" "))  
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="25", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("Type name of model to plot"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text=" (Should be name of existing model object, e.g, omn1)  "), sticky="w")
  onOK <- function(){
    closeDialog()
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="") 
    doItAndPrint(paste("qqnorm(", stat, " )", sep=""))
  }
  OKCancelHelp(helpSubject="qqnorm")
  
  tkgrid(labelsFrame, sticky="w") 
  tkgrid.configure(helpButton, sticky="e")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=5, columns=3)
}

radialcmd <- function(){
  initializeDialog(title=gettextRcmdr("Radial plot"))
  labelsFrame <- tkframe(top)
  statVar <- tclVar(gettextRcmdr(" "))  
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="25", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("Type name of model to plot"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text=" (Should be name of existing model object, e.g, omn1)  "), sticky="w")
  onOK <- function(){
    closeDialog()
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="") 
    doItAndPrint(paste("radial(", stat, " )", sep=""))
  }
  OKCancelHelp(helpSubject="radial")
  
  tkgrid(labelsFrame, sticky="w") 
  tkgrid.configure(helpButton, sticky="e")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=5, columns=3)
}

##==== Graphics ====##

#Forest Plot

ForestPlotcmd <- function(){
  initializeDialog(title=gettextRcmdr("Forest Plot"))
  labelsFrame <- tkframe(top)
  statVar <- tclVar(gettextRcmdr(" "))  
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="25", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("Type name of model to plot"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text="(Should be name of existing model object, e.g, omn1)"), sticky="w")
  onOK <- function(){
    closeDialog()
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="") 
    doItAndPrint(paste("forest(", stat, " )", sep=""))
  }
  OKCancelHelp(helpSubject="forest", model=TRUE)
  
  tkgrid(labelsFrame, sticky="w") 
  tkgrid.configure(helpButton, sticky="e")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=5, columns=3)
}


#ForestPlotcmd <- function(){
#  initializeDialog(title=gettextRcmdr("Forest Plot"))
#  variablesFrame <- tkframe(top)
#  UpdateModelNumber()
#  modelName <- tclVar(paste("ForestPlotGraph.", getRcmdr("modelNumber"), sep=""))
#  labelsFrame <- tkframe(top)
#  titleVar <- tclVar(gettextRcmdr("<auto>"))
#  #ylabVar <- tclVar(gettextRcmdr("<auto>"))
#  titleFrame <- tkframe(labelsFrame)
#  titleEntry <- ttkentry(titleFrame, width="25", textvariable=titleVar)
#  titleScroll <- ttkscrollbar(titleFrame, orient="horizontal",
#  command=function(...) tkxview(titleEntry, ...))
#  tkconfigure(titleEntry, xscrollcommand=function(...) tkset(titleScroll, ...))
#  tkgrid(labelRcmdr(titleFrame, text=gettextRcmdr("title"), fg="blue"), sticky="w")
#  tkgrid(titleEntry, sticky="w")
#  tkgrid(titleScroll, sticky="ew")
#  tkgrid(titleFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
#  onOK <- function(){
#    closeDialog()
#    .activeDataSet <- ActiveDataSet()
#    title <- trim.blanks(tclvalue(titleVar))
#    title <- if(title == gettextRcmdr("<auto>")) "" else paste(', title="', title, '"', sep="")
#    modelFR <- as.character(tclvalue(modelFRVariable)) 
#    doItAndPrint(paste("ForestPlot(", .activeDataSet, ",  method='" ,modelFR, "'",title,")", sep="")) 
#    activateMenus()
#    tkfocus(CommanderWindow())
#  }
#  OKCancelHelp(helpSubject="ForestPlot")
#  radioButtons(name="modelFR", buttons=c("Fixed", "Random"), 
#               values=c("fixed", "random"),
#               labels=gettextRcmdr(c("fixed", "random")), 
#               title=gettextRcmdr("Model"))
#  tkgrid(modelFRFrame, sticky="w")
#  tkgrid(labelsFrame, sticky="w")
#  tkgrid(labelRcmdr(top, text=" "))
#  tkgrid(buttonsFrame, columnspan=2, sticky="w")
#  dialogSuffix(rows=8, columns=2)
#}


#Funnel Plot

FunnelPlotcmd <- function(){
  initializeDialog(title=gettextRcmdr("Funnel Plot"))
  labelsFrame <- tkframe(top)
  statVar <- tclVar(gettextRcmdr(" "))  
  statFrame <- tkframe(labelsFrame)
  statEntry <- ttkentry(statFrame, width="25", textvariable=statVar)
  tkgrid(labelRcmdr(statFrame, text=gettextRcmdr("Type name of model to plot"), fg="blue"), sticky="w")
  tkgrid(statEntry, sticky="w")
  tkgrid(statFrame, labelRcmdr(labelsFrame, text=" (Should be name of existing model object, e.g, omn1)  "), sticky="w")
  onOK <- function(){
    closeDialog()
    stat <- trim.blanks(tclvalue(statVar))
    stat <- paste (' ', stat, '', sep="") 
    doItAndPrint(paste("funnel(", stat, " )", sep=""))
  }
  OKCancelHelp(helpSubject="funnel", model=TRUE)
  
  tkgrid(labelsFrame, sticky="w") 
  tkgrid.configure(helpButton, sticky="e")
  tkgrid(buttonsFrame, stick="w")
  dialogSuffix(rows=5, columns=3)
}


#FunnelPlotcmd <- function(){
#  initializeDialog(title=gettextRcmdr("Funnel Plot"))
#  variablesFrame <- tkframe(top)
#  UpdateModelNumber()
#  modelName <- tclVar(paste("FunnelGraph.", getRcmdr("modelNumber"), sep=""))
#  labelsFrame <- tkframe(top)
#  titleVar <- tclVar(gettextRcmdr("<auto>"))
#  titleFrame <- tkframe(labelsFrame)
#  titleEntry <- ttkentry(titleFrame, width="25", textvariable=titleVar)
#  titleScroll <- ttkscrollbar(titleFrame, orient="horizontal",
#  command=function(...) tkxview(titleEntry, ...))
#  tkconfigure(titleEntry, xscrollcommand=function(...) tkset(titleScroll, ...))
#  tkgrid(labelRcmdr(titleFrame, text=gettextRcmdr("title"), fg="blue"), sticky="w")
#  tkgrid(titleEntry, sticky="w")
#  tkgrid(titleScroll, sticky="ew")
#  tkgrid(titleFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
#  onOK <- function(){
#    closeDialog()
#    .activeDataSet <- ActiveDataSet()
#    title <- trim.blanks(tclvalue(titleVar))
#    title <- if(title == gettextRcmdr("<auto>")) "" else paste(', title="', title, '"', sep="")
#    activateMenus()
#    tkfocus(CommanderWindow())
#  }
#  OKCancelHelp(helpSubject="FunnelPlot")
#  radioButtons(name="modelFR", buttons=c("Fixed", "Random"), 
#               values=c("fixed", "random"),
#               labels=gettextRcmdr(c("fixed", "random")), 
#               title=gettextRcmdr("Model"))
#  tkgrid(modelFRFrame, sticky="w")
#  tkgrid(labelsFrame, sticky="w")
#  tkgrid(labelRcmdr(top, text=" "))
#  tkgrid(buttonsFrame, columnspan=2, sticky="w")
#  dialogSuffix(rows=8, columns=2)
#}


#MultiModGraph

MultiModGraphcmd <- function(){
  initializeDialog(title=gettextRcmdr("Multi-Predictor Moderator Graph"))
  variablesFrame <- tkframe(top)
  UpdateModelNumber()
  modelName <- tclVar(paste("MultiModGraph.", getRcmdr("modelNumber"), sep=""))
  .factor <- Factors()
  variablesFrame <- tkframe(top)
  xBox <- variableListBox(variablesFrame, .factor, title=gettextRcmdr("categorical moderator (pick one)"))
  .numeric <- Numeric()
  yBox <- variableListBox(variablesFrame, .numeric, title=gettextRcmdr("continuous moderator (pick one)"))
  labelsFrame <- tkframe(top)
  titleVar <- tclVar(gettextRcmdr("<auto>"))
  conmod.nameVar <- tclVar(gettextRcmdr("<auto>"))
  titleFrame <- tkframe(labelsFrame)
  titleEntry <- ttkentry(titleFrame, width="25", textvariable=titleVar)
  titleScroll <- ttkscrollbar(titleFrame, orient="horizontal",
                              command=function(...) tkxview(titleEntry, ...))
  tkconfigure(titleEntry, xscrollcommand=function(...) tkset(titleScroll, ...))
  tkgrid(labelRcmdr(titleFrame, text=gettextRcmdr("title"), fg="blue"), sticky="w")
  tkgrid(titleEntry, sticky="w")
  tkgrid(titleScroll, sticky="ew")
  conmod.nameFrame <- tkframe(labelsFrame)
  conmod.nameEntry <- ttkentry(conmod.nameFrame, width="25", textvariable=conmod.nameVar)
  conmod.nameScroll <- ttkscrollbar(conmod.nameFrame, orient="horizontal",
                                    command=function(...) tkxview(conmod.nameEntry, ...))
  tkconfigure(conmod.nameEntry, xscrollcommand=function(...) tkset(conmod.nameScroll, ...))
  tkgrid(labelRcmdr(conmod.nameFrame, text=gettextRcmdr("continuous moderator label"), fg="blue"), sticky="w")
  tkgrid(conmod.nameEntry, sticky="w")
  tkgrid(conmod.nameScroll, sticky="ew")
  tkgrid(titleFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  tkgrid(conmod.nameFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  onOK <- function(){
    x <- getSelection(xBox)
    conx <- getSelection(yBox)
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    title <- trim.blanks(tclvalue(titleVar))
    title <- if(title == gettextRcmdr("<auto>")) "" else paste(', title="', title, '"', sep="")
    conmod.name <- trim.blanks(tclvalue(conmod.nameVar))
    conmod.name <- if(conmod.name == gettextRcmdr("<auto>")) "" else paste(', conmod.name="', conmod.name, '"', sep="")
    modelFR <- as.character(tclvalue(modelFRVariable)) 
    doItAndPrint(paste("MultiModGraph(", .activeDataSet, ", ",.activeDataSet, "$", conx, ",  ",.activeDataSet, 
                       "$", x, ", method='" ,modelFR, "'",title, conmod.name,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="MultiModGraph")
  radioButtons(name="modelFR", buttons=c("Fixed", "Random"), 
               values=c("fixed", "random"),
               labels=gettextRcmdr(c("fixed", "random")), 
               title=gettextRcmdr("Model"))
  tkgrid(modelFRFrame, sticky="w")
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(xBox), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  dialogSuffix(rows=8, columns=2)
}

#MAregGraph

MAregGraphcmd <- function(){
  initializeDialog(title=gettextRcmdr("Meta-Regression Graph"))
  variablesFrame <- tkframe(top)
  UpdateModelNumber()
  modelName <- tclVar(paste("MAregGraph.", getRcmdr("modelNumber"), sep=""))
  .variables <- Variables()
  variablesFrame <- tkframe(top)
  xBox <- variableListBox(variablesFrame, .variables, title=gettextRcmdr("moderator variable (pick one)"))
  labelsFrame <- tkframe(top)
  titleVar <- tclVar(gettextRcmdr("<auto>"))
  titleFrame <- tkframe(labelsFrame)
  titleEntry <- ttkentry(titleFrame, width="25", textvariable=titleVar)
  titleScroll <- ttkscrollbar(titleFrame, orient="horizontal",
                              command=function(...) tkxview(titleEntry, ...))
  tkconfigure(titleEntry, xscrollcommand=function(...) tkset(titleScroll, ...))
  tkgrid(labelRcmdr(titleFrame, text=gettextRcmdr("title"), fg="blue"), sticky="w")
  tkgrid(titleEntry, sticky="w")
  tkgrid(titleScroll, sticky="ew")
  tkgrid(titleFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  modnameVar <- tclVar(gettextRcmdr("<auto>"))
  modnameFrame <- tkframe(labelsFrame)
  modnameEntry <- ttkentry(modnameFrame, width="25", textvariable=modnameVar)
  modnameScroll <- ttkscrollbar(modnameFrame, orient="horizontal",
                                command=function(...) tkxview(modnameEntry, ...))
  tkconfigure(modnameEntry, xscrollcommand=function(...) tkset(modnameScroll, ...))
  tkgrid(labelRcmdr(modnameFrame, text=gettextRcmdr("moderator label on plot"), fg="blue"), sticky="w")
  tkgrid(modnameEntry, sticky="w")
  tkgrid(modnameScroll, sticky="ew")
  tkgrid(modnameFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  ylimVar <- tclVar(gettextRcmdr("c(0,1)"))
  ylimFrame <- tkframe(labelsFrame)
  ylimEntry <- ttkentry(ylimFrame, width="25", textvariable=ylimVar)
  ylimScroll <- ttkscrollbar(ylimFrame, orient="horizontal",
                             command=function(...) tkxview(ylimEntry, ...))
  tkconfigure(ylimEntry, xscrollcommand=function(...) tkset(ylimScroll, ...))
  tkgrid(labelRcmdr(ylimFrame, text=gettextRcmdr("y-axis limits"), fg="blue"), sticky="w")
  tkgrid(ylimEntry, sticky="w")
  tkgrid(ylimScroll, sticky="ew")
  tkgrid(ylimFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  onOK <- function(){
    x <- getSelection(xBox)
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    title <- trim.blanks(tclvalue(titleVar))
    title <- if(title == gettextRcmdr("<auto>")) "" else paste(', title="', title, '"', sep="")
    modname <- trim.blanks(tclvalue(modnameVar))
    modname <- if(modname == gettextRcmdr("<auto>")) "" else paste(', modname="', modname, '"', sep="")
    ylim <- trim.blanks(tclvalue(ylimVar))
    ylim <- if(modname == gettextRcmdr("<auto>")) "" else paste(', ylim=', ylim, '', sep="")                               ### will use this format for the comp r section
    modelFR <- as.character(tclvalue(modelFRVariable)) 
    doItAndPrint(paste("MAregGraph(", .activeDataSet, ",  ",.activeDataSet, "$", x, ", method='" ,modelFR, "'",title, modname, ylim,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="MAregGraph")
  radioButtons(name="modelFR", buttons=c("Fixed", "Random"), 
               values=c("fixed", "random"),
               labels=gettextRcmdr(c("fixed", "random")), 
               title=gettextRcmdr("Model"))
  
  tkgrid(modelFRFrame, sticky="w")
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  tkgrid(labelRcmdr(variablesFrame, text="    "), getFrame(xBox), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  dialogSuffix(rows=8, columns=2)
}

#CatModGraph

CatModGraphcmd <- function(){
  initializeDialog(title=gettextRcmdr("Categorical Moderator Graph"))
  variablesFrame <- tkframe(top)
  UpdateModelNumber()
  modelName <- tclVar(paste("CatGraph.", getRcmdr("modelNumber"), sep=""))
  .factor <- Factors()
  variablesFrame <- tkframe(top)
  xBox <- variableListBox(variablesFrame, .factor, title=gettextRcmdr("moderator variable (pick one)"))
  labelsFrame <- tkframe(top)
  titleVar <- tclVar(gettextRcmdr("<auto>"))
  modnameVar <- tclVar(gettextRcmdr("<auto>"))
  titleFrame <- tkframe(labelsFrame)
  titleEntry <- ttkentry(titleFrame, width="25", textvariable=titleVar)
  titleScroll <- ttkscrollbar(titleFrame, orient="horizontal",
                              command=function(...) tkxview(titleEntry, ...))
  tkconfigure(titleEntry, xscrollcommand=function(...) tkset(titleScroll, ...))
  tkgrid(labelRcmdr(titleFrame, text=gettextRcmdr("title"), fg="blue"), sticky="w")
  tkgrid(titleEntry, sticky="w")
  tkgrid(titleScroll, sticky="ew")
  modnameFrame <- tkframe(labelsFrame)
  modnameEntry <- ttkentry(modnameFrame, width="25", textvariable=modnameVar)
  modnameScroll <- ttkscrollbar(modnameFrame, orient="horizontal",
                                command=function(...) tkxview(modnameEntry, ...))
  tkconfigure(modnameEntry, xscrollcommand=function(...) tkset(modnameScroll, ...))
  tkgrid(labelRcmdr(modnameFrame, text=gettextRcmdr("moderator label on plot"), fg="blue"), sticky="w")
  tkgrid(modnameEntry, sticky="w")
  tkgrid(modnameScroll, sticky="ew")
  tkgrid(titleFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  tkgrid(modnameFrame, labelRcmdr(labelsFrame, text="     "), sticky="w")
  onOK <- function(){
    x <- getSelection(xBox)
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    title <- trim.blanks(tclvalue(titleVar))
    title <- if(title == gettextRcmdr("<auto>")) "" else paste(', title="', title, '"', sep="")
    modname <- trim.blanks(tclvalue(modnameVar))
    modname <- if(modname == gettextRcmdr("<auto>")) "" else paste(', modname="', modname, '"', sep="")
    modelFR <- as.character(tclvalue(modelFRVariable)) 
    doItAndPrint(paste("CatModGraph(", .activeDataSet, ",  ",.activeDataSet, "$", x, ", method='" ,modelFR, "'",title, modname,")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="CatModGraph")
  radioButtons(name="modelFR", buttons=c("Fixed", "Random"), 
               values=c("fixed", "random"),
               labels=gettextRcmdr(c("fixed", "random")), 
               title=gettextRcmdr("Model"))
  
  tkgrid(modelFRFrame, sticky="w")
  tkgrid(getFrame(xBox), sticky="nw")
  tkgrid(labelsFrame, sticky="w")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  tkgrid(labelRcmdr(variablesFrame, text="    "), getFrame(xBox), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  dialogSuffix(rows=8, columns=2)
}

##==== Other ====##

# CorAtten

CorAttencmd <- function(){
  initializeDialog(title=gettextRcmdr("Correction for Attenuation"))
  variablesFrame <- tkframe(top)
  UpdateModelNumber()
  modelName <- tclVar(paste("CorAtt.", getRcmdr("modelNumber"), sep=""))
  .variables <-Variables()
  variablesFrame <- tkframe(top)
  xBox <- variableListBox(variablesFrame, .variables, title=gettextRcmdr("reliability of predictor variable"))
  yBox <- variableListBox(variablesFrame, .variables, title=gettextRcmdr("reliability of outcome variable"))
  labelsFrame <- tkframe(top)
  onOK <- function(){
    xx <- getSelection(xBox)
    yy <- getSelection(yBox)
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    doItAndPrint(paste("atten(", .activeDataSet, ", ",.activeDataSet, "$", xx, ",  ",.activeDataSet, 
                       "$", yy, ")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="CorAtten")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(xBox), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  dialogSuffix(rows=8, columns=2)
}

# Kappa

Kappacmd <- function(){
  initializeDialog(title=gettextRcmdr("Reliability: Kappa"))
  variablesFrame <- tkframe(top)
  UpdateModelNumber()
  modelName <- tclVar(paste("Kappa.", getRcmdr("modelNumber"), sep=""))
  .variables <-Variables()
  variablesFrame <- tkframe(top)
  xBox <- variableListBox(variablesFrame, .variables, title=gettextRcmdr("first rater of categorical variable"))
  yBox <- variableListBox(variablesFrame, .variables, title=gettextRcmdr("second rater of categorical variable"))
  labelsFrame <- tkframe(top)
  onOK <- function(){
    one <- getSelection(xBox)
    two <- getSelection(yBox)
    closeDialog()
    .activeDataSet <- ActiveDataSet()
    doItAndPrint(paste("Kappa(",.activeDataSet, "$", one, ",  ",.activeDataSet, 
                       "$", two, ")", sep="")) 
    activateMenus()
    tkfocus(CommanderWindow())
  }
  OKCancelHelp(helpSubject="Kappa")
  tkgrid(labelRcmdr(top, text=" "))
  tkgrid(buttonsFrame, columnspan=2, sticky="w")
  tkgrid(getFrame(yBox), labelRcmdr(variablesFrame, text="    "), getFrame(xBox), sticky="nw")
  tkgrid(variablesFrame, sticky="w")
  dialogSuffix(rows=8, columns=2)
}


