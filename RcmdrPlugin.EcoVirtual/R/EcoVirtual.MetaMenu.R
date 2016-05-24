
########### Meta Menu #################
intCol <-function () 
{
initializeDialog(title = gettextRcmdr("Internal Colonization"))
#### Salva dados
dsname <- tclVar("Do_Not_Save")
entryDsname <- tkentry(top, width="20", textvariable=dsname)
####
tfVar <- tclVar("100")
tfEntry <- tkentry(top, width = "4", textvariable = tfVar)
clVar <- tclVar("20")
lnVar <- tclVar("20")
clEntry <- tkentry(top, width = "4", textvariable = clVar)
lnEntry <- tkentry(top, width = "4", textvariable = lnVar)
### controle de barra
fiVar <- tclVar(0.4)
fiVarSlider <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=fiVar, resolution=0.01, orient="horizontal")
#####
iVar <- tclVar("0.1") ## nclassVar ->iVar
iEntry <- tkentry(top, width = "6", textvariable = iVar)
peVar <- tclVar(0.05) 
peEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=peVar, resolution=0.01, orient="horizontal")
	onOK <- function() 
	{
	#closeDialog()
   tf <- round(as.numeric(tclvalue(tfVar)))
        if (is.na(tf) || tf <= 0) 
        {
        errorCondition(message = "Number of simulations must be a positive integer")
        return()
        }
  	cl <- round(as.numeric(tclvalue(clVar)))
        if (is.na(cl) || cl <= 0) 
        {
        errorCondition(message = "Number of columns on the simulated arena must be a positive integer.")
        return()
        }
	ln <- round(as.numeric(tclvalue(lnVar)))
        if (is.na(ln) || ln <= 0) 
        {
        errorCondition("Number of lines on the simulated arena must be a positive integer.")
        return()
        }
   i <- as.numeric(tclvalue(iVar))
        if (i<0 || i > 10) 
        {
        errorCondition(message = "Colonization constant must be between 0 and 10")
        return()
        }
   fi <- as.numeric(tclvalue(fiVar))
   pe <- as.numeric(tclvalue(peVar))
############ Data name
   dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue=="") 
        {
        	command <- paste("metaCi(tmax = ",tf, ", cl = ", cl,", fi = ", fi,", rw =", ln,", ci = ", i,", pe = ", pe, ")", sep = "")
        }
        else  
		  {
		  command <- paste(dsnameValue,"<- metaCi(tmax = ",tf, ", rw = ", cl,", fi = ", fi,", ln =", ln,", ci = ", i,", pe = ", pe, ")", sep = "")
		  }
########
	doItAndPrint(command)
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "metapopulation")
# data name
tkgrid(tklabel(top, text="Enter name for data set: "), entryDsname, sticky="e")
tkgrid(tklabel(top, text="Simulation Arena Conditions : ", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Maximum time  "), tfEntry, sticky = "e")
tkgrid(tklabel(top, text = "Columns  "), clEntry, sticky = "e")
tkgrid(tklabel(top, text = "Rows  "), lnEntry, sticky = "e")
#
tkgrid(tklabel(top, text="Species parameters :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text="Initial occupancy  "), fiVarSlider, sticky="se")
tkgrid(tklabel(top, text = "Colonization coefficient  "), iEntry, sticky = "e")
tkgrid(tklabel(top, text = "Prob. Extinction  "), peEntry, sticky = "se")
tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
tkgrid.configure(entryDsname, sticky = "w")
tkgrid.configure(tfEntry, sticky = "w")
tkgrid.configure(clEntry, sticky = "w")
tkgrid.configure(lnEntry, sticky = "w")
tkgrid.configure(fiVarSlider, sticky = "w")
tkgrid.configure(iEntry, sticky = "w")
tkgrid.configure(peEntry, sticky = "w")
dialogSuffix(rows = 8, columns = 2, focus = tfEntry)
}


######################################
##########################################
propRain <-function() 
{
initializeDialog(title = gettextRcmdr("Propagulus Rain"))
####
dsname <- tclVar("Do_Not_Save")
entryDsname <- tkentry(top, width="20", textvariable=dsname)
###
tfVar <- tclVar("100")
tfEntry <- tkentry(top, width = "4", textvariable = tfVar)
clVar <- tclVar("20")
lnVar <- tclVar("20")
clEntry <- tkentry(top, width = "4", textvariable = clVar)
lnEntry <- tkentry(top, width = "4", textvariable = lnVar)
fiVar <- tclVar(0.25) ## nclassVar ->fiVar
fiEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=fiVar, resolution=0.01, orient="horizontal")
pcVar <- tclVar("0.1") ## nclassVar ->fiVar
pcEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=pcVar, resolution=0.01, orient="horizontal")
peVar <- tclVar("0.05") ## nclassVar ->fiVar
peEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=peVar, resolution=0.01, orient="horizontal")
	onOK <- function() 
	{
        #closeDialog()
        tf <- round(as.numeric(tclvalue(tfVar)))
        if (is.na(tf) || tf <= 0) 
        {
            errorCondition(message = "Number of simulations must be a positive integer")
            return()
        }
        cl <- round(as.numeric(tclvalue(clVar)))
        if (is.na(cl) || cl <= 0) 
        {
            errorCondition(message = "Number of columns on the simulated arena must be a positive integer.")
            return()
        }
        ln <- round(as.numeric(tclvalue(lnVar)))
        if (is.na(ln) || ln <= 0) 
        {
            errorCondition("Number of lines on the simulated arena must be a positive integer.")
            return()
        }
        fi <- as.numeric(tclvalue(fiVar))
        pc <- as.numeric(tclvalue(pcVar))
        pe <- as.numeric(tclvalue(peVar))
############ Data name
   dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue=="") 
        {
        	command <- paste("metaPop(tmax = ",tf, ", cl = ", cl,", fi = ", fi,", rw =", ln,", pc = ", pc,", pe = ", pe, ")", sep = "")
        }
        else  
		  {
		  command <- paste(dsnameValue,"<-metaPop(tmax = ",tf, ", cl = ", cl,", fi = ", fi,", rw =", ln,", pc = ", pc,", pe = ", pe, ")", sep = "")
		  }
########
#	command <- paste("metaPop(tf = ",tf, ", cl = ", cl,", fi = ", fi,", ln =", ln,", pc = ", pc,", pe = ", pe, ")", sep = "")
	doItAndPrint(command)
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "metapopulation")
tkgrid(tklabel(top, text="Enter name for data set: "), entryDsname, sticky="e")
tkgrid(tklabel(top, text="Simulation Arena Conditions : ", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Maximum time  "), tfEntry, sticky = "e")
tkgrid(tklabel(top, text = "Columns  "), clEntry, sticky = "e")
tkgrid(tklabel(top, text = "Rows  "), lnEntry, sticky = "e")
tkgrid(tklabel(top, text="Species parameters : ", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Initial occupancy  "), fiEntry, sticky = "se")
tkgrid(tklabel(top, text = "Colonization probability "), pcEntry, sticky = "se")
tkgrid(tklabel(top, text = "Extinction probability "), peEntry, sticky = "se")
tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
tkgrid.configure(entryDsname, sticky = "w")
tkgrid.configure(tfEntry, sticky = "w")
tkgrid.configure(clEntry, sticky = "w")
tkgrid.configure(lnEntry, sticky = "w")
tkgrid.configure(fiEntry, sticky = "w")
tkgrid.configure(pcEntry, sticky = "w")
tkgrid.configure(peEntry, sticky = "w")
dialogSuffix(rows = 7, columns = 2, focus = tfEntry)
}
###############################################################
##################################################
#metaCiEr <-function(tf,cl,ln,fi,i,e)
#metaEr <-function(tf,cl,ln,fi,pc,e)
resEff <-function () 
{
initializeDialog(title = gettextRcmdr("Rescue Effect"))
####
dsname <- tclVar("Do_Not_Save")
entryDsname <- tkentry(top, width="20", textvariable=dsname)
###
tfVar <- tclVar("100")
tfEntry <- tkentry(top, width = "4", textvariable = tfVar)
clVar <- tclVar("20")
lnVar <- tclVar("20")
clEntry <- tkentry(top, width = "4", textvariable = clVar)
lnEntry <- tkentry(top, width = "4", textvariable = lnVar)
fiVar <- tclVar(0.25) ## nclassVar ->fiVar
fiEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=fiVar, resolution=0.01, orient="horizontal")
pcVar <- tclVar("0.1") ## nclassVar ->fiVar
pcEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=pcVar, resolution=0.01, orient="horizontal")
eVar <- tclVar("0.05") ## nclassVar ->fiVar
eEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable = eVar, resolution=0.01, orient="horizontal")
	onOK <- function() 
	{
        #closeDialog()
        tf <- round(as.numeric(tclvalue(tfVar)))
        if (is.na(tf) || tf <= 0) 
        {
            errorCondition(message = "Number of simulations must be a positive integer")
            return()
        }
        cl <- round(as.numeric(tclvalue(clVar)))
        if (is.na(cl) || cl <= 0) 
        {
            errorCondition(message = "Number of columns on the simulated arena must be a positive integer.")
            return()
        }
        ln <- round(as.numeric(tclvalue(lnVar)))
        if (is.na(ln) || ln <= 0) 
        {
            errorCondition("Number of lines on the simulated arena must be a positive integer.")
            return()
        }
        e <- as.numeric(tclvalue(eVar))
        if (e<0 || e > 1) 
        {
            errorCondition(message = "Extinction constant must be between 0 and 1")
            return()
        }
        fi <- as.numeric(tclvalue(fiVar))
        pc <- as.numeric(tclvalue(pcVar))
############ Data name
   dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue=="") 
        {
        	command <- paste("metaEr(tmax = ",tf, ", cl = ", cl,", fi = ", fi,", rw =", ln,", pc = ", pc,", ce = ", e, ")", sep = "")
        }
        else  
		  {
		  command <- paste(dsnameValue,"<-metaEr(tmax = ",tf, ", cl = ", cl,", fi = ", fi,", rw =", ln,", pc = ", pc,", ce = ", e, ")", sep = "")
		  }
########   
#	command <- paste("metaEr(tf = ",tf, ", cl = ", cl,", fi = ", fi,", ln =", ln,", pc = ", pc,", e = ", e, ")", sep = "")
	doItAndPrint(command)
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "metapopulation")
tkgrid(tklabel(top, text="Enter data set name: "), entryDsname, sticky="e")
tkgrid(tklabel(top, text="Simulation Arena Conditions : ", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Maximum time  "), tfEntry, sticky = "e")
tkgrid(tklabel(top, text = "Columns  "), clEntry, sticky = "e")
tkgrid(tklabel(top, text = "Rows  "), lnEntry, sticky = "e")
tkgrid(tklabel(top, text="Species parameters : ", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Initial occupancy "), fiEntry, sticky = "e")
tkgrid(tklabel(top, text = "Colonization probability  "), pcEntry, sticky = "e")
tkgrid(tklabel(top, text = "Extinction coefficient  "), eEntry, sticky = "e")
tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
tkgrid.configure(entryDsname, sticky = "w")
tkgrid.configure(tfEntry, sticky = "w")
tkgrid.configure(clEntry, sticky = "w")
tkgrid.configure(lnEntry, sticky = "w")
tkgrid.configure(fiEntry, sticky = "w")
tkgrid.configure(pcEntry, sticky = "w")
tkgrid.configure(eEntry, sticky = "w")
dialogSuffix(rows = 8, columns = 2, focus = tfEntry)
}
################################################
#metaCiEr <-function(tf,cl,ln,fi,i,e)
#metaEr <-function(tf,cl,ln,fi,pc,e)


resEffcol <-function () 
{
initializeDialog(title = gettextRcmdr("Rescue and Internal Colonization"))
####
dsname <- tclVar("Do_Not_Save")
entryDsname <- tkentry(top, width="20", textvariable=dsname)
###
tfVar <- tclVar("100")
tfEntry <- tkentry(top, width = "4", textvariable = tfVar)
clVar <- tclVar("20")
lnVar <- tclVar("20")
clEntry <- tkentry(top, width = "4", textvariable = clVar)
lnEntry <- tkentry(top, width = "4", textvariable = lnVar)
fiVar <- tclVar(0.25) ## nclassVar ->fiVar
fiEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=fiVar, resolution=0.01, orient="horizontal")
iVar <- tclVar("0.1") ## nclassVar ->fiVar
iEntry <- tkentry(top, width = "6", textvariable = iVar)
eVar <- tclVar("0.05") ## nclassVar ->fiVar
eEntry <- tkentry(top, width = "6", textvariable = eVar)
	onOK <- function() 
	{
        #closeDialog()
        tf <- round(as.numeric(tclvalue(tfVar)))
        if (is.na(tf) || tf <= 0) 
        {
            errorCondition(message = "Number of simulations must be a positive integer")
            return()
        }
        cl <- round(as.numeric(tclvalue(clVar)))
        if (is.na(cl) || cl <= 0) 
        {
            errorCondition(message = "Number of columns on the simulated arena must be a positive integer.")
            return()
        }
        ln <- round(as.numeric(tclvalue(lnVar)))
        if (is.na(ln) || ln <= 0) 
        {
            errorCondition("Number of lines on the simulated arena must be a positive integer.")
            return()
        }
        e <- as.numeric(tclvalue(eVar))
        if (e<0 || e > 100) 
        {
            errorCondition(message = "Extinction constant must be positive ")
            return()
        }
        i <- as.numeric(tclvalue(iVar))
        if (i<0 || i > 1000) 
        {
            errorCondition(message = "Extinction constant must be positive")
            return()
        }
        fi <- as.numeric(tclvalue(fiVar))
############ Data name
   dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue=="") 
        {
        	command <- paste("metaCiEr(tmax = ",tf, ", cl = ", cl,", fi = ", fi,", rw =", ln,", ci = ", i,", ce = ", e, ")", sep = "")
        }
        else  
		  {
		  command <- paste(dsnameValue,"<-metaCiEr(tmax = ",tf, ", rw = ", cl,", fi = ", fi,", ln =", ln,", ci = ", i,", ce = ", e, ")", sep = "")
		  }
########   
#	command <- paste("metaEr(tf = ",tf, ", cl = ", cl,", fi = ", fi,", ln =", ln,", pc = ", pc,", e = ", e, ")", sep = "")
	doItAndPrint(command)
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "metapopulation")
tkgrid(tklabel(top, text="Enter data set name:"), entryDsname, sticky="e")
tkgrid(tklabel(top, text="Simulation Arena Conditions :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Maximum time  "), tfEntry, sticky = "e")
tkgrid(tklabel(top, text = "Columns  "), clEntry, sticky = "e")
tkgrid(tklabel(top, text = "Rows  "), lnEntry, sticky = "e")
tkgrid(tklabel(top, text="Species parameters  :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Initial occupancy  "), fiEntry, sticky = "e")
tkgrid(tklabel(top, text = "Colonization coefficient  "), iEntry, sticky = "e")
tkgrid(tklabel(top, text = "Extinction coefficient  "), eEntry, sticky = "e")
tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
tkgrid.configure(entryDsname, sticky = "w")
tkgrid.configure(tfEntry, sticky = "w")
tkgrid.configure(clEntry, sticky = "w")
tkgrid.configure(lnEntry, sticky = "w")
tkgrid.configure(fiEntry, sticky = "w")
tkgrid.configure(iEntry, sticky = "w")
tkgrid.configure(eEntry, sticky = "w")
dialogSuffix(rows = 8, columns = 2, focus = tfEntry)
}


####################################
########################################################
#metaComp(tmax=100,cl=100,ln=100,fi1=0.1,fi2=0.4,i1=0.4,i2=0.5,pe=0.25)
metacompDb <-function () 
{
initializeDialog(title = gettextRcmdr("Meta Competition"))
#### Salva dados
dsname <- tclVar("Do_Not_Save")
entryDsname <- tkentry(top, width="20", textvariable=dsname)
####
tmaxVar <- tclVar("100")
tmaxEntry <- tkentry(top, width = "4", textvariable = tmaxVar)
clVar <- tclVar("20")
lnVar <- tclVar("20")
clEntry <- tkentry(top, width = "4", textvariable = clVar)
lnEntry <- tkentry(top, width = "4", textvariable = lnVar)
### controle de barra

fi1Var <- tclVar(0.4)
fi1VarSlider <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=fi1Var, resolution=0.01, orient="horizontal")
fi2Var <- tclVar(0.4)
fi2VarSlider <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=fi2Var, resolution=0.01, orient="horizontal")
#####
i1Var <- tclVar("0.1") ## nclassVar ->iVar
i1Entry <- tkentry(top, width = "6", textvariable = i1Var)
i2Var <- tclVar("0.1") ## nclassVar ->iVar
i2Entry <- tkentry(top, width = "6", textvariable = i2Var)
peVar <- tclVar(0.05) 
peEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=peVar, resolution=0.01, orient="horizontal")
distrVar <- tclVar(0.00) 
distrEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=distrVar, resolution=0.01, orient="horizontal")
cantoVar <- tclVar("1")
cantoBox <- tkcheckbutton(top, variable = cantoVar)


	onOK <- function() 
	{
	#closeDialog()
   tmax <- round(as.numeric(tclvalue(tmaxVar)))
        if (is.na(tmax) || tmax <= 0) 
        {
        errorCondition(message = "Number of simulations must be a positive integer")
        return()
        }
  	cl <- round(as.numeric(tclvalue(clVar)))
        if (is.na(cl) || cl <= 0) 
        {
        errorCondition(message = "Number of columns on the simulated arena must be a positive integer.")
        return()
        }
	ln <- round(as.numeric(tclvalue(lnVar)))
        if (is.na(ln) || ln <= 0) 
        {
        errorCondition("Number of lines on the simulated arena must be a positive integer.")
        return()
        }
   i1 <- as.numeric(tclvalue(i1Var))
        if (i1<0 || i1 > 100) 
        {
        errorCondition(message = "Colonization constant must be between 0 and 10")
        return()
        }
   i2 <- as.numeric(tclvalue(i2Var))
        if (i2<0 || i2 > 100) 
        {
        errorCondition(message = "Colonization constant must be between 0 and 10")
        return()
        }
   fi1 <- as.numeric(tclvalue(fi1Var))
   fi2 <- as.numeric(tclvalue(fi2Var))
   pe <- as.numeric(tclvalue(peVar))
   distr <- as.numeric(tclvalue(distrVar))
   cantoVF <- as.logical(as.numeric(tclvalue(cantoVar)))
############ Data name
   dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue=="") 
        {
        	command <- paste("metaComp(tmax = ",tmax, ", cl = ", cl,", fi1 = ", fi1,", fi2 = ", fi2,", rw =", ln,", i1 = ", i1,", i2 = ", i2,", pe = ", pe,", D = ", distr,", anima = ", cantoVF, ")", sep = "")
        }
        else  
		  {
		  command <- paste(dsnameValue,"<- metaComp(tmax = ",tmax, ", cl = ", cl,", fi1 = ", fi1,", fi2 = ", fi2,", rw =", ln,", i1 = ", i1,", i2 = ", i2,", pe = ", pe,", D = ", distr,", anima = ", cantoVF, ")", sep = "")
		  }
########
#metaComp(tmax=100,cl=100,ln=100,fi1=0.1,fi2=0.4,i1=0.4,i2=0.5,pe=0.25)
	doItAndPrint(command)
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "metaComp")
# data name
tkgrid(tklabel(top, text="Enter name for data set:"), entryDsname, sticky="e")
##
tkgrid(tklabel(top, text="Simulation Arena Conditions :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Maximum time"), tmaxEntry, sticky = "e")
tkgrid(tklabel(top, text = "Columns"), clEntry, sticky = "e")
tkgrid(tklabel(top, text = "Rows"), lnEntry, sticky = "e")
#
tkgrid(tklabel(top, text="Best Competitor Species :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text="Initial occupancy "), fi1VarSlider, sticky="se")
tkgrid(tklabel(top, text = "Colonization coefficient "), i1Entry, sticky = "e")
tkgrid(tklabel(top, text="Inferior Competitor Species :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text="Initial occupancy "), fi2VarSlider, sticky="se")
tkgrid(tklabel(top, text = "Colonization coefficient  "), i2Entry, sticky = "e")
tkgrid(tklabel(top, text="Both Species :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Prob. Extinction  "), peEntry, sticky = "se")
tkgrid(tklabel(top, text = "Habitat Destruction  "), distrEntry, sticky = "se")
tkgrid(tklabel(top, text = "Show simulation frames"), cantoBox, sticky = "e")

tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
tkgrid.configure(entryDsname, sticky = "w")
tkgrid.configure(tmaxEntry, sticky = "w")
tkgrid.configure(clEntry, sticky = "w")
tkgrid.configure(lnEntry, sticky = "w")
tkgrid.configure(fi1VarSlider, sticky = "w")
tkgrid.configure(i1Entry, sticky = "w")
tkgrid.configure(fi2VarSlider, sticky = "w")
tkgrid.configure(i2Entry, sticky = "w")
tkgrid.configure(peEntry, sticky = "w")
tkgrid.configure(distrEntry, sticky = "w")
tkgrid.configure(cantoBox, sticky = "w")

dialogSuffix(rows = 12, columns = 2, focus = entryDsname)
}
############################################################

