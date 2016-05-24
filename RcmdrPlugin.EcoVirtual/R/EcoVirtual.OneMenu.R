#### Dialogs boxes for RcmdrPlugin.EcoVirtual package
### Alexandre Adalardo de Oliveira 17 fevereiro 2010
####################################################
## por padrao .First.lib nao existe mais
# .First.lib <- function(libname, pkgname){
#     if (!interactive()) return()
#     Rcmdr <- options()$Rcmdr
#     plugins <- Rcmdr$plugins
#     if ((!pkgname %in% plugins) && !getRcmdr("autoRestart")) {
#         Rcmdr$plugins <- c(plugins, pkgname)
#         options(Rcmdr=Rcmdr)
#         closeCommander(ask=FALSE, ask.save=TRUE)
#         Commander()
#         }
#     }
####################################################
## usando a funcao do RcmdrPlugin.TeachingDemo (Jonh Fox)
# Note: the following function (with contributions from Richard Heiberger and Milan Bouchet-Valat) 
# can be included in any Rcmdr plug-in package to cause the package to load
# the Rcmdr if it is not already loaded
.onAttach <- function(libname, pkgname){
    if (!interactive()) return()
    putRcmdr("slider.env", new.env())    
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


########################################################
estDemDb <-function () 
{
initializeDialog(title = gettextRcmdr("Demographic Stochasticity"))
#### Salva dados
dsname <- tclVar("Do_Not_Save")
entryDsname <- tkentry(top, width="20", textvariable=dsname)
####
noVar <- tclVar("100")
noEntry <- tkentry(top, width = "4", textvariable = noVar)
bVar <- tclVar(0.2)
bEntry <- tkentry(top, width = "6", textvariable = bVar)
dVar <- tclVar(0.2)
dEntry <- tkentry(top, width = "6", textvariable = dVar)
tmaxVar <- tclVar(50)
tmaxEntry <- tkentry(top, width = "4", textvariable = tmaxVar)
npopVar <-tclVar("20")
npopEntry <- tkentry(top, width = "3", textvariable = npopVar)
	onOK <- function() 
	{
        #closeDialog()
        N0 <- round(as.numeric(tclvalue(noVar)))
        if (is.na(N0) || N0 <= 0) 
        {
            errorCondition(message = "Number of individuals at the simulation start must be a positive integer")
            return()
        }
        tmax <- round(as.numeric(tclvalue(tmaxVar)))
        if (is.na(tmax) || tmax <= 0) 
        {
            errorCondition("Number of simulations must be a positive integer")
            return()
        }
        b=as.numeric(tclvalue(bVar))
        d=as.numeric(tclvalue(dVar))
        npop=as.numeric(tclvalue(npopVar))
############ Data name
##estDem(tmax=10, n=0.2, m=0.2, N0=10, nsim=20, ciclo=5000)
   dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue=="") 
        {
        	command <- paste("estDem(tmax= ",tmax, ", b = ", b,", d = ", d,", N0 =", N0, ", nsim = ", npop, ")", sep = "")
        }
        else  
		  {
		  command <- paste(dsnameValue, " <- estDem(tmax= ",tmax, ", b = ", b,", d = ", d,", N0 =", N0, ", nsim = ", npop, ")", sep = "")
		  }
	doItAndPrint(command)
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "dynPop")
tkgrid(tklabel(top, text="Enter name for data set: "), entryDsname, sticky="e")
#tkgrid(tklabel(top, text="Simulation Conditions :  ", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Maximum time "), tmaxEntry, sticky = "e")
tkgrid(tklabel(top, text = "Number of simulations "), npopEntry, sticky = "e")
tkgrid(tklabel(top, text="Population parameters :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Initial size  "), noEntry, sticky = "e")
tkgrid(tklabel(top, text = "Birth rate  "), bEntry, sticky = "se")
tkgrid(tklabel(top, text = "Death rate  "), dEntry, sticky = "e")
tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#tkgrid(tklabel(top, text="Enter name for data set:"), entryDsname, sticky="e")
tkgrid.configure(entryDsname, sticky = "w")
tkgrid.configure(tmaxEntry, sticky = "w")
tkgrid.configure(npopEntry, sticky = "w")
tkgrid.configure(noEntry, sticky = "w")
tkgrid.configure(bEntry, sticky = "w")
tkgrid.configure(dEntry, sticky = "w")
dialogSuffix(rows = 6, columns = 2, focus = tmaxEntry)
}
#########################################################
#########################################################
popExpDb<-function()
{
initializeDialog(title = gettextRcmdr("Exponential Growth"))
#### Salva dados
dsname <- tclVar("Do_Not_Save")
entryDsname <- tkentry(top, width="20", textvariable=dsname)
####
noVar <- tclVar("10")
noEntry <- tkentry(top, width = "4", textvariable = noVar)
lambVar <- tclVar(1.05)
tmaxVar <- tclVar(10)
tmaxEntry <- tkentry(top, width = "4", textvariable = tmaxVar)
inttVar<-tclVar(1.00)
#########################
	set.gr=function(...)
	{
#extVF <- as.logical(as.numeric(tclvalue(extVar)))
	command=paste("popExp(N0 = ", as.numeric(tclvalue(noVar)), ", lamb = ", as.numeric(tclvalue(lambVar)), ", tmax = ", as.numeric(tclvalue(tmaxVar)),", intt = ", as.numeric(tclvalue(inttVar)),")", sep="")
	doItAndPrint(command)
#	tkfocus(CommanderWindow())
	}
##########################
lambEntry<-tkscale(top, from=0, to=5, showvalue=TRUE, variable=lambVar, resolution=0.01, orient="horizontal", 
command=set.gr)
inttEntry <- tkscale(top, from=0.01, to=1, showvalue=TRUE, variable=inttVar, resolution=0.01, orient="horizontal",command=set.gr)
################################################
	onOK <- function() 
	{	
############ Data name
   dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue=="") 
        {
        	command <- paste("popExp(N0 = ", as.numeric(tclvalue(noVar)), ", lamb = ", as.numeric(tclvalue(lambVar)), ", tmax = ", as.numeric(tclvalue(tmaxVar)),", intt = ", as.numeric(tclvalue(inttVar)),")", sep="")
        }
        else  
		  {
		  command <- paste(dsnameValue, " <- popExp(N0 = ", as.numeric(tclvalue(noVar)), ", lamb = ", as.numeric(tclvalue(lambVar)), ", tmax = ", as.numeric(tclvalue(tmaxVar)),", intt = ", as.numeric(tclvalue(inttVar)),")", sep="")
		  }
	doItAndPrint(command)
	tkfocus(CommanderWindow())
	}
#popExp(N0,lamb,tmax, intt= 1) 
OKCancelHelp(helpSubject = "dynPop")
tkgrid(tklabel(top, text="Enter name for last simulation data set: "), entryDsname, sticky="e")
#tkgrid(tklabel(top, text="Simulation Arena Conditions :  ", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Maximum time"), tmaxEntry, sticky = "e")
tkgrid(tklabel(top, text = "Interval time size "), inttEntry, sticky = "e")
tkgrid(tklabel(top, text="Species parameters :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Initial population size  "), noEntry, sticky = "e")
tkgrid(tklabel(top, text = "Population growth rate (lambda)  "), lambEntry, sticky = "e")

tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#tkgrid(tklabel(top, text="Enter name for data set:"), entryDsname, sticky="e")
tkgrid.configure(entryDsname, sticky = "w")
tkgrid.configure(tmaxEntry, sticky = "w")
tkgrid.configure(inttEntry, sticky = "w")
tkgrid.configure(noEntry, sticky = "w")
tkgrid.configure(lambEntry, sticky = "w")
dialogSuffix(rows = 6, columns = 2, focus = tmaxEntry)
}
############################################################
estEnvDb <-function () 
{
initializeDialog(title = gettextRcmdr("Environmental Sthocasticity"))
#### Salva dados
dsname <- tclVar("Do_Not_Save")
entryDsname <- tkentry(top, width="20", textvariable=dsname)
####
noVar <- tclVar("10")
noEntry <- tkentry(top, width = "4", textvariable = noVar)
lambVar <- tclVar(1.05)
varrVar <- tclVar(0.02)
tmaxVar <- tclVar(100)
tmaxEntry <- tkentry(top, width = "4", textvariable = tmaxVar)
npopVar<-tclVar("20")
npopEntry <- tkentry(top, width = "4", textvariable = npopVar)
extVar <- tclVar("0")
extBox <- tkcheckbutton(top, variable = extVar)
#################################################################
	set.gr=function(...)
	{
#extVF <- as.logical(as.numeric(tclvalue(extVar)))
	command=paste("estEnv(N0 = ", as.numeric(tclvalue(noVar)), ", lamb = ", as.numeric(tclvalue(lambVar)), ", varr = ", as.numeric(tclvalue(varrVar)),", tmax = ", as.numeric(tclvalue(tmaxVar)),", npop = ", as.numeric(tclvalue(npopVar)), ", ext =", as.logical(as.numeric(tclvalue(extVar))), ")", sep="")
	doItAndPrint(command)
#	tkfocus(CommanderWindow())
	}
##############################################################
lambEntry<-tkscale(top, from=0, to=5, showvalue=TRUE, variable=lambVar, resolution=0.01, orient="horizontal", 
command=set.gr)
varrEntry <- tkscale(top, from=0, to=15, showvalue=TRUE, variable=varrVar, resolution=0.01, orient="horizontal", 
command=set.gr)
################################################
	onOK <- function() 
	{	
############ Data name
   dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue=="") 
        {
        	command <- paste("estEnv(N0 = ", as.numeric(tclvalue(noVar)), ", lamb = ", as.numeric(tclvalue(lambVar)), ", varr = ", as.numeric(tclvalue(varrVar)),", tmax = ", as.numeric(tclvalue(tmaxVar)),", npop = ", as.numeric(tclvalue(npopVar)), ", ext =", as.logical(as.numeric(tclvalue(extVar))), ")", sep = "")
        }
        else  
		  {
		  command <- paste(dsnameValue, " <- estEnv(N0 = ", as.numeric(tclvalue(noVar)), ", lamb = ", as.numeric(tclvalue(lambVar)), ", varr = ", as.numeric(tclvalue(varrVar)),", tmax = ", as.numeric(tclvalue(tmaxVar)),", npop = ", as.numeric(tclvalue(npopVar)), ", ext =", as.logical(as.numeric(tclvalue(extVar))), ")", sep = "")
		  }
	doItAndPrint(command)
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "dynPop")
tkgrid(tklabel(top, text="Enter name for last simulation data set: "), entryDsname, sticky="e")
#tkgrid(tklabel(top, text="Simulation Arena Conditions :  ", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Maximum time"), tmaxEntry, sticky = "e")
tkgrid(tklabel(top, text = "Number of stochastic simulations "), npopEntry, sticky = "e")
tkgrid(tklabel(top, text="Species parameters :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Initial population size  "), noEntry, sticky = "e")
tkgrid(tklabel(top, text = "Population growth rate (lambda) "), lambEntry, sticky = "e")
tkgrid(tklabel(top, text="Environment stochasticity: ", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "lambda variance "), varrEntry, sticky = "e")
tkgrid(tklabel(top, text = "Population Extinction"), extBox, sticky = "e")

tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
#tkgrid(tklabel(top, text="Enter name for data set:"), entryDsname, sticky="e")
tkgrid.configure(entryDsname, sticky = "w")
tkgrid.configure(tmaxEntry, sticky = "w")
tkgrid.configure(npopEntry, sticky = "w")
tkgrid.configure(noEntry, sticky = "w")
tkgrid.configure(varrEntry, sticky = "w")
tkgrid.configure(lambEntry, sticky = "w")
tkgrid.configure(extBox, sticky = "w")
dialogSuffix(rows = 6, columns = 2, focus = tmaxEntry)
}
#########################################################
##########################################################
#crescLog(N0=10, r=0.05, K=80, tmax=100)
popLogDb<-function () 
{
initializeDialog(title = gettextRcmdr("Logistic Growth"))
#### Salva dados
dsname <- tclVar("Do_Not_Save")
entryDsname <- tkentry(top, width="20", textvariable=dsname)
####
noVar <- tclVar("10")
rVar <- tclVar(0.05)
kVar <- tclVar("20")
tmaxVar <- tclVar(100)
noEntry <- tkentry(top, width = "4", textvariable = noVar)
tmaxEntry <- tkentry(top, width = "4", textvariable = tmaxVar)
kEntry<-tkentry(top, width = "4", textvariable = kVar)
extVar <- tclVar("0")
extBox <- tkcheckbutton(top, variable = extVar)
#################################################################
	set.gr=function(...)
	{
	command <- paste("popLog(N0= ", as.numeric(tclvalue(noVar)), ", r = ", as.numeric(tclvalue(rVar)),", K = ", as.numeric(tclvalue(kVar)),", tmax =", round(as.numeric(tclvalue(tmaxVar))), ", ext =", as.logical(as.numeric(tclvalue(extVar))),")", sep = "")
	doItAndPrint(command)
#	tkfocus(CommanderWindow())
	}
##############################################################
rEntry <-tkscale(top, from=-5, to=5, showvalue=TRUE, variable=rVar, resolution=0.01, orient="horizontal", command=set.gr)
#kEntry <- tkscale(top, from=0, to= 10*as.numeric(tclvalue(noVar)), showvalue=TRUE, variable=kVar, resolution=round(0.1*as.numeric(tclvalue(noVar))), orient="horizontal", command=set.gr)
##############################################################
	onOK <- function() 
	{
N0 <- round(as.numeric(tclvalue(noVar)))
#kEntry <- tkscale(top, from=0, to= 10*N0, showvalue=TRUE, variable=kVar, resolution=round(0.1*N0), orient="horizontal")
tmax <- round(as.numeric(tclvalue(tmaxVar)))
K <- as.numeric(tclvalue(kVar))
r=as.numeric(tclvalue(rVar))
#closeDialog()
############ Data name
   dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue=="") 
        {
        	command <- paste("popLog(N0= ",N0, ", r = ", r,", K = ", K,", tmax =", tmax, ", ext =", as.logical(as.numeric(tclvalue(extVar))),")", sep = "")
        }
        else  
		  {
		  command <- paste(dsnameValue,"<-popLog(N0= ",N0, ", r = ", r,", K = ", K,", tmax =", tmax, ", ext =", as.logical(as.numeric(tclvalue(extVar))),")", sep = "")
		  }
########
	doItAndPrint(command)
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "dynPop")
tkgrid(tklabel(top, text="Enter name for data set: "), entryDsname, sticky="e")
#tkgrid(tklabel(top, text="Simulation Arena Conditions : ", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Maximum time "), tmaxEntry, sticky = "e")
tkgrid(tklabel(top, text="Species parameters: ", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Initial population size "), noEntry, sticky = "e")
tkgrid(tklabel(top, text = "Carrying capacity (K) "), kEntry, sticky = "e")
tkgrid(tklabel(top, text = "Intrinsic growth rate (r) "), rEntry, sticky = "se")
tkgrid(tklabel(top, text = "Population Extinction"), extBox, sticky = "e")
tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
tkgrid.configure(entryDsname, sticky = "w")
tkgrid.configure(tmaxEntry, sticky = "w")
tkgrid.configure(noEntry, sticky = "w")
tkgrid.configure(kEntry, sticky = "w")
tkgrid.configure(rEntry, sticky = "w")
tkgrid.configure(extBox, sticky = "w")
dialogSuffix(rows = 6, columns = 2, focus = tmaxEntry)
}
##########################################
#compLV(n01=10, n02=10,r1=0.05, r2=0.03, k1=80, k2=50, alfa=1.2, beta=0.5, tmax=200)
##########################################
logBifDb<-function () 
{
initializeDialog(title = gettextRcmdr("Logistic Bifurcation"))
#### Salva dados
dsname <- tclVar("Do_Not_Save")
entryDsname <- tkentry(top, width="20", textvariable=dsname)
####
noVar <- tclVar("10")
rdminVar <- tclVar(1)
rdmaxVar <- tclVar(3)
kVar <- tclVar("20")
nrdVar <- tclVar("500")
tmaxVar <- tclVar(200)
noEntry <- tkentry(top, width = "4", textvariable = noVar)
tmaxEntry <- tkentry(top, width = "4", textvariable = tmaxVar)
kEntry<-tkentry(top, width = "4", textvariable = kVar)
nrdEntry<-tkentry(top, width = "4", textvariable = nrdVar)
#extVar <- tclVar("0")
#extBox <- tkcheckbutton(top, variable = extVar)
##################################################################
#	set.gr=function(...)
#	{
#	command <- paste("popLog(N0= ", as.numeric(tclvalue(noVar)), ", r = ", as.numeric(tclvalue(rVar)),", K = ", as.numeric(tclvalue(kVar)),", tmax =", round(as.numeric(tclvalue(tmaxVar))), ", ext =", as.logical(as.numeric(tclvalue(extVar))),")", sep = "")
#	doItAndPrint(command)
##	tkfocus(CommanderWindow())
#	}
##############################################################
rdminEntry <-tkscale(top, from= 1, to=2.999, showvalue=TRUE, variable=rdminVar, resolution=0.001, orient="horizontal")#, command=set.gr
rdmaxEntry <-tkscale(top, from= 2, to=3, showvalue=TRUE, variable=rdmaxVar, resolution=0.001, orient="horizontal")
nrdEntry<-tkscale(top, from= 100, to=1000, showvalue=TRUE, variable=nrdVar, resolution=100, orient="horizontal")
#kEntry <- tkscale(top, from=0, to= 10*as.numeric(tclvalue(noVar)), showvalue=TRUE, variable=kVar, resolution=round(0.1*as.numeric(tclvalue(noVar))), orient="horizontal", command=set.gr)
##############################################################
	onOK <- function() 
	{
N0 <- round(as.numeric(tclvalue(noVar)))
#kEntry <- tkscale(top, from=0, to= 10*N0, showvalue=TRUE, variable=kVar, resolution=round(0.1*N0), orient="horizontal")
tmax <- round(as.numeric(tclvalue(tmaxVar)))
K <- as.numeric(tclvalue(kVar))
   if (K <= N0) 
     {
     errorCondition("K should be larger than N0")
     return()
     }
nrd <- as.numeric(tclvalue(nrdVar))
rdmin=as.numeric(tclvalue(rdminVar))
rdmax=as.numeric(tclvalue(rdmaxVar))
#bifAttr(N0, K, tmax, nrd,maxrd=3, minrd=1)
############ Data name
   dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue=="") 
        {
        	command <- paste("bifAttr(N0= ",N0, ", K = ", K, ", tmax =", tmax, ", nrd = ", nrd, ", minrd = ", rdmin,", maxrd = ", rdmax,")", sep = "")
        }
        else  
		  {
		  command <- paste(dsnameValue,"<-bifAttr(N0= ",N0, ", K = ", K, ", tmax =", tmax, ", nrd = ", nrd, ", minrd = ", rdmin,", maxrd = ", rdmax,")", sep = "")
		  }
########
	doItAndPrint(command)
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "dynPop")
tkgrid(tklabel(top, text="Enter name for data set:  "), entryDsname, sticky="e")
#tkgrid(tklabel(top, text="Simulation Arena Conditions : ", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Time to convergence  "), tmaxEntry, sticky = "e")
tkgrid(tklabel(top, text="Species parameters :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Initial population size  "), noEntry, sticky = "e")
tkgrid(tklabel(top, text = "Carrying capacity (K) "), kEntry, sticky = "e")
tkgrid(tklabel(top, text = "Minimum discrete growth rate (rd) "), rdminEntry, sticky = "se")
tkgrid(tklabel(top, text = "Maximum discrete growth rate (rd) "), rdmaxEntry, sticky = "se")
tkgrid(tklabel(top, text = "Number of rds"), nrdEntry, sticky = "se")
#tkgrid(tklabel(top, text = "Population Extinction"), extBox, sticky = "e")
tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
tkgrid.configure(entryDsname, sticky = "w")
tkgrid.configure(tmaxEntry, sticky = "w")
tkgrid.configure(noEntry, sticky = "w")
tkgrid.configure(kEntry, sticky = "w")
tkgrid.configure(rdminEntry, sticky = "w")
tkgrid.configure(rdmaxEntry, sticky = "w")
tkgrid.configure(nrdEntry, sticky = "w")
#tkgrid.configure(extBox, sticky = "w")
dialogSuffix(rows = 7, columns = 2, focus = tmaxEntry)
}
##########################################
##########################################
compDb<-function () 
{
initializeDialog(title = gettextRcmdr("Competition LV Model"))
#### Salva dados
dsname <- tclVar("Do_Not_Save")
entryDsname <- tkentry(top, width="20", textvariable=dsname)
####
n01Var <- tclVar("10")
n01Entry <- tkentry(top, width = "4", textvariable = n01Var)
n02Var <- tclVar("10")
n02Entry <- tkentry(top, width = "4", textvariable = n02Var)
r1Var <- tclVar(0.05)
r1Entry=tkscale(top, from=-5, to=5, showvalue=TRUE, variable=r1Var, resolution=0.01, orient="horizontal")
r2Var <- tclVar(0.05)
r2Entry=tkscale(top, from=-5, to=5, showvalue=TRUE, variable=r2Var, resolution=0.01, orient="horizontal")
k1Var <- tclVar("20")
k1Entry <- tkentry(top, width = "4", textvariable = k1Var)
k2Var <- tclVar("20")
k2Entry <- tkentry(top, width = "4", textvariable = k2Var)
alfaVar <- tclVar("1.2")
alfaEntry <- tkentry(top, width = "6", textvariable = alfaVar)
betaVar <- tclVar("0.5")
betaEntry <- tkentry(top, width = "6", textvariable = betaVar)
tmaxVar <- tclVar(100)
tmaxEntry <- tkentry(top, width = "4", textvariable = tmaxVar)
	onOK <- function() 
	{
        closeDialog()
        n01 <- round(as.numeric(tclvalue(n01Var)))
        if (is.na(n01) || n01 <= 0) 
        {
            errorCondition(message = "Number of individuals at the simulation start must be a positive integer")
            return()
        }
        n02 <- round(as.numeric(tclvalue(n02Var)))
        if (is.na(n02) || n02 <= 0) 
        {
            errorCondition(message = "Number of individuals at the simulation start must be a positive integer")
            return()
        }
        tmax <- round(as.numeric(tclvalue(tmaxVar)))
        if (is.na(tmax) || tmax <= 0) 
        {
            errorCondition("Number of simulations must be a positive integer")
            return()
        }
        k1 <- as.numeric(tclvalue(k1Var))
        if (is.na(k1) || k1 <= 0)
        {
            errorCondition(message = "Carrying Capacity (K) must be a positive integer")
            return()
        }
        k2 <- as.numeric(tclvalue(k2Var))
        if (is.na(k2) || k2 <= 0)
        {
            errorCondition(message = "Carrying Capacity (K) must be a positive integer")
            return()
        }
        r1=as.numeric(tclvalue(r1Var))
        r2=as.numeric(tclvalue(r2Var))
        alfa=as.numeric(tclvalue(alfaVar))
        beta=as.numeric(tclvalue(betaVar))
############ Data name
   dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue=="") 
        {
        	command <- paste("compLV(n01= ",n01, ",n02= ",n02, ", r1 = ", r1,", r2 = ", r2,", k1 = ", k1,", k2 = ", k2,", alfa = ", alfa,", beta = ", beta,", tmax =", tmax,")", sep = "")
        }
        else  
		  {
		  command <- paste("compLV(n01= ",n01, "n02= ",n02, ", r1 = ", r1,", r2 = ", r2,", k1 = ", k1,", k2 = ", k2,", alfa = ", alfa,", beta = ", beta,", tmax =", tmax,")", sep = "")
		  }
#compLV(n01=10, n02=10,r1=0.05, r2=0.03, k1=80, k2=50, alfa=1.2, beta=0.5, tmax=200)
########
	doItAndPrint(command)
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "compLV")
tkgrid(tklabel(top, text="Enter name for data set:"), entryDsname, sticky="e")
tkgrid(tklabel(top, text = "Number of simulations  "), tmaxEntry, sticky = "e")
tkgrid(tklabel(top, text="Best competitor species parameters : ", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Initial population  "), n01Entry, sticky = "e")
tkgrid(tklabel(top, text = "Carrying capacity (K)  "), k1Entry, sticky = "e")
tkgrid(tklabel(top, text = "Intrinsic growth rate  "), r1Entry, sticky = "se")
tkgrid(tklabel(top, text = "Alpha coefficiente"), alfaEntry, sticky = "e")
tkgrid(tklabel(top, text="Worse Competitor Species :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Initial population  "), n02Entry, sticky = "e")
tkgrid(tklabel(top, text = "Carrying capacity (K)  "), k2Entry, sticky = "e")
tkgrid(tklabel(top, text = "Intrinsic growth rate  "), r2Entry, sticky = "se")
tkgrid(tklabel(top, text = "Beta coefficiente"), betaEntry, sticky = "e")
tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
tkgrid.configure(entryDsname, sticky = "w")
tkgrid.configure(tmaxEntry, sticky = "w")
tkgrid.configure(n01Entry, sticky = "w")
tkgrid.configure(n02Entry, sticky = "w")
tkgrid.configure(k1Entry, sticky = "w")
tkgrid.configure(k2Entry, sticky = "w")
tkgrid.configure(r1Entry, sticky = "w")
tkgrid.configure(r2Entry, sticky = "w")
tkgrid.configure(alfaEntry, sticky = "w")
tkgrid.configure(betaEntry, sticky = "w")
dialogSuffix(rows = 11, columns = 2, focus = tmaxEntry)
}
##########################################
######################################
popStrDb<-function () 
{
initializeDialog(title = gettextRcmdr("Structured Population"))
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
sjVar <- tclVar(0.4)
sjEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=sjVar, resolution=0.01, orient="horizontal")
jjVar <- tclVar(0.4)
jjEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=jjVar, resolution=0.01, orient="horizontal")
jaVar <- tclVar(0.4)
jaEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=jaVar, resolution=0.01, orient="horizontal")
aaVar <- tclVar(0.4)
aaEntry <- tkscale(top, from=0, to=1, showvalue=TRUE, variable=aaVar, resolution=0.01, orient="horizontal")
#####
fecVar <- tclVar("1.2") 
fecEntry <- tkentry(top, width = "6", textvariable = fecVar)
nsVar <- tclVar("50") 
nsEntry <- tkentry(top, width = "6", textvariable = nsVar)
njVar <- tclVar("20") 
njEntry <- tkentry(top, width = "6", textvariable = njVar)
naVar <- tclVar("10") 
naEntry <- tkentry(top, width = "6", textvariable = naVar)
	onOK <- function() 
	{
	closeDialog()
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
        errorCondition("Number of rows on the simulated arena must be a positive integer.")
        return()
        }
   ns <- as.numeric(tclvalue(nsVar))
        if (ns<0 ) 
        {
        errorCondition(message = "Mean number of propagules must be positive")
        return()
        }
   nj <- as.numeric(tclvalue(njVar))
        if (nj<0 | nj > (cl*ln) ) 
        {
        errorCondition(message = "Number of juveniles must be positive and less than number of patches")
        return()        
        }
   na <- as.numeric(tclvalue(naVar))
        if (na<0 | (na+nj) > (cl*ln) ) 
        {
        errorCondition(message = "Number of adults must be positive and adults plus juvenils must be less than number of patches")
        return()
        }
   p.sj <- as.numeric(tclvalue(sjVar))
   p.jj <- as.numeric(tclvalue(jjVar))
   p.ja <- as.numeric(tclvalue(jaVar))
   p.aa <- as.numeric(tclvalue(aaVar))
   fec<- as.numeric(tclvalue(fecVar))
############ Data name
#popStr(p.sj, p.jj, p.ja, p.aa, fec, ns,nj,na, ln, cl, tmax)
   dsnameValue <- trim.blanks(tclvalue(dsname))
        if (dsnameValue == "Do_Not_Save" | dsnameValue=="") 
        {
        	command <- paste("popStr(tmax = ",tmax, ", cl = ", cl,", rw = ", ln, ", p.sj = ", p.sj,", p.jj = ", p.jj,", p.ja =", p.ja,", p.aa = ", p.aa,", fec = ", fec,", ns = ", ns,", nj = ", nj,", na = ", na, ")", sep = "")
        }
        #popStr=function(p.sj, p.jj, p.ja, p.aa, fec, ns,nj,na, ln, cl, tmax)
        else  
		  {
		  command <- paste(dsnameValue,"<- popStr(tmax = ",tmax, ", cl = ", cl,", rw = ", ln, ", p.sj = ", p.sj,", p.jj = ", p.jj,", p.ja =", p.ja,", p.aa = ", p.aa,", fec = ", fec,", ns = ", ns,", nj = ", nj,", na = ", na, ")", sep = "")
		  }
########
##popStr(p.sj, p.jj, p.ja, p.aa, fec, ns,nj,na, ln, cl, tmax)
	doItAndPrint(command)
	tkfocus(CommanderWindow())
	}
OKCancelHelp(helpSubject = "dynPop")
# data name
tkgrid(tklabel(top, text="Enter name for data set:"), entryDsname, sticky="e")
##
tkgrid(tklabel(top, text="Simulation Arena Conditions :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Maximum time"), tmaxEntry, sticky = "e")
tkgrid(tklabel(top, text = "Columns"), clEntry, sticky = "e")
tkgrid(tklabel(top, text = "Rows"), lnEntry, sticky = "e")
#
tkgrid(tklabel(top, text="Initial Population Size :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text = "Seeds "), nsEntry, sticky = "e")
tkgrid(tklabel(top, text = "Saplings "), njEntry, sticky = "e")
tkgrid(tklabel(top, text = "Adults "), naEntry, sticky = "e")
tkgrid(tklabel(top, text="Transitions Probabilities :", fg="blue"), sticky="w")
tkgrid(tklabel(top, text="Seed to sapling "), sjEntry, sticky="se")
tkgrid(tklabel(top, text="Sapling persistence "), jjEntry, sticky="se")
tkgrid(tklabel(top, text="Sapling to adult "), jaEntry, sticky="se")
tkgrid(tklabel(top, text="Adult survival "), aaEntry, sticky="se")
tkgrid(tklabel(top, text="Fecundity "), fecEntry, sticky="e")
##
tkgrid(buttonsFrame, sticky = "w", columnspan = 2)
##
tkgrid.configure(entryDsname, sticky = "w")
tkgrid.configure(tmaxEntry, sticky = "w")
tkgrid.configure(clEntry, sticky = "w")
tkgrid.configure(lnEntry, sticky = "w")
tkgrid.configure(sjEntry, sticky = "w")
tkgrid.configure(jjEntry, sticky = "w")
tkgrid.configure(jaEntry, sticky = "w")
tkgrid.configure(aaEntry, sticky = "w")
tkgrid.configure(fecEntry, sticky = "w")
tkgrid.configure(nsEntry, sticky = "w")
tkgrid.configure(njEntry, sticky = "w")
tkgrid.configure(naEntry, sticky = "w")
dialogSuffix(rows = 14, columns = 2, focus = tmaxEntry)
}
############################################################

