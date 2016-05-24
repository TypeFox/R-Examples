# R MTC Jags
# Copyright (C) 2014. Marcelo Goulart Correia

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# Marcelo Goulart Correia
# Rua das Laranjeiras, 374 - 5o. Andar
# Laranjeiras - Rio de Janeiro - RJ
# Zip code: 22240-005
# marcelo.goulart@inc.saude.gov.br
# mgoulart.inc@gmail.com

# last modified: 2015-09-08 by M.G. Correia

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

RunMTCModel <- function(){
    Library('runjags')
    Library('rmeta')
    Library('coda')
    defaults <- list(ArmsVar="", TreatVar="", StudVar="", ChainsVar="4", SimsVar="10000", BurnInVar="4000", models="fem")
    dialog.values <- getDialog("RunMTCModel", defaults)
    initializeDialog(title=gettextRcmdr("Setting up MTC model"))
    
    .activeModel <- ActiveModel()
    currentModel <- if (!is.null(.activeModel))
				class(get(.activeModel, envir=.GlobalEnv))[1] == "runjags.model"
			else FALSE
	if (currentModel) {
		currentFields <- formulaFields(get(.activeModel, envir=.GlobalEnv), hasLhs=TRUE)
		if (currentFields$data != ActiveDataSet()) currentModel <- FALSE
			  }

    UpdateModelNumber()
    modelName <- tclVar(paste("MTCModel.", getRcmdr("modelNumber"), sep = ""))
    modelFrame <- tkframe(top)
    model <- ttkentry(modelFrame, width = "20", textvariable = modelName) 
    ArmsVar <- tclVar(dialog.values$ArmsVar)
    ArmsEntry <- tkentry(top, width="6", textvariable=ArmsVar)
    TreatVar <- tclVar(dialog.values$TreatVar)
    TreatEntry <- tkentry(top, width="6", textvariable=TreatVar)
    StudVar <- tclVar(dialog.values$StudVar)
    StudEntry <- tkentry(top, width="6", textvariable=StudVar)
    ChainsVar <- tclVar(dialog.values$ChainsVar)
    ChainsEntry <- tkentry(top, width="6", textvariable=ChainsVar)
    SimsVar <- tclVar(dialog.values$SimsVar)
    SimsEntry <- tkentry(top, width="6", textvariable=SimsVar)  
    BurnInVar <- tclVar(dialog.values$BurnInVar)
    BurnInEntry <- tkentry(top, width="6", textvariable=BurnInVar)
    radioButtons(name="models", buttons=c("fem", "rem", "threerem", "marem"), labels=gettextRcmdr(c("FE Model", "RE Model ignoring multi-arms effect", "RE Model for 2- and 3-arms trial", "RE Model for multi-arms trial")), title=gettextRcmdr("Models"))   
    
    onOK <- function(){
        modelValue <- trim.blanks(tclvalue(modelName))
        closeDialog()      
	if (!is.valid.name(modelValue)){
      	    errorCondition(recall=RunMTCModel, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
      	    return()
    	    }  
        Arms <- as.numeric(tclvalue(ArmsVar))
        if (is.na(Arms) || Arms <= 0){
            errorCondition(recall=RunMTCModel, message="Number of arms must be a positive number.")
            return()
            }
        Treat <- as.numeric(tclvalue(TreatVar))
        if (is.na(Treat) || Treat <= 0){
            errorCondition(recall=RunMTCModel, 
                message="Number of treatments must be a positive number.")
            return()
            }
        Stud <- round(as.numeric(tclvalue(StudVar)))
        if (is.na(Stud) || Stud <= 0){
            errorCondition(recall=RunMTCModel, message="Number of studies must be a positive integer.")
            return()
            }
        Chain <- round(as.numeric(tclvalue(ChainsVar)))
        if (is.na(Chain) || Chain <= 0){
            errorCondition(recall=RunMTCModel, message="Number of chains must be a positive integer.")
            return()
            }
        Sims <- as.numeric(tclvalue(SimsVar))
        if (is.na(Sims) || Sims <= 0){
            errorCondition(recall=RunMTCModel, message="Number of simulations must be a positive integer.")
            return()
            }
	BurnIn <- as.numeric(tclvalue(BurnInVar))
        if (is.na(BurnIn) || BurnIn <= 0){
            errorCondition(recall=RunMTCModel, message="Number of Burn-in sims must be a positive integer.")
            return()
            }
	modelValue <- trim.blanks(tclvalue(modelName))
		if (!is.valid.name(modelValue)) {
			UpdateModelNumber(-1)
			errorCondition(recall = linearRegressionModel, message = sprintf(gettextRcmdr("\"%s\" is not a valid name."), 
							modelValue))
			return()
		}
		if (is.element(modelValue, listLinearModels())) {
			if ("no" == tclvalue(checkReplace(modelValue, type = gettextRcmdr("Model")))) {
				UpdateModelNumber(-1)
				linearRegressionModel()
				return()
			}
		}
	putDialog("RunMTCModel", lapply(list(ArmsVar=Arms, TreatVar=Treat, StudVar=Stud, ChainsVar=Chain, 
                                                      SimsVar=Sims, BurnInVar=BurnIn), as.character))
	models <- tclvalue(modelsVariable)

	if(models=="marem")
		{
		justDoIt(paste("newdata <- subset(", ActiveDataSet(),",", "t..1. == 1" ,")", sep = ""))
		doItAndPrint(paste("newdata$Prop <- with(newdata,ifelse(r..1. > 0, r..1./n..1., 'NA'))", sep =""))
		doItAndPrint(paste("newdata$Prop <- as.numeric(paste(newdata$Prop))",sep=""))
		doItAndPrint(paste("newdata <- newdata[complete.cases(newdata$Prop),]",sep=""))
		doItAndPrint(paste("newdata$LN_Prop <- with(newdata,ifelse(r..1. > 0, log(r..1./n..1.), 'NA'))",sep=""))
		doItAndPrint(paste("newdata$LN_Prop_SE <- with(newdata,ifelse(r..1. > 0, log(sqrt(Prop*(1-Prop)/n..1.)), 'NA'))",sep=""))
		doItAndPrint(paste("meta_mA <- with(newdata,meta.summaries(LN_Prop, LN_Prop_SE,method=c('random')), logeffect=T)",sep=""))
		doItAndPrint(paste("mean_mA <- meta_mA$summary",sep=""))
		doItAndPrint(paste("prec_mA <- 1/with(newdata,var(LN_Prop,na.rm=T))",sep=""))
		doItAndPrint(paste("output_mA <- capture.output(cat('mA ~ dnorm(',mean_mA,',',prec_mA,')'))",sep = ""))
		}
	if(models== "fem" | models== "rem" | models== "threerem")
	        {
		justDoIt(paste("newdata <- subset(", ActiveDataSet(),",", "t == 1" ,")", sep = ""))
		doItAndPrint(paste("newdata$Prop <- with(newdata,ifelse(r > 0, r/n, 'NA'))", sep =""))
		doItAndPrint(paste("newdata$Prop <- as.numeric(newdata$Prop)", sep = ""))
		doItAndPrint(paste("newdata <- newdata[complete.cases(newdata),]", sep = ""))
		doItAndPrint(paste("newdata$LN_Prop <- with(newdata,ifelse(r > 0, log(r/n), 'NA'))", sep = ""))
		doItAndPrint(paste("newdata$LN_Prop_SE <- with(newdata,ifelse(r > 0, log(sqrt(Prop*(1-Prop)/n)), 'NA'))", sep = ""))
		doItAndPrint(paste("meta_mA <- with(newdata,meta.summaries(LN_Prop, LN_Prop_SE,method=c('random')), logeffect=T)", sep = ""))
		doItAndPrint(paste("mean_mA <- meta_mA$summary",sep=""))
		doItAndPrint(paste("prec_mA <- 1/with(newdata,var(LN_Prop,na.rm=T))",sep=""))
		doItAndPrint(paste("output_mA <- capture.output(cat('mA ~ dnorm(',mean_mA,',',prec_mA,')'))", sep = ""))
	        }

	if(models=="fem")
		{
		Teste <- paste("Model <- 'model{", "for(i in 1:N){", "logit(p[s[i],t[i]]) <- mu[s[i]] + d[t[i]] - d[b[i]]",
		"r[i]~dbin(p[s[i],t[i]],n[i])}", "for(j in 1:NS){", "mu[j]~dnorm(0,.0001)}", "for (k in 2:NT){", 
		"d[k] ~ dnorm(0,.0001)}", "d[1] <- 0", output_mA, "for (l in 1:NT) { logit(T[l]) <- mA + d[l] }}'", sep = "\n")
		doItAndPrint(Teste)
		command <- paste(modelValue," <- run.jags(model=Model, monitor=c('d','T'), n.chains=", Chain, ", method='rjags', data = list('s'=",ActiveDataSet(),"$s, 't'=",ActiveDataSet(),"$t, 'r'=",ActiveDataSet(),"$r, 'n'=",ActiveDataSet(),"$n, 'b'=",ActiveDataSet(),"$b,'N'=",Arms,",'NS'=",Stud,",'NT'=",Treat,"), burnin =", BurnIn, ",sample =", Sims,",silent.jags=T",")", sep="")
		command_aux <- paste(modelValue,"_aux <- run.jags(model=Model, monitor=c('T'), n.chains=", Chain, ", method='rjags', data = list('s'=",ActiveDataSet(),"$s, 't'=",ActiveDataSet(),"$t, 'r'=",ActiveDataSet(),"$r, 'n'=",ActiveDataSet(),"$n, 'b'=",ActiveDataSet(),"$b,'N'=",Arms,",'NS'=",Stud,",'NT'=",Treat,"), burnin =", BurnIn, ",sample =", Sims,",silent.jags=T",")", sep="")
		command_aux2 <- paste(modelValue,"_aux2 <- run.jags(model=Model, monitor=c('d'), n.chains=", Chain, ", method='rjags', data = list('s'=",ActiveDataSet(),"$s, 't'=",ActiveDataSet(),"$t, 'r'=",ActiveDataSet(),"$r, 'n'=",ActiveDataSet(),"$n, 'b'=",ActiveDataSet(),"$b,'N'=",Arms,",'NS'=",Stud,",'NT'=",Treat,"), burnin =", BurnIn, ",sample =", Sims,",silent.jags=T",")", sep="")
		}

	if(models=="rem")
		{
		Teste <- paste("Model <- 'model{", "for(i in 1:N){", "logit(p[i])<-mu[s[i]]+delta[i]  * (1-equals(t[i],b[i]))",
	        "r[i]~dbin(p[i],n[i])", "delta[i] ~ dnorm(md[i],tau)","md[i] <-	d[t[i]] - d[b[i]]}","for(j in 1:NS){ mu[j]~dnorm(0,.0001) }",
		"d[1]<-0","for (k in 2:NT)  {d[k] ~ dnorm(0,.0001) }","sd~dunif(0,2)","tau<-1/pow(sd,2)",output_mA,"for (k in 1:NT) { logit(T[k])<- mA +d[k] }}'", sep ="\n")
		doItAndPrint(Teste)
		command <- paste(modelValue," <- run.jags(model=Model, monitor=c('d','T'), n.chains=", Chain, ", method='rjags', data = list('s'=",ActiveDataSet(),"$s, 't'=",ActiveDataSet(),"$t, 'r'=",ActiveDataSet(),"$r, 'n'=",ActiveDataSet(),"$n, 'b'=",ActiveDataSet(),"$b,'N'=",Arms,",'NS'=",Stud,",'NT'=",Treat,"), burnin =", BurnIn, ",sample =", Sims,",silent.jags=T",")", sep="")
		command_aux <- paste(modelValue,"_aux <- run.jags(model=Model, monitor=c('T'), n.chains=", Chain, ", method='rjags', data = list('s'=",ActiveDataSet(),"$s, 't'=",ActiveDataSet(),"$t, 'r'=",ActiveDataSet(),"$r, 'n'=",ActiveDataSet(),"$n, 'b'=",ActiveDataSet(),"$b,'N'=",Arms,",'NS'=",Stud,",'NT'=",Treat,"), burnin =", BurnIn, ",sample =", Sims,",silent.jags=T",")", sep="")
		command_aux2 <- paste(modelValue,"_aux2 <- run.jags(model=Model, monitor=c('d'), n.chains=", Chain, ", method='rjags', data = list('s'=",ActiveDataSet(),"$s, 't'=",ActiveDataSet(),"$t, 'r'=",ActiveDataSet(),"$r, 'n'=",ActiveDataSet(),"$n, 'b'=",ActiveDataSet(),"$b,'N'=",Arms,",'NS'=",Stud,",'NT'=",Treat,"), burnin =", BurnIn, ",sample =", Sims,",silent.jags=T",")", sep="")
		}

	if(models=="threerem")
		{
		Teste <- paste("Model <- 'model{","sw[1] <- 0","for(i in 1:N)  {","logit(p[i])<-mu[s[i]]+ delta[i]  * (1-equals(t[i],b[i]))",
		"r[i]~dbin(p[i],n[i])","delta[i] ~ dnorm(md[i],taud[i])","taud[i] <- tau * (1 + equals(m[i],3) /3)","md[i] <- d[t[i]] - d[b[i]]  +  equals(m[i],3) * sw[i]}",
		"for (i in 2:N)  {sw[i] <- (delta[i-1] -  d[t[i-1]] + d[b[i-1]])/2}","for(j in 1:NS){ mu[j]~dnorm(0,.0001) }","d[1]<-0",
		"for (k in 2:NT)  {d[k] ~ dnorm(0,.0001) }","sd~dunif(0,2)","tau<-1/pow(sd,2)",output_mA,"for (k in 1:NT)   { logit(T[k])<- mA  +d[k] }}'", sep="\n")
		doItAndPrint(Teste)
		command <- paste(modelValue," <- run.jags(model=Model, monitor=c('d','T'), n.chains=", Chain, ", method='rjags', data = list('s'=",ActiveDataSet(),"$s, 't'=",ActiveDataSet(),"$t, 'r'=",ActiveDataSet(),"$r, 'n'=",ActiveDataSet(),"$n, 'b'=",ActiveDataSet(),"$b, 'm'=",ActiveDataSet(),"$m, 'N'=",Arms,",'NS'=",Stud,",'NT'=",Treat,"), burnin =", BurnIn, ",sample =", Sims,",silent.jags=T",")", sep="")
		command_aux <- paste(modelValue,"_aux <- run.jags(model=Model, monitor=c('T'), n.chains=", Chain, ", method='rjags', data = list('s'=",ActiveDataSet(),"$s, 't'=",ActiveDataSet(),"$t, 'r'=",ActiveDataSet(),"$r, 'n'=",ActiveDataSet(),"$n, 'b'=",ActiveDataSet(),"$b, 'm'=",ActiveDataSet(),"$m, 'N'=",Arms,",'NS'=",Stud,",'NT'=",Treat,"), burnin =", BurnIn, ",sample =", Sims,",silent.jags=T",")", sep="")
		command_aux2 <- paste(modelValue,"_aux2 <- run.jags(model=Model, monitor=c('d'), n.chains=", Chain, ", method='rjags', data = list('s'=",ActiveDataSet(),"$s, 't'=",ActiveDataSet(),"$t, 'r'=",ActiveDataSet(),"$r, 'n'=",ActiveDataSet(),"$n, 'b'=",ActiveDataSet(),"$b, 'm'=",ActiveDataSet(),"$m, 'N'=",Arms,",'NS'=",Stud,",'NT'=",Treat,"), burnin =", BurnIn, ",sample =", Sims,",silent.jags=T",")", sep="")
		}

	if(models=="marem")
		{
		justDoIt(paste('n <- as.matrix(',ActiveDataSet(),'[,substr(colnames(',ActiveDataSet(),'), 1, 2) == "n."])', sep = ""))
		justDoIt(paste('t <- as.matrix(',ActiveDataSet(),'[,substr(colnames(',ActiveDataSet(),'), 1, 2) == "t."])', sep = ""))
		justDoIt(paste('r <- as.matrix(',ActiveDataSet(),'[,substr(colnames(',ActiveDataSet(),'), 1, 2) == "r."])', sep = ""))
		Teste <- paste("Model <- 'model{","for(i in 1:NS){","w[i,1] <-0","delta[i,t[i,1]]<-0","mu[i] ~ dnorm(0,.0001)","for (k in 1:na[i])  {",
	        "r[i,k] ~ dbin(p[i,t[i,k]],n[i,k])","logit(p[i,t[i,k]])<-mu[i] + delta[i,t[i,k]] }","for (k in 2:na[i]) {",
                "delta[i,t[i,k]] ~ dnorm(md[i,t[i,k]],taud[i,t[i,k]])","md[i,t[i,k]] <-  d[t[i,k]] - d[t[i,1]]  + sw[i,k]",
                "taud[i,t[i,k]] <- tau *2*(k-1)/k","w[i,k] <- (delta[i,t[i,k]]  - d[t[i,k]] + d[t[i,1]])","sw[i,k] <-sum(w[i,1:(k-1)])/(k-1) }}", 					 
                "d[1]<-0","for (k in 2:NT){d[k] ~ dnorm(0,.0001) }","sd~dunif(0,2)","tau<-1/pow(sd,2)",output_mA,
		"for (k in 1:NT)   { logit(T[k])<- mA  +d[k] }}'", sep="\n")
		doItAndPrint(Teste)
		command <- paste(modelValue," <- run.jags(model=Model, monitor=c('d','T'), n.chains=", Chain, ", method='rjags', data = list('t'=t, 'r'=r, 'n'=n, 'na'=",ActiveDataSet(),"$na, 'N'=",Arms,",'NS'=",Stud,",'NT'=",Treat,"), burnin =", BurnIn, ",sample =", Sims,",silent.jags=T",")", sep="")
		command_aux <- paste(modelValue,"_aux <- run.jags(model=Model, monitor=c('T'), n.chains=", Chain, ", method='rjags', data = list('t'=t, 'r'=r, 'n'=n, 'na'=",ActiveDataSet(),"$na, 'N'=",Arms,",'NS'=",Stud,",'NT'=",Treat,"), burnin =", BurnIn, ",sample =", Sims,",silent.jags=T",")", sep="")
		command_aux2 <- paste(modelValue,"_aux2 <- run.jags(model=Model, monitor=c('d'), n.chains=", Chain, ", method='rjags', data = list('t'=t, 'r'=r, 'n'=n, 'na'=",ActiveDataSet(),"$na, 'N'=",Arms,",'NS'=",Stud,",'NT'=",Treat,"), burnin =", BurnIn, ",sample =", Sims,",silent.jags=T",")", sep="")
		}

        doItAndPrint(command)
	justDoIt(command_aux)
	justDoIt(command_aux2)
        tkfocus(CommanderWindow())
        	}

    OKCancelHelp(helpSubject="run.jags", reset="RunMTCModel")
    tkgrid(labelRcmdr(modelFrame, text = gettextRcmdr("Enter name for model:")), model, sticky = "w")
    tkgrid(modelFrame, sticky = "w")
    modelFormula()
    tkgrid(tklabel(top, text="Number of arms"), ArmsEntry, sticky="e")
    tkgrid(tklabel(top, text="Number of treatments"), TreatEntry, sticky="e")
    tkgrid(tklabel(top, text="Number of studies"), StudEntry, sticky="e")
    tkgrid(tklabel(top, text="Number of chains"), ChainsEntry, sticky="e")
    tkgrid(tklabel(top, text="Number of simulations"), SimsEntry, sticky="e")
    tkgrid(tklabel(top, text="Number of burn-in sims"), BurnInEntry, sticky="e")
    tkgrid(modelsFrame, sticky="w")    
    tkgrid.configure(ArmsEntry, sticky="w")
    tkgrid.configure(TreatEntry, sticky="w")
    tkgrid.configure(StudEntry, sticky="w")
    tkgrid.configure(ChainsEntry, sticky="w")
    tkgrid.configure(SimsEntry, sticky="w")
    tkgrid.configure(BurnInEntry, sticky="w")    
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()

    }
    
#.activeModel <- ActiveModel()

#------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------

SummaryModel <- function(){
    
    initializeDialog(title=gettextRcmdr("Model summary"))
    modelName <- tclVar(paste("MTCModel.", "1", sep = ""))
    modelFrame <- tkframe(top)    
    model <- ttkentry(modelFrame, width = "20", textvariable = modelName) 
    UpdateModelNumber()

    onOK <- function(){
        modelValue <- trim.blanks(tclvalue(modelName))
        closeDialog()      
	if (!is.valid.name(modelValue)){
      	    errorCondition(recall=SummaryModel, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
      	    return()
    	    } 
	command <- paste('summary(',modelValue,')', sep="")
        doItAndPrint(command)
	command2 <- modelValue
        doItAndPrint(command2)
	tkfocus(CommanderWindow())
	               }
    OKCancelHelp(helpSubject="run.jags", reset="SummaryModel")
    tkgrid(labelRcmdr(modelFrame, text = gettextRcmdr("Enter a name for existing model:")), model, sticky = "w")
    tkgrid(modelFrame, sticky = "w")    
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()
	                   }

#------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------

BestStrat <- function(){

    initializeDialog(title=gettextRcmdr("Check for best strategy"))
    modelName <- tclVar(paste("MTCModel.", "1", sep = ""))
    modelFrame <- tkframe(top)    
    model <- ttkentry(modelFrame, width = "20", textvariable = modelName)
    #defaults <- list(TreatVar="")
    #detailVariable <- tclVar("1")
    #detailCheckBox <- ttkcheckbutton(optionsFrame, variable=detailVariable)
    #TreatVar <- tclVar(dialog.values$TreatVar)
    #TreatEntry <- tkentry(top, width="6", textvariable=TreatVar)
    radioButtons(name="outcome", buttons=c("favor","unfavor"), labels=gettextRcmdr(c("Favorable", "Unfavoravle")), title=gettextRcmdr("Outcome"))     

    onOK <- function(){
        modelValue <- trim.blanks(tclvalue(modelName))
        closeDialog()      
	if (!is.valid.name(modelValue)){
      	    errorCondition(recall=BestStrat, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
      	    return()
    	    }
	outcome <- tclvalue(outcomeVariable)
	
        command_stack <- paste('base <- do.call(rbind.data.frame, ',modelValue,'_aux$mcmc)', sep = "")
	justDoIt(command_stack)

	if (outcome == 'favor'){
	justDoIt(paste('max_base <- apply(base,1,max)',sep = ""))
	justDoIt(paste('p_best <- ifelse(base == max_base, 1, 0)',sep = ""))
	}
	
	if (outcome == 'unfavor'){
	justDoIt(paste('min_base <- apply(base,1,min)',sep = ""))
	justDoIt(paste('p_best <- ifelse(base == min_base, 1, 0)',sep = ""))
	}
	command <- paste('summary(p_best)', sep="")
        doItAndPrint(command)

	#if (tclvalue(detailVariable) == "1")
	#{
	#
	#if (outcome == 'favor'){
	#
	#for (n in 1:Treat-1){
	#for (y in 2:Treat){
	#justDoIt(paste('Diffs <- (teste[,',n,'] - teste[,',y,'])',sep = ""))
	#justDoIt(paste('table(do.call(rbind,lapply(Diffs,function(x) x<=0)))',sep = ""))
	#}
	#}
	#
	#for (y in 1:Treat){
	#	justDoIt(paste('exp(quantile(do.call(rbind,MTCModel.1$mcmc)[,',y,'],c(0.1, 0.25, 0.5, 0.75, 0.9)))', sep = ""))
	#	justDoIt(paste('table(do.call(rbind,lapply(MTCModel.1$mcmc,function(x) x<=0))[,',y,'])', sep = ""))
	#		  }
	#
	#}
	#
	#
	#
	#if (outcome == 'unfavor'){
	#	for (n in 1:Treat-1){
	#		for (y in 2:Treat){
	#			justDoIt(paste('Diffs <- (teste[,',n,'] - teste[,',y,'])',sep = ""))
	#			justDoIt(paste('table(do.call(rbind,lapply(Diffs,function(x) x>=0)))',sep = ""))
	#				  }
	#			    }
	#			 }
	#
	#
	#
	#}
	
	
	tkfocus(CommanderWindow())

	               }
    OKCancelHelp(helpSubject="run.jags", reset="BestStrat")
    tkgrid(labelRcmdr(modelFrame, text = gettextRcmdr("Enter a name for existing model:")), model, sticky = "w")
    tkgrid(labelRcmdr(optionsFrame, text=gettextRcmdr("Detailed analysis:")), detailCheckBox, sticky="w")
    #tkgrid(tklabel(top, text="Number of treatments"), TreatEntry, sticky="e")	
    tkgrid(modelFrame, sticky = "w") 
    tkgrid(outcomeFrame, sticky="w")  
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()
	                }

#------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------

DensityPlots <- function(){
    
    initializeDialog(title=gettextRcmdr("Density plots")) 
    modelName <- tclVar(paste("MTCModel.", "1", sep = ""))
    modelFrame <- tkframe(top)    
    model <- ttkentry(modelFrame, width = "20", textvariable = modelName) 
    UpdateModelNumber()

    onOK <- function(){
        modelValue <- trim.blanks(tclvalue(modelName))
        closeDialog()      
	if (!is.valid.name(modelValue)){
      	    errorCondition(recall=DensityPlots, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
      	    return()
    	    } 
	command <- paste('plot(',modelValue,',layout=c(2,1))', sep="")
        doItAndPrint(command)
	tkfocus(CommanderWindow())
	               }
    OKCancelHelp(helpSubject="run.jags", reset="DensityPlots")
    tkgrid(labelRcmdr(modelFrame, text = gettextRcmdr("Enter a name for existing model:")), model, sticky = "w")
    tkgrid(modelFrame, sticky = "w")    
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()
			}

#------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------

ForestPlot <- function(){
    
    initializeDialog(title=gettextRcmdr("Forest Plot")) 
    modelName <- tclVar(paste("MTCModel.", "1", sep = ""))
    modelFrame <- tkframe(top)    
    model <- ttkentry(modelFrame, width = "20", textvariable = modelName) 
    UpdateModelNumber()

    onOK <- function(){
        modelValue <- trim.blanks(tclvalue(modelName))
        closeDialog()      
	if (!is.valid.name(modelValue)){
      	    errorCondition(recall=ForestPlot, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
      	    return()
    	    } 
	command <- paste('Teste <- summary(',modelValue,'_aux2)', sep="")
        doItAndPrint(command)
	justDoIt(paste('lod <- Teste[,"Mean"]', sep=""))
	justDoIt(paste('sd_lod <- Teste[,"SD"]', sep=""))
	justDoIt(paste('N_Trat <- length(Teste[,"Mean"])', sep=""))

	justDoIt(paste('matriz_lor <- outer(lod,lod,"-")', sep=""))
	justDoIt(paste('LOR <- matriz_lor[lower.tri(matriz_lor)]', sep=""))
	justDoIt(paste('matriz_sd_lor <- outer(sd_lod,sd_lod,function(a,b) {sqrt(a^2+b^2)})', sep=""))
	justDoIt(paste('SD_LOR <- matriz_sd_lor[lower.tri(matriz_sd_lor)]', sep=""))
	justDoIt(paste('OR <- exp(LOR)', sep=""))

	justDoIt(paste('nomes <- paste("OR",col(matriz_lor),row(matriz_lor), sep = "_")', sep=""))
	justDoIt(paste('nomes <- matrix(nomes,ncol=N_Trat,byrow=T,dimnames=list(letters[1:N_Trat],letters[1:N_Trat]))', sep=""))
	justDoIt(paste('nomes <- nomes[lower.tri(nomes)]', sep=""))

	justDoIt(paste('LClr_OR <- exp(LOR - 1.96 * (SD_LOR))', sep=""))
	justDoIt(paste('UClr_OR <- exp(LOR + 1.96 * (SD_LOR))', sep=""))
	justDoIt(paste('summary_OR <- t(rbind(OR, LClr_OR, UClr_OR))', sep=""))

	justDoIt(paste('m <- c(OR)', sep=""))
	justDoIt(paste('l <- c(LClr_OR)', sep=""))
	justDoIt(paste('u <- c(UClr_OR)', sep=""))

	justDoIt(paste('tabletext<-cbind(c("Comparisons",nomes),c("OR",format(m,digits=2)),c("LClr_OR",format(l,digits=2)),c("UClr_OR",format(u,digits=2)))', sep=""))
	
	command <- paste('print(summary_OR)', sep="")
        doItAndPrint(command)
	command_fplot <- paste('forestplot(tabletext,c(NA,OR),c(NA,LClr_OR),c(NA,UClr_OR),zero=1,,clip=c(0,5))', sep = "")
	doItAndPrint(command_fplot)

	tkfocus(CommanderWindow())
	               }
    OKCancelHelp(helpSubject="forestplot", reset="ForestPlot")
    tkgrid(labelRcmdr(modelFrame, text = gettextRcmdr("Enter a name for existing model:")), model, sticky = "w")
    tkgrid(modelFrame, sticky = "w")    
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()
			}

#------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------

Relationship <- function(){
    Library('igraph')
    initializeDialog(title=gettextRcmdr("Relationship graph"))
    defaults <- list(models="fem")
    dialog.values <- getDialog("Relationship", defaults)
    radioButtons(name="models", buttons=c("fem", "rem", "threerem", "marem"), labels=gettextRcmdr(c("FE Model", "RE Model ignoring multi-arms effect", "RE Model for 2- and 3-arms trial", "RE Model for multi-arms trial")), title=gettextRcmdr("Models"))

	onOK <- function(){
        closeDialog()

	models <- tclvalue(modelsVariable)
	justDoIt(paste('Size <- ncol(t)', sep = ""))

	if(models=="marem")
	{
	final <- 0
	for (x in 2:Size-1)
	{
	justDoIt(paste('temp <- as.data.frame(cbind(t[,',x,'],t[,',x+1,']))', sep = ""))
	final <- rbind(final, temp)
	}
	plota <- cbind(final$V1, final$V2) 
	plota <- as.data.frame(plota) 
	plota <- plota[plota$V1 != plota$V2,] 
	plota <- plota[complete.cases(plota$V1),] 
	g <- graph.data.frame(plota, directed=FALSE) 
	summary(g) 
	g$layout <- layout.fruchterman.reingold(g) 
	plot(g)        
        }

	if(models== "fem" | models== "rem" | models== "threerem")
	{
	command <- paste("plota <- cbind(",ActiveDataSet(),"$b, ",ActiveDataSet(),"$t)", sep="")
	doItAndPrint(command)
	plota <- as.data.frame(plota) 
	plota <- plota[plota$V1 != plota$V2,] 
	g <- graph.data.frame(plota, directed=FALSE) 
	summary(g)
	g$layout <- layout.fruchterman.reingold(g) 
	plot(g)
        }
	tkfocus(CommanderWindow())
	}

    OKCancelHelp(helpSubject="layout.fruchterman.reingold", reset="Relationship")
    tkgrid(modelsFrame, sticky="w")
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()
    }

#------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------

GelmanRubinPlot <- function(){
    
    initializeDialog(title=gettextRcmdr("Gelman-Rubin Plot")) 
    modelName <- tclVar(paste("MTCModel.", "1", sep = ""))
    modelFrame <- tkframe(top)    
    model <- ttkentry(modelFrame, width = "20", textvariable = modelName) 
    UpdateModelNumber()

    onOK <- function(){
        modelValue <- trim.blanks(tclvalue(modelName))
        closeDialog()      
	if (!is.valid.name(modelValue)){
      	    errorCondition(recall=GelmanRubinPlot, message=sprintf(gettextRcmdr('"%s" is not a valid name.'), modelValue), model=TRUE)
      	    return()
    	    } 
	command <- paste('Teste <- summary(',modelValue,'_aux)', sep="")
        justDoIt(command)
	justDoIt(paste('N_Trat <- length(Teste[,"Mean"])', sep=""))
	justDoIt(paste('Teste_d <- ',modelValue,'_aux2$mcmc', sep=""))
	doItAndPrint(paste('gelman.plot(Teste_d[,2:N_Trat])', sep=""))
	justDoIt(paste('x11()', sep=""))
	justDoIt(paste('Teste_T <- ',modelValue,'_aux$mcmc', sep=""))
	doItAndPrint(paste('gelman.plot(Teste_T[,1:N_Trat])', sep=""))

	tkfocus(CommanderWindow())
	               }
    OKCancelHelp(helpSubject="gelman.plot", reset="GelmanRubinPlot")
    tkgrid(labelRcmdr(modelFrame, text = gettextRcmdr("Enter a name for existing model:")), model, sticky = "w")
    tkgrid(modelFrame, sticky = "w")    
    tkgrid(buttonsFrame, sticky="w")
    dialogSuffix()
			}

#------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------//------------
