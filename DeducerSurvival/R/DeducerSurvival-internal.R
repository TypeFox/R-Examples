
######################################################################################
######################################################################################
#
# DeducerSurvival survival curves and tests for Deducer Ockendon & Cool
# Current Version 1.3
# 07.07.2012
#
######################################################################################
######################################################################################


.KMCheckFunction <- function(state){
	#make sure at least one variables is selected for each of FollowUp and Event
	if(length(state$FollowUp)<1){
		return("Please select a 'Follow Up' variable")}
	if(length(state$FollowUp)>1){
		return("Please select only one 'Follow Up' variable")}
	if(length(state$Event)<1){
		return("Please select an 'Event' variable")}
	if(length(state$Event)>1){
		return("Please select only one 'Event' variable")}
	if(length(state$Group)>1){
		return("Please select only one 'Group' variable")}
	if(state$Confidence == "Custom" & state$Custom_Confidence >1){
		return("Please enter a confidence level between 0 and 1")}
	if(state$Confidence == "Custom" & state$Custom_Confidence <0){
		return("Please enter a confidence level between 0 and 1")}
	return("")
}



.KMRunFunction <- function(state){

#make formula for Survival Object
cmd <- paste ("mySurvival = with(",state$data)
	cmd <- paste(cmd,",")
	cmd <- paste(cmd,"Surv(")
	cmd <- paste(cmd,state$FollowUp)
	cmd <- paste(cmd,",")
	cmd <- paste(cmd,state$Event)
	cmd <- paste(cmd,"))\n")

#present summary of survival object
	cmd <- paste(cmd,"summary(mySurvival)\n")
	
#fit a and plot KM curve 
# enter here the colour selection from the dialogue for con interval and colour
	conf <- state$Confidence
	if (conf == "Custom"){conf <- state$Custom_Confidence}
	colour1 <- state$Colour1
	colour2 <- state$Colour2
	colour3 <- state$Colour3


if(length(state$Group)<1){
	cmd <- paste(cmd,"mycurve=survfit(mySurvival~1,conf.int=")
	
	cmd <- paste(cmd,conf)
	
	cmd <- paste(cmd,")\n")

}else if (length(state$Group)>0){
	cmd <- paste(cmd,"mycurve=with(",state$data)
	
	cmd <- paste(cmd,",")
	
	cmd <- paste(cmd,"survfit(mySurvival~")
	cmd <- paste(cmd,state$Group)
	cmd <- paste(cmd,",conf.int=")
	cmd <- paste(cmd,conf)
	
	cmd <- paste(cmd,"))\n")
}
	
	cmd <- paste(cmd,"summary(mycurve)\n")
	
	cmd <- paste(cmd,"plot(mycurve, col=c(")
	dblquote <- "\""
	cmd <- paste(cmd,dblquote)
	cmd <- paste(cmd,colour1)
	cmd <- paste(cmd,dblquote)
	cmd <- paste(cmd,",")
	cmd <- paste(cmd,dblquote)
	cmd <- paste(cmd,colour2)
	cmd <- paste(cmd,dblquote)
	cmd <- paste(cmd,",")
	cmd <- paste(cmd,dblquote)
	cmd <- paste(cmd,colour3)
	cmd <- paste(cmd,dblquote)
	cmd <- paste(cmd,"))\n")

	title <- "\"Kaplan Meier Curve\""
	xlab <- "\"Follow Up Duration\""
	ylab <- "\"Probability\""
	cmd <- paste(cmd,"title(main=")
	cmd <- paste(cmd,title)
	cmd <- paste(cmd,",xlab=")
	cmd <- paste(cmd,xlab)
	cmd <- paste(cmd,",ylab=")
	cmd <- paste(cmd,ylab,")\n")
	


if(length(state$Group)>0){
		
	cmd <- paste(cmd,"with(",state$data)
	
	cmd <- paste(cmd,",")
	
	cmd <- paste(cmd,"legend(0,0.3,unique(")
	
	cmd <- paste(cmd, state$Group)
	
	cmd <- paste(cmd, "),lty=c(1,1,1),col=c(")
	
	cmd <- paste(cmd,dblquote)
	cmd <- paste(cmd,colour1)
	cmd <- paste(cmd,dblquote)
	cmd <- paste(cmd,",")
	cmd <- paste(cmd,dblquote)
	cmd <- paste(cmd,colour2)
	cmd <- paste(cmd,dblquote)
	cmd <- paste(cmd,",")
	cmd <- paste(cmd,dblquote)
	cmd <- paste(cmd,colour3)
	cmd <- paste(cmd,dblquote)
	
	cmd <- paste(cmd,")))\n")}

execute(cmd)
}



makeKMDialog <- function (){

KMDialog <- new(SimpleRDialog)
KMDialog$setTitle("KM Survival")
KMDialog$setSize(600L,500L)
#add variable selector
	variableSelector <- new(VariableSelectorWidget)
	variableSelector$setTitle("data")
	addComponent(KMDialog,variableSelector,10,200,850,10)
#add a list for the FollowUp variable
	FollowUpList<- new(SingleVariableWidget,variableSelector)
	FollowUpList$setTitle("FollowUp",TRUE)
	addComponent(KMDialog, FollowUpList,10,500,260, 200)
#add a list for the Event variable
	EventList<- new(SingleVariableWidget,variableSelector)
	EventList$setTitle("Event",TRUE)
	addComponent(KMDialog, EventList,290,500,540, 200)
#add a list for the Group variable
	GroupList<- new(SingleVariableWidget,variableSelector)
	GroupList$setTitle("Group",TRUE)
	addComponent(KMDialog, GroupList,570,500,820, 200)
#add three text boxes to set the colour of up to three survival curves
	colour1Area <- new(TextAreaWidget,"Colour1")	
	addComponent(KMDialog, colour1Area,10,990,120, 600)
	colour1Area$setDefaultModel(c("black"))
	colour2Area <- new(TextAreaWidget,"Colour2")	
	addComponent(KMDialog, colour2Area,150,990,260, 600)
	colour2Area$setDefaultModel(c("black"))
	colour3Area <- new(TextAreaWidget,"Colour3")
	addComponent(KMDialog, colour3Area,290,990,400, 600)
	colour3Area$setDefaultModel(c("black"))
#add a variable for the Confidence interval
	confBoxes <- new(ButtonGroupWidget,"Confidence",c("0","0.95","0.99","Custom"))
	addComponent(KMDialog, confBoxes,430,990,680, 600)
	confBoxes$setDefaultModel(c("0.95"))
	confCustom <- new(TextFieldWidget,"Custom_Confidence")
	addComponent(KMDialog, confCustom,710,990,800,600)
#output options
	KMDialog$setCheckFunction(toJava(.KMCheckFunction))
	KMDialog$setRunFunction(toJava(.KMRunFunction))
	
	return(KMDialog)
}




.CompareCheckFunction <- function(state){
	#make sure at least one variables is selected for each: FollowUp, Event and Group
	if(length(state$FollowUp)<1){
		return("Please select a 'Follow Up' variable")}
	if(length(state$FollowUp)>1){
		return("Please select only one 'Follow Up' variable")}
	if(length(state$Event)<1){
		return("Please select an 'Event' variable")}
	
	if(length(state$Event)>1){
		return("Please select only one 'Event' variable")}
	
	if(length(state$Group)<1){
		return("Please select a 'Group' variable")}
	if(length(state$Group)>1){
		return("Please select only one 'Group' variable")}

	return("")
}


.CompareRunFunction <- function(state){
#make formula for Survival Object
	cmd <- paste ("mySurvival = with(",state$data)
	cmd <- paste(cmd,",")
	cmd <- paste(cmd,"Surv(")
	cmd <- paste(cmd,state$FollowUp)
	cmd <- paste(cmd,",")
	cmd <- paste(cmd,state$Event)
	cmd <- paste(cmd,"))\n")
	

# set parameters for formula below:
	paratest <- state$Parametric
	
	dblquote <- "\""


	nonparaoption <- state$NonParametric

if (nonparaoption == "Run"){

# Perform the logrank test on the mySurvival object:
	cmd <- paste(cmd,"mylogrank=with(",state$data)
	
	cmd <- paste(cmd,",")
	
	cmd <- paste(cmd,"survdiff(mySurvival~")
	cmd <- paste(cmd,state$Group)
	
	cmd <- paste(cmd,"))\n")
	
	cmd <- paste(cmd,"mylogrank\n")



# Perform the cox proportional hazards on the mySurvival object:

	cmd <- paste(cmd,"mycoxph=with(",state$data)
	
	cmd <- paste(cmd,",")
	
	cmd <- paste(cmd,"coxph(mySurvival~")
	cmd <- paste(cmd,state$Group)
	
	cmd <- paste(cmd,"))\n")
	
	cmd <- paste(cmd,"mycoxph\n")
	
}


# Perform parametric regression analysis on the mySurvival object:
	cmd <- paste(cmd,"myparametric=with(",state$data)
	
	cmd <- paste(cmd,",")
	
	cmd <- paste(cmd,"survreg(mySurvival~")
	cmd <- paste(cmd,state$Group)
	
	cmd <- paste(cmd,",")
	
	cmd <- paste(cmd,"dist =")
	
	cmd <- paste(cmd,dblquote)
	
	cmd <- paste(cmd,paratest,sep='')
	
	cmd <- paste(cmd,dblquote,sep='')
	
	cmd <- paste(cmd,"))\n")
	
	cmd <- paste(cmd,"myparametric\n")
execute(cmd)
}




makeCompareDialog <- function (){

Compare <- new(SimpleRDialog)
Compare$setTitle("Compare Survival Curves")
Compare$setSize(600L,500L)
#add variable selector
	variableSelector <- new(VariableSelectorWidget)
	variableSelector$setTitle("data")
	addComponent(Compare,variableSelector,10,200,850,10)
#add a list for the FollowUp variables
	FollowUpList<- new(SingleVariableWidget,variableSelector)
	FollowUpList$setTitle("FollowUp",TRUE)
	addComponent(Compare, FollowUpList,10,500,260, 200)
#add a list for the Event variables
	EventList<- new(SingleVariableWidget,variableSelector)
	EventList$setTitle("Event",TRUE)
	addComponent(Compare, EventList,290,500,540, 200)

#add a list for the Group variable
	GroupList<- new(SingleVariableWidget,variableSelector)
	GroupList$setTitle("Group",TRUE)
	addComponent(Compare, GroupList,570,500,820, 200)

#add a selection parameter for the parametric regression analysis
	paraBoxes <- new(ButtonGroupWidget,"Parametric",c("extreme","logistic","gaussian","weibull","exponential","rayleigh","loggaussian","lognormal","loglogistic","t"))
	addComponent(Compare, paraBoxes,10,990,680, 600)
	paraBoxes$setDefaultModel(c("weibull"))

#add option for non-parametric tests
	NonParaBoxes <-new(ButtonGroupWidget,"NonParametric",c("Run","None"))
	addComponent(Compare, NonParaBoxes,700 ,990,900 ,600)
	#NonParaBoxes$setDefaultModel(c("Run"))

#output options
	Compare$setCheckFunction(toJava(.CompareCheckFunction))
	Compare$setRunFunction(toJava(.CompareRunFunction))
	
	return(Compare)
}


.onLoad <- function(libname, pkgname) {

#if deducer gui is not running, do minimal load
	deducerLoaded <- try(.deducer == .jnull(),silent=TRUE)
	if(inherits(deducerLoaded,"try-error") || deducerLoaded)
		return(NULL)

.jpackage(pkgname,lib.loc=libname)

.registerDialog("Kaplan Meier", makeKMDialog)
.registerDialog("Compare", makeCompareDialog)


deducer.addMenu("Survival")
deducer.addMenuItem("Kaplan Meier",,".getDialog('Kaplan Meier')$run()","Survival")
deducer.addMenuItem("Compare",,".getDialog('Compare')$run()","Survival")


if(.windowsGUI){
	winMenuAdd("Survival")
	winMenuAddItem("Survival", "Kaplan Meier", "deducer('Kaplan Meier')")
	winMenuAddItem("Survival", "Compare", "deducer('Compare')")

}else if(.jgr){
	#jgr.addMenu("Survival")
	jgr.insertMenu("Survival",7)
	jgr.addMenuItem("Survival", "Kaplan Meier", "deducer('Kaplan Meier')")
	jgr.addMenuItem("Survival", "Compare", "deducer('Compare')")
}
	

}

