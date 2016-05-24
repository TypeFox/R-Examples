mutossGUIVars <- local({
			outputCon <- NULL
			errorCon <- NULL
			list(getOutput=function() outputCon, setOutput=function(value) outputCon <<- value,
					getError=function() errorCon, setError=function(value) errorCon <<- value)
		})

startRecording <- function() {
	mutossGUIVars$setOutput(textConnection(".mutossGUIoutput", open="w"))
	mutossGUIVars$setError(textConnection(".mutossGUIerrorMsg", open="w"))
	sink(mutossGUIVars$getOutput())
	sink(mutossGUIVars$getError(), type="message")
}

stopRecording <- function() {
	sink()
	sink(type="message")
	#withCallingHandlers( {			
	close(mutossGUIVars$getOutput())
	close(mutossGUIVars$getError())
	#}, error = function(e) recover())
	return(list(output=getOutput(), errorMsg=getErrorMsg()))
}

myContrMat <- function(type,l,df,group) {
	n <- table(df[,group])[as.numeric(factor(l,levels=levels(df[,group])))]
	x <- contrMat(n=n,type=type)
}

myContrMat <- function(type,l,df,group) {
	n <- table(df[,group])[as.numeric(factor(l,levels=levels(df[,group])))]
	x <- contrMat(n=n,type=type)
}

getOutput <- function() {
	#return(textConnectionValue(mutossGUI.vars$outputCon))
	return(get(".mutossGUIoutput"))
}

getErrorMsg <- function() {
	#return(textConnectionValue(mutossGUI.vars$errorCon))
	return(get(".mutossGUIerrorMsg"))
}

showRejected <- function(obj) {
	obj <- get(obj)
	print(which(obj@rejected))
}