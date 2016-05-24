########################## PREVISIONES.R ############################

.TSTutorial=function(series,student,report,contRep){

	if(missing(student)) student=T
	if(missing(report)) report=list()
	if(missing(contRep)) contRep=list()
   	report  <- do.call(reportCntrl, report)
	contRep  <- do.call(contRepCntrl, contRep)
	if(is.null(contRep$name))	contRep$name="Series"

	
	if(class(series)!="ts"){
		stop("'series' has to be an object of class 'ts'")
	}else{
		series=ts(data=as.numeric(unlist(series)),start=attr(series,"tsp")[1],end=attr(series,"tsp")[2],frequency=attr(series,"tsp")[3])
	}
	
	session=initializer(object=series,student=student,report=report,contRep=contRep)
	while(getcami(session@cami)!="salir"){
		session=activate(session)
	}
	graphics.off()
	if(session["report"]["report"]){
		ender(session)
	}	
	return(invisible())
}
setGeneric("TSTutorial", function(series, student, report, contRep){ standardGeneric("TSTutorial")})
setMethod(f="TSTutorial",signature=c("ts","ANY","ANY","ANY"),definition=.TSTutorial)