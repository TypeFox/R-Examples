	setGeneric("initModelClass", function(model)
	standardGeneric("initModelClass"))
	
	setGeneric("residPart", function(model, group, 
	multimodel, thetalist, clpindepX, finished, returnX, rawtheta) 
	standardGeneric("residPart"))

	setGeneric("getClpindepX", function(model, multimodel, theta,
        returnX, rawtheta, dind)
        standardGeneric("getClpindepX"))
	
	setGeneric("plotter", function(model, multimodel, multitheta,  
	plotoptions)
        standardGeneric("plotter"))

