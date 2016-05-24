  setMethod("initModelClass", signature(model="mass"), 
    function (model) 
    {
	model@clpdep <- (model@weight || model@lclp0 || model@lclpequ)
	
	model@ncomp <- length(model@peakpar) + if(model@extracomp) 1 else 0
	model@ncolc <- model@ncomp	
	model
   }) 

