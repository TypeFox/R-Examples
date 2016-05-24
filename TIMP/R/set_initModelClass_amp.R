  setMethod("initModelClass", signature(model="amp"), 
    function (model) 
    {
      model@ncomp <- length(model@amps)
      if(length(model@nl)==0) {
        model@ncolc <- array(model@ncomp, 1)
	} else {
	model@ncolc <- array(model@ncomp, model@nl)
	}    
      model
   }) 

