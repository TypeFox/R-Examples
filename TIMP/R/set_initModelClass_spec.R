  setMethod("initModelClass", signature(model="spec"), 
    function (model) 
    {
	model@specdisp <- length(model@specdisppar) != 0 
	model@timedep <- model@weight || model@lclp0 || model@specdisp
	model@clpdep <- model@timedep
        model@ncomp <- length(model@specpar)
	if(length(model@nt)==0) {
        model@ncole <- array(model@ncomp, 1)
	} else {
	model@ncole <- array(model@ncomp, model@nt)
	}
	model
   }) 
