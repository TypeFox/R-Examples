  setMethod("getClpindepX", signature(model="spec"), 
    function(model, multimodel, theta, returnX, rawtheta, dind) 
    {
       if(returnX) 
		 theta <- getThetaCl(rawtheta, multimodel)[[dind]]
       x <- specModel(theta@specpar, model) 
       
        if(returnX) 
		    x <- as.vector(x) 
       x
   }) 
