plotOne = function(mat, main=NULL) {
	matplot(mat, lty=1, type="l",  ylab=main)
}

checkChain = function(chain, parameters=NULL, oneFigure=TRUE) {

if(is.array(chain)) {
  chain = list(beta=chain)
}

if(is.null(parameters)) {
	parameters = grep("^sd|range|intercept|betas",
			names(chain),value=TRUE,ignore.case=TRUE)
}
  # find out the number of parameters
  	scalars = unlist(lapply(chain[parameters], is.matrix))
	notScalars = parameters[!scalars]
	scalars = parameters[scalars]
	betas = NULL
	if(length(notScalars) & any(names(chain)=="betas"))  {
		betas = notScalars[
				notScalars %in% dimnames(chain$betas)[[3]] 
						]
	}
	parArray = notScalars[
			unlist(lapply(chain[notScalars], 
					function(qq) length(dim(qq))==3)
		)
	]
	

	
	if(oneFigure) {
		Nplots = length(scalars) + length(betas)

		for(Darray in parArray)
			Nplots = Nplots + dim(chain[[Darray]])[3]
	
		if(Nplots==0) warning("Nothing to plot")
		if(Nplots > 16) warning("Creating", Nplots, "plots, this might not work")
		
	  par(mfrow=c(ceiling(Nplots/4),min(c(4, Nplots))),
			  mar=c(1.5,2.5,0,0), mgp=c(1.5, 0.5, 0))
  }
  
  for(D in scalars)
	  plotOne(chain[[D]], D)

  for(D in betas)
	  plotOne(chain$betas[,,D], D)
  
  for(D in parArray){
	  Dhere = chain[[D]]
	  Shere = dimnames(Dhere)[[3]]
	  if(is.null(Shere)) Shere = seq(1, dim(Dhere)[3])
	  for(D2 in Shere) {
		  plotOne(Dhere[,,D2], paste(D, ":", D2,sep=""))		 
	  }
  }
	  
	  
  invisible() 
}

