`exclude.influence` <-  
function(model, grouping=NULL, level=NULL, obs=NULL, gf="single", delete=TRUE)
{
	data.adapted <- model.frame(model)
	added.variables <- character()
	ranef <- NA
	rm(ranef)

	####
	# Code kindly provided by Jennifer Bufford
	if("(weights)" %in% names(data.adapted)) {
    names(data.adapted)[names(data.adapted)=="(weights)"] <-
    as.character(model@call$weights)}
    if("(offset)" %in% names(data.adapted)) {
    names(data.adapted)[names(data.adapted)=="(offset)"] <-
    as.character(model@call$offset)}
    if(sum(grepl("offset", names(data.adapted)))>0) {
    names(data.adapted)[grep("offset", names(data.adapted))] <-
    gsub('offset\\(|\\)',"",names(data.adapted)[grep("offset", names(data.adapted))])}
	####
	
	
	if(!is.null(obs))
	{
    
    if(!is.null(grouping) | !is.null(level))
    {
      warning("Specification of the 'obs' parameter overrules specification of the 'grouping' and 'level' parameters.")
    }
    
	  data.adapted <- data.adapted[-obs,]
	  # For some reason, update() can not find the data.adapted within the function without the line below. 
	  data.update <- data.adapted		
	  # It can find data.update
	  
	  model.updated <- update(model, data=data.update)
	  return(model.updated)
	}
  
  
  
  
  
  
	if(delete==TRUE)
	{
    
    ## Only works when length(level) == 1, this needs to be enhanced
	 	group.var <- which(names(data.adapted) == grouping)
	 	
		for (i in 1:length(level))
			{
			data.adapted <- subset(data.adapted, data.adapted[,group.var]!=level[i])
     }
    
	 	# For some reason, update() can not find the data.adapted within the function without the line below. 
	 	data.update <- data.adapted		
	 	# It can find data.update
	 	
	 	model.updated <- update(model, data=data.update)
	 	
		return(model.updated)
    
	 }
	
  
  
  
  
  
	if(names(data.adapted)[2] != "intercept.alt")
		{

		data.adapted$intercept.alt <- ifelse(model@flist[[grouping]]==level[1], 0, 1)

		data.adapted[, ncol(data.adapted)+1] <- 
			ifelse(model@flist[[grouping]]==level[1], 1, 0)
		
		added.variables <- make.names(paste("estex.", as.character(level[1]), sep=""))
		colnames(data.adapted)[ncol(data.adapted)] <- added.variables
		

		if(length(level) > 1)
			{
			for (i in 2:length(level))
  				{
				
				data.adapted$intercept.alt[model@flist[[grouping]]==level[i]] <- 0
				
				data.adapted[, ncol(data.adapted)+1] <- 
					ifelse(model@flist[[grouping]]==level[i], 1, 0)
		
				added.variables <- append(added.variables, values = make.names(paste("estex.", as.character(level[i]), sep="")))
				
				colnames(data.adapted)[ncol(data.adapted)] <- added.variables[length(added.variables)]

				}
			}
		
			if(gf=="single")
				{
				# grnr refers to "grouping number"
				grnr <- which(names(ranef(model))==grouping)
		
				if (length(names(ranef(model)[[grnr]])) == 1)
					{
						model.updated <- update(model, 
						   formula = as.formula(paste(". ~ 0 + intercept.alt +", 
						   paste(added.variables, collapse="+"), 	
						   "+ .",
						   "- (1 |", grouping, ") + (0 + intercept.alt |", grouping, ")")),
							data = data.adapted)

					}
		
				if (length(names(ranef(model)[[grnr]])) > 1)
					{
						model.updated <- update(model, 
						   formula = as.formula(paste(". ~ 0 + intercept.alt + ", 
						   paste(added.variables, collapse="+"),
						   " + .",
						   paste(" - (", paste(names(ranef(model)[[grnr]])[-1], collapse="+"), "|", grouping, ")"),
					   		" + (0 + intercept.alt +", paste(names(ranef(model)[[grnr]])[-1], collapse="+"), "|", grouping, ")")),
							data = data.adapted)
					}
				}
			
			if(gf=="all")
				{
					delete.gf <- vector()
				  	for (i in 1:length(ranef(model)))
						{
						if(length(names(ranef(model)[[i]])) > 1)
							{
							delete.gf[i] <- paste(
								"- (", 
								paste(names(ranef(model)[[i]][-1]), collapse="+"), 
								"|", 
								names(ranef(model))[i],
								")")
							}
							
						if(length(names(ranef(model)[[i]])) == 1)
							{
							delete.gf[i] <- paste(
								"- ( 1 |", 
								names(ranef(model))[i],
								")")
							}	
					  	}  
				  	delete.gf <- paste(delete.gf, collapse=" ")
				  
				  	new.gf <- vector()
				  	for (i in 1:length(ranef(model)))
					  {
						if(length(names(ranef(model)[[i]])) > 1)
							{
						  	new.gf[i] <- paste(
							  	"+ (0 + intercept.alt +", 
								paste(names(ranef(model)[[i]][-1]), collapse="+"), 
								"|", 
								names(ranef(model))[i],
								")")
							}
							
						if(length(names(ranef(model)[[i]])) == 1)
							{
						  	new.gf[i] <- paste(
							  	"+ (0 + intercept.alt |", 
								names(ranef(model))[i],
								")")
							}
					  } 
				  new.gf <- paste(new.gf, collapse=" ")
				  
				model.updated <- update(model,
					formula = as.formula(
						paste(
							". ~ 0 + intercept.alt + ",
							paste(added.variables, collapse="+"),
							"+ . ",
							delete.gf,
							new.gf)),
					data=data.adapted)
				
				}
		}
		
		
	if(names(data.adapted)[2] == "intercept.alt")
		{
  		
  		for (i in 1:length(level))
  			{
				data.adapted$intercept.alt[model@flist[[grouping]]==level[i]] <- 0
				
				data.adapted[, ncol(data.adapted)+1] <- 
					ifelse(model@flist[[grouping]]==level[i], 1, 0)
		
				added.variables <- append(added.variables, values = make.names(paste("estex.", as.character(level[i]), sep="")))
				
				colnames(data.adapted)[ncol(data.adapted)] <- added.variables[length(added.variables)]
			}

		model.updated <- update(model, 
			formula = as.formula(paste(
				". ~ 0 + intercept.alt + ", 
				paste(added.variables, collapse="+"),
				"+ .")),
			data = data.adapted)
		}		

	return(model.updated)
	
}