`startingFunction` <-
function(startingValues, file="getInits.R") {
    

    indent = encodeString(" ", width=7)
    
    
    sink(file)
    cat("getInits = function() { \n\n")
    cat("scale = 1.5\n")
    cat("SDscale = 4\n\n")
    
        
#    cat("startingValues =", as.character(sys.calls()[[1]])[2], "\n\n")
   
    fixed = names(startingValues)
    fixed = fixed[! fixed %in% c("vars","phi")]

    random = paste("R", names(startingValues$vars),sep="")
    fixed = fixed[! fixed %in% random]

    cat("result = list()\n\n")
    for(Dbeta in fixed) {
     cat("result[[\"", Dbeta, "\"]] = sign(startingValues[[\"",
      Dbeta, "\" ]]) *\n",  encodeString(" ", width=4), 
      "runif(length(startingValues[[\"",
      Dbeta, "\" ]]),\n", indent, "abs(startingValues[[\"", Dbeta, 
      "\"]])/scale,\n", indent, "scale * abs(startingValues[[\"", 
      Dbeta, "\"]]))\n\n", sep="")
    }
    cat("\n")    

     
  for(Dvar in names(startingValues$vars)) {
    cat("result[[\"", paste("SD", Dvar, sep=""), 
      "\"]] = sqrt(runif(1,\n", indent, "startingValues$vars[[\"", Dvar, 
        "\"]]/scale,\n", indent, "startingValues$vars[[\"",
          Dvar, "\"]]*scale))\n\n", sep="")

     theR =  paste("R", Dvar, sep="")
     theSV = paste("startingValues[[\"", theR, "\"]]", sep="")
    cat("result[[\"", theR, 
      "\"]] = rnorm(length(", theSV, "),\n ", indent, 
      theSV, ", startingValues$vars[[\"", Dvar, "\"]]/SDscale)\n\n",sep="")
   
	# if it's a spatial random effect, make it sum to zero
	if(length(grep("Spatial$", Dvar))) {
		cat("
		if(all(abs(result[[ \"", theR, "\" ]])<0.1)){
			result[[ \"", theR, "\" ]] = rep(0, 
					length(result[[ \"", theR, "\" ]] ))
		} else {	
			result[[ \"", theR, "\" ]] = result[[ \"", theR, "\" ]] -
					mean(result[[ \"", theR, "\" ]])
		}",sep="")
	}

  }
 
  # range parameters for geostatiatical models
  for(Dvar in names(startingValues$phi)) {
      cat("result[[\"", paste("phi", Dvar, sep=""), 
        	  "\"]] = runif(1,\n", indent, "startingValues$phi[[\"", Dvar, 
        	  "\"]]/scale,\n", indent, "startingValues$phi[[\"",
              Dvar, "\"]]*scale)\n\n", sep="")
  } 
  
  
  
  cat("\nreturn(result)\n")
  cat("\n}")
  sink()

}

