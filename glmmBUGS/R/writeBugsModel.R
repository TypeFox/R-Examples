`writeBugsModel` <-
function(file, effects, covariates, observations, 
  family=c("bernoulli", "binomial", "poisson", "normal",  "other"),
  spatial=NULL, geostat=FALSE, prefix="", reparam=NULL, brugs=TRUE, priors=NULL) {

# spatial is a character string of names of random effects
 if(!length(reparam)) {
 interceptString = "intercept"
} else {
interceptString = "interceptUnparam"
}

if (brugs) inprod="inprod" else inprod="inprod2"


  if(!is.character(family)) {
    warning("family must be a character string, ie \"poisson\" or \"binomial\" ")
  }
  family = family[1]
  if(length(spatial)) {
    if(!all(spatial %in% effects))
      warning("spatial effects used which are not specified as random effects")
  }
  
  offset = observations[-1]
  observations = observations[1]
  
  if(family=="poisson") {
    link = "log"
    ddist = "dpois"
  } else if(family=="bernoulli") {
    link="logit"
    ddist="dbern"
  } else if(family=="binomial") {
    link="logit"
    ddist="dbin"
  } else if(family %in% c("normal", "gaussian")) {
    link=""
    ddist="dnorm"
  } else {
    ddist = paste("d", family, sep="")
    link=""
  }
  if(link != "") {
    link = paste(link, "(", sep="")
    endlink = ")"
  } else {
    endlink=""
  }
  
  
if(!is.null(file)) {
  sink(file)
  
  cat("model{\n\n")
}  

  # the first effect
  Deffect=1
  indent = encodeString(" ", width=2*Deffect)
  theE = paste(effects[Deffect], sep="")
  theD =  paste( "D", theE, sep="") 
  cat("for(", theD, " in 1:N", theE, ") {\n\n", sep="")
  if(!geostat | (!(effects[Deffect] %in% spatial)) ) {
    cat(indent, "R", effects[Deffect], "[", theD, "] ~ dnorm(mean", 
     theE, "[", theD, "], T", theE, ")\n",sep="")
  }
 cat(indent, "mean", theE, "[", theD, "] <- ", interceptString, prefix,  sep="")
  
  # the covariates
  # check to see if there's more than one
  if(length(covariates[[theE]])==1) {
    cat(" + beta", theE, " * X", theE, "[", theD, "]", sep="")
  } else if (length(covariates[[theE]]) > 1) {
    cat(" +",  inprod, "(beta", theE, "[] , X", theE, "[", theD, ",])", sep="")
    
  }   
  # spatial
  if(theE %in% spatial & !geostat) {
     cat("+ R", theE, "Spatial[Sspatial", theE, "[", theD, "]]", sep="")
  }
  cat("\n")
  
  # subsequent effects, if any
  if(length(effects)>1) {
    for(Deffect in seq(2, length(effects))) {
      theE = effects[Deffect]
      theD = paste("D", theE, sep="")
      thePastD = paste("D", effects[Deffect-1], sep="")
      thePastS = paste("S", effects[Deffect-1], sep="")
      cat("\n", indent, "for(", theD, " in ", thePastS, "[",
        thePastD, "]:(", thePastS, "[", thePastD, "+1]-1)){\n\n", sep="")
      indent = encodeString(" ", width=2*Deffect)
      if(!geostat | (!(effects[Deffect] %in% spatial)) ) {
        cat(indent, "R", theE, "[", theD, "] ~ dnorm(mean", 
          theE, "[", theD, "], T", theE, ")\n",sep="")
	  }
      cat(indent, "mean", theE, "[", theD, "] <- R", 
        effects[Deffect-1], "[",
        thePastD, "]", sep="")
  
  # the covariates
  # check to see if there's more than one
      if(length(covariates[[theE]])==1) {
        cat(" + beta", theE, " * X", theE, "[", theD, "]", sep="")
      } else if (length(covariates[[theE]]) > 1) {
        cat(" +", inprod, "(beta", theE, "[] , X", theE, "[", theD, ",])", sep="")
    
      }   
  # spatial
  if(theE %in% spatial & !geostat) {
        	cat("+ R", theE, "Spatial[", theD, "]", sep="")
  }

  cat("\n")
        
    }
    }

  # the observations
      theE = paste(prefix, "observations", sep="")
      theD = paste("D", theE, sep="")
      thePastD = paste("D", effects[length(effects)], sep="")
      thePastS = paste("S", effects[length(effects)], sep="")
      
      # the loop
      cat("\n", indent, "for(", theD, " in ", thePastS, "[",
        thePastD, "]:(", thePastS, "[", thePastD, "+1]-1)){\n\n", sep="")
      indent = encodeString(" ", width=2*length(effects)+2)
      # distribution of observations
 cat(indent, prefix, observations,  "[", theD, "] ~ ", ddist, "(mean", 
        theE, "[", theD, "]", sep="")
      # if binomial, add offsets
      if( family=="binomial" & length(offset)) {
        cat(", ")
        cat(toString(paste( prefix, offset,"[", theD, "]", sep="")))   
      } else if(family %in% c("normal", "gaussian"))
        cat(", Tobservations")  
      cat(")\n",sep="")
      # mean of observations  
      cat(indent, link, "mean", theE, "[", theD, "]", endlink, " <- R",  
        effects[length(effects)], "[", thePastD, "]", sep="")

      if(length(covariates[[theE]])==1) {
         cat(" + beta", theE, " * X", theE, "[", theD, "]", sep="")
      } else if (length(covariates[[theE]]) > 1) {
        cat(" +", inprod,"(beta", theE, "[] , X", theE, "[", theD, ",])", sep="")
      }   

    
#      if(length(covariates[[theE]])==1) {
#        cat(" + betaobservations * Xobservations[", theD, "]", sep="")
#      } else if (length(covariates[["observations"]]) > 1) {
#        cat(" + inprod2(betaobservations[] , Xobservations[", theD, ",])", sep="")
#      }
       if(family!="binomial" & length(offset)) {
           cat("+", gsub(",", "+", toString(paste(prefix,offset,  "[", theD, "]", sep="")))) 
      }   
  cat("\n\n")
    
    
    # the closing brackets
   effects = c(effects, "observations")
    
    for(Deffect in seq(length(effects),1)) {
       cat(encodeString("", width=2*Deffect-2), "}#", effects[Deffect],"\n",sep="")
    }

    cat("\n\n")

  # the spatial distributions
	if(geostat){
	    for(Deffect in spatial) {  
        	cat("R", Deffect, "[1:N", Deffect, "Spatial] ~ spatial.exp(",sep="")
			cat("mean", Deffect, "[1:N", Deffect, "Spatial], xSpatial", Deffect, 
				"[1:N", Deffect, "Spatial],\n    ySpatial",
				Deffect, "[1:N", Deffect, "Spatial],",sep="")
			if(brugs) {
				cat("V",prefix,  Deffect,", ",
					"scale", prefix, Deffect, ", ", 
					"1)\n", sep="")
				cat("V",prefix,  Deffect," <- ",
					"pow(SD",prefix,  Deffect,",2)\n", sep="") 					
			} else {
			cat("T", prefix, Deffect,
					", scale",prefix,  Deffect,", 1)\n", sep="")
		}
		# prior on phi
		parName = paste("phi", prefix, Deffect, sep="")
		if(parName %in% names(priors))
			cat(parName, " ~ ", unlist(priors[parName]), "\n")
		else
		 cat(parName, " ~ dgamma(0.01, 0.01)")
	 	# transform phi to range
		cat("scale", prefix, Deffect, " <- 2*", 
				"SD", prefix,Deffect,"/",parName,sep="")
		}
	} else { # a BYM model
	  for(Deffect in spatial) {  
        cat("R", Deffect, "Spatial[1:N", Deffect, "Spatial] ~ car.normal(adj", Deffect, 
          "[], weights", Deffect, "[], num", Deffect, "[], T", Deffect, "Spatial)\n", sep="")  
      }  
}

  # the priors
  cat("\n\n# priors\n\n")  
  
  intName = paste("intercept", prefix, sep="")
  if(intName %in% names(priors)) {
	  cat(intName, " ~ ", unlist(priors[intName]), "\n")
  } else {
   cat(intName, " ~ dunif(-1000,1000)\n")
  }
  for(Deffect in names(covariates)) {
    thiscov = covariates[[Deffect]]
    if(length(thiscov)==1) {
		betaName = paste("beta", prefix, Deffect, sep="")
	  if(betaName %in% names(priors)) {
		 cat(betaName, " ~ ", unlist(priors[betaName]), "\n")
	  } else {
        	cat(betaName, " ~ dunif(-1000,1000)\n", sep="")
	  }
    } else if(length(thiscov)>1) {
      for(Dpar in 1:length(covariates[[Deffect]])) {
          parName = paste("beta", prefix, Deffect, "[", Dpar, "]", sep="")
          if(any(names(priors)==parName)) {
           cat(parName, " ~ ", unlist(priors[parName]), "\n")
          } else {     
          cat(parName,  " ~ dunif(-1000,1000)\n",sep="")
        }
    }
  }
}

  
if(length(reparam)){ 
 cat("interceptUnparam", prefix, "<- intercept", prefix, sep="")
 for(Deffect in reparam) { #names(covariates)){
   if(length(covariates[[Deffect]])==1) {
   # add prefix here
     cat("- beta", Deffect, " * X", Deffect, "reparam", sep="")
   } else {
     if (length(covariates[[Deffect]])>1){
          cat("-", inprod, "(beta", prefix, Deffect, "[]," , "X", Deffect, "reparam[])", sep="")
     }
   }
 }
}
 
cat("\n")
  if(! family %in% c("normal", "gaussian"))
    effects = effects[-length(effects)]
  for(Deffect in effects) {
	  parName = paste("SD", prefix, Deffect, sep="")
    cat("T", prefix, Deffect, " <- pow(", parName, ", -2)\n", sep="")
     if(any(names(priors) == parName)){
         cat(parName, "~", unlist(priors[parName]), "\n", sep="")
     }else{
     cat(parName, " ~ dunif(0, 100)\n", sep="")
     }
  }
 
	# if a BYM model, add prior for standard deviation of spatial effect
  if (length(spatial) & !geostat) {
  for(Deffect in spatial) {
	  parName = paste("SD", prefix, Deffect, "Spatial", sep="")
	  cat("T", prefix, Deffect, "Spatial <- pow(", 
			   parName, ", -2)\n", sep="")
       if(any(names(priors) == parName)){
       cat(parName, "~", priors[parName], "\n", sep="")
       }else{
       cat(parName, "~ dunif(0, 100)\n", sep="")    
       }
  }
}
  
if(!is.null(file)) {
  cat("\n} # model\n") 

  sink()
}
 return(invisible(scan(file, what='a', sep="\n",quiet=TRUE)))
}
