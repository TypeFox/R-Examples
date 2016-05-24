`restoreParams` <-
function(bugsResult, ragged=NULL,extraX=NULL) {
                              
thearray = bugsResult$sims.array
parnames = dimnames(thearray)[[3]]
# vector valued parameters
vecPars = grep("\\[[[:digit:]]+\\]$", parnames, value=TRUE)
# matrix valued parameters
matPars =  grep("[[:digit:]+],[[:digit:]]+\\]$", parnames, value=TRUE)
# scalar parameters
scPars = parnames[! parnames %in% c(vecPars, matPars)]
scPars = scPars[grep("^beta",scPars,invert=TRUE)]


result = list()


# if there are any random slopes (matrix parameters), put the slopes in their
# own element of the result list
if(length(matPars)) {
# find the precision matrix
precisionIndex = grep("^T", matPars)
precisions = matPars[precisionIndex]
matPars = matPars[-precisionIndex]
result$precision = thearray[,,precisions]

# convert precisions to variances
precisions = unique(gsub("\\[[[:digit:]]+,[[:digit:]]+\\]", "", precisions))


for(D in precisions) {
  result[[paste("var", D, sep="")]] = cholInvArray(result$precision, D)
}

# find the last column, which is the intercept
colno = substr(matPars, regexpr(",[[:digit:]]+\\]",  matPars)+1, 10000)
colno = substr(colno, 1, nchar(colno)-1)
maxcol = max(as.integer(colno))
if(is.na(maxcol))
  warning("can't find max col number")
# put the slopes in result
interceptcols = grep(paste(",", maxcol, "\\]", sep=""), matPars)
slopecols = matPars[-interceptcols]
result$slopes = thearray[,,slopecols]
# turn the intercept columns into vector parameters
matPars = which(parnames %in% matPars)
parnames[matPars] = gsub(paste(",", maxcol, "\\]", sep=""), "\\]", parnames[matPars])
dimnames(thearray)[[3]] = parnames
# update vecPars to reflect these intercepts
vecPars = grep("\\[[[:digit:]]+\\]$", parnames, value=TRUE)
} #end matPars


# if it's a geostatistical model, convert phi to range
	thephi = grep("^phi", scPars,value=TRUE)
	for(D in thephi) {
		thesd = gsub("^phi", "SD", D)
		if(thesd %in% scPars){
			# get rid of phi
			scPars = grep(D, scPars, invert=TRUE,value=TRUE)
			# and add range
			therange = gsub("^phi", "range", D)
			result[[therange]] = thearray[,,D] / thearray[,,thesd]			
		}
		
	}

for(D in scPars)
  result[[D]] = thearray[,,D]

if(!length(scPars))
	warning("no parameter names")

  fixedEffects = grep("^X", names(ragged), value=TRUE)

  betas=NULL
 
  
  # extract betas 
for(D in fixedEffects) {
	tobind = thearray[,,
			grep(gsub("^X","beta", D), dimnames(thearray)[[3]]),drop=F]
	if(!dim(tobind)[3])
		warning("can't find fixed effect parameters for ", D)
	
	newnames=dimnames(ragged[[D]])[[2]]

	if(length(newnames)== (dim(tobind)[3]) )
		dimnames(tobind)[[3]] = newnames

	betas = abind::abind(betas, tobind, along=3)
	
}  

result$betas = betas
  


#the random effects
# find grouping variables, all the variables with one dimensional indices
groups = unique(gsub("\\[[[:digit:]]+\\]$", "", vecPars))
thebetas = grep("^beta", groups)
if(length(thebetas))
	groups = groups[-thebetas]

for(D in groups) {
  thisGroup = grep(paste("^", D, "\\[", sep=""), vecPars, value=TRUE)
  result[[D]] = thearray[,,thisGroup]
}
# for all random effects other than spatial ones, create fitted Values
theSpatialGroups = grep("Spatial$",groups)
if(length(theSpatialGroups)) {
  notSpatial = groups[-theSpatialGroups]
} else {
  notSpatial = groups
}
for(D in notSpatial) {
  result[[paste("Fitted",D, sep="")]] = result[[D]]
}

if(is.null(ragged)) {
  return(result)
}

  # if a ragged option is given,
  # undo the reparametrisation and add better names to parameters

  groups = paste("S", substr(groups, 2, nchar(groups)), sep="")

  randomEffects = groups[groups %in% names(ragged)]
  randomEffects = names(sort(unlist(lapply(ragged[randomEffects], length))))
  randomEffects = substr(randomEffects, 2, nchar(randomEffects))
 
  if(!length(randomEffects)) {
     warning(paste(toString(groups), ":cannot find random effects"))
    return(result)
  }
  # the reparametrised mean of the random effects (just the intercept for now)
  theMeanOld = array(result[["intercept"]], 
       c(dim(result[["intercept"]]), 1))
  Nchain = dim(theMeanOld)[2]
  # "torep" is how to expand out the means, to assign fitted values from level D
  #   to random effects at level D+1
  #   at the first level, torep assigns the intercept to each group   
  torep = rep(1, length(ragged[[paste("S", randomEffects[1], sep="")]])-1)

  betanames = dimnames(result$betas)[[3]]
  
  for(D in randomEffects) {  
        theR = paste("R", D, sep="")

        theS = ragged[[paste("S", D, sep="")]]
        thenames = names(theS)[-length(theS)]
   
        if(length(thenames) != (dim(result[[theR]])[3]) )
          warning(D, "different dimensions in bugsResult and ragged")

       dimnames(result[[theR]])[[3]] = thenames   
       dimnames(result[[paste("Fitted",theR,sep="")]])[[3]] = thenames   
       
       
       DX =  paste("X", D, sep="")
       Dbeta = paste("beta", D, sep="")

      themean = theMeanOld[,,torep]
 
      
       # expand the previous mean vector to the number of current random effects
#       return(list(themean=themean, rag = ragged[[DX]], res = result[[Dbeta]][,,]))
      if(!is.null(ragged[[DX]])) {
       # if there are covariates at this level

        
	theX = t(ragged[[DX]])
	# if only one covariate at this level
	if(Dbeta %in% betanames){
		theseBetas =  result$betas[,,Dbeta,drop=F]		
	} else { # more than one covariate, have matrices
	    if(all(rownames(theX) %in% dimnames(result$betas)[[3]]) ) {
			theseBetas = result$betas[,,
					rownames(theX),drop=FALSE]
		} else {
			warning("cannot find ",D," , ", DX)
		}		
	} # end else more than one covariate
	
    for(Dchain in 1:Nchain) {
        themean[,Dchain,] = themean[,Dchain,] + 
            	abind::adrop(
					theseBetas[,Dchain,,drop=FALSE], drop=2
				)%*% theX
    } 
	}	# end there are covariates here
	
	      # get ready for the next random effect
      # the random effects (reparameterised) are the means of the next level
      theMeanOld = result[[theR]]
      # see torep above
      torep = diff(theS)
      torep = rep(1:length(torep), torep)

      # un-re-parametrise
      result[[theR]] = result[[theR]] - themean
     }         # end for D in randomEffects
     
     # add names to the spatial bit
     spatialEffects = paste("R", randomEffects, "Spatial", sep="")
     spatialEffects = spatialEffects[spatialEffects %in% names(result)]
     for(D in spatialEffects) {
        DsubR = gsub("Spatial$", "", D)
        Dsub = gsub("^R", "", DsubR)
        Dfitted = paste("Fitted",DsubR, sep="")

        #names to the spatial bit
       thenames = names(ragged[[paste("num", Dsub, sep="")]])
       if(is.null(thenames)) {
          # if there were no names in the adjacency bit, take them from Sspatial
          thenames = paste("noname", 
            1:ragged[[paste("N", Dsub, "Spatial",sep="")]], sep="")
            
          thenames[ ragged[[paste("Sspatial", Dsub, sep="")]] ] = 
            names(ragged[[paste("Sspatial", Dsub, sep="")]])
       }
        
       theID = dimnames(result[[D]])[[3]]
       theID = gsub("[[:graph:]]+\\[", "", theID)
       theID = gsub("\\]$", "", theID)
       dimnames(result[[D]])[[3]] = thenames[as.integer(theID)]

# if there are any regions without data, they'll not have an RCSDUID component
# so add them in.  Note covariates aren't added, so in effect 
# regions without data are assumed to have baseline covariates
       # regions which dont have Rstuff
       regionsNoV = thenames[!thenames %in% dimnames(result[[DsubR]])[[3]] ]

       if(length(regionsNoV)) {
        # dimension of array to hold the realisations for regions without Rstuff
        dimNoV = c(dim(result$intercept),length(regionsNoV) ) 

       # standard deviationsf or realisations for these regions
       # as an array  
       sdBig = array(result[[paste("SD",Dsub,sep="")]], dimNoV)

       # realisations of spatially independent effect for regions without Rstuff
       VfornoV =  rnorm(prod(dim(sdBig)), 0, sdBig)
       # convert to an array       
       VfornoV = array(VfornoV, dimNoV) 
       # add names
       dimnames(VfornoV)[[3]] = regionsNoV

       # find regions with a spatial random effect  
       DsubSpatial=paste(DsubR, "Spatial", sep="")      
       withSpatial = regionsNoV[regionsNoV %in% 
        dimnames(result[[DsubSpatial]])[[3]] ]
       # and add it to the realised non-spatail effect       
       VfornoV[,,withSpatial] =   VfornoV[,,withSpatial] +
            result[[DsubSpatial]][,,withSpatial]
        
      # put these results into the posterior simulations array
      result[[DsubR]] = abind(result[[DsubR]], VfornoV, along=3)
      fittedForNoV = VfornoV + array(result$intercept, dimNoV)
# add covariates if we have them
if(!is.null(extraX)) {
haveExtraX = rownames(extraX)[rownames(extraX) %in% VfornoV]
theBeta = result[[paste("beta", Dsub, sep="")]]
haveBeta = colnames(extraX)[colnames(extraX) %in% rownames(theBeta)]
  if(length(haveExtraX) & length(haveBeta)){
    fittedForNoV = fittedForNoV + 
      extraX[haveExtraX,haveBeta] %*% theBeta[haveBeta,]
    }
} #end extraX loop

      result[[Dfitted]] = abind(result[[Dfitted]], fittedForNoV, along=3)
  
      result[[DsubR]] = result[[DsubR]][,,thenames]
      result[[Dfitted]] = result[[Dfitted]][,,thenames]
  }     # end if regions with no V
       
}   #end for D in spatialEffects
     


  return(result)
}

