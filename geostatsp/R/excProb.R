excFunQQ = function(themat, threshold) {


	
	if(length(themat)) {
	over = themat[,"x"]>threshold
	
	if(any(is.na(over))) {
		return(NA)
	}
	if(all(over)){
		prob = 1
	} else if(any(over)) {
		
		toInt = rbind(c(threshold, approx(themat[,"x"], themat[,"y"], threshold)$y),
				themat[over,]
		)
		if(requireNamespace('pracma', quietly=TRUE)){
	  	prob = pracma::trapz(toInt[,"x"], toInt[,"y"])
    } else {
      theDiff = diff(toInt[,'x'])/2
      theDiff = c(0,theDiff) + c(theDiff,0)
      prob = sum(theDiff * toInt[,"y"])
    }
	} else {
		prob = 0
	} 
	} else { # matrix has length zero
		prob = NA
	}
	
	prob
} 

excProb = function(x, threshold=0, random=FALSE, template=NULL, 
		templateIdCol=NULL, nuggetInPrediction=TRUE) {
	

	
elementsColumnwise = FALSE # default is indexes for x[[index]] refer to 
# cells numbered columnwise
	lowertail=FALSE	
	
if(is.list(x))	{
 # model is from lgm
	if(length(grep("^predict|^param|^summary", names(x)))>2){
		template = raster(x$predict)
		if(!random) {
			# check for boxcox
			if(any(names(x$predict)=="predict.boxcox")) {
				resp = x$predict[["predict.boxcox"]]	
				if(x$param["boxcox"]!=0) {
					threshold = (threshold^x$param["boxcox"] - 1) /
						x$param["boxcox"]
#					lowertail = x$param["boxcox"] < 0
				} else {
					threshold = log(threshold)
				}
			} else if(any(names(x$predict)=="predict.log")){ 
				resp = log(x$predict[["predict"]])
				threshold = log(threshold)
			} else {
				resp = x$predict[["predict"]]
			}
		}	else { # end if !random
			resp = x$predict[["random"]]
		}
		
		theSD = values(x$predict[["krigeSd"]])
		if(nuggetInPrediction) {
			theSD = sqrt(theSD^2 + x$param["nugget"])
		} 
		
		
		excProbAll = pnorm( threshold , mean=values(resp), sd= theSD,
				lower.tail=lowertail)

		
	} else { # isn't from lgm but is a list
	
	# check if it's from glgm
	 if(all(c("inla", "raster","parameters") %in% names(x))) {
		template = x$raster[['space']]
		if(!random ) {
			x = x$inla$marginals.lincomb.derived
		} else {
			x = x$inla$marginals.random$space
		}
	
	 }
	 excProbAll = unlist(lapply(x, excFunQQ, threshold=threshold))
	 names(excProbAll) = names(x)
	 
	 
	 } # end not from lgm	
		
	} else { # not a list, must be a matrix or data frame with two columns.
		excProbAll = excFunQQ(x, threshold)
		
	}		


	# make sure probabilities are between zero and 1
excProbAll = pmax(pmin(excProbAll, 1),0)

if(length(grep("^Raster", class(template)))) {
	# fill in NA's for cells with no predictions
	if(any(names(template)=='space')){
		# names of excProbAll should refer to space ID's
		allNames = paste(
				gsub("[[:digit:]]+$", "", 
						names(excProbAll)[1]),
				values(template[['space']]), 
				sep='')
		excProbAll = excProbAll[allNames]
	}
	template = raster(template)
	if(elementsColumnwise) {
		values(template) = matrix(excProbAll, 
							nrow=nrow(template),ncol=ncol(template),
							byrow=FALSE)
	} else {
		values(template) = excProbAll
	}
	excProbAll = template
	names(excProbAll) = paste("exc","random"[random], threshold, sep="")
	
} 

if(length(grep("(SpatialPolygonsDataFrame|SpatialPointsDataFrame)", class(template)))) {
	newcol=paste("exc", "random"[random],threshold, sep="")
	if(is.null(templateIdCol)) {
		template[[newcol]] = excProbAll
	} else {
	template[[newcol]] = 
			excProbAll[template@data[templateIdCol]]
	}
	excProbAll = template
} 


excProbAll
} 