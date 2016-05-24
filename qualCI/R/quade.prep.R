quade.prep <- function(data, set, treatment, withinRank, unit=NULL, betweenRank){
    # save call
    call <- match.call()
    
	# check for presence of needed objects
    ind <- match(c("data","set","treatment","withinRank","betweenRank"),names(call),nomatch=0)
    if(ind[1]==0){
    	stop("A 'data' argument is required.",call.=FALSE)
    }
    if(ind[2]==0){
    	stop("A 'set' argument is required.",call.=FALSE)
    }
    if(ind[3]==0){
    	stop("You must specify the treatment indicator using 'treatment'.",call.=FALSE)
    }
    if(ind[4]==0){
    	stop("You must specify the within-set rank using 'withinRank'.",call.=FALSE)
    }
    if(ind[5]==0){
    	stop("A 'betweenRank' argument is required.",call.=FALSE)
    }
    
    # check that data is dataframe
	if(class(data)!="data.frame"){
		stop("'data' object must be a data frame.")
	}
	
	# check that variables in dataframe
   	if (!(set %in% colnames(data))) {
        stop("Please make sure that 'set' is in the data object.",call.=FALSE)
    }
   	if (!(treatment %in% colnames(data))) {
        stop("Please make sure that 'treatment' is in the data object.",call.=FALSE)
    }
   	if (!(withinRank %in% colnames(data))) {
        stop("Please make sure that 'withinRank' is in the data object.",call.=FALSE)
    } 

	# check that no missing data
	dsub <- data[,c(set,treatment,withinRank)]
	if(any(is.na(dsub))){
		stop("Please make sure that there is no missing data in set identifiers, treatment indicators, or the within-set rank.")	
	}
	
	# check that treatment is binary
    if(any(data[,treatment] %in% c(0,1)==FALSE)){
    	stop("Treatment indicator must be binary.")
    }
    


    if(!is.numeric(data[,withinRank])){
    	stop("Within-set rank must be numeric / integer.")
    }
  
    # get all unique sets
    sets <- as.character(as.vector(unique(data[,set])))
	
	# turn variables into lists by set
	treatList <- by(data, data[,set],function(x) x[, treatment])
	withinRankList <- by(data, data[,set],function(x) x[, withinRank])
	if(!is.null(unit)){
		units <- by(data[,unit], data[,set],function(x) as.character(x))
	}
	
	# check that each set has treated and control
    if(!all(sapply(treatList, function(x) all(c(0,1) %in% x)))){
    		stop("Each set must contain at least one treated and one control unit.")
    	}
	
	# units per set
	units.per.set <- sapply(withinRankList, function(x) length(x))
	
  	# check that ranks don't exceed number of units per set (by set)
    if(!all(sapply(withinRankList, function(x) all(sort(x)==1:length(x)) | all(sort(x)==0:(length(x)-1))))) {
    		stop("Within-set ranks must be sequential integers from 1 to number of units in the set (or from 0 to number of units in the set minus one). Ties are not allowed.")
    }
    
    # check that betweenRank is a vector
    if(!is.vector(betweenRank) | !is.numeric(betweenRank)){
    	stop('betweenRank should be a numeric vector of between-set ranks.')
    }
    
    # if no names for betweenRank, set to integers
    if(is.null(names(betweenRank))){
    	names(betweenRank) <- as.character(1:length(betweenRank))
    }
    
    # check that sets are contained in betweenRanks
 	if(any(sets %in% names(betweenRank)==FALSE)){
 		stop("'betweenRank' object is missing between-set rank for at least one set.")
 	}   
    
    # calculate possible treatments	
	out <- lapply(treatList,function(x) possible.treatments(x))
	names(out) <- as.character(1:length(out))
	for(s in 1:length(out)){
		out[[s]]$obsTreat <- treatList[[s]]	
		out[[s]]$withinRank <- withinRankList[[s]]
		bw <- betweenRank[names(out)[s]]
		names(bw) <- NULL
		out[[s]]$rank <- bw
		if(!is.null(unit)){
			set.labels <- units[[s]]
		} else{
			set.labels <- as.character(1:length(withinRankList[[s]]))	
			}
		colnames(out[[s]]$possibleTreat) <- set.labels
		names(out[[s]]$obsTreat) <- set.labels
		names(out[[s]]$withinRank) <- set.labels
		}
	
	attr(out,"unitNames") <- ifelse(!is.null(unit), TRUE, FALSE)
	attr(out,"pairs") <- ifelse(all(units.per.set==2), TRUE, FALSE)
	class(out) <- "matchedSets"
	return(out)
}