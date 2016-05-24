clean <-
function(mast, set, v.avg.min=0.4, v.avg.max=50, dir.clean=TRUE, turb.clean=4, icing=FALSE, rep=NULL, n.rep=5) {
### cleaning faulty values of mast, set or specified set of mast
	
	r <- NULL
		
	if(missing(mast) && missing(set)) stop("No data to clean - please specify mast and/or set")
	if(!is.null(v.avg.min)) if(!is.numeric(v.avg.min)) stop("'v.avg.min' must be numeric or NULL")
	if(!is.null(v.avg.max)) if(!is.numeric(v.avg.max)) stop("'v.avg.max' must be numeric or NULL")
	if(is.null(dir.clean)) stop("'dir.clean' must be specified as TRUE or FALSE")
	if(!is.null(turb.clean)) if(!is.numeric(turb.clean)) stop("'turb.clean' must be numeric or NULL")
	if(is.null(icing)) stop("'icing' must be specified as TRUE or FALSE")
	if(!is.null(rep)) if(any(is.character(rep)==FALSE)) stop("'rep' must be a vector of characters or NULL")
	if(!is.null(rep) && is.null(n.rep)) stop("Please specify 'n.rep'")
	if(!is.null(rep)) if(!is.null(n.rep)) if(!is.numeric(n.rep)) stop("'n.rep' must be numeric or NULL") 
	if(missing(mast) && !missing(set)) { # set
		if(class(set)!="set") stop(substitute(set), " is no set object")
		set$data <- clean.int(set$data, v.avg.min, v.avg.max, dir.clean, turb.clean, icing, rep, n.rep+1)
		r <- set
	} else if(!missing(mast) && missing(set)) { # mast
		if(class(mast)!="mast") stop(substitute(mast), " is no mast object")
		num.sets <- length(mast$sets)
		for(s in 1:num.sets) {
			message("Cleaning set ", s, "...")
			mast$sets[[s]]$data <- clean.int(mast$sets[[s]]$data, v.avg.min, v.avg.max, dir.clean, turb.clean, icing, rep, n.rep+1)
		}
		r <- mast
	} else if(!is.null(mast) && !is.null(set)) { # set of mast
		if(class(mast)!="mast") stop(substitute(mast), " is no mast object")
		num.sets <- length(mast$sets)
		if(!is.numeric(set)) set <- match(set, names(mast$sets))
		if(is.na(set)) stop("Set not found")
		if(set<0 || set>num.sets) stop("Set not found")
		message("Cleaning set ", set, "...")
		mast$sets[[set]]$data <- clean.int(mast$sets[[set]]$data, v.avg.min, v.avg.max, dir.clean, turb.clean, icing, rep, n.rep+1)
		r <- mast
	}
	
	return(r)
}
