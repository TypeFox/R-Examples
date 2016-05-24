"reco.bulk" <- function(formula, data, INDEX, window = 1, hook = "mean", remove.outliers = FALSE, fall.back = TRUE, ...){
	## order data
	data <- data[order(INDEX),]
	INDEX <- sort(INDEX)
	## extract terms
	tt <- terms(formula)
	## perpare data when moving window is wanted..
	if(window != 1){
		ind <- unique(INDEX)
		tmp <- embed(ind, window)
		nms <- apply(tmp, 1, paste, collapse=".")
		dat <- apply(tmp, 1, function(x) data[INDEX %in% x,])
	}
	## ..or just prepare by chopping
	else{
		dat <- split(data, INDEX)
		nms <- unique(INDEX)
	}
	## extract the relevant pieces of data needed for modelling
	dat.fm <- lapply(dat, function(x) get_all_vars(formula, x))
	## model approach 1 (one temperature)
	if(length(attr(tt, "variables")) == 3){
		mods <- lapply(dat.fm, function(x) reco(x[,1], x[,2], ...))
		which.Temp <- sapply(c(1:length(mods)), function(x) as.character(attr(tt, "variables"))[-c(1:2)])
	}
	## model approach 2 (more than one temperature)
	else{
		mods <- lapply(dat.fm, function(x) lapply(c(2:ncol(x)), function(y) reco(x[,1], x[,y], ...)))
		## handle complete NA attempts
		all.na <- sapply(mods, function(x) sum(sapply(x, is.na))!=length(x))
		which.not <- names(all.na)[!all.na]
		if(sum(!all.na)>0){
			warning(paste("model for", paste(which.not, collapse=", "), "could not be build"))
			## extract the selector for the best model
			mods.sel <- lapply(mods[all.na], function(x) which.min(lapply(x, function(y) ifelse(class(y[[1]]) == "try-error", NA, AIC(y[[1]])))))
			## get the best model each
			mods[all.na] <- lapply(c(1:length(mods[all.na])), function(x) mods[all.na][[x]][[mods.sel[[x]]]])
			## which temperature did the best model?
			which.Temp <- rep(NA, length(mods))
			which.Temp[all.na] <- sapply(c(1:length(mods[all.na])), function(x) as.character(attr(tt, "variables"))[-c(1:2)][[mods.sel[[x]]]])
		}
		else{
			mods.sel <- lapply(mods, function(x) which.min(lapply(x, function(y) ifelse(class(y[[1]]) == "try-error", NA, AIC(y[[1]])))))
			## get the best model each
			mods <- lapply(c(1:length(mods)), function(x) mods[all.na][[x]][[mods.sel[[x]]]])
			## which temperature did the best model?
			which.Temp <- sapply(c(1:length(mods)), function(x) as.character(attr(tt, "variables"))[-c(1:2)][[mods.sel[[x]]]])
		}	
	}
	## outlier handling
	## first the function is defined on a per model basis
	if(remove.outliers){
		get.out <- function(mod){
			if(length(mod[[1]])==1){mod <- mod}
			if(length(mod[[1]])>1){
				outlier <- boxplot.stats(resid(mod[[1]]))$out
        # define invalid function
        invalid <- function(x){
          if (missing(x) || is.null(x) || length(x) == 0) 
            return(TRUE)
          if (is.list(x)) 
            return(all(sapply(x, invalid)))
          else if (is.vector(x)) 
            return(all(is.na(x)))
          else return(FALSE)
        }
				if(invalid(outlier)){mod <- mod}
				else{
					sel2 <- match(outlier, resid(mod[[1]]))
					nd <- data.frame(R = mod[[1]]$model$R, Temp = mod[[1]]$model$Temp)[-sel2,]
					mod.red <- reco(nd$R, nd$Temp, ...)
					# fall back to unreduced version if reduced version fails
					if(class(mod.red[[1]])=="try-error"){mod <- mod[[1]]}
					else{
						mod <- mod.red[1]
						if(length(sel2)>0){mod <- c(mod, n.out = length(sel2))}
					}
				}
			}
		}
		mods <- lapply(mods, function(x) get.out(x))
		mods <- lapply(mods, function(x) get.out(x))
	}
	
	### handling NA and non-sufficient models by falling back to a mean model
	## handling NA models
	sel.na <- sapply(mods, function(x) class(x[[1]]))!="nls"
	mods[sel.na] <- lapply(dat.fm[sel.na], function(x) list(linear = lm(R ~ Temp, data = data.frame(R = rep(mean(x[,1]), 5), Temp = 1:5))))
	which.Temp[sel.na] <- unlist(attr(tt, "term.labels")[1])
	## handling models with negative R ~ Temp relationships
	if(fall.back){
		# define function
		cure <- function(mod){
			if(coef(lm(R ~ Temp, mod[[1]]$model))[2] < 0){
				data.art <- data.frame(R = rep(mean(mod[[1]]$model$R), 30), Temp = 1:30)
				neu <- lm(R ~ Temp, data.art)
				alt <- lm(R ~ Temp, mod[[1]]$model)
				list(neu=neu, alt=alt)
			}
			else{mod}
		}
		# run it
		mods <- lapply(mods, cure)
	}
	### end handling NA and non-sufficient models
	
	## cure missing class attribute for later plotting and handling
	for(i in c(1:length(mods))){
		if(!is.na(mods[[i]][1])){class(mods[[i]][[1]]) <- c("reco", class(mods[[i]][[1]]))}
	}
	
	## extract timestamp
	extract.ts <- function(x, type){
		switch(type, mean = mean(x, na.rm=TRUE), median = median(x, na.rm=TRUE), min = min(x, na.rm=TRUE), max = max(x, na.rm=TRUE))
	}
	ts <- lapply(dat.fm, function(x) extract.ts(x[,ncol(x)], hook))
	
	## compile result parts for output
	out <- lapply(c(1:length(mods)), function(x) list(ts = ts[[x]], mod = mods[[x]], which.Temp = which.Temp[[x]]))
	names(out) <- nms
	class(out) <- "breco"
	return(out)
}