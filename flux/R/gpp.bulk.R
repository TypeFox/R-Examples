"gpp.bulk" <- function(formula, data, INDEX, window = 1, hook = "mean", oot.id = c("D", "T"), min.dp = 5, Reco.m = NULL, ts.Reco = NULL, fall.back = TRUE, ...){
	## order data
	data <- data[order(INDEX),]
	INDEX <- sort(INDEX)
	## prepare data when moving window is wanted..
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
	
	## version when Reco fluxes are taken from data as the measurements that are closest to NEE measurements
	if(is.null(Reco.m)){
		## check for sufficient number of flux measurement in both transparent and opaque measuerments
		enough.opaque <- sapply(dat.fm, function(x) sum(x[ncol(x)]==oot.id[1]))>=min.dp
    	enough.transparent <- sapply(dat.fm, function(x) sum(x[ncol(x)]==oot.id[2]))>=min.dp
		enough <- (enough.opaque + enough.transparent) > 1
		if(sum(!enough)>0){message(paste(sum(!enough), "campaign(s) did not contain enough data points: ", paste(names(dat.fm)[!enough], collapse=", ")))}
		## model gpp
		mods <- lapply(dat.fm[enough], function(x) gpp2(x[,1], x[,2], x[,3], x[,4], oot.id=oot.id, ...))
	}
	
	## version when Reco models or Reco fluxes are provided additionally
	else{
		enough <- sapply(dat.fm, function(x) nrow(x)) >= min.dp
		if(sum(!enough)>0){message(paste(sum(!enough), "campaign(s) did not contain enough data points: ", paste(names(dat.fm)[!enough], collapse=", ")))}
		## model gpp
		# when several (or only one) Reco model(s) are (is) provided
		if(is.null(ts.Reco)){
			mods <- lapply(dat.fm[enough], function(x) gpp(x[,1], x[,2], x[,3], x[,4:ncol(x)], Reco.m = Reco.m, ...))
		}
		# when Reco data are provided
		else{
			mods <- lapply(dat.fm[enough], function(x) gpp(x[,1], x[,2], x[,3], Reco.m = Reco.m, ts.Reco = ts.Reco, ...))
		}
	}
	
	### handling NA and non-sufficient models by falling back to a mean model
	## handling NA models
	if(fall.back){
		# define function
		cure <- function(x){
			if(class(x[[1]])=="try-error"){
				x$data$dat$GPP <- x$data$dat$GPP + x$data$offset
				x$data$offset <- 0
				GPP <- rnorm(50, mean = mean(x$data$dat$GPP), sd = sd(x$data$dat$GPP))
				PAR <- seq(0, 2000, length.out=50)
				x$mg <- lm(GPP ~ PAR)
				return(x)
			}
			else if(coef(x$mg)[2] > 0){
				x$data$dat$GPP <- x$data$dat$GPP + x$data$offset
				x$data$offset <- 0
				GPP <- rnorm(50, mean = mean(x$data$dat$GPP), sd = sd(x$data$dat$GPP))
				PAR <- seq(0, 2000, length.out=50)
				x$mg <- lm(GPP ~ PAR)
				return(x)
			}
			else{return(x)}
		}
		# run it
		mods <- lapply(mods, cure)
	}
	### end handling NA and non-sufficient models
	
	## extract timestamps
	extract.ts <- function(x, type){
		switch(type, mean = mean(x, na.rm=TRUE), median = median(x, na.rm=TRUE), min = min(x, na.rm=TRUE), max = max(x, na.rm=TRUE))
	}
	ts <- lapply(dat.fm[enough], function(x) extract.ts(x[,3], hook))
	
	## compile result parts for output
	out <- lapply(c(1:length(mods)), function(x) list(ts = ts[[x]], mod = mods[[x]]))
	names(out) <- nms[enough]
	class(out) <- "bgpp"
	return(out)
}