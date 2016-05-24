"budget.reco" <- function(models, new.data, set.back = NULL, time.unit = "extract", adjust = TRUE, correct = list(thresh = "get", cvrm = TRUE, wndw = 12, intvl = 0.95), return.models = FALSE){
	cl <- match.call()
	# when a set.back matrix is given to set fluxes at times back to certain values
	# (incorporating cuts or smooth starts and endings)
	if(!is.null(set.back)){
		# define set.back function (makes artifiical models)
		sb <- function(timestamp, R, models){
			if(R <= (-999)){
				modsel <- which.min(abs(sapply(models, function(x) julian(x$ts)) - as.numeric(julian(timestamp))))
				out <- models[[modsel]]
				out$ts <- timestamp
			}
			else{
				R <- rep(R, 10)
				Temp <- c(1:10)
				mod <- vector("list", 1)
				mod[[1]] <- lm(R ~ Temp)
				names(mod) <- "linear"
				out <- list(ts = timestamp, mod = mod, which.Temp = NA)
			}
			return(out)
		}
		# make them models
		set.back.models <- lapply(c(1:nrow(set.back)), function(x) sb(set.back[x,1], set.back[x,2], models))
		nms <- rep("cut", nrow(set.back))
		nms[set.back[,2] == (-999)] <- "start"
		nms[set.back[,2] == (-9999)] <- "end"
		names(set.back.models) <- nms
		# put them together with original models
		models <- c(models, set.back.models)
		# order them according to timestamp
		timeline <- sapply(models, function(x) julian(x$ts))
		keep <- (timeline >= timeline[names(timeline)=="start"]) & (timeline <= timeline[names(timeline)=="end"])
		timeline <- timeline[keep]
		models <- models[keep]
		models <- models[order(timeline)]
	}
	cut.models <- models[names(models)=="cut"]
	models <- models[names(models)!="cut"]
	# predict the whole time series in new.data with all models (quite bold when resources are critical)
	mmh <- lapply(models, function(x) predict(x$mod[[1]], newdata=data.frame(Temp = new.data[,x$which.Temp])))
	# get standard errors of the models
	seom <- lapply(models, function(x) summary(x$mod[[1]])$sigma)
	# correct
	if(adjust){
		cormod <- function(x){
			preds <- predict(x)
			meas <- x$model$R
			coef(lm(preds~meas+0))
		}
		cors <- lapply(models, function(x) cormod(x$mod[[1]]))
		mmh <- lapply(seq(length(mmh)), function(x) mmh[[x]] / cors[[x]])
	}
	# get timestamps of the models
	times <- lapply(models, function(x) x$ts)
	times1 <- times[-length(times)]
	times2 <- times[-1]
	# extract modeled values for time spans between models (take the ones valid for 1st model in span)
	mmh1 <- lapply(c(1:length(times1)), function(x) mmh[[x]][(new.data$timestamp >= times1[[x]]) & (new.data$timestamp <= times2[[x]])])
	# extract modeled values for time spans between models (take the ones valid for 2nd model in span)
	mmh2 <- lapply(c(1:length(times1)), function(x) mmh[[x+1]][(new.data$timestamp >= times1[[x]]) & (new.data$timestamp <= times2[[x]])])
	# get standard errors of the models for time spans between models (the ones valid for 1st model in span)
	seom1 <- lapply(c(1:length(times1)), function(x) rep(seom[[x]], sum((new.data$timestamp >= times1[[x]]) & (new.data$timestamp <= times2[[x]]))))
	# get standard errors of the models for time spans between models (the ones valid for 2nd model in span)
	seom2 <- lapply(c(1:length(times1)), function(x) rep(seom[[x+1]], sum((new.data$timestamp >= times1[[x]]) & (new.data$timestamp <= times2[[x]]))))
	# same for weights
	weights2 <- lapply(mmh1, function(x) c(1:length(x)))
	weights1 <- lapply(mmh2, function(x) c(length(x):1))
	# do weighted means
	reco.fluxes <- lapply(c(1:length(mmh1)), function(x) apply(cbind(mmh1[[x]], mmh2[[x]], weights1[[x]], weights2[[x]]), 1, function(y) weighted.mean(y[1:2], y[3:4])))
	seoms <- lapply(c(1:length(seom1)), function(x) apply(cbind(seom1[[x]], seom2[[x]], weights1[[x]], weights2[[x]]), 1, function(y) weighted.mean(y[1:2], y[3:4])))
	ids <- lapply(c(1:length(seoms)), function(x) rep(x, length(seoms[[x]])))
	# prepare to cut to length of whole period (from first timestamp in models to last one)
	times.spanned <- (new.data$timestamp >= times[[1]]) & (new.data$timestamp <= times[[length(times)]])
	# cut and correct time unit
	reco.fluxes <- unlist(reco.fluxes)[c(1:sum(times.spanned))]
	ses <- unlist(seoms)[c(1:sum(times.spanned))]
	ids <- unlist(ids)[c(1:sum(times.spanned))]
	# compile results
	out <- data.frame(timestamp = new.data$timestamp[times.spanned], reco.flux = reco.fluxes, reco.se = ses, reco.id = ids)
	# final outlier tests (correct modelled values that are way off limits)
	if(!is.null(correct)){
		# get simple threshold above which fluxes should be omitted
		if(correct$thresh=="get"){
			thresh <- ceiling(max(unlist(lapply(models, function(x) x$mod[[1]]$model$R))))
		}
		# if not default ("get") function expects a number which is a maximum flux allowed
		sel <- out$reco.flux < thresh
		# check for more correction
		if(correct$cvrm){
			cv <- unlist(apply(embed(out$reco.flux, correct$wndw), 1, function(x) sd(x)/mean(x)))
			sel2 <- cv < quantile(cv, correct$intvl)
			sel2[1:correct$wndw] <- TRUE
			sel <- (sel + c(sel2, rep(TRUE, correct$wndw-1))) == 2
		}
		out <- out[sel,]
		if(time.unit == "extract"){
			tu <- new.data$timestamp[-1] - new.data$timestamp[-length(new.data$timestamp)]
			t.unit <- attributes(tu)$units
			time.units <- summary(factor(tu))
			time.unit <- as.numeric(names(time.units)[which.max(time.units)])
			time.unit <- switch(t.unit, "secs" = time.unit, "mins" = time.unit*60, "hours" = time.unit*60*60, "days" = time.unit*60*60*24)
		}
		tmp <- lips(out$timestamp, out$reco.flux, time.unit)
		names(tmp) <- c("timestamp", "reco.flux")
		tmp$reco.se <- lips(out$timestamp, out$reco.se, time.unit)[,2]
		tmp$reco.id <- lips(out$timestamp, as.factor(out$reco.id), time.unit)[,2]
		out <- tmp
	}
	out <- out[,c(2:4,1)]
	out$reco.id <- as.integer(out$reco.id)
	# return models if wanted
	if(return.models){
		models <- c(models, cut.models)
		timeline <- sapply(models, function(x) julian(x$ts))
		models <- models[order(timeline)]
		class(models) <- "breco"
		out <- list(tbl = out, models = models, call=cl)
	}
	return(out)
}