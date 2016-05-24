
select.mpt <- function(mpt.results, output = c("standard", "full"), round.digit = 6, dataset) {
	if(!is.list(mpt.results)) stop("mpt.results need to be a list.")
	if(length(mpt.results)< 2) stop("length(mpt.results) needs to be >= 2.") 
	n.models <- length(mpt.results)
	class.data <- vapply(mpt.results, function(x) class(x[["data"]][["observed"]]), "")
	if (!all(class.data == class.data[1])) stop("observed data from results differ (i.e., some seem to be individual data others not)")
	observed.data <- lapply(mpt.results, function(x) return(x[["data"]][["observed"]]))
	equal.data <- sapply(observed.data, function(x) identical(observed.data[[1]],x))
	if(!all(equal.data)) stop("observed data for the models differ.")
	if (is.matrix(observed.data[[1]])) {
		n.data <- dim(observed.data[[1]])[1]
	} else {
		if (is.list(observed.data[[1]])) {
			if (!missing(dataset)) {
					if (length(dataset) == 1) {
						n.data <- 1
						red.res <- lapply(mpt.results, "[", i = 1:3)
						mpt.results <- lapply(red.res, function(x) lapply(x, function(x) x[[1]][dataset,]))
					} else {
						# stop("'length(dataset > 1)' currently not implemented!")
						red.res <- lapply(mpt.results, "[", i = 1:3)
						mpt.results <- (lapply(red.res, function(x) lapply(x, function(x) {
							 tmp <- x
							 tmp[[1]] <- tmp[[1]][dataset,]
							 #browser()
							 if (isTRUE(names(tmp)[2] == "sum")) {
								tmp[[2]] <- as.data.frame(t(apply(tmp[[1]], 2, sum)))
								if ("p.value" %in% colnames(tmp[[2]])) {
									tmp[[2]][1,"p.value"] <- 1-pchisq(tmp[[2]][1,"G.Squared"], tmp[[2]][1,"df"])
								}
							 }
							 tmp
							 })))
            #browser()
            observed.data <- lapply(observed.data, function(x) list(individual = x[["individual"]][dataset,]))
						n.data <- length(dataset)
						c.aggregated <- FALSE
					}
			} else {
				n.data <- dim(observed.data[[1]][[1]])[1]
				if (any(vapply(mpt.results, function(x) length(x[[1]]), 0) != 3)) c.aggregated <- FALSE
				else c.aggregated <- TRUE
			}
		} else stop("problem with data object of mpt.results")
	}
	#browser()
  mc <- match.call()
  
  if (!is.null(names(mpt.results))) m.names <- names(mpt.results)
	else m.names <- all.vars(mc[[2]])
	
	if (n.data == 1) {
		c.fia <- sapply(mpt.results, function(x) any(grepl("^FIA$", colnames(x[["information.criteria"]]))))
		if (any(c.fia)) if (all(!c.fia)) warning(paste("FIA not available for model(s):", paste(m.names[which(c.fia == FALSE)], collapse = ", ")))
		n.parameters <- vapply(mpt.results, function(x) x[["model.info"]][["n.parameters"]], 0)
		G.Squared <- vapply(mpt.results, function(x) x[["goodness.of.fit"]][1,"G.Squared"], 0)
		df <- vapply(mpt.results, function(x) x[["goodness.of.fit"]][1,"df"], 0)
		p <- vapply(mpt.results, function(x) x[["goodness.of.fit"]][1,"p.value"], 0)
		if (any(c.fia)) {
			FIA <- vapply(mpt.results, function(x) tryCatch(x[["information.criteria"]][["FIA"]], error = function(e) NA), 0)
			delta.FIA <- FIA - min(FIA, na.rm = TRUE)
			FIA.penalty <- vapply(mpt.results, function(x) tryCatch(x[["information.criteria"]][["FIA.penalty"]], error = function(e) NA), 0)
		}
		AIC <- sapply(mpt.results, function(x) x[["information.criteria"]][["AIC"]])
		BIC <- sapply(mpt.results, function(x) x[["information.criteria"]][["BIC"]])
		delta.AIC <- AIC - min(AIC)
		denom.wAIC <- sum(exp(-0.5*(delta.AIC)))
		wAIC <- sapply(delta.AIC, function(x) exp(-0.5*(x))/denom.wAIC)
		delta.BIC <- BIC - min(BIC)
		denom.wBIC <- sum(exp(-0.5*(delta.BIC)))
		wBIC <- sapply(delta.BIC, function(x) exp(-0.5*(x))/denom.wBIC)
		df.out <- data.frame(model = m.names, n.parameters, G.Squared, df, p)
		if (any(c.fia)) {
			df.out <- cbind(df.out, FIA.penalty, delta.FIA)
			if (output[1] != "standard") df.out <- cbind(df.out, FIA)
		}
		if (output[1] == "standard") {
			df.out <- cbind(df.out, delta.AIC, wAIC, delta.BIC, wBIC)   
		} else {
			df.out <- cbind(df.out, delta.AIC, wAIC, AIC, delta.BIC, wBIC, BIC)
		}
	} else {
		#browser()
		c.fia <- sapply(mpt.results, function(x) any(grepl("^FIA$", colnames(x[["information.criteria"]][["individual"]]))))
		if (any(c.fia)) if (all(!c.fia)) warning(paste("FIA not available for model(s):", paste(m.names[which(c.fia == FALSE)], collapse = ", ")))
		n.parameters <- vapply(mpt.results, function(x) x[["model.info"]][["individual"]][1,"n.parameters"], 0)
		G.Squared.sum <- vapply(mpt.results, function(x) x[["goodness.of.fit"]][["sum"]][1,"G.Squared"], 0)
		df.sum <- vapply(mpt.results, function(x) x[["goodness.of.fit"]][["sum"]][1,"df"], 0)
		p.sum <- vapply(mpt.results, function(x) x[["goodness.of.fit"]][["sum"]][1,"p.value"], 0)
		p.smaller.05 <- vapply(mpt.results, function(x) sum(x[["goodness.of.fit"]][["individual"]][,"p.value"] < .05), 0)
		if (any(c.fia)) {
			FIA.sum <- vapply(mpt.results, function(x) {if (any(grepl("^FIA$", colnames(x[["information.criteria"]][["sum"]])))) x[["information.criteria"]][["sum"]][["FIA"]] else NA}, 0)
			FIA.penalty.sum <- vapply(mpt.results, function(x) {if (any(grepl("^FIA.penalty$", colnames(x[["information.criteria"]][["sum"]])))) x[["information.criteria"]][["sum"]][["FIA.penalty"]] else NA}, 0)
			delta.FIA.sum <- FIA.sum - min(FIA.sum, na.rm = TRUE)
			if (c.aggregated) {
				FIA.aggregated <- vapply(mpt.results, function(x) {if (any(grepl("^FIA$", colnames(x[["information.criteria"]][["aggregated"]])))) x[["information.criteria"]][["aggregated"]][["FIA"]] else NA}, 0)
				FIA.penalty.aggregated <- vapply(mpt.results, function(x) {if (any(grepl("^FIA.penalty$", colnames(x[["information.criteria"]][["aggregated"]])))) x[["information.criteria"]][["aggregated"]][["FIA.penalty"]] else NA}, 0)
				delta.FIA.aggregated <- FIA.aggregated - min(FIA.aggregated, na.rm = TRUE)
			}
			FIAs <- vapply(mpt.results, function(x) {if (any(grepl("^FIA$", colnames(x[["information.criteria"]][["individual"]])))) x[["information.criteria"]][["individual"]][["FIA"]] else rep(NA, n.data)}, rep(0, n.data))
			FIA.best <- rowSums(apply(FIAs, 1, function(x) round(x, round.digit) == min(round(x, round.digit), na.rm = TRUE)))
		}
		AIC.sum <- sapply(mpt.results, function(x) x[["information.criteria"]][["sum"]][["AIC"]])
		BIC.sum <- sapply(mpt.results, function(x) x[["goodness.of.fit"]][["sum"]][["G.Squared"]] + sum(x$model.info$individual$n.parameters) * log(sum(observed.data[[1]]$individual)))
		if (c.aggregated) {
            G.Squared.aggregated <- vapply(mpt.results, function(x) x[["goodness.of.fit"]][["aggregated"]][1,"G.Squared"], 0)
            df.aggregated <- vapply(mpt.results, function(x) x[["goodness.of.fit"]][["aggregated"]][1,"df"], 0)
            p.aggregated <- vapply(mpt.results, function(x) x[["goodness.of.fit"]][["aggregated"]][1,"p.value"], 0)
			AIC.aggregated <- sapply(mpt.results, function(x) x[["information.criteria"]][["aggregated"]][["AIC"]])
			BIC.aggregated <- sapply(mpt.results, function(x) x[["information.criteria"]][["aggregated"]][["BIC"]])
			delta.AIC.aggregated <- AIC.aggregated - min(AIC.aggregated)
			denom.wAIC.aggregated <- sum(exp(-0.5*(delta.AIC.aggregated)))
			wAIC.aggregated <- sapply(delta.AIC.aggregated, function(x) exp(-0.5*(x))/denom.wAIC.aggregated)
			delta.BIC.aggregated <- BIC.aggregated - min(BIC.aggregated)
			denom.wBIC.aggregated <- sum(exp(-0.5*(delta.BIC.aggregated)))
			wBIC.aggregated <- sapply(delta.BIC.aggregated, function(x) exp(-0.5*(x))/denom.wBIC.aggregated)
		}
		delta.AIC.sum <- AIC.sum - min(AIC.sum)
		denom.wAIC.sum <- sum(exp(-0.5*(delta.AIC.sum)))
		wAIC.sum <- sapply(delta.AIC.sum, function(x) exp(-0.5*(x))/denom.wAIC.sum)
		delta.BIC.sum <- BIC.sum - min(BIC.sum)
		denom.wBIC.sum <- sum(exp(-0.5*(delta.BIC.sum)))
		wBIC.sum <- sapply(delta.BIC.sum, function(x) exp(-0.5*(x))/denom.wBIC.sum)
		AICs <- sapply(mpt.results, function(x) x[["information.criteria"]][["individual"]][["AIC"]])
		AIC.best <- rowSums(apply(AICs, 1, function(x) round(x, round.digit) == min(round(x, round.digit))))
		BICs <- sapply(mpt.results, function(x) x[["information.criteria"]][["individual"]][["BIC"]])
		BIC.best <- rowSums(apply(BICs, 1, function(x) round(x, round.digit) == min(round(x, round.digit))))
		df.out <- data.frame(model = m.names, n.parameters, G.Squared.sum, df.sum, p.sum, p.smaller.05)
        if (c.aggregated & (output[1] != "standard")) {
            df.out <- cbind(df.out, G.Squared.aggregated, df.aggregated, p.aggregated)
        }
		if (any(c.fia)) {
			df.out <- cbind(df.out, FIA.penalty.sum, delta.FIA.sum, FIA.best)
			if (output[1] != "standard") {
				df.out <- cbind(df.out, FIA.sum)
				if (c.aggregated) df.out <- cbind(df.out, FIA.penalty.aggregated, delta.FIA.aggregated, FIA.aggregated)
			}
		}
		if (output[1] == "standard") {
			df.out <- cbind(df.out, delta.AIC.sum, wAIC.sum, AIC.best, delta.BIC.sum, wBIC.sum, BIC.best)   
		} else {
			df.out <- cbind(df.out, delta.AIC.sum, wAIC.sum, AIC.best, AIC.sum)
			if (c.aggregated) df.out <- cbind(df.out, delta.AIC.aggregated, wAIC.aggregated, AIC.aggregated)
			df.out <- cbind(df.out, delta.BIC.sum, wBIC.sum, BIC.best, BIC.sum)
			if (c.aggregated) df.out <- cbind(df.out, delta.BIC.aggregated, wBIC.aggregated, BIC.aggregated)
		}
	}
	
	rownames(df.out) <- NULL
	
	df.out[,-1] <- round(df.out[,-1], round.digit)
	df.out
	
}

