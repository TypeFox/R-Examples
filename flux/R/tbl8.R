tbl8 <- function(models){
	smc <- class(models)
	if(smc=="breco"){
		tbl <- data.frame(t(sapply(models, function(x) coef(x$mod[[1]]))))
		tbl$R2 <- sapply(models, function(x) sum(resid(x$mod[[1]])^2) / sum((x$mod[[1]]$model$R-mean(x$mod[[1]]$model$R))^2))
		tbl$n <- sapply(models, function(x) length(x$mod[[1]]$model$R))
		tbl$ts <- sapply(models, function(x) as.character(x$ts))
		tbl$which.Temp <- sapply(models, function(x) as.character(x$which.Temp))
	}
	else if(smc=="bgpp"){
		tbl <- data.frame(nms = names(models))
		tbl$te <- sapply(models, function(x) class(x$mod[[1]]))
		tbl$GPmax <- NA
		tbl$GPmax[tbl$te=="nls"] <- sapply(models[tbl$te=="nls"], function(x) coef(x$mod[[1]])[1])
		tbl$alpha <- NA
		tbl$alpha[tbl$te=="nls"] <- sapply(models[tbl$te=="nls"], function(x) coef(x$mod[[1]])[2])
		tbl$a <- NA
		tbl$a[tbl$te=="lm"] <- sapply(models[tbl$te=="lm"], function(x) coef(x$mod[[1]])[1])
		tbl$b <- NA
		tbl$b[tbl$te=="lm"] <- sapply(models[tbl$te=="lm"], function(x) coef(x$mod[[1]])[2])
		tbl$offset <- sapply(models, function(x) x$mod$data$offset)
		tbl$R2 <- sapply(models, function(x) sum(resid(x$mod[[1]])^2) / sum((x$mod[[1]]$model$GPP-mean(x$mod[[1]]$model$GPP))^2))
		tbl$n <- sapply(models, function(x) length(x$mod[[1]]$model$GPP))
		tbl$ts <- sapply(models, function(x) as.character(x$ts))
	}
	else{stop("no recognised class", call. = FALSE)}
	return(tbl)
} 