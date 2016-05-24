ipwtm <- function(
	exposure,
	family,
	link,
	numerator = NULL,
	denominator,
	id,
	tstart,
	timevar,
	type,
	data,
	corstr = "ar1",
	trunc = NULL,
	...)
	{
		#save input
			tempcall <- match.call()
		#some basic input checks
			if (!("exposure" %in% names(tempcall))) stop("No exposure variable specified")
			if (!("family" %in% names(tempcall)) | ("family" %in% names(tempcall) & !(tempcall$family %in% c("binomial", "survival", "multinomial", "ordinal", "gaussian")))) stop("No valid family specified (\"binomial\", \"survival\", \"multinomial\", \"ordinal\", \"gaussian\")")
			if (tempcall$family == "binomial") {if(!(tempcall$link %in% c("logit", "probit", "cauchit", "log", "cloglog"))) stop("No valid link function specified for family = binomial (\"logit\", \"probit\", \"cauchit\", \"log\", \"cloglog\")")}
			if (tempcall$family == "ordinal" ) {if(!(tempcall$link %in% c("logit", "probit", "cauchit", "cloglog"))) stop("No valid link function specified for family = binomial (\"logit\", \"probit\", \"cauchit\", \"cloglog\")")}
			if (!("denominator" %in% names(tempcall))) stop("No denominator model specified")
			if (!is.null(tempcall$numerator) & !is(eval(tempcall$numerator), "formula")) stop("Invalid numerator formula specified")
			if (!is.null(tempcall$denominator) & !is(eval(tempcall$denominator), "formula")) stop("Invalid denominator formula specified")
			if (!("id" %in% names(tempcall))) stop("No patient id specified")
			if (tempcall$family == "survival" & !("tstart" %in% names(tempcall))) stop("No tstart specified, is necessary for family = \"survival\"")
			if (!("timevar" %in% names(tempcall))) stop("No timevar specified")
			if (!("type" %in% names(tempcall))) stop("No type specified (\"first\" or \"all\")")
			if (!(tempcall$type %in% c("first", "all"))) stop("No type specified (\"first\" or \"all\")")
			if (tempcall$family %in% c("survival", "multinomial", "ordinal") & tempcall$type == "all") stop(paste("Type \"all\" not yet implemented for family = ", deparse(tempcall$family, width.cutoff = 500), sep = ""))
			if (tempcall$family %in% c("gaussian") & tempcall$type == "first") stop(paste("Type \"first\" not implemented for family = ", deparse(tempcall$family, width.cutoff = 500), sep = ""))
			if (tempcall$family %in% c("gaussian") & !("numerator" %in% names(tempcall))) stop("Numerator necessary for family = \"gaussian\"")
			if (!("data" %in% names(tempcall))) stop("No data specified")
			if (!is.null(tempcall$trunc)) {if(tempcall$trunc < 0 | tempcall$trunc > 0.5) stop("Invalid truncation percentage specified (0-0.5)")}
		#record original order of dataframe so that the output can be returned in the same order
			order.orig <- 1:nrow(data)
			order.orig <- order.orig[order(
				eval(parse(text = paste("data$", deparse(tempcall$id, width.cutoff = 500), sep = ""))),
				eval(parse(text = paste("data$", deparse(tempcall$timevar, width.cutoff = 500), sep = "")))
				)] #sort as below
		#sort dataframe on follow-up time within each individual, necessary for cumulative products below
			data <- data[order(
				eval(parse(text = paste("data$", deparse(tempcall$id, width.cutoff = 500), sep = ""))),
				eval(parse(text = paste("data$", deparse(tempcall$timevar, width.cutoff = 500), sep = "")))
				),]
		#make new dataframe for newly computed variables, to prevent variable name conflicts
			tempdat <- data.frame(
				id = data[,as.character(tempcall$id)],
				timevar = data[,as.character(tempcall$timevar)],
				exposure = data[,as.character(tempcall$exposure)]
			)
		#make selection variable, time points up to first switch from lowest value, or all time points
			if (type == "first" & (family == "binomial" | family == "survival"))
				{tempdat$selvar <- do.call("c", lapply(split(tempdat$exposure, tempdat$id),function(x)if (!is.na(match(1, x))) return(c(rep(1,match(1, x)),rep(0,length(x)-match(1, x)))) else return(rep(1,length(x)))))}
			if (type == "first" & (family == "multinomial" | family == "ordinal")){
					z <- unique(tempdat$exposure)[unique(tempdat$exposure) != sort(unique(tempdat$exposure))[1]]
					min2 <- function(x)ifelse(min(is.na(unique(x))) == 1, NA, min(x, na.rm = TRUE))
					tempdat$selvar <- do.call("c", lapply(split(tempdat$exposure, tempdat$id),function(x)if (!is.na(min2(match(z, x)))) return(c(rep(1,min2(match(z, x))),rep(0,length(x)-min2(match(z, x))))) else return(rep(1,length(x)))))
			}
			if (type == "all")
				{tempdat$selvar <- rep(1, nrow(tempdat))}
		#weights binomial, type "first"
			if (tempcall$family == "binomial" & tempcall$type == "first") {
				if(tempcall$link == "logit") lf <- binomial(link = logit)
				if(tempcall$link == "probit") lf  <- binomial(link = probit)
				if(tempcall$link == "cauchit") lf  <- binomial(link = cauchit)
				if(tempcall$link == "log") lf  <- binomial(link = log)
				if(tempcall$link == "cloglog") lf  <- binomial(link = cloglog)
				if (is.null(tempcall$numerator)) tempdat$w.numerator <- 1
				else {
					mod1 <- glm(
						formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
						family = lf,
						data = data,
						subset = tempdat$selvar == 1,
						na.action = na.fail,
						...)
					tempdat$p.numerator <- vector("numeric", nrow(tempdat))
					tempdat$p.numerator[tempdat$exposure == 0 & tempdat$selvar == 1] <- 1 - predict.glm(mod1, type = "response")[tempdat$exposure[tempdat$selvar == 1] == 0]
					tempdat$p.numerator[tempdat$exposure == 1 & tempdat$selvar == 1] <- predict.glm(mod1, type = "response")[tempdat$exposure[tempdat$selvar == 1] == 1]
					tempdat$p.numerator[tempdat$selvar == 0] <- 1
					tempdat$w.numerator <- unlist(lapply(split(tempdat$p.numerator, tempdat$id), function(x)cumprod(x)))
						mod1$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
						mod1$call$family <- tempcall$link
						mod1$call$data <- tempcall$data
						mod1$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
				}
				mod2 <- glm(
					formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
					family = lf,
					data = data,
					subset = tempdat$selvar == 1,
					na.action = na.fail,
					...)
				tempdat$p.denominator <- vector("numeric", nrow(tempdat))
				tempdat$p.denominator[tempdat$exposure == 0 & tempdat$selvar == 1] <- 1 - predict.glm(mod2, type = "response")[tempdat$exposure[tempdat$selvar == 1] == 0]
				tempdat$p.denominator[tempdat$exposure == 1 & tempdat$selvar == 1] <- predict.glm(mod2, type = "response")[tempdat$exposure[tempdat$selvar == 1] == 1]
				tempdat$p.denominator[tempdat$selvar == 0] <- 1
				tempdat$w.denominator <- unlist(lapply(split(tempdat$p.denominator, tempdat$id), function(x)cumprod(x)))
				mod2$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
				mod2$call$family <- tempcall$link
				mod2$call$data <- tempcall$data
				mod2$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
				tempdat$ipw.weights <- tempdat$w.numerator/tempdat$w.denominator
			}
		#weights binomial, type "all"
			if (tempcall$family == "binomial" & tempcall$type == "all") {
				if(tempcall$link == "logit") lf <- binomial(link = logit)
				if(tempcall$link == "probit") lf  <- binomial(link = probit)
				if(tempcall$link == "cauchit") lf  <- binomial(link = cauchit)
				if(tempcall$link == "log") lf  <- binomial(link = log)
				if(tempcall$link == "cloglog") lf  <- binomial(link = cloglog)
				if (is.null(tempcall$numerator)) tempdat$w.numerator <- 1
				else {
					mod1 <- glm(
						formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
						family = lf,
						data = data,
						na.action = na.fail,
						...)
					tempdat$p.numerator <- vector("numeric", nrow(tempdat))
					tempdat$p.numerator[tempdat$exposure == 0] <- 1 - predict.glm(mod1, type = "response")[tempdat$exposure == 0]
					tempdat$p.numerator[tempdat$exposure == 1] <- predict.glm(mod1, type = "response")[tempdat$exposure == 1]
					tempdat$w.numerator <- unlist(lapply(split(tempdat$p.numerator, tempdat$id), function(x)cumprod(x)))
						mod1$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
						mod1$call$family <- tempcall$link
						mod1$call$data <- tempcall$data
				}
				mod2 <- glm(
					formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
					family = lf,
					data = data,
					na.action = na.fail,
					...)
				tempdat$p.denominator <- vector("numeric", nrow(tempdat))
				tempdat$p.denominator[tempdat$exposure == 0] <- 1 - predict.glm(mod2, type = "response")[tempdat$exposure == 0]
				tempdat$p.denominator[tempdat$exposure == 1] <- predict.glm(mod2, type = "response")[tempdat$exposure == 1]
				tempdat$w.denominator <- unlist(lapply(split(tempdat$p.denominator, tempdat$id), function(x)cumprod(x)))
				mod2$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
				mod2$call$family <- tempcall$link
				mod2$call$data <- tempcall$data
				tempdat$ipw.weights <- tempdat$w.numerator/tempdat$w.denominator
			}
		#weights Cox
			if (tempcall$family == "survival") {
				if (is.null(tempcall$numerator)) tempdat$w.numerator <- 1
					else {
					mod1 <- coxph(
						formula = eval(parse(text = paste("Surv(", deparse(tempcall$tstart), ", ", deparse(tempcall$timevar, width.cutoff = 500), ", ", deparse(tempcall$exposure, width.cutoff = 500), ") ", deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
						data = data,
						subset = tempdat$selvar == 1,
						na.action = na.fail,
						method = "efron",
						...)
					bh <- basehaz(mod1, centered = TRUE)
					temp <- data.frame(timevar = sort(unique(tempdat$timevar)))
					temp <- merge(temp, data.frame(timevar = bh$time, bashaz.cum.numerator = bh$hazard), by = "timevar", all.x = TRUE);rm(bh)
					if (is.na(temp$bashaz.cum.numerator[temp$timevar == min(unique(tempdat$timevar))])) temp$bashaz.cum.numerator[temp$timevar == min(unique(tempdat$timevar))] <- 0
					temp$bashaz.cum.numerator <- approx(x = temp$timevar, y = temp$bashaz.cum.numerator, xout = temp$timevar, method = "constant", rule = 2)$y
					temp$bashaz.numerator[1] <- temp$bashaz.cum.numerator[1]
					temp$bashaz.numerator[2:nrow(temp)] <- diff(temp$bashaz.cum.numerator, 1)
					temp$bashaz.cum.numerator <- NULL
					tempdat <- merge(tempdat, temp, by = "timevar", all.x = TRUE);rm(temp)
					tempdat <- tempdat[order(tempdat$id, tempdat$timevar),]
					tempdat$risk.numerator[tempdat$selvar == 1] <- predict(mod1, type="risk", centered = TRUE)
					tempdat$hazard.numerator[tempdat$selvar == 1] <- with(tempdat[tempdat$selvar == 1,], bashaz.numerator*risk.numerator)
					tempdat$p.numerator[with(tempdat, selvar == 1 & exposure == 0)] <- with(tempdat[with(tempdat, selvar == 1 & exposure == 0),], exp(-1*bashaz.numerator*risk.numerator))
					tempdat$p.numerator[with(tempdat, selvar == 1 & exposure == 1)] <- 1 - with(tempdat[with(tempdat, selvar == 1 & exposure == 1),], exp(-1*bashaz.numerator*risk.numerator))
					tempdat$p.numerator[tempdat$selvar == 0] <- 1
					tempdat$w.numerator <- unsplit(lapply(split(tempdat$p.numerator, tempdat$id), function(x) cumprod(x)), tempdat$id)
					mod1$call$formula <- eval(parse(text = paste("Surv(", deparse(tempcall$tstart), ", ", deparse(tempcall$timevar, width.cutoff = 500), ", ", deparse(tempcall$exposure, width.cutoff = 500), ") ", deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
					mod1$call$data <- tempcall$data
				}
				mod2 <- coxph(
					formula = eval(parse(text = paste("Surv(", deparse(tempcall$tstart), ", ", deparse(tempcall$timevar, width.cutoff = 500), ", ", deparse(tempcall$exposure, width.cutoff = 500), ") ", deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
					data = data,
					subset = tempdat$selvar == 1,
					na.action = na.fail,
					method = "efron",
					...)
				bh <- basehaz(mod2, centered = TRUE)
				temp <- data.frame(timevar = sort(unique(tempdat$timevar)))
				temp <- merge(temp, data.frame(timevar = bh$time, bashaz.cum.denominator = bh$hazard), by = "timevar", all.x = TRUE);rm(bh)
				if (is.na(temp$bashaz.cum.denominator[temp$timevar == min(unique(tempdat$timevar))])) temp$bashaz.cum.denominator[temp$timevar == min(unique(tempdat$timevar))] <- 0
				temp$bashaz.cum.denominator <- approx(x = temp$timevar, y = temp$bashaz.cum.denominator, xout = temp$timevar, method = "constant", rule = 2)$y
				temp$bashaz.denominator[1] <- temp$bashaz.cum.denominator[1]
				temp$bashaz.denominator[2:nrow(temp)] <- diff(temp$bashaz.cum.denominator, 1)
				temp$bashaz.cum.denominator <- NULL
				tempdat <- merge(tempdat, temp, by = "timevar", all.x = TRUE);rm(temp)
				tempdat <- tempdat[order(tempdat$id, tempdat$timevar),]
				tempdat$risk.denominator[tempdat$selvar == 1] <- predict(mod2, type="risk", centered = TRUE)
				tempdat$hazard.denominator[tempdat$selvar == 1] <- with(tempdat[tempdat$selvar == 1,], bashaz.denominator*risk.denominator)
				tempdat$p.denominator[with(tempdat, selvar == 1 & exposure == 0)] <- with(tempdat[with(tempdat, selvar == 1 & exposure == 0),], exp(-1*bashaz.denominator*risk.denominator))
				tempdat$p.denominator[with(tempdat, selvar == 1 & exposure == 1)] <- 1 - with(tempdat[with(tempdat, selvar == 1 & exposure == 1),], exp(-1*bashaz.denominator*risk.denominator))
				tempdat$p.denominator[tempdat$selvar == 0] <- 1
				tempdat$w.denominator <- unsplit(lapply(split(tempdat$p.denominator, tempdat$id), function(x)cumprod(x)), tempdat$id)
				mod2$call$formula <- eval(parse(text = paste("Surv(", deparse(tempcall$tstart), ", ", deparse(tempcall$timevar, width.cutoff = 500), ", ", deparse(tempcall$exposure, width.cutoff = 500), ") ", deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
				mod2$call$data <- tempcall$data
				mod2$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
				tempdat$ipw.weights <- tempdat$w.numerator/tempdat$w.denominator
			}
		#weights multinomial
			if (tempcall$family == "multinomial") {
				if (is.null(tempcall$numerator)) tempdat$p.numerator <- 1
				else {
					mod1 <- multinom(
						formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
						data = data,
						subset = tempdat$selvar == 1,
						na.action = na.fail,
						...)
					pred1 <- as.data.frame(predict(mod1, type = "probs"))
					tempdat$p.numerator[tempdat$selvar == 0] <- 1
					for (i in 1:length(unique(tempdat$exposure)))tempdat$p.numerator[with(tempdat, tempdat$selvar == 1 & exposure == sort(unique(tempdat$exposure))[i])] <- pred1[tempdat$exposure[tempdat$selvar == 1] == sort(unique(tempdat$exposure))[i],i]
					mod1$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
					mod1$call$data <- tempcall$data
					mod1$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
				}
				mod2 <- multinom(
					formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
					data = data,
					subset = tempdat$selvar == 1,
					na.action = na.fail,
					...)
				pred2 <- as.data.frame(predict(mod2, type = "probs"))
				tempdat$p.denominator[tempdat$selvar == 0] <- 1
				for (i in 1:length(unique(tempdat$exposure)))tempdat$p.denominator[with(tempdat, tempdat$selvar == 1 & exposure == sort(unique(tempdat$exposure))[i])] <- pred2[tempdat$exposure[tempdat$selvar == 1] == sort(unique(tempdat$exposure))[i],i]
				mod2$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
				mod2$call$data <- tempcall$data
				mod2$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
				tempdat$ipw.weights <- unsplit(lapply(split(with(tempdat, p.numerator/p.denominator), tempdat$id), function(x)cumprod(x)), tempdat$id)
			}
		#weights ordinal
			if (tempcall$family == "ordinal") {
				if(tempcall$link == "logit") m <- "logistic"
				if(tempcall$link == "probit") m  <- "probit"
				if(tempcall$link == "cloglog") m  <- "cloglog"
				if(tempcall$link == "cauchit") m  <- "cauchit"
				if (is.null(tempcall$numerator)) tempdat$p.numerator <- 1
				else {
					mod1 <- polr(
						formula = eval(parse(text = paste("as.factor(", deparse(tempcall$exposure, width.cutoff = 500), ")", deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
						data = data,
						method = m,
						subset = tempdat$selvar == 1,
						na.action = na.fail,
						...)
					pred1 <- as.data.frame(predict(mod1, type = "probs"))
					tempdat$p.numerator[tempdat$selvar == 0] <- 1
					for (i in 1:length(unique(tempdat$exposure)))tempdat$p.numerator[with(tempdat, tempdat$selvar == 1 & exposure == sort(unique(tempdat$exposure))[i])] <- pred1[tempdat$exposure[tempdat$selvar == 1] == sort(unique(tempdat$exposure))[i],i]
					mod1$call$formula <- eval(parse(text = paste("as.factor(", deparse(tempcall$exposure, width.cutoff = 500), ")", deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
					mod1$call$data <- tempcall$data
					mod1$call$method <- m
					mod1$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
				}
				mod2 <- polr(
					formula = eval(parse(text = paste("as.factor(", deparse(tempcall$exposure, width.cutoff = 500), ")", deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
					data = data,
					method = m,
					subset = tempdat$selvar == 1,
					na.action = na.fail,
					...)
				pred2 <- as.data.frame(predict(mod2, type = "probs"))
				tempdat$p.denominator[tempdat$selvar == 0] <- 1
				for (i in 1:length(unique(tempdat$exposure)))tempdat$p.denominator[with(tempdat, tempdat$selvar == 1 & exposure == sort(unique(tempdat$exposure))[i])] <- pred2[tempdat$exposure[tempdat$selvar == 1] == sort(unique(tempdat$exposure))[i],i]
				mod2$call$formula <- eval(parse(text = paste("as.factor(", deparse(tempcall$exposure, width.cutoff = 500), ")", deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
				mod2$call$data <- tempcall$data
				mod2$call$method <- m
				mod2$call$subset <- paste("up to first instance of ", deparse(tempcall$exposure, width.cutoff = 500), " = 1 (selvar == 1)", sep = "")
				tempdat$ipw.weights <- unsplit(lapply(split(with(tempdat, p.numerator/p.denominator), tempdat$id), function(x)cumprod(x)), tempdat$id)
			}
		#weights gaussian
			if (tempcall$family == "gaussian") {
				mod1 <- geeglm(
					formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = ""))),
					data = data,
					id = tempdat$id,
					corstr = tempcall$corstr,
					waves = tempdat$timevar,
					...)
				mod1$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$numerator, width.cutoff = 500), sep = "")))
				mod1$call$data <- tempcall$data
				mod1$call$id <- tempcall$id
				mod1$call$corstr <- tempcall$corstr
				mod1$call$waves <- tempcall$waves
				mod2 <- geeglm(
					formula = eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = ""))),
					data = data,
					id = tempdat$id,
					corstr = tempcall$corstr,
					waves = tempdat$timevar,
					...)
				mod2$call$formula <- eval(parse(text = paste(deparse(tempcall$exposure, width.cutoff = 500), deparse(tempcall$denominator, width.cutoff = 500), sep = "")))
				mod2$call$data <- tempcall$data
				mod2$call$id <- tempcall$id
				mod2$call$corstr <- tempcall$corstr
				mod2$call$waves <- tempcall$waves
				tempdat$kdens1 <- dnorm(tempdat$exposure, predict(mod1), as.numeric(sqrt(summary(mod1)$dispersion)[1]))
				tempdat$kdens2 <- dnorm(tempdat$exposure, predict(mod2), as.numeric(sqrt(summary(mod2)$dispersion)[1]))
				tempdat$ipw.weights <- unsplit(lapply(split(with(tempdat, kdens1/kdens2), tempdat$id), function(x)cumprod(x)), tempdat$id)
			}
		#check for NA's in weights
			if (sum(is.na(tempdat$ipw.weights)) > 0) stop ("NA's in weights!")
		#truncate weights, when trunc value is specified (0-0.5)
			if (!(is.null(tempcall$trunc))){
				tempdat$weights.trunc <- tempdat$ipw.weights
				tempdat$weights.trunc[tempdat$ipw.weights <= quantile(tempdat$ipw.weights, 0+trunc)] <- quantile(tempdat$ipw.weights, 0+trunc)
				tempdat$weights.trunc[tempdat$ipw.weights >  quantile(tempdat$ipw.weights, 1-trunc)] <- quantile(tempdat$ipw.weights, 1-trunc)
			}
		#return results in the same order as the original input dataframe
			if (is.null(tempcall$trunc)){
				if (is.null(tempcall$numerator)) return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], call = tempcall, selvar = tempdat$selvar[order(order.orig)], den.mod = mod2))
				else return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], call = tempcall, selvar = tempdat$selvar[order(order.orig)], num.mod = mod1, den.mod = mod2))
			}
			else{
				if (is.null(tempcall$numerator)) return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], weights.trunc = tempdat$weights.trunc[order(order.orig)], call = tempcall, selvar = tempdat$selvar[order(order.orig)], den.mod = mod2))
				else return(list(ipw.weights = tempdat$ipw.weights[order(order.orig)], weights.trunc = tempdat$weights.trunc[order(order.orig)], call = tempcall, selvar = tempdat$selvar[order(order.orig)], num.mod = mod1, den.mod = mod2))
			}
}
