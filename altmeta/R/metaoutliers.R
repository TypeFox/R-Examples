metaoutliers <- function(y, s2, model){
	if(length(y) != length(s2) | any(s2 < 0)) stop("error in the input data.")
	w <- 1/s2
	y.p <- sum(y*w)/sum(w)
	n <- length(y)

	if(missing(model)){
		hetmeasure <- metahet.base(y, s2)
		Ir2 <- hetmeasure$Ir2
		if(Ir2 < 0.3){
			model <- "FE"
			cat("This function uses fixed-effect meta-analysis because Ir2 < 30%.\n")
		}else{
			model <- "RE"
			cat("This function uses random-effects meta-analysis because Ir2 >= 30%.\n")
		}
	}

	if(!is.element(model, c("FE", "RE"))) stop("wrong input for the argument model.")

	y.p.i <- res <- std.res <- numeric(n)
	if(model == "FE"){
		for(i in 1:n){
			w.temp <- w[-i]
			y.temp <- y[-i]
			y.p.i[i] <- sum(y.temp*w.temp)/sum(w.temp)
			res[i] <- y[i] - y.p.i[i]
			var.res.i <- 1/sum(w.temp) + s2[i]
			std.res[i] <- res[i]/sqrt(var.res.i)
		}
	}else{
		for(i in 1:n){
			s2.temp <- s2[-i]
			y.temp <- y[-i]
			tau2.temp <- metahet.base(y.temp, s2.temp)$tau2.DL
			w.temp <- 1/(s2.temp + tau2.temp)
			y.p.i[i] <- sum(y.temp*w.temp)/sum(w.temp)
			res[i] <- y[i] - y.p.i[i]
			var.res.i <- 1/sum(w.temp) + s2[i] + tau2.temp
			std.res[i] <- res[i]/sqrt(var.res.i)
		}
	}

	outliers <- which(abs(std.res) >= 3)
	if(length(outliers) == 0) outliers <- "All the standardized residuals are smaller than 3"

	out <- NULL
	out$model <- model
	out$std.res <- std.res
	out$outliers <- outliers

	class(out) <- "metaoutliers"
	return(out)
}