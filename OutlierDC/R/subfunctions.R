.onAttach <- function(libname, pkgname) {
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage("")
    packageStartupMessage("Package ", pkgname, " (",ver,") loaded.")
}

##################################################
setGeneric("show")
setMethod("show", "OutlierDC", function(object){

		cat("\n     Outlier Detection for Censored Data\n\n")
		cat(" Call: ")
		print(object@call)
		cat(" Algorithm: ")
		switch(object@method, 
			score = cat("Scoring algorithm"),
			boxplot = cat("Boxplot algorithm"),
			residual = cat("Residual-based algorithm")
		)
		cat(paste(" (",object@method,")", sep =""),"\n")

		cat(" Model: ")
		switch(object@rq.model, 
			PengHuang = cat("Peng and Huang estimator"),
			Powell = cat("Powell estimator"),
			Portnoy = cat("Portnoy estimator"),
			Wang = cat("Locally weighted censored quantile regression")
		)
		cat(paste(" (",object@rq.model,")", sep = ""),"\n")

		mf1 <- model.frame(object@formula, data = object@raw.data)
		resp <- model.response(mf1)
		times = resp[ ,1]
		delta = resp[ ,2]
		X = model.matrix(object@formula, data = object@raw.data)
		n <- length(times)

		if(object@method == "residual") cat(" Value for cut-off k_r: ",object@k_r,"\n")
		else if(object@method == "boxplot") cat(" Value for cut-off k_b: ",object@k_b,"\n")
		else if(object@method == "score") cat(" Value for cut-off k_s: ",object@k_s,"\n")

		cat(" # of outliers detected: ", object@n.outliers, "\n")
			
		if(object@method == "residual"){
			print.data <- cbind(times, delta, X, residual = object@score, sigma = object@k_r * object@cutoff)
			order.row <- order(times)
			print.data <- zapsmall(print.data, digits = 3)
			print.data <- as.data.frame(print.data[object@outliers, , drop = FALSE])
			order.row <- order(print.data$times, decreasing = TRUE)
			print.data <- cbind(print.data, Outlier = "*")
			#order <- order(object@score, decreasing = TRUE)
			#print.data <- zapsmall(print.data, digits = 3)
			#Signif <- ifelse(object@outliers, "*", "")
			#print.data <- cbind(print.data, Outlier = Signif)
			#print.data <- as.data.frame(print.data)
			#print.data <- print.data[order, , drop = FALSE]

			if(object@bound %in% c("both", "UB")){
				cat("\n Outliers detected:\n")
				print.head <- head(print.data)
				print(print.head)
				cat("\n", nrow(print.head) ,"of all", object@n.outliers, "outliers were displayed. \n")
			}

			if(object@bound %in% c("both", "LB")){
				cat("\n Outliers detected by lower fence:\n")
				#print(print.data[n:(n-5),], digits = 3)
				print(tail(print.data))
			}
		}
		else if(object@method ==  "score"){
			print.data <- cbind(times, delta, X, score = object@score)
			order <- order(object@score, decreasing = TRUE)
			print.data <- zapsmall(print.data, digits = 3)
			Signif <- ifelse(object@outliers, "*", "")
			print.data <- cbind(print.data, Outlier = Signif)
			print.data <- as.data.frame(print.data)
			print.data <- print.data[order, , drop = FALSE]

			if(object@bound %in% c("both", "UB")){
				cat("\n Top 6 outlying scores:\n")
				print(head(print.data))
			}
			
			if(object@bound %in% c("both", "LB")){
				cat("\n Bottom 6 outlying scores:\n")
				#print(print.data[n:(n-5),], digits = 3)
				print(tail(print.data))
			}						
		}
		else if(object@method == "boxplot"){
			cat("\n Outliers detected:\n")
			if(object@bound == "UB"){
				print.data <- cbind(times, delta, X, UB = object@upper)
			}
			else if(object@bound == "LB"){
				print.data <- cbind(times, delta, X, LB = object@lower)				
			}
			else if(object@bound == "both"){
				print.data <- cbind(times, delta, X, LB = object@lower, UB = object@upper)				
			}

			order.row <- order(times)
			print.data <- zapsmall(print.data, digits = 3)
			print.data <- as.data.frame(print.data[object@outliers, , drop = FALSE])
			order.row <- order(print.data$times, decreasing = TRUE)
			print.data <- cbind(print.data, Outlier = "*")
			print.order <- print.data[order.row, ]
			print(print.order)
			cat("\n", nrow(print.order) ,"of all", object@n.outliers, "outliers were displayed. \n")
		}
	}
)

####################################################################
setGeneric("plot")
setMethod("plot", "OutlierDC", function(x, y = NA, ...){
		mf1 <- model.frame(x@formula, x@raw.data)
		resp <- model.response(mf1)
		Times <- resp[ ,1]
		status <- resp[ ,2]

		if(x@method == "residual"){
			Residuals = x@score
			Fitted.values = x@fitted.mat[ ,3]
			limit.y <- max(Residuals, abs(x@k_r) * x@cutoff)
			plot(Fitted.values, Residuals,pch = c(1,3)[status+1], ylim = c(-1 * limit.y, limit.y), ...)
			grid()
			points(Fitted.values[x@outliers], Residuals[x@outliers], pch = c(1,3)[status[x@outliers]+1], col = "blue")
			if(x@bound %in% c("both", "UB")) abline(h = x@k_r * x@cutoff, col = "blue", lty = 2, lwd = 1.5)
			if(x@bound %in% c("both", "LB")) abline(h = -1 * x@k_r * x@cutoff, col = "blue", lty = 2, lwd = 1.5)
			legend("bottomleft",c("Censored","Event"), cex=1, pch=c(1,3), bty = "n")
		}
		else if(x@method == "score"){
			Scores = x@score
			tmp <- qqnorm(Scores, main = "Q-Q plot of outlying scores", pch = c(1,3)[status+1])
			qqline(Scores, col = "tomato", lwd = 1.5)
			grid()
			if(!is.logical(x@upper)) abline(h = x@upper, col = "blue", lwd = 2, lty = 2)
			if(!is.logical(x@lower)) abline(h = x@lower, col = "blue", lwd = 2, lty = 2)

			points(tmp$x[x@outliers], tmp$y[x@outliers], col = "blue", pch = c(1,3)[status[x@outliers]+1])
			legend("bottomright",c("Censored","Event"), cex=1, pch=c(1,3), bty = "n")
		}
		else if(x@method == "boxplot"){
			cov.x <- model.matrix(x@formula, data = x@raw.data)
			n <- ncol(cov.x) - 1
			for(i in 1:n){
				covariate <- cov.x[ ,i+1]
				order <- order(covariate)
				plot(covariate, Times, pch = c(1,3)[status+1], axes=F, main = paste("Scatter plot for covariate", i, sep = ""), ...)
				points(covariate[x@outliers], Times[x@outliers], pch =  c(1,3)[status[x@outliers]+1], col = "blue")
				
				if(x@bound %in% c("both", "UB")){
				lines(covariate[order], x@fitted.mat[order,4], lwd = 1.5,lty = 3)
				lines(covariate[order], x@upper[order], col = "blue", lwd = 1.5,lty = 2)
				}
				lines(covariate[order], x@fitted.mat[order,3], lwd = 1.5)
				if(x@bound %in% c("both", "LB")){
				lines(covariate[order], x@fitted.mat[order,2], lwd = 1.5,lty = 3)
				lines(covariate[order], x@lower[order], col = "blue", lwd = 1.5,lty = 2)
				}
				axis(1)
				axis(2, at= round(quantile(Times, probs = 0:5/5),1), labels= round(quantile(Times, probs = 0:5/5),1))
				box()
				rug(jitter(cov.x[order], amount = 0.01), ticksize = 0.01)
				legend("topright",c("Censored","Event"), cex=1, pch=c(1,3), bty = "n")
			}
		}
	}
)

###################################################################
setGeneric("coef")
setMethod("coef", "OutlierDC", function(object) round(object@coefficients,3) )

####################################################################
setGeneric("summary")
setMethod("summary","OutlierDC", function(object, taus = c(.1, .25, .5, .75, .9)){
		fit <- crq(object@formula, data = object@raw.data, method = object@method)
		summary(fit, taus = taus)
	}
)

####################################################################
setGeneric("update")
setMethod("update","OutlierDC", function(object, k_s = NA, LB = NA){
# object: OutlierDC object
# UB, LB: sample quantiles
# This function is designed for the scoring algorithm
		Scores = object@score
		UB = k_s
		if(is.na(UB) & is.na(LB)) stop("Please, update the object using the argument UB and LB")
		else if(!is.na(UB) & is.na(LB)){
			if(!is.logical(object@lower)) object@outliers <- Scores > UB | object@outliers
			else object@outliers <- Scores > UB
			object@k_s <- UB
			object@upper <- UB
		}
		else if(is.na(UB) & !is.na(LB)){
			if(!is.logical(object@upper)) object@outliers <- Scores < LB | object@outliers	
			else object@outliers <-  Scores < LB
			
			object@lower <- LB
		}
		else{
			mf1 <- model.frame(object@formula, data = object@raw.data)
			resp <- model.response(mf1)
			times = resp[,1]
		
			object@outliers <- ifelse(times > object@fitted.mat[,3],
					Scores > UB, 
					Scores < LB)
			object@lower <- LB
			object@upper <- UB
			object@k_s <- UB
		}
		object@n.outliers <- sum(object@outliers)
		object@refined.data <- object@raw.data[-object@outliers,, drop = FALSE]
		return(object)
	}
)
# End @ Feb 2013 by Soo-Heang Eo
