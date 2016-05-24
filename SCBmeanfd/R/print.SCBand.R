print.SCBand <-
function(x, ...) 
{
	object <- x
	if (is.null(object$par) && (!is.matrix(object$nonpar)) ) {		
		cat("\nMean function estimation\n")
		cat("Bandwidth:", round(object$bandwidth,4), "\n")
		cat("SCB type:", switch(object$scbtype, no = "no SCB", normal = "normal", 
			bootstrap = "boostrap", both = "normal and bootstrap"),"\n")
		if (object$scbtype != "no") {
				cat("Confidence level:", object$level, "\n")
				cat("Quantile used for SCB:\n")		
				statresult <- if (object$scbtype == "normal") { data.frame(object$qnorm)
								} else if (object$scbtype == "boot") { data.frame(object$qboot)
								} else data.frame(object$qnorm, object$qboot)
				names(statresult) <- c(switch(object$scbtype, normal = "normal", 
				bootstrap = "bootstrap", both = c("normal", "bootstrap")))
				print(statresult, print.gap = 2L, right = FALSE, digits = 4L, row.names = FALSE)
		}

	} 
	else {
		
		Pnorm <- if (is.null(object$pnorm)) { 
			NULL
			} else if (object$pnorm == 0) { 
			"<1e-16"
		 	} else if (object$pnorm < 1e-4) { 
			paste0("<1e-", as.integer(ceiling(log(object$pnorm, 10))))  
		   	} else object$pnorm

		Pboot <- if (is.null(object$pboot)) { 
			NULL
	   		} else if (object$pboot == 0) { 
			"<1e-16"
	   		} else if (object$pboot < 1e-4) { 
			paste0("<1e-", as.integer(ceiling(log(object$pboot, 10))))  
	   		} else object$pboot
		
		if (object$scbtype == "normal") {
			statresult <- data.frame(object$teststat, Pnorm)
			names(statresult) <- c("stat", "p")
		} 
		else if (object$scbtype == "bootstrap") {
			statresult <- data.frame(object$teststat, Pboot)
			names(statresult) <- c("stat", "p")
		} 
		else { 
			statresult <- data.frame(object$teststat, Pnorm, Pboot)
			names(statresult) <- c("stat", "normal p", "bootstrap p")	
		}
		
		if (length(object$model)) {		
			cat("\nGoodness-of-fit test\n") 	
			cat ("Model for the mean function: ")
			if (length(object$model) == 1) {
				if (object$model == 0) { cat ("zero\n") 
				} else if (object$model == 1) { cat("linear\n")
				} else cat("polynomial of degree <=", object$model,"\n")
			} 	
			else cat("function space of dimension", ncol(object$model),"\n")
			cat("Bandwidth:", round(object$bandwidth, 4), "\n")
			cat("SCB type:", switch(object$scbtype, no = "no SCB", normal = "normal", 
				bootstrap = "boostrap", both = "normal and bootstrap"),"\n")
			if (object$scbtype != "no") {
			cat("Significance level:", object$level, "\n")
			cat("Test statistic and p value\n")
			print(statresult, right = FALSE, print.gap = 2L, digits = 4L, row.names = FALSE)
			} 
	
		} 
		else {
			cat("\nEquality test for mean functions\n\n")
			cat("Bandwidths:", round(object$bandwidth,4),"\n")
			cat("SCB type:", switch(object$scbtype, no = "no SCB", normal = "normal", 
				bootstrap = "boostrap", both = "normal and bootstrap"),"\n")
			if (object$scbtype != "no") {
				cat("Significance level:", object$level, "\n")
				cat("Test statistic and p value\n")
				print(statresult, right = FALSE, print.gap = 2L, digits = 4L, row.names = FALSE)
			} 
		} 
	}
}
