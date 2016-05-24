summary.SCBand <-
function(object, ...) 
{
	cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

	cat("Data:\n")
	if (is.list(object$x)) {
	df <- data.frame( c(paste0("[", round(min(object$x[[1]]),3),", ",round(max(object$x[[1]]),3),"]"),		
			if (!is.null(object$y)) {
			c(paste0("[",round(min(object$y[[1]]),3),", ",round(max(object$y[[1]]),3),"]"), 
			nrow(object$y[[1]])) }), 
			c(paste0("[", round(min(object$x[[1]]),3),", ",round(max(object$x[[1]]),3),"]"),		
			if (!is.null(object$y)) {
			c(paste0("[",round(min(object$y[[2]]),3),", ",round(max(object$y[[2]]),3),"]"), 
			nrow(object$y[[2]]))
		}))
	colnames(df) <- paste("Data",1:2)
	rownames(df) <- c("Range.x", if(!is.null(object$y)) c("Range.y", "Sample size"))
	print(df, print.gap = 2L, right = FALSE, digits = 4L, row.names = TRUE)
	}
	else {
	cat("Range(x): [", round(min(object$x),3),", ",round(max(object$x),3),"]\n", sep="")
		if (!is.null(object$y)) {
			cat("Range(y): [",round(min(object$y),3),", ",round(max(object$y),3),"]\n", sep="")
			cat("Sample size:", nrow(object$y),"\n")
		} 
	}
		
	cat("\nAnalysis:\n")
	if (is.null(object$par) && (!is.matrix(object$nonpar)) ) {		
		cat("Mean function estimation\n")
		cat("Bandwidth:", round(object$bandwidth,4), "\n")
		cat("Grid size:", object$gridsize, "\n")
		cat("SCB type:", switch(object$scbtype, no = "no SCB", normal = "normal", 
			bootstrap = "boostrap", both = "normal and bootstrap"),"\n")
		if (object$scbtype != "no") {
				cat("Confidence level:", object$level, "\n")
				cat("Replicates:", switch(object$scbtype, normal = paste(object$nrep,"(normal)\n"), 
					bootstrap = paste(object$nboot,"(bootstrap)\n"), both = paste(object$nrep,"(normal)", 
					object$nboot,"(bootstrap)\n")))
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
			names(statresult) <- c("supnorm", "p")
		} 
		else if (object$scbtype == "bootstrap") {
			statresult <- data.frame(object$teststat, Pboot)
			names(statresult) <- c("supnorm", "p")
		} 
		else { 
			statresult <- data.frame(object$teststat, Pnorm, Pboot)
			names(statresult) <- c("supnorm", "normal p", "bootstrap p")	
		}
		
		if (length(object$model)) {		
			cat("Goodness-of-fit test\n") 	
			if (object$model == 0) {
				cat ("Model: zero mean function\n") 
			} else if (object$model == 1) {
				cat("Model: linear mean function\n")
			} else if (length(object$model) == 1) {
				cat("Model: polynomial mean function of degree <=", object$model,"\n")
			} else cat("Model: mean function in function space of dimension", ncol(object$model),"\n")	
			cat("Bandwidth:", round(object$bandwidth, 4), "\n")
			cat("Grid size:", object$gridsize, "\n")
			cat("SCB type:", switch(object$scbtype, no = "no SCB", normal = "normal", 
				bootstrap = "boostrap", both = "normal and bootstrap"),"\n")
			if (object$scbtype != "no") {
				cat("Significance level:", object$level, "\n")
				cat("Replicates:", switch(object$scbtype, normal = paste(object$nrep,"(normal)\n"), 
					bootstrap = paste(object$nboot,"(bootstrap)\n"), both = paste(object$nrep,"(normal)", 
					object$nboot,"(bootstrap)\n")))
				cat("Test statistic and p value(s)\n")
				print(statresult, right = FALSE, print.gap = 2L, digits = 4L, row.names = FALSE)
			} 
	
		} 
		else {
			cat("Equality test for mean functions\n")
			cat("Bandwidths:", round(object$bandwidth,4),"\n")
			cat("Grid size:", object$gridsize, "\n")
			cat("SCB type:", switch(object$scbtype, no = "no SCB", normal = "normal", 
				bootstrap = "boostrap", both = "normal and bootstrap"),"\n")
			if (object$scbtype != "no") {
				cat("Significance level:", object$level, "\n")
				cat("Replicates:", switch(object$scbtype, normal = paste(object$nrep,"(normal)\n"), 
					bootstrap = paste(object$nboot,"(bootstrap)\n"), both = paste(object$nrep,"(normal)", 
					object$nboot,"(bootstrap)\n")))
				cat("Test statistic and p value(s)\n")
				print(statresult, right = FALSE, print.gap = 2L, digits = 4L, row.names = FALSE)
			} 
		}
	}

}
