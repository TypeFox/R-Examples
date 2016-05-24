#' Summary of the results obtained from Tariff algorithm
#'
#' This function prints the summary message of the fitted results.
#'
#' @param object fitted object from \code{\link{tariff}}
#' @param top number of top CSMF to show
#' @param id the ID of a specific death to show
#' @param ... not used
#' @keywords Tariff
#' @examples
#' 
#' \donttest{
#' data("RandomVA3")
#' test <- RandomVA3[1:200, ]
#' train <- RandomVA3[201:400, ]
#' allcauses <- unique(train$cause)
#' fit <- tariff(causes.train = "cause", symps.train = train, 
#' 			symps.test = test, causes.table = allcauses)
#' correct <- which(fit$causes.test[,2] == test$cause)
#' accuracy <- length(correct) / dim(test)[1]
#' summary(fit)
#' summary(fit, top = 10)
#' summary(fit, id = "p849", top = 3)
#' }
#' 
summary.tariff <- function(object, top = 5, id = NULL, ...){

	out <- NULL
	out$top <- top
	out$id.toprint <- id
	out$N <- dim(object$score)[1]
	if(!is.null(id)){
		index <- which(rownames(object$score) == id)
		if(length(index) == 0){
			stop("Error: provided ID not found")
		}else{
			out$indiv.top <- 1:top
			out$indiv.top <- cbind(out$indiv.top, names(sort(object$score[index, ], decreasing = TRUE)[1:top]))
			out$indiv.top <- data.frame(out$indiv.top)
			colnames(out$indiv.top) <- c("Rank", "Cause")
		}
	}else{
		out$csmf.ordered <- sort(object$csmf, decreasing = TRUE)
	}


	class(out) <- "tariff_summary"
	return(out)
}

#' Print method for the summary of the results obtained from Tariff algorithm
#'
#' This function prints the summary message of the fitted results.
#'
#' @param x summary object for Tariff fit
#' @param ... not used
#' @keywords Tariff

print.tariff_summary <- function(x, ...){
	# print single death summary
	if(!is.null(x$id.toprint)){
		cat(paste0("Tariff fitted top ", x$top, " causes for death ID: ", x$id.toprint, "\n\n"))		
		print(x$indiv.top, row.names = FALSE, right = FALSE)

	# print population summary
	}else{
		cat(paste("Tariff fitted on", x$N, "deaths\n"))
		cat("\n")
		cat(paste("Top", x$top,  "CSMFs:\n"))
		csmf.out.ordered <- x$csmf.ordered
		csmf.out.ordered <- as.matrix(csmf.out.ordered[1:x$top])
		colnames(csmf.out.ordered) <- "CSMF"
	    print(csmf.out.ordered)
	}	
}