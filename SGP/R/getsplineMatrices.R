`getsplineMatrices` <- 
function(
	my.matrices,
	my.matrix.content_area.progression,
	my.matrix.grade.progression,
	my.matrix.time.progression,
	my.matrix.time.progression.lags,
	my.exact.grade.progression.sequence=FALSE,
	return.highest.order.matrix=FALSE,
	return.multiple.matrices=FALSE,
	my.matrix.order=NULL,
	my.matrix.highest.order=NULL,
	my.matrix.time.dependency=NULL,
	what.to.return="MATRICES") {

	Matrix_TF <- Order <- Grade <- NULL
	if (is.null(my.exact.grade.progression.sequence)) my.exact.grade.progression.sequence <- FALSE

	### Utility functions

	splineMatrix_equality <- function(my.matrix, my.order=NULL, my.exact.grade.progression.sequence=FALSE) {
		tmp.df <- data.frame()
		if (is.null(my.order)) my.order <- (2:length(my.matrix.time.progression))-1
		if (my.exact.grade.progression.sequence) my.order <- length(my.matrix.time.progression)-1
		for (i in seq_along(my.order)) {
			tmp.df[i,1] <- identical(my.matrix@Content_Areas[[1]], tail(my.matrix.content_area.progression, my.order[i]+1)) & 
					identical(my.matrix@Grade_Progression[[1]], as.character(tail(my.matrix.grade.progression, my.order[i]+1))) & 
					identical(my.matrix@Time[[1]], as.character(tail(my.matrix.time.progression, my.order[i]+1))) &
					identical(all.equal(as.numeric(my.matrix@Time_Lags[[1]]), as.numeric(tail(my.matrix.time.progression.lags, my.order[i]))), TRUE) &
					identical(names(my.matrix@Version[['Matrix_Information']][['SGPt']][['VARIABLES']]), names(my.matrix.time.dependency))
			tmp.df[i,2] <-	my.order[i]
			tmp.df[i,3] <- tail(my.matrix@Grade_Progression[[1]], 1)
		}
		names(tmp.df) <- c("Matrix_TF", "Order", "Grade")
		return(tmp.df)
	} ### END splineMatrix_equality

	getsplineMatrix <- function(
				my.exact.grade.progression.sequence=FALSE,
				my.matrix.order=NULL,
				what.to.return="MATRICES") {

		if (what.to.return=="ORDERS") {
			tmp.list <- lapply(my.matrices, splineMatrix_equality, my.exact.grade.progression.sequence=my.exact.grade.progression.sequence)
			tmp.orders <-  as.numeric(unlist(sapply(tmp.list[sapply(tmp.list, function(x) any(x[['Matrix_TF']]))], subset, Matrix_TF==TRUE, select=Order)))
			return(sort(tmp.orders))
		}
		if (what.to.return=="GRADES") {
			tmp.list <- lapply(my.matrices, splineMatrix_equality, my.exact.grade.progression.sequence=my.exact.grade.progression.sequence)
			tmp.grades <-  as.numeric(unlist(sapply(tmp.list[sapply(tmp.list, function(x) any(x[['Matrix_TF']]))], subset, Matrix_TF==TRUE, select=Grade)))
			return(sort(tmp.grades))
		}
		if (what.to.return=="MATRICES") {
			if (is.null(my.matrix.order)) my.matrix.order <- length(my.matrix.time.progression.lags)
			my.tmp.index <- which(sapply(lapply(my.matrices, splineMatrix_equality, my.order=my.matrix.order), function(x) x[['Matrix_TF']]))
			if (length(my.tmp.index)==0) {
				stop(paste("\tNOTE: No splineMatrix exists with designated content_area.progression:", paste(my.matrix.content_area.progression, collapse=", "), "year.progression:", 
					paste(my.matrix.time.progression, collapse=", "), "and grade.progression", paste(my.matrix.grade.progression, collapse=", ")))
			}
			if (length(my.tmp.index)>1) {
				if (!return.multiple.matrices) {
					stop(paste("\tNOTE: Multiple splineMatrix objects exists with designated content_area.progression:", paste(my.matrix.content_area.progression, collapse=", "), 
					"year.progression:", paste(my.matrix.time.progression, collapse=", "), "grade.progression:", paste(my.matrix.grade.progression, collapse=", "), 
					"time.progression.lags:", paste(my.matrix.time.progression.lags, collapse=", ")))
				} else return(my.matrices[my.tmp.index])
			} else {
				return(my.matrices[[my.tmp.index]])
			}
		}
	} ### END getsplineMatrix


	### Call getsplineMatrix

	if (is.null(my.matrix.order)) {
		tmp.orders <- getsplineMatrix(my.exact.grade.progression.sequence=my.exact.grade.progression.sequence, what.to.return="ORDERS")
	} else {
		tmp.orders <- my.matrix.order
	}

	if (what.to.return=="ORDERS")  return(tmp.orders)

	if (!is.null(my.matrix.highest.order)) {
		tmp.orders <- tmp.orders[tmp.orders <= my.matrix.highest.order]
	}

	if (what.to.return=="MATRICES") {
		if (return.highest.order.matrix) tmp.orders <- tail(tmp.orders, 1)
		return(lapply(tmp.orders, function(x) getsplineMatrix(my.matrix.order=x)))
	}
} ### END getsplineMatrices
