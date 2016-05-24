`prepare.data` <-
function(dataset) {

	# rename columns of dataset
	names(dataset) <- c("id", "class", "transition", "performance")

	# convert class & transition to factors
	dataset$class <- as.factor(dataset$class)
	dataset$transition <- as.factor(dataset$transition)

	# remove missing values
	dataset <- dataset[!is.na(dataset$class),]
	dataset <- dataset[!is.na(dataset$transition),]
	dataset <- dataset[!is.na(dataset$performance),]

	invisible(dataset)
}

