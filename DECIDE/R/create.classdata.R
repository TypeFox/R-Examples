`create.classdata` <-
function(dataset) {

	# create separate datasets for each class

	data_class <- list()
	no_classes <- length(levels(dataset$class))

	for (i in 1:no_classes) {
		data_class[[i]] <- dataset[which(dataset$class==i),]
		}
	for (i in 1:no_classes) {
		names(data_class)[i] <- paste("class", i, sep="")
		}

	return(data_class)
}

