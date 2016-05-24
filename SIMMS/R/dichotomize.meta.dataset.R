dichotomize.meta.dataset <- function(feature.name, expression.data, survival.data, other.data = NULL, data.type.ordinal = FALSE) {

	# we'll return the overall groups and survival data to the caller
	groups   <- vector();
	survival <- vector();
	status   <- vector();

	# loop over all datasets
	for (i in 1:length(expression.data)) {

		# localize the current dataset
		expression.values <- expression.data[[i]];
		expression.values <- as.vector( unlist( expression.values[feature.name,] ) );
		survival.object   <- survival.data[[i]];

		# dichotomize this dataset if needed
		if (!data.type.ordinal) {
			dichotomized.results <- SIMMS::dichotomize.dataset(expression.values);
			}
		else {
			dichotomized.results <- expression.values;
			}

		# add the returned data to return-objects
		groups   <- c(groups,   dichotomized.results);
		survival <- c(survival, survival.object[,1]);
		status   <- c(status,   survival.object[,2]);

		}

	# create an object of data to return
	to.return <- list(
		groups = groups,
		survtime = survival,
		survstat = status
		);

	# if there is an extra data unlist it and pass it back
	if (!is.null(other.data) & class(other.data) == 'list') {

		# loop over all elements
		for (i in 1:length(other.data)) {

			# create a temporary vector to store the data
			temp <- vector();

			# loop over all datasets and concatenate the data
			for (j in 1:length(other.data[[i]])) {

				# check that the feature exists in this dataset
				expression.values <- expression.data[[j]];
				expression.values <- as.vector( unlist( expression.values[feature.name,] ) );
				if ( all( is.na(expression.values) ) ) { next; }

				# if it does exist in this dataset keep the annotation data
				temp <- c(temp, as.character(other.data[[i]][[j]]));

				}

			# ensure data is right length
			if (length(temp) != length(to.return$groups)) { next; }

			# add the data to the return object
			to.return[[3+i]] <- temp;
			names(to.return)[[3+i]] <- names(other.data)[[i]];

			}

		}

	return(to.return);

	}
