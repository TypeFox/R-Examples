fit.LTR <- function(dataset1, dataset2) {

	# try to treat both datasets as matrices
	dataset1 <- as.matrix(dataset1);
	dataset2 <- as.matrix(dataset2);

	# do some error checking on the input parameters
	if (class(dataset1) != 'matrix') { stop('dataset1 is not a matrix'); }
	if (class(dataset2) != 'matrix') { stop('dataset2 is not a matrix'); }
	if (ncol(dataset1) != ncol(dataset2)) { stop('dataset1 and dataset2 have different numbers of columns'); }
	if (nrow(dataset1) != nrow(dataset2)) { stop('dataset1 and dataset2 have different numbers of rows'); }
	if (!all(rownames(dataset1) == rownames(dataset2))) { warning('non-matching rownames for dataset1 and dataset2'); }

	# keep track of fit characteristics
	to.return <- list(
		rownames = rownames(dataset1),
		intercepts = vector(mode = "numeric", length = nrow(dataset1)),
		slopes = vector(mode = "numeric", length = nrow(dataset1)),
		r.squared = vector(mode = "numeric", length = nrow(dataset1)),
		residuals = matrix(
			nrow = nrow(dataset1),
			ncol = ncol(dataset1)
			)
		)

	# loop over all ProbeSets and perform fit
	for (i in 1:nrow(dataset1)) {

		# fit the linear model
		fit <- lm(
			formula = dataset2 ~ dataset1,
			data = data.frame(
				dataset2 = dataset2[i,],
				dataset1 = dataset1[i,]
				)
			);

		# save model-fit characteristics
		to.return$intercepts[i] <- as.numeric( fit$coefficients[1] );
		to.return$slopes[i]     <- as.numeric( fit$coefficients[2] );
		to.return$r.squared[i]  <- summary(fit)$r.squared;
		to.return$residuals[i,] <- as.numeric( fit$residuals );

		}

	# return a list of all fit-characteristics
	return(to.return);

	}
