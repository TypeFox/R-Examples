apply.LTR <- function(x, model.fit) {

	# try to treat the input parameters as appropriate classes
	x <- as.matrix(x);
	model.fit <- as.list(model.fit)

	# ensure that the input parameters actually have the correct classes
	if ( class(x) != 'matrix' ) { stop('x is not a matrix'); }
	if ( class(model.fit) != 'list' ) { stop('model.fit is not a list'); }
	if (! 'intercepts' %in% names(model.fit) ) { stop('model.fit does not contain intercepts'); }
	if (! 'slopes' %in% names(model.fit) ) { stop('model.fit does not contain slopes'); }
	if (! 'r.squared' %in% names(model.fit) ) { stop('model.fit does not contain r.squared'); }
	if (! 'residuals' %in% names(model.fit) ) { stop('model.fit does not contain residuals'); }

	# ensure that the input parameters have the correct dimensions
	model.fit$residuals <- as.matrix(model.fit$residuals);
	if ( class(model.fit$residuals) != 'matrix' ) { stop('model.fit$residuals is not a matrix'); }
	if ( length(model.fit$slopes) != length(model.fit$intercepts) ) { stop('invalid model: intercepts and slopes of different lengths'); }
	if ( length(model.fit$slopes) != nrow(model.fit$residuals) ) { stop('invalid model: slopes and residuals of different lengths'); }
	if ( length(model.fit$slopes) < nrow(x) ) { stop('invalid input data: contains probes not in the model'); }
	if ( length(model.fit$slopes) < nrow(x) ) { stop('invalid input data: missing probes in the model'); }
	
	# create a matrix to store the results
	to.return <- matrix( 
		nrow = nrow(x),
		ncol = ncol(x)
		);
	rownames(to.return) <- rownames(x);
	colnames(to.return) <- colnames(x);

	# loop through all ProbeSets applying the fit
	for (i in 1:nrow(x)) {
		to.return[i,] <- model.fit$intercepts[i] + model.fit$slopes[i] * x[i,];
		}

	# return the fitted expression matrix
	return(to.return);

	}
