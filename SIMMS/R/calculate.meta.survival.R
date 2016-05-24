calculate.meta.survival <- function(feature.name, expression.data, survival.data, rounding = 3, other.data = NULL, data.type.ordinal = FALSE) {

	# verify that we got appropriate input data
	to.abort <- FALSE;
	if ("list" != class(expression.data)) { to.abort <- TRUE; }
	if ("list" != class(survival.data)) { to.abort <- TRUE; }
	if (length(expression.data) != length(survival.data)) { to.abort <- TRUE; }

	# stop processing if we have bad data
	if (to.abort) {
		warning('data failed sanity-checking');
		return(
			list(
				"cox.stats" = rep(NA,5),
			 	"cox.obj" = NA
				)
			);
		}

	# dichotomize meta
	dichotomized.data <- SIMMS::dichotomize.meta.dataset(
		feature.name = feature.name,
		expression.data = expression.data,
		survival.data = survival.data,
		other.data = other.data,
		data.type.ordinal = data.type.ordinal
		);

	# handle all-NA values (i.e. feature not in the dataset)
	if (0 == length(dichotomized.data$groups) || 
		0 == length(dichotomized.data$survtime) || 
		0 == length(dichotomized.data$survstat) || 
		all(is.na(dichotomized.data$groups))
		) {
		cat('\n\tfeature not in the dataset: ', feature.name);
		return( 
			list(
				"cox.stats" = rep(NA,5),
			 	"cox.obj" = NA
				)
			);
		}

	# fit coxph model
	return(
		SIMMS::fit.coxmodel(
			groups = dichotomized.data$groups,
			survobj = Surv(dichotomized.data$survtime, dichotomized.data$survstat),
			rounding = rounding,
			other.data = other.data,
			data.type.ordinal = data.type.ordinal
			)
		);
	}
