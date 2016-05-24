create.survobj <- function(annotation = NULL, truncate.survival = 100) {

	# localize the survival data
	survtime <- annotation$survtime;
	survstat <- annotation$survstat;

	# check the annotation for each patient
	if (length(annotation$survtime.unit) == 0 | length(annotation$survtime) == 0) {
		warning('patient annotation missing');
		return(NA);
		}
	
	for (j in 1:nrow(annotation)) {
		if (is.na(annotation$survtime.unit[j]) | is.na(survtime[j])) { survtime[j] <- NA; }
		else if (annotation$survtime.unit[j] == "days")   { survtime[j] <- survtime[j] / 365.25; }
		else if (annotation$survtime.unit[j] == "weeks")  { survtime[j] <- survtime[j] / 52.18; }
		else if (annotation$survtime.unit[j] == "months") { survtime[j] <- survtime[j] / 12; }
		else if (annotation$survtime.unit[j] == "years")  { } # do nothing, this is the default}
		}

	# create the Surv object
	if (!is.numeric(survtime) || !is.numeric(survstat)) {
		survobj <- rep(NA, length(survtime));
		}
	else {
		# truncate survival
		survstat[survtime > truncate.survival] <- 0;
		survtime[survtime > truncate.survival] <- truncate.survival;

		# create survival object
		survobj <- Surv(survtime, survstat);
		}

	return(survobj);

	}
