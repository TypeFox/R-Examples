fit.interaction.model <- function(feature1, feature2, expression.data, survival.data, data.type.ordinal = FALSE) {

	groups1 <- SIMMS::dichotomize.meta.dataset(
		expression.data = expression.data,
		survival.data = survival.data,
		feature.name = feature1,
		other.data = NULL,
		data.type.ordinal = data.type.ordinal
		);

	groups2 <- SIMMS::dichotomize.meta.dataset(
		expression.data = expression.data,
		survival.data = survival.data,
		feature.name = feature2,
		other.data = NULL,
		data.type.ordinal = data.type.ordinal
		);

	# fit the interaction Cox model
	interaction.HR <- interaction.P <- NA;
	coxmodel <- NULL;
	survobj <- Surv(groups1$survtime, groups1$survstat);

	# handle the all-missing case smoothly
	if (all(is.na(groups1$groups * groups2$groups)) || 
		length(groups1$survtime) != length(groups1$groups) ||
		length(unique(groups1$groups[which(!is.na(groups1$groups))])) < 2 ||
		length(unique(groups2$groups[which(!is.na(groups2$groups))])) < 2
		) {
		cat('\nskipping interaction model, as all-missing from one or more feature groups or improperly formatted data (see below)');
		cat('\n\tunique levels in feature ', feature1, ":", length(unique(groups1$groups[which(!is.na(groups1$groups))])));
		cat('\n\tunique levels in feature ', feature2, ":", length(unique(groups2$groups[which(!is.na(groups2$groups))])));
		}
	else {

		# levels not specified as the beta of interaction terms remains unchanged regardless of 
		# which groups is used as base-line for groups1 and groups2
		tryCatch(
			expr = {
				coxmodel <- coxph(
					survobj ~ factor(groups1$groups) + factor(groups2$groups) + as.numeric(groups1$groups == groups2$groups)
					)
				},
			error = function(ex) {
				cat("\nInteraction model failed to converge (a known coxph issue) for features: ", feature1, " and ", feature2);
				}
			);

		# extract summary statistics from cox model
		coxmodel <- summary(coxmodel);

		# check fail to converge case
		if (c("Class", "Mode") %in% names(coxmodel) 
			&& coxmodel[["Class"]] == "NULL" && coxmodel[["Mode"]] == "NULL") {
			# no extra warning needed at this stage, but keep this placeholders for unforeseen coxph cases
			}
		else {
			interaction.HR <- as.numeric(coxmodel$coefficients[3,2]);
			interaction.P  <- as.numeric(coxmodel$coefficients[3,5]);
			}
		}

	# fit the two univariate models
	survival1 <- SIMMS::calculate.meta.survival(
		feature.name = feature1,
		expression.data = expression.data,
		survival.data = survival.data,
		data.type.ordinal = data.type.ordinal
		);

	survival2 <- SIMMS::calculate.meta.survival(
		feature.name = feature2,
		expression.data = expression.data,
		survival.data = survival.data,
		data.type.ordinal = data.type.ordinal
		);

	# return the survival statistics
	return(
		list(
			"cox.uv.1" = survival1,
			"cox.uv.2" = survival2,
			"cox.int" = c(interaction.HR, interaction.P)
			)
		);

	}
