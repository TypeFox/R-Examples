calculate.sensitivity.stats <- function(all.data = NULL) {

	# to store output of this method
	sen.spe <- rep(NA, 7);

	names(sen.spe) <- c("TP", "FP", "TN", "FN", "Sensitivity", "Specificity", "Accuracy");

	# TP
	sen.spe["TP"] <- nrow(
		all.data[
			all.data[, "riskgroup"] ==  "1"
			&
			all.data[, "real.survival.group"] ==  "high",
			]
		);
		
	# FP
	sen.spe["FP"] <- nrow(
		all.data[
			all.data[, "riskgroup"] ==  "1"
			&
			all.data[, "real.survival.group"] ==  "low",
			]
		);

	# TN
	sen.spe["TN"] <- nrow(
		all.data[
			all.data[, "riskgroup"] ==  "0"
			&
			all.data[, "real.survival.group"] ==  "low",
			]
		);

	# FN
	sen.spe["FN"] <- nrow(
		all.data[
			all.data[, "riskgroup"] ==  "0"
			&
			all.data[, "real.survival.group"] ==  "high",
			]
		);

	# Sensitivity
	sen.spe["Sensitivity"] <- sen.spe["TP"] / (sen.spe["TP"] + sen.spe["FN"]) * 100;

	# Specificity
	sen.spe["Specificity"] <- sen.spe["TN"] / (sen.spe["TN"] + sen.spe["FP"]) * 100;

	# Overall accuracy
	sen.spe["Accuracy"] <- (sen.spe["TP"] + sen.spe["TN"]) / (sen.spe["TP"] + sen.spe["TN"] + sen.spe["FP"] + sen.spe["FN"]) * 100;
	
	return (sen.spe);
	}