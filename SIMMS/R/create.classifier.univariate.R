create.classifier.univariate <- function(data.directory = ".", output.directory = ".", feature.selection.datasets = NULL, feature.selection.p.threshold = 0.05, training.datasets = NULL, validation.datasets = NULL, top.n.features = 25, models = c("1", "2", "3")) {

	# output directories
	out.dir <- paste(output.directory, "/output/", sep = "");
	graphs.dir <- paste(output.directory, "/graphs/", sep = "");

	# program data files and initialise variables
	all.feature.selection.names <- paste(sort(feature.selection.datasets), collapse="_");
	all.training.names <- paste(sort(training.datasets), collapse="_");
	top.subnets <- list();
	all.subnet.scores <- list();
	coxph.header <- c("HR", "95l", "95u", "P", "n");

	# lets read top subnets for the give feature.selection.datasets
	for (model in models) {

		# read in all subnets scores
		top.subnets[[model]] <- as.matrix(
			read.table(
				file = paste(out.dir, "/top_subnets_score__TRAINING_", all.feature.selection.names, "__model_", model, "__PV_", feature.selection.p.threshold, ".txt", sep = ""),
				header = TRUE,
				row.names = 1,
				sep = "\t"
				)
			);

		# read training and validation datasets
		for (dataset in c(training.datasets, validation.datasets)) {

			all.subnet.scores[[dataset]][[model]] <- as.matrix(
				read.table(
					file = paste(out.dir, "patientwise_subnets_score__", dataset, "__TRAINING_", all.feature.selection.names, "__model_", model, ".txt", sep = ""),
					header = TRUE,
					row.names = 1,
					sep = "\t"
					)
				);
			}
		}

	for (model in models) {

		coxph.uv <- list();
		risk.groups.uv <- list();
		risk.scores.uv <- list();
		training.medians <- vector();
		feature.names <- rownames(top.subnets[[model]])[1:top.n.features];

		# go over all the features one by one
		for (feature.name in feature.names) {

			# to store combined risk score of all TRAINING datasets
			all.training.risk.scores <- NULL;
			all.training.groups <- NULL;
			all.training.survtime <- NULL;
			all.training.survstat <- NULL;

			# lets process the training datasets
			for (dataset in training.datasets) {
				risk.scores <- all.subnet.scores[[dataset]][[model]][feature.name, ];
				risk.groups <- SIMMS::dichotomize.dataset(risk.scores);
				survtime <- all.subnet.scores[[dataset]][[model]]["survtime", ];
				survstat <- all.subnet.scores[[dataset]][[model]]["survstat", ];

				risk.groups.uv[[dataset]] <- cbind(risk.groups.uv[[dataset]], risk.groups);
				colnames(risk.groups.uv[[dataset]])[ncol(risk.groups.uv[[dataset]])] <- feature.name;

				risk.scores.uv[[dataset]] <- cbind(risk.scores.uv[[dataset]], risk.scores);
				colnames(risk.scores.uv[[dataset]])[ncol(risk.scores.uv[[dataset]])] <- feature.name;

				coxph.res <- SIMMS::fit.coxmodel(risk.groups, Surv(survtime, survstat))$cox.stats;
				names(coxph.res) <- coxph.header;
				coxph.uv[[dataset]] <- rbind(coxph.uv[[dataset]], coxph.res);
				rownames(coxph.uv[[dataset]])[nrow(coxph.uv[[dataset]])] <- feature.name;

				all.training.risk.scores <- c(all.training.risk.scores, risk.scores);
				all.training.survtime <- c(all.training.survtime, survtime);
				all.training.survstat <- c(all.training.survstat, survstat);
				}

			# combined training set
			all.training.groups <- SIMMS::dichotomize.dataset(all.training.risk.scores);

			risk.groups.uv[["all_training"]] <- cbind(risk.groups.uv[["all_training"]], all.training.groups);
			colnames(risk.groups.uv[["all_training"]])[ncol(risk.groups.uv[["all_training"]])] <- feature.name;

			risk.scores.uv[["all_training"]] <- cbind(risk.scores.uv[["all_training"]], all.training.risk.scores);
			colnames(risk.scores.uv[["all_training"]])[ncol(risk.scores.uv[["all_training"]])] <- feature.name;

			coxph.res <- SIMMS::fit.coxmodel(all.training.groups, Surv(all.training.survtime, all.training.survstat))$cox.stats;
			names(coxph.res) <- coxph.header;
			coxph.uv[["all_training"]] <- rbind(coxph.uv[["all_training"]], coxph.res);
			rownames(coxph.uv[["all_training"]])[nrow(coxph.uv[["all_training"]])] <- feature.name;

			training.medians <- c(training.medians, median(all.training.risk.scores));
			names(training.medians)[length(training.medians)] <- feature.name;
			}

		# write training results to file system
		for (dataset in c(training.datasets, "all_training")) {
			survtime <- all.subnet.scores[[dataset]][[model]]["survtime", ];
			survstat <- all.subnet.scores[[dataset]][[model]]["survstat", ];

			if (dataset == "all_training") {
				survtime <- all.training.survtime;
				survstat <- all.training.survstat;
				}

			write.table(
				x = cbind("survtime" = survtime, "survstat" = survstat, risk.groups.uv[[dataset]]),
				file = paste(out.dir, "riskgroups_uv__", dataset, "__TRAINING_", all.training.names, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);

			write.table(
				x = cbind("survtime" = survtime, "survstat" = survstat, risk.scores.uv[[dataset]]),
				file = paste(out.dir, "riskscores_uv__", dataset, "__TRAINING_", all.training.names, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);

			write.table(
				x = cbind(coxph.uv[[dataset]], "Q" = p.adjust(coxph.uv[[dataset]][, "P"], method = "BH")),
				file = paste(out.dir, "coxph_uv__", dataset, "__TRAINING_", all.training.names, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);
			}

		# save training set medians of top n features, just in case
		write.table(
			x = training.medians,
			file =  paste(out.dir, "riskscore_uv_median__", "all_training", "__TRAINING_", all.training.names, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
			row.names = TRUE,
			col.names = NA,
			sep = "\t"
			);

		# get risk score of VALIDATION datasets
		coxph.uv <- list();
		risk.groups.uv <- list();
		risk.scores.uv <- list();

		# go over all the features one by one (not doing it in the 
		# for loop above as users might wanna mix training & validation datasets)
		for (feature.name in feature.names) {

			# to store combined risk score of all VALIDATION datasets
			all.validation.risk.scores <- NULL;
			all.validation.groups <- NULL;
			all.validation.survtime <- NULL;
			all.validation.survstat <- NULL;

			# let's process the validation datasets
			for (dataset in validation.datasets) {
				risk.scores <- all.subnet.scores[[dataset]][[model]][feature.name, ];
				risk.groups <- SIMMS::dichotomize.dataset(risk.scores, split.at = training.medians[feature.name]);
				survtime <- all.subnet.scores[[dataset]][[model]]["survtime", ];
				survstat <- all.subnet.scores[[dataset]][[model]]["survstat", ];

				risk.groups.uv[[dataset]] <- cbind(risk.groups.uv[[dataset]], risk.groups);
				colnames(risk.groups.uv[[dataset]])[ncol(risk.groups.uv[[dataset]])] <- feature.name;

				risk.scores.uv[[dataset]] <- cbind(risk.scores.uv[[dataset]], risk.scores);
				colnames(risk.scores.uv[[dataset]])[ncol(risk.scores.uv[[dataset]])] <- feature.name;

				coxph.res <- SIMMS::fit.coxmodel(risk.groups, Surv(survtime, survstat))$cox.stats;
				names(coxph.res) <- coxph.header;
				coxph.uv[[dataset]] <- rbind(coxph.uv[[dataset]], coxph.res);
				rownames(coxph.uv[[dataset]])[nrow(coxph.uv[[dataset]])] <- feature.name;

				all.validation.risk.scores <- c(all.validation.risk.scores, risk.scores);
				all.validation.survtime <- c(all.validation.survtime, survtime);
				all.validation.survstat <- c(all.validation.survstat, survstat);

				}

			# combined validation set
			all.validation.groups <- SIMMS::dichotomize.dataset(all.validation.risk.scores, split.at = training.medians[feature.name]);

			risk.groups.uv[["all_validation"]] <- cbind(risk.groups.uv[["all_validation"]], all.validation.groups);
			colnames(risk.groups.uv[["all_validation"]])[ncol(risk.groups.uv[["all_validation"]])] <- feature.name;

			risk.scores.uv[["all_validation"]] <- cbind(risk.scores.uv[["all_validation"]], all.validation.risk.scores);
			colnames(risk.scores.uv[["all_validation"]])[ncol(risk.scores.uv[["all_validation"]])] <- feature.name;

			coxph.res <- SIMMS::fit.coxmodel(all.validation.groups, Surv(all.validation.survtime, all.validation.survstat))$cox.stats;
			names(coxph.res) <- coxph.header;
			coxph.uv[["all_validation"]] <- rbind(coxph.uv[["all_validation"]], coxph.res);
			rownames(coxph.uv[["all_validation"]])[nrow(coxph.uv[["all_validation"]])] <- feature.name;
			}

		# write validation results to file system
		for (dataset in c(validation.datasets, "all_validation")) {
			survtime <- all.subnet.scores[[dataset]][[model]]["survtime", ];
			survstat <- all.subnet.scores[[dataset]][[model]]["survstat", ];

			if (dataset == "all_validation") {
				survtime <- all.validation.survtime;
				survstat <- all.validation.survstat;
				}

			write.table(
				x = cbind("survtime" = survtime, "survstat" = survstat, risk.groups.uv[[dataset]]),
				file = paste(out.dir, "riskgroups_uv__", dataset, "__TRAINING_", all.training.names, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);

			write.table(
				x = cbind("survtime" = survtime, "survstat" = survstat, risk.scores.uv[[dataset]]),
				file = paste(out.dir, "riskscores_uv__", dataset, "__TRAINING_", all.training.names, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);

			write.table(
				x = cbind(coxph.uv[[dataset]], "Q" = p.adjust(coxph.uv[[dataset]][, "P"], method = "BH")),
				file = paste(out.dir, "coxph_uv__", dataset, "__TRAINING_", all.training.names, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);
			}

		}
	}