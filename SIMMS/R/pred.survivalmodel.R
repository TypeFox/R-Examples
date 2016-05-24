pred.survivalmodel <- function(data.directory = ".", output.directory = ".", feature.selection.datasets = NULL, feature.selection.p.threshold = 0.05, training.datasets = NULL, validation.datasets = NULL, top.n.features = 25, models = c("1", "2", "3"), write.risk.data = TRUE) {

	# output directories
	out.dir <- paste(output.directory, "/output/", sep = "");
	graphs.dir <- paste(output.directory, "/graphs/", sep = "");

	# program data files and initialise variables
	all.feature.selection.names <- paste(sort(feature.selection.datasets), collapse="_");
	all.training.names <- paste(sort(training.datasets), collapse="_");
	top.subnets <- list();
	stepAIC.betas <- NULL;
	all.subnet.scores <- list();
	all.training.median = NULL;
	coxph.header <- c("HR", "95l", "95u", "P", "n");

	# lets read in the training & validation datasets
	for (dataset in c(training.datasets, validation.datasets)) {
		for (model in models) {
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

	# lets read in the selected subnets and their betas
	for (model in models) {

		for (direction in c("forward", "backward")) {
			stepAIC.betas <- as.matrix(
				read.table(
					file = paste(out.dir, "beta__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep = ""),
					header = TRUE,
					row.names = 1,
					sep = "\t"
					)
				);
			selected.variables <- rownames(stepAIC.betas);

			# compute risk score of TRAINING datasets
			all.training.risk.scores <- NULL;
			all.training.groups <- NULL;
			all.training.survtime <- NULL;
			all.training.survstat <- NULL;

			for (dataset in training.datasets) {

				# multiply betas of this variable with corresponding per-patient scores
				risk.scores <- stepAIC.betas[selected.variables, "beta"] *
					all.subnet.scores[[dataset]][[model]][selected.variables, ];

				# if variables > 1, then sum up all the scores once its times corresponding beta
				if(length(selected.variables) > 1) {risk.scores <- apply(risk.scores, 2, sum)};

				risk.groups <- SIMMS::dichotomize.dataset(risk.scores);
				names(risk.groups) <- colnames(all.subnet.scores[[dataset]][[model]]);
				survtime <- all.subnet.scores[[dataset]][[model]]["survtime", ];
				survstat <- all.subnet.scores[[dataset]][[model]]["survstat", ];

				coxph.res <- SIMMS::fit.coxmodel(risk.groups, Surv(survtime, survstat))$cox.stats;
				if (!is.na(coxph.res[4])) {
					survdiff.res <- get.chisq.stats(risk.groups, Surv(survtime, survstat));
					}
				names(coxph.res) <- coxph.header;

				# save the results to file system for later visualizations (either BL or ordinary graphics)
				if (write.risk.data == TRUE) {
					write.table(
						x = cbind("riskgroup" = risk.groups, "survtime" = survtime, "survstat" = survstat),
						file = paste(out.dir, "riskgroups__", dataset, "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
						row.names = TRUE,
						col.names = NA,
						sep = "\t"
						);

					write.table(
						x = cbind("riskscore" = risk.scores, "survtime" = survtime, "survstat" = survstat),
						file = paste(out.dir, "riskscores__", dataset, "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
						row.names = TRUE,
						col.names = NA,
						sep = "\t"
						);
					}

				write.table(
					x = t(coxph.res),
					file = paste(out.dir, "coxph__", dataset, "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);

				if (!is.na(coxph.res[4])) {
					write.table(
						x = t(survdiff.res),,
						file = paste(out.dir, "survdiff__", dataset, "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
						row.names = TRUE,
						col.names = NA,
						sep = "\t"
						);
					}

				all.training.risk.scores <- c(all.training.risk.scores, risk.scores);
				all.training.survtime <- c(all.training.survtime, survtime);
				all.training.survstat <- c(all.training.survstat, survstat);
				}

			# combined training set
			all.training.groups <- SIMMS::dichotomize.dataset(all.training.risk.scores);
			coxph.res <- SIMMS::fit.coxmodel(all.training.groups, Surv(all.training.survtime, all.training.survstat))$cox.stats;
			if (!is.na(coxph.res[4])) {
				survdiff.res <- get.chisq.stats(all.training.groups, Surv(all.training.survtime, all.training.survstat));
				}
			names(coxph.res) <- coxph.header;
	
			if (write.risk.data == TRUE) {
				write.table(
					x = cbind("riskgroup" = all.training.groups, "survtime" = all.training.survtime, "survstat" = all.training.survstat),
					file = paste(out.dir, "riskgroups__", "all_training", "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);

				write.table(
					x = cbind("riskscore" = all.training.risk.scores, "survtime" = all.training.survtime, "survstat" = all.training.survstat),
					file = paste(out.dir, "riskscores__", "all_training", "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);
				}

			write.table(
				x = t(coxph.res),
				file = paste(out.dir, "coxph__", "all_training", "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);

			if (!is.na(coxph.res[4])) {
				write.table(
					x = t(survdiff.res),
					file = paste(out.dir, "survdiff__", "all_training", "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);
				}

			all.training.median <- median(all.training.risk.scores);

			# store all_training set median score
			write(
				x = all.training.median,
				file =  paste(out.dir, "riskscore_median__", "all_training", "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep="")
				);
			
			# compute risk score of VALIDATION datasets
			all.validation.risk.scores <- NULL;
			all.validation.groups <- NULL;
			all.validation.survtime <- NULL;
			all.validation.survstat <- NULL;

			for (dataset in validation.datasets) {

				# multiply betas of this variable with corresponding per-patient scores
				risk.scores <- stepAIC.betas[selected.variables, "beta"] *
					all.subnet.scores[[dataset]][[model]][selected.variables, ];
				if(length(selected.variables) > 1) {risk.scores <- apply(risk.scores, 2, sum)};

				risk.groups <- SIMMS::dichotomize.dataset(risk.scores, split.at = all.training.median);
				names(risk.groups) <- colnames(all.subnet.scores[[dataset]][[model]]); 
				survtime <- all.subnet.scores[[dataset]][[model]]["survtime", ];
				survstat <- all.subnet.scores[[dataset]][[model]]["survstat", ];


				coxph.res <- SIMMS::fit.coxmodel(risk.groups, Surv(survtime, survstat))$cox.stats;
				if (!is.na(coxph.res[4])) {
					survdiff.res <- get.chisq.stats(risk.groups, Surv(survtime, survstat));
					}
				names(coxph.res) <- coxph.header;

				# save the results to file system for later visualizations (either BL or ordinary graphics)
				if (write.risk.data == TRUE) {
					write.table(
						x = cbind("riskgroup" = risk.groups, "survtime" = survtime, "survstat" = survstat),
						file = paste(out.dir, "riskgroups__", dataset, "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
						row.names = TRUE,
						col.names = NA,
						sep = "\t"
						);

					write.table(
						x = cbind("riskscore" = risk.scores, "survtime" = survtime, "survstat" = survstat),
						file = paste(out.dir, "riskscores__", dataset, "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
						row.names = TRUE,
						col.names = NA,
						sep = "\t"
						);
					}

				write.table(
					x = t(coxph.res),
					file = paste(out.dir, "coxph__", dataset, "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);
				
				if (!is.na(coxph.res[4])) {
					write.table(
						x = t(survdiff.res),
						file = paste(out.dir, "survdiff__", dataset, "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
						row.names = TRUE,
						col.names = NA,
						sep = "\t"
						);
					}
				
				all.validation.risk.scores <- c(all.validation.risk.scores, risk.scores);
				all.validation.survtime <- c(all.validation.survtime, survtime);
				all.validation.survstat <- c(all.validation.survstat, survstat);
				}

			# combined validation set
			all.validation.groups <- SIMMS::dichotomize.dataset(all.validation.risk.scores, split.at = all.training.median);
			coxph.res <- SIMMS::fit.coxmodel(all.validation.groups, Surv(all.validation.survtime, all.validation.survstat))$cox.stats;
			if (!is.na(coxph.res[4])) {
				survdiff.res <- get.chisq.stats(all.validation.groups, Surv(all.validation.survtime, all.validation.survstat));
				}
			names(coxph.res) <- coxph.header;

			if (write.risk.data == TRUE) {
				write.table(
					x = cbind("riskgroup" = all.validation.groups, "survtime" = all.validation.survtime, "survstat" = all.validation.survstat),
					file = paste(out.dir, "riskgroups__", "all_validation", "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);

				write.table(
					x = cbind("riskscore" = all.validation.risk.scores, "survtime" = all.validation.survtime, "survstat" = all.validation.survstat),
					file = paste(out.dir, "riskscores__", "all_validation", "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);
				}

			write.table(
				x = t(coxph.res),
				file = paste(out.dir, "coxph__", "all_validation", "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t"
				);

			if (!is.na(coxph.res[4])) {
				write.table(
					x = t(survdiff.res),
					file = paste(out.dir, "survdiff__", "all_validation", "__TRAINING_", all.training.names, "__", direction, "__model_", model, "__top_", top.n.features, ".txt", sep=""),
					row.names = TRUE,
					col.names = NA,
					sep = "\t"
					);
				}

			}
		}

	# let's apply univariate cox model to each of the top features, independently, for all models
	# let's read top subnets for the give feature.selection.datasets
	for (model in models) {
		top.subnets[[model]] <- as.matrix(
			read.table(
				file = paste(out.dir, "top_subnets_score__TRAINING_", all.feature.selection.names, "__model_", model, "__PV_", feature.selection.p.threshold, ".txt", sep = ""),
				header = TRUE,
				row.names = 1,
				sep = "\t"
				)
			);
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
