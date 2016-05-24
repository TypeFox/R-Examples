create.sensitivity.plot <- function(riskscore = NULL, riskgroup = NULL, survtime =  NULL, survstat = NULL, survtime.cutoffs = c(seq(5,10,1)), output.directory = ".", file.stem = NULL, main.title = "", resolution = 100) {

	all.data <- data.frame("riskscore" = riskscore, "riskgroup" = riskgroup, "survtime" = survtime, "survstat" =  survstat);

	# output directories
	out.dir <- paste(output.directory, "/output/", sep = "");
	graphs.dir <- paste(output.directory, "/graphs/", sep = "");

	# add a column for hypothetical class label based on survtime based dichotomization
	all.data <- cbind(all.data, "real.survival.group" = rep(NA, nrow(all.data)));

	# to store results of all cutoff time points
	all.sen.spe <- list();
	survtime.cutoffs <- unique(survtime.cutoffs);

	# set the graphics driver
	current.type <- getOption("bitmapType");
	options(bitmapType = "cairo");

	png(
		filename = paste(graphs.dir, file.stem, ".png", sep = ""),
		height = 5,
		width = 5,
		units = "in",
		res = resolution
		#,compression = "lzw"
		);
	par(
		mfrow = c(1,1),
		tcl = -0.1, 
		las = 1,
		font.lab = 2,
		font.axis = 2,
		mar = c(2.5, 2.3, 1, 0.2),
		mgp = c(1.4, 0.3, 0),
		oma = c(0, 0, 0, 0)
		);
	# colour scheme for each ROC curve
	colour <- 1;

	for (cutoff in survtime.cutoffs) {

		all.data.trimmed <- NULL;

		# remove patients that have survtime < cutoff and no event (i-e censored before cutoff survtime)
		all.data.trimmed <- all.data[
			!(as.numeric(all.data[, "survtime"]) <= cutoff
			&
			all.data[, "survstat"] ==  "0")
			, ];

		# assign risk groups based on survival time
		all.data.trimmed[
			as.numeric(all.data.trimmed[, "survtime"]) <= cutoff
			&
			all.data.trimmed[, "survstat"] == "1"
			, "real.survival.group"] <- "high"; # 1

		all.data.trimmed[
			as.numeric(all.data.trimmed[, "survtime"]) > cutoff
			, "real.survival.group"] <- "low";  # 0

		# PCB: this whole section (above and below) might be easier to read if you abstracted all the indexing into separate lines
		# high.risk.patients <- as.numeric(all.data.trimmed[, "survtime"]) <= cutoff & all.data.trimmed[, "survstat"] == "1"
		# all.data.trimmed[high.risk.patients, "real.survival.group"] <- 'high';

		# lets compute the sensitivity stats for predefined (median based) risk groups
		all.sen.spe[[as.character(cutoff)]] <- calculate.sensitivity.stats(all.data.trimmed);

		# lets compute ROC
		x.FPR <- vector();
		y.TPR <- vector();

		# first sort the risk scores in ascending order, for sliding ROC dichotomization
		all.data.trimmed <- all.data.trimmed[order(all.data.trimmed[, "riskscore"]), ];

		# set all patients to high risk group
		all.data.trimmed[, "riskgroup"] <- 1;
		for (i.roc in 1:nrow(all.data.trimmed) ) {

			# slide the dichotomization window
			all.data.trimmed[i.roc, "riskgroup"] <- "0";

			sensitivity.stats <- calculate.sensitivity.stats(all.data.trimmed);

			x.FPR <- c(
				x.FPR,
				( 1 - (sensitivity.stats["Specificity"] / 100) )
				);

			y.TPR <- c(
				y.TPR,
				(sensitivity.stats["Sensitivity"] / 100)
				);
			}

		# lets plot the curve
		if (colour == 1) {
			plot(x.FPR, y.TPR, xlab = "False positive rate", ylab = "True positive rate", type = "l", col = colour, lwd = 2);
			}
		else {
			lines(x.FPR, y.TPR, col = colour, lwd = 2);
			}

		colour <- colour + 1;

		}

	title(main.title, cex.main = 1, font.main = 2);

	# add figure legends
	par(mar = c(0, 0, 0, 0), font = 2);
	legend("bottomright", lty = 1, cex = 0.8, legend = paste(survtime.cutoffs, "yr."), col = 1:(colour - 1), bty = "n", lwd = 2);

	dev.off();
	options(bitmapType = current.type);

	# store the results
	all.sen.spe <- t(as.data.frame(all.sen.spe));
	rownames(all.sen.spe) <- paste("cutoff_", survtime.cutoffs, sep = "");

	write.table(
		all.sen.spe,
		file = paste(out.dir, file.stem, ".txt", sep=""),
		row.names = TRUE,
		col.names = NA,
		sep = "\t"
		);
	}
