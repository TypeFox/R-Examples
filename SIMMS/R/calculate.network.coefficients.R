calculate.network.coefficients <- function(data.directory = ".", output.directory = ".", training.datasets = NULL, data.types = c("mRNA"), data.types.ordinal = c("cnv"), subnets.file.flattened = NULL, truncate.survival = 100, subset = NULL) {

	all.training.names <- paste(sort(training.datasets), collapse="_");
	gene.pairs <- read.table(
		file = subnets.file.flattened,
		header = TRUE,
		row.names = NULL,
		sep = "\t",
		as.is = TRUE
		);

	# focus only on unique records
	gene.pairs <- unique(gene.pairs);
	subnets.genes <- unique(c(gene.pairs$GeneID1, gene.pairs$GeneID2));
	subnets.genes <- paste(subnets.genes, "_at", sep = "");
	# PCB: does the above line this only works with Affymetrix data or data with an _at trailing and a unique ID before it?

	cancer.data <- load.cancer.datasets(
		truncate.survival = truncate.survival,
		datasets.to.load = training.datasets, 
		data.types = data.types, 
		data.directory = data.directory,
		subset = subset
		);

	# DATA PROCESSING
	scaled.data <- cancer.data;

	# SCALE DATA
	for (data.type in data.types) {
		for (i in 1:length(scaled.data[["all.data"]][[data.type]]) ) {
			if (data.type %in% data.types.ordinal) {
				scaled.data[["all.data"]][[data.type]][[i]] <- as.data.frame( 
					scaled.data[["all.data"]][[data.type]][[i]]
					);
				}
			else {
				scaled.data[["all.data"]][[data.type]][[i]] <- as.data.frame( 
					t( scale( t(scaled.data[["all.data"]][[data.type]][[i]]) ) ) 
					);
				}
			}
		}

	# ANALYZE EACH SAMPLE
	# make the ProbeSet IDs
	gene.pairs$ProbeID1 <- paste(as.character(gene.pairs$GeneID1), "_at", sep = "");
	gene.pairs$ProbeID2 <- paste(as.character(gene.pairs$GeneID2), "_at", sep = "");

	# initialise return oject
	nodes.edges.stats <- list();

	for (data.type in data.types) {

		data.type.ordinal <- FALSE;
		if (data.type %in% data.types.ordinal) {
			data.type.ordinal <- TRUE;
			}

		coxph.nodes <- matrix(
			data = NA,
			nrow = length(subnets.genes),
			ncol = 2,
			dimnames = list(
				subnets.genes,
				c("coef", "P")
				)
			);

		coxph.nodes.obj <- list();

		coxph.edges.coef <- matrix(
			data = 1,
			nrow = length(subnets.genes),
			ncol = length(subnets.genes),
			dimnames = list(
				subnets.genes,
				subnets.genes
				)
			);

		# just copy the structure to get the p-value matrix
		coxph.edges.P <- coxph.edges.coef;

		# fill in the empty objects
		for (i in 1:nrow(gene.pairs)) {

			results <- SIMMS::fit.interaction.model(
				feature1 = gene.pairs$ProbeID1[i],
				feature2 = gene.pairs$ProbeID2[i],
				expression.data = scaled.data[["all.data"]][[data.type]],
				survival.data = scaled.data$all.survobj,
				data.type.ordinal = data.type.ordinal
				);

			#cat("\nNew pair\t", gene.pairs$ProbeID1[i], "\t", gene.pairs$ProbeID2[i]);
			#print(results);

			# isolate results of g1 and g2
			if (is.list(results)) {

				if (!data.type.ordinal) {
					results[["cox.uv.1"]][["cox.stats"]][1] <- log2(
						results[["cox.uv.1"]][["cox.stats"]][1]
						);
					results[["cox.uv.2"]][["cox.stats"]][1] <- log2(
						results[["cox.uv.2"]][["cox.stats"]][1]
						);
					results.g1 <- results[["cox.uv.1"]][["cox.stats"]];
					results.g2 <- results[["cox.uv.2"]][["cox.stats"]];
					}
				else {
					results.gx <- list();
					for (cox.uv.i in c("cox.uv.1", "cox.uv.2")) {

						# check if cox model failed
						if (length(results[[cox.uv.i]][["cox.obj"]]) > 1 && !is.na(results[[cox.uv.i]][["cox.obj"]])) {
							results[[cox.uv.i]][["cox.stats"]][, "HR"] <- log2(
								results[[cox.uv.i]][["cox.stats"]][, "HR"]
								);
							colnames(results[[cox.uv.i]][["cox.stats"]])[1] <- "coef";

							# switch signs for levels less that zero because in cox model, zero was baseline
							# and switching of signs is needed as e.g deletions are smaller in number than zero
							neg.groups <- which(as.numeric(rownames(results[[cox.uv.i]][["cox.stats"]])) < 0);
							if (length(neg.groups) > 0) { 
								results[[cox.uv.i]][["cox.stats"]][neg.groups, "coef"] <- 
									-results[[cox.uv.i]][["cox.stats"]][neg.groups, "coef"];
								}

							# if ordinal data type - pick the smallest P
							min.p <- which(
								results[[cox.uv.i]][["cox.stats"]][, "p"] ==
									min(results[[cox.uv.i]][["cox.stats"]][, "p"], na.rm = T)
								);
							if (length(min.p) > 0) {
								results.gx[[cox.uv.i]] <- results[[cox.uv.i]][["cox.stats"]][min.p, ];
								}
							else {
								results.gx[[cox.uv.i]] <- NA;
								}
							}
						else {
							results.gx[[cox.uv.i]] <- NA;
							}
						}
					results.g1 <- results.gx[["cox.uv.1"]];
					results.g2 <- results.gx[["cox.uv.2"]];
					}

				# sometimes same gene is in multiple interactions and hence return all NULL row
				# and with other interactions, return numeric HR and we would like to keep numeric
				# HRs not NA when numeric HR is available
				if (!is.na(results.g1[1])) { 

					# store HR and P
					coxph.nodes[gene.pairs$ProbeID1[i], "coef"] <- results.g1[1];
					coxph.nodes[gene.pairs$ProbeID1[i], "P"]  <- results.g1[4];

					# store cox fit objects
					coxph.nodes.obj[[gene.pairs$ProbeID1[i]]] <- results[["cox.uv.1"]][["cox.stats"]];
					}

				if (!is.na(results.g2[1])) { 

					# store HR and P
					coxph.nodes[gene.pairs$ProbeID2[i], "coef"] <- results.g2[1];
					coxph.nodes[gene.pairs$ProbeID2[i], "P"]  <- results.g2[4];

					# store cox fit objects
					coxph.nodes.obj[[gene.pairs$ProbeID2[i]]] <- results[["cox.uv.2"]][["cox.stats"]];				
					}

				# store the edges HR and P in a matrix - for faster lookups
				coxph.edges.coef[gene.pairs$ProbeID1[i], gene.pairs$ProbeID2[i]] <-
					log2(results[["cox.int"]][1]);
				coxph.edges.coef[gene.pairs$ProbeID2[i], gene.pairs$ProbeID1[i]] <-
					log2(results[["cox.int"]][1]);
				coxph.edges.P[ gene.pairs$ProbeID1[i], gene.pairs$ProbeID2[i]] <- 
					results[["cox.int"]][2];
				coxph.edges.P[ gene.pairs$ProbeID2[i], gene.pairs$ProbeID1[i]] <- 
					results[["cox.int"]][2];

				}
			}

		# populate return object
		nodes.edges.stats[[data.type]][["nodes.coef"]] <- coxph.nodes;
		nodes.edges.stats[[data.type]][["edges.coef"]] <- coxph.edges.coef;
		nodes.edges.stats[[data.type]][["edges.P"]] <- coxph.edges.P;
		nodes.edges.stats[[data.type]][["cox.uv"]] <- coxph.nodes.obj;
		}

	return(nodes.edges.stats);

	}
