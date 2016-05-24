load.cancer.datasets <- function(tumour.only = TRUE, with.survival.only = TRUE, truncate.survival = 100, datasets.to.load = 'all', data.types = c('mRNA'), datasets.file = 'datasets.txt', data.directory = '.', verbose = FALSE, subset = NULL) {

	# read in the listing of all datasets
	datasets <- read.table(
		file = paste(data.directory, "/", datasets.file, sep = ''),
		header = TRUE,
		sep = "\t",
		row.names = NULL,
		as.is = TRUE
		);

	# handle the user requesting all datasets
	if ('all' %in% datasets.to.load) { datasets.to.load <- datasets$dataset; }

	# make sure all requested datasets exist
	if ( !all( datasets.to.load %in% datasets$dataset ) ) { stop('Non-existant dataset requested'); }

	# initialize data objects
	all.data <- list();
	all.survobj <- list();
	all.probesets <- vector();

	# load all requested datasets one-by-one
	for (dataset in datasets.to.load) {

		# give the user some status info if requested
		if (verbose) { cat('\nreading annotations for dataset: ', dataset); }

		# set up the path to the dataset's files
		dataset.directory <- paste(data.directory, dataset, '/', sep = '');

		annotation <- read.table(
			paste(dataset.directory, datasets$annotation[datasets$dataset == dataset], sep = ''),
			header = TRUE,
			row.names = 1,
			sep = "\t"
			);

		# only extract tumour samples
		if (tumour.only) {
			annotation <- annotation[annotation$Tumour == "Yes",];
			}

		# if a subset of samples has been specified, discard any samples not belonging to the subset
		if (!is.null(subset)) {
			annotation <- annotation[annotation[[subset$Field]] == subset$Entry, ];
			}

		# select the appropriate survtime and survstat variable for each dataset
		if (!("survstat" %in% colnames(annotation))) {
			annotation$survstat <- annotation[,datasets$survstat[datasets$dataset == dataset]];
			}
		if (!("survtime" %in% colnames(annotation))) {
			annotation$survtime <- annotation[,datasets$survtime[datasets$dataset == dataset]];
			}
		if (!("survtime.unit" %in% colnames(annotation))) {
			annotation$survtime.unit <- annotation[,datasets$survtime.unit[datasets$dataset == dataset]];
			}

		# handle survtime <= 0
		annotation$survtime[annotation$survtime <= 0] <- 1e-05;

		# only keep samples with survival data
		if (with.survival.only) {
			annotation <- annotation[!is.na(annotation$survtime) & !is.na(annotation$survstat),];
			}

		# make R like rownames
		rownames(annotation) <- make.names(as.character(rownames(annotation)));

		common.samples <- rownames(annotation);
		for (data.type in data.types) {
			# ensure the data type has a corresponding column with exact name in datasets.txt
			tryCatch(
				expr = {
					all.data[[data.type]][[dataset]] <- read.table(
						paste(dataset.directory, datasets[[data.type]][datasets$dataset == dataset], sep = ""),
						header = TRUE,
						row.names = 1,
						sep = "\t"
						)
					},
				error = function(ex) {
					stop("\nWell... you asked for data.type: [", data.type, "] but i cant find a column named as this data.type in datasets.txt - hence dieing");
					}
				);

			# remove samples (columns) which have NA for all genes
			all.data[[data.type]][[dataset]] <- all.data[[data.type]][[dataset]][, 
				apply(
					all.data[[data.type]][[dataset]],
					2,
					function(x) any(!is.na(x))
					)
				];
			
			# compile a list of samples that are common across all molecular features & survival
			common.samples <- intersect(
				common.samples, 
				colnames(all.data[[data.type]][[dataset]])
				);
			}

		# limit to common samples only
		annotation <- annotation[common.samples,];

		# only keep the dataset if it has at least one sample
		if (nrow(annotation) >= 1) {
			for (data.type in data.types) {
				# ensure equivalent sorting of the annotation and data objects
				tryCatch(
					expr = { all.data[[data.type]][[dataset]] <- all.data[[data.type]][[dataset]][, common.samples] }, 
					error = function(ex) {
						cat("\n\nSome columns in the data file do not match the rows in the annotation	 file\n\n");
						}
					);

				all.probesets <- unique( c(all.probesets, rownames(all.data[[data.type]][[dataset]])) );
				}

			# save survival objects
			all.survobj[[dataset]] <- SIMMS::create.survobj(annotation = annotation, truncate.survival = truncate.survival);
			}
		else {
			stop("\n\nDataset: [", dataset, "] does not have any valid data to process. Please remove from feature selection, training and validation dataset vectors\n\n");
			}
		}
	
	# return the final object
	return(
		list(
			all.data = all.data,
			all.survobj = all.survobj,
			all.probesets = all.probesets
			)
		);

	}
