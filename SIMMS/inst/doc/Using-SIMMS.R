## ----xtable, echo=FALSE, results="asis"----------------------------------
library(xtable);

x <- matrix(
	data = c(
		"{ *x ∈ ℝ* : *x ≥ 0* }", "0, 2, 1.37, 7.04, 9.68", "log2 mRNA or miRNA abundance", 
		"{ *x ∈ ℝ* : *x ≥ 0* }", "0.1, 0.38, 0.78, 0.22, 0.98", "DNA methylation beta values",
		"{ *x ∈ ℤ* }", "-2, -1, 0, 1, 2", "copy-number calls. Baseline = 0 (Neutral/Diploid)",
		"{ *x ∈ ℤ* : *x ≥ 0* }", "0, 1, 2, 3", "mutation data. Baseline = 0 (Wildtype)",
		"{ *x ∈ ℤ* : *x ≠ 0* }", "-2, -1, 1, 2", "Unsupported due to missing baseline 0",
		"{ *x ∉ ℝ* }", "WT, Mutant, Gain, Deleted", "Unsupported due to alphabets"
		), 
	nrow = 6, ncol = 3, byrow = TRUE,
	dimnames = list(
		NULL,
		c("Data type", "Example data", "Example profile")
		)
	);

tab <- xtable(x);

print(tab, type = "html", include.rownames = FALSE);

## ---- results = "hide", message = FALSE, eval = TRUE---------------------
options("warn" = -1);

# load SIMMS library
library("SIMMS");

# path of the data directory containing Breastdata1/ and Breastdata2/ subdirectories
data.directory <- get.program.defaults(networks.database = "test")[["test.data.dir"]];

# path of the directory where results will be stored
output.directory <- tempdir();

# molecular profiles to be used
data.types <- c("mRNA");

# feature selection datasets
# name of the dataset directory containing mRNA abundance and annotation profiles of training dataset/s
feature.selection.datasets <- c("Breastdata1");

# model training datasets, ideally same as feature selection datasets
# name of the dataset directory containing mRNA abundance and annotation profiles of training dataset/s
training.datasets <- feature.selection.datasets;

# validation datasets
# name of the dataset directory containing mRNA abundance and annotation profiles of validation dataset/s
validation.datasets <- c("Breastdata2");

# all the possible P value thresholds that one may consider applying to feature selection process.
# its the P value of univariate (genewise) Cox model statistics
feature.selection.p.thresholds <- c(0.5);

# one of the P values above, to be used for subsequent analysis. Not a vector for performance reasons
feature.selection.p.threshold <- 0.5;

# names of the learning algorithms to be used for the final multivarite model
learning.algorithms <- c("backward", "forward", "glm");

# top features to be used for model selection (Backwards elimination, Forward selection, GLM)
# you can try a number of different model selection runs by specifying a vector of top n features
top.n.features <- c(5);

# truncate survival
truncate.survival <- 10;

# calculate per feature univariate coefficients in training sets
derive.network.features(
	data.directory = data.directory,
	output.directory = output.directory,
	data.types = data.types,
	feature.selection.datasets = feature.selection.datasets,
	feature.selection.p.thresholds = feature.selection.p.thresholds,
	networks.database = "test", # or "default" for Reactome/BioCarta/NCI-PID
	truncate.survival = truncate.survival
	);

# calculate per-subnetwork scores in both training and validation sets
prepare.training.validation.datasets(
	data.directory = data.directory,
	output.directory = output.directory,
	data.types = data.types,
	feature.selection.datasets = feature.selection.datasets,
	datasets = c(training.datasets, validation.datasets),
	networks.database = "test", # or "default" for Reactome/BioCarta/NCI-PID
	truncate.survival = truncate.survival
	);

# iterate over varying top n features, identify and validate survival models
for (top.n in top.n.features) {

	# create classifier assessing univariate prognostic power of subnetwork modules (Train and Validate)
	ret <- create.classifier.univariate(
		data.directory = data.directory,
		output.directory = output.directory,
		feature.selection.datasets = feature.selection.datasets,
		feature.selection.p.threshold = feature.selection.p.threshold,
		training.datasets = training.datasets,
		validation.datasets = validation.datasets,
		top.n.features = top.n
		);

	# create a multivariate classifier (Train and Validate)
	ret <- create.classifier.multivariate(
		data.directory = data.directory,
		output.directory = output.directory,
		feature.selection.datasets = feature.selection.datasets,
		feature.selection.p.threshold = feature.selection.p.threshold,
		training.datasets = training.datasets,
		validation.datasets = validation.datasets,
		learning.algorithms = learning.algorithms,
		top.n.features = top.n
		);

	# perform Kaplan-Meier analysis and generate plots
	create.survivalplots(
		data.directory = data.directory,
		output.directory = output.directory,
		training.datasets = training.datasets,
		validation.datasets = validation.datasets,
		top.n.features = top.n,
		learning.algorithms = learning.algorithms,
		truncate.survival = truncate.survival,
		survtime.cutoffs = c(5),
		main.title = FALSE,
		KM.plotting.fun = "create.KM.plot",
		resolution = 100
		);
	}

