ALF <- function(data, ...) UseMethod("ALF")
# This function implements the workflow used by (1) to select a suitable model for absolute label-free quantification.

# 1.	Ludwig, C., Claassen, M., Schmidt, A. \& Aebersold, R. Estimation of Absolute Protein Quantities of Unlabeled Samples by Selected Reaction Monitoring Mass Spectrometry. Molecular \& Cellular Proteomics 11, M111.013987-M111.013987 (2012).

ALF.default <- function(data, report_filename="ALF_report.pdf", prediction_filename="ALF_prediction.csv", peptide_methods = c("top"), peptide_topx = c(1,2,3), peptide_strictness = "loose", peptide_summary = "mean", transition_topx = c(1,2,3), transition_strictness = "loose", transition_summary = "sum", fasta = NA, apex_model = NA, combine_precursors = FALSE, combine_peptide_sequences = FALSE, consensus_proteins = TRUE, consensus_peptides = TRUE, consensus_transitions = TRUE, scampi_method = "LSE", scampi_iterations = 10, scampi_outliers = FALSE, scampi_outliers_iterations = 2, scampi_outliers_threshold = 2, cval_method = "boot", cval_mcx = 1000, ...) {
	pdf(file=report_filename)
	
	# nr_peptides nr_transitions tuning
	data.tune <- tune.ALF(data, peptide_methods = peptide_methods, peptide_topx = peptide_topx, peptide_strictness = peptide_strictness, peptide_summary = peptide_summary, transition_topx = transition_topx, transition_strictness = transition_strictness, transition_summary = transition_summary, fasta = fasta, apex_model = apex_model, combine_precursors = combine_precursors, combine_peptide_sequences = combine_peptide_sequences, consensus_proteins = consensus_proteins, consensus_peptides = consensus_peptides, consensus_transitions = consensus_transitions, scampi_method = scampi_method, scampi_iterations = scampi_iterations, scampi_outliers = scampi_outliers, scampi_outliers_iterations = scampi_outliers_iterations, scampi_outliers_threshold = scampi_outliers_threshold, cval_method = cval_method, cval_mcx = cval_mcx)

	transition_topx_min <- as.numeric(rownames(data.tune)[which(data.tune == min(data.tune), arr.ind = TRUE)[1]])
	peptide_method_min <- colnames(data.tune)[which(data.tune == min(data.tune), arr.ind = TRUE)[2]]

	if (!(peptide_method_min %in% c("all","iBAQ","APEX","NSAF","SCAMPI"))) {
		peptide_method_min = "top"
		peptide_topx_min <- as.numeric(strsplit(colnames(data.tune)[which(data.tune == min(data.tune), arr.ind = TRUE)[2]],"top")[[1]][2])
	}
	else {
		peptide_topx_min <- NA
	}

	performanceplot.ALF(data.tune)
		
	# calculate optimal model
	optimal.ProteinInference <- ProteinInference(data, peptide_method = peptide_method_min, peptide_topx = peptide_topx_min, peptide_strictness = peptide_strictness, peptide_summary = peptide_summary, transition_topx = transition_topx_min, transition_strictness = transition_strictness, transition_summary = transition_summary, fasta = fasta, apex_model = apex_model, combine_precursors = combine_precursors, combine_peptide_sequences = combine_peptide_sequences, consensus_proteins = consensus_proteins, consensus_peptides = consensus_peptides, consensus_transitions = consensus_transitions, scampi_method = scampi_method, scampi_iterations = scampi_iterations, scampi_outliers = scampi_outliers, scampi_outliers_iterations = scampi_outliers_iterations, scampi_outliers_threshold = scampi_outliers_threshold)
	optimal.AbsoluteQuantification <- AbsoluteQuantification(optimal.ProteinInference)
	optimal.AbsoluteQuantification <- predict(optimal.AbsoluteQuantification)
	plot(optimal.AbsoluteQuantification)
	
	optimal.AbsoluteQuantification.cval <- cval.AbsoluteQuantification(optimal.AbsoluteQuantification,method=cval_method, mcx = cval_mcx)
	plot(optimal.AbsoluteQuantification.cval)
	hist.AbsoluteQuantification(optimal.AbsoluteQuantification.cval)
	
	dev.off()
	
	export.AbsoluteQuantification(optimal.AbsoluteQuantification, file = prediction_filename)

	print(paste("PDF Report written to:",normalizePath(report_filename)))
	print(paste("CSV Report written to:",normalizePath(prediction_filename)))
}

tune.ALF <- function(data, peptide_methods = "top", peptide_topx = 2, peptide_strictness = "strict", peptide_summary = "mean", transition_topx = 3, transition_strictness = "strict", transition_summary = "sum", fasta = NA, apex_model = NA, combine_precursors = FALSE, combine_peptide_sequences = FALSE, consensus_proteins = TRUE, consensus_peptides = TRUE, consensus_transitions = TRUE, scampi_method = "LSE", scampi_iterations = 10, scampi_outliers = FALSE, scampi_outliers_iterations = 2, scampi_outliers_threshold = 2, cval_method = cval_method, cval_mcx = cval_mcx, ...) {
	if ("top" %in% peptide_methods) {
		cvmfe.mx <- matrix(nrow = length(transition_topx), ncol = length(peptide_topx)+length(peptide_methods)-1)
		peptide_topx_number <- length(peptide_topx)

	}
	else {
		cvmfe.mx <- matrix(nrow = length(transition_topx), ncol = length(peptide_methods))
		peptide_topx_number <- 0

	}
	peptide_methods_names <- c()
	rownames(cvmfe.mx) <- transition_topx

	if ("top" %in% peptide_methods) {
		i <- 1
		while (i <= length(peptide_topx)) {
			j <- 1
			while (j <= length(transition_topx)) {
				cvmfe.mx[j,i] <- cval.AbsoluteQuantification(AbsoluteQuantification(ProteinInference(data, peptide_method = "top", peptide_topx = peptide_topx[i], peptide_strictness = peptide_strictness, peptide_summary = peptide_summary, transition_topx = transition_topx[j], transition_strictness = transition_strictness, transition_summary = transition_summary, fasta = fasta, apex_model = apex_model, combine_precursors = combine_precursors, combine_peptide_sequences = combine_peptide_sequences, consensus_proteins = consensus_proteins, consensus_peptides = consensus_peptides, consensus_transitions = consensus_transitions, scampi_method = scampi_method, scampi_iterations = scampi_iterations, scampi_outliers = scampi_outliers, scampi_outliers_iterations = scampi_outliers_iterations, scampi_outliers_threshold = scampi_outliers_threshold)), method = cval_method, mcx = cval_mcx)$cv$mfe
				j <- j + 1
			}
			peptide_methods_names<-c(peptide_methods_names,paste("top",peptide_topx[i],sep=""))
			i <- i + 1
		}
		peptide_methods <- peptide_methods[-which(peptide_methods=="top")]
	}

	if (length(peptide_methods) > 0) {
		i <- 1
		while (i <= length(peptide_methods)) {
			j <- 1
			while (j <= length(transition_topx)) {
				cvmfe.mx[j,i+peptide_topx_number] <- cval.AbsoluteQuantification(AbsoluteQuantification(ProteinInference(data, peptide_method = peptide_methods[i], peptide_topx = NA, peptide_strictness = NA, peptide_summary = peptide_summary, transition_topx = transition_topx[j], transition_strictness = transition_strictness, transition_summary = transition_summary, fasta = fasta, apex_model = apex_model, combine_precursors = combine_precursors, combine_peptide_sequences = combine_peptide_sequences, consensus_proteins = consensus_proteins, consensus_peptides = consensus_peptides, consensus_transitions = consensus_transitions, scampi_method = scampi_method, scampi_iterations = scampi_iterations, scampi_outliers = scampi_outliers, scampi_outliers_iterations = scampi_outliers_iterations, scampi_outliers_threshold = scampi_outliers_threshold)), method = cval_method, mcx = cval_mcx)$cv$mfe
				j <- j + 1
			}
			peptide_methods_names<-c(peptide_methods_names,peptide_methods[i])
			i <- i + 1
		}
	}

	colnames(cvmfe.mx) <- peptide_methods_names
		
	return(cvmfe.mx)
}

performanceplot.ALF <- function(x, ...) {
	print(levelplot(t(x),xlab="Peptides",ylab="Transitions",main=paste("Optimal model: ",rownames(x)[which(x == min(x), arr.ind = TRUE)[1]]," Transitions, ",colnames(x)[which(x == min(x), arr.ind = TRUE)[2]]," Peptides",sep="")))
}
