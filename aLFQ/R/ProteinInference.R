# Protein inference for aLFQ import data frame
ProteinInference <- function(data, ...)  UseMethod("ProteinInference")

ProteinInference.default <- function(data, peptide_method = "top", peptide_topx = 2, peptide_strictness = "strict", peptide_summary = "mean", transition_topx = 3, transition_strictness = "strict", transition_summary = "sum", fasta = NA, apex_model = NA, combine_precursors = FALSE, combine_peptide_sequences = FALSE, consensus_proteins = TRUE, consensus_peptides = TRUE, consensus_transitions = TRUE, scampi_method = "LSE", scampi_iterations = 10, scampi_outliers = FALSE, scampi_outliers_iterations = 2, scampi_outliers_threshold = 2, ...) {
	peptide_sequence <- response <- concentration <- peptide_intensity <- NULL

	data <- data.table(data)

	# if the data is on the transition level
	if ("transition_intensity" %in% names(data)) {
		peptide <- PeptideInference(data, transition_topx = transition_topx, transition_strictness = transition_strictness, transition_summary = transition_summary, consensus_proteins = consensus_proteins, consensus_transitions = consensus_transitions)
	}
	else {
		peptide <- data
	}

	# if the aLFQ import data frame contains any anchor peptides (not proteins!), the concentrations of the according endogenous peptides (and proteins) is inferred
	if (dim(unique(peptide[,c("run_id","protein_id","concentration"), with = FALSE]))[1]!=dim(unique(peptide[,c("run_id","protein_id"), with = FALSE]))[1]) {
		# first step: all peptides with the same sequence as the anchor peptides are selected
	
		setkeyv(peptide,c("run_id","protein_id","peptide_id"))

		conc_peptides<-subset(peptide,peptide_sequence %in% subset(peptide,concentration != "?")$peptide_sequence & peptide_intensity > 0)
		setkeyv(conc_peptides,c("run_id","protein_id","peptide_sequence"))

		# second step: the intensity of the endogenous peptides is divided by the anchor peptide intensity and multiplied by the concentration of the anchor peptide
		conc_peptides<-conc_peptides[, list("concentration"=(peptide_intensity[which(concentration=="?")]/peptide_intensity[which(concentration!="?")])*as.numeric(concentration[which(concentration!="?")])), by=key(conc_peptides)]

		# third step: the endogenous peptides with know assigned concentrations are averaged to compute a mean protein concentration
		conc_proteins<-conc_peptides
		setkeyv(conc_peptides,c("run_id","protein_id"))

		conc_proteins<-conc_proteins[, list("concentration"=mean(concentration)), by=key(conc_proteins)]

		peptide<-merge(peptide[concentration=="?"][,concentration:=NULL],conc_proteins, by=c("run_id","protein_id"), all.x=TRUE)
		peptide$concentration[ is.na(peptide$concentration) ] <- "?"
	}

	# if the data is on the peptide level
	if ("peptide_intensity" %in% names(peptide)) {
		protein <- protein_inference.ProteinInference(peptide, peptide_method = peptide_method, peptide_topx = peptide_topx, peptide_strictness = peptide_strictness, peptide_summary = peptide_summary, fasta = fasta, apex_model = apex_model, combine_precursors = combine_precursors, combine_peptide_sequences = combine_peptide_sequences, consensus_proteins = consensus_proteins, consensus_peptides = consensus_peptides, scampi_method = scampi_method)
	}
	else {
		protein <- peptide
	}

	# if the data is on the protein level
	if ("protein_intensity" %in% names(protein)) {
		result <- protein
		setnames(result,"protein_intensity","response")
	}
	else {
		result <- protein
	}

	# return(data.frame(subset(result, response > 0)))
	return(data.frame(result))
}

protein_inference.ProteinInference <- function(data.dt, peptide_method = "top", peptide_topx = 1, peptide_strictness = "loose", peptide_summary = "mean", fasta = NA, apex_model = NA, combine_precursors = FALSE, combine_peptide_sequences = FALSE, consensus_proteins = TRUE, consensus_peptides = TRUE, scampi_method = "LSE", scampi_iterations = 10, scampi_outliers = FALSE, scampi_outliers_iterations = 2, scampi_outliers_threshold = 2) {
	run_id <- protein_id <- peptide_id <- peptide_sequence <- precursor_charge <- omni <- omni_reference <- peptide_intensity <- mean_peptide_intensity <- min_mean_peptide_intensity <- protein_sequence <- protein_sequence_length <- concentration <- response <- NULL

	if (!is.na(fasta)) {
		peptide_sequences.fasta <- trypsin(fasta)
	}

	# consensus filter
	if (consensus_peptides){
		# only use peptide_ids that occur in all run_ids
		setkeyv(data.dt,c("protein_id","peptide_id","precursor_charge"))

		omni.dt <- data.dt[, list("omni" = length(run_id)), by = key(data.dt)]

		data.dt<-data.dt[omni.dt]

		if (consensus_proteins) {
			data.dt<-subset(data.dt, omni == length(unique(data.dt$run_id)))
		}
		else {
			omni_reference.dt <- data.dt[, list("omni_reference" = length(unique(run_id))), by = protein_id]
			data.dt<-data.dt[omni_reference.dt][omni == omni_reference,]
		}
	}

	# should precursors be summed?
	if (combine_precursors) {
		setkeyv(data.dt,c("run_id","protein_id","peptide_id","peptide_sequence","concentration"))
		data.dt<-data.dt[, list("precursor_charge"=0, "peptide_intensity"=sum(peptide_intensity)), by=key(data.dt)]
	}

	# should peptide sequences be summed?
	if (combine_peptide_sequences) {
		setkeyv(data.dt,c("run_id","protein_id","peptide_sequence","concentration"))
		data.dt<-data.dt[, list("peptide_id"=peptide_sequence[1], "precursor_charge"=0, "peptide_intensity"=sum(peptide_intensity)), by=key(data.dt)]
	}
			
	if (peptide_method == "top") {
		# consensus filter
		if (consensus_peptides){
			setkeyv(data.dt,c("protein_id","peptide_id","precursor_charge"))

			# calculate mean of per peptide_id
			mean_peptides.dt <- data.dt[, list("mean_peptide_intensity"=mean(peptide_intensity)), by=key(data.dt)]

			data.dt<-data.dt[mean_peptides.dt]
			setkeyv(data.dt,c("run_id","protein_id","concentration"))

			# select top consensus peptides
			data.dt<-data.dt[data.dt[, list("min_mean_peptide_intensity"=min(strictnessfilter.ProteinInference(sort(mean_peptide_intensity,decreasing=TRUE)[1:peptide_topx],strictness=peptide_strictness)),na.rm=T), by=key(data.dt)]]

			data.dt<-subset(data.dt, mean_peptide_intensity >= min_mean_peptide_intensity)
		}

		setkeyv(data.dt,c("run_id","protein_id","concentration"))

		if (peptide_summary == "mean") {
			data.dt<-data.dt[, list("response"=mean(strictnessfilter.ProteinInference(sort(peptide_intensity,decreasing=TRUE)[1:peptide_topx],strictness=peptide_strictness))), by=key(data.dt)]

		} else if (peptide_summary == "median") {
			data.dt<-data.dt[, list("response"=median(strictnessfilter.ProteinInference(sort(peptide_intensity,decreasing=TRUE)[1:peptide_topx],strictness=peptide_strictness))), by=key(data.dt)]
		} else if (peptide_summary == "sum") {
			data.dt<-data.dt[, list("response"=sum(strictnessfilter.ProteinInference(sort(peptide_intensity,decreasing=TRUE)[1:peptide_topx],strictness=peptide_strictness))), by=key(data.dt)]
		}
		data.dt<-subset(data.dt,!is.na(response))
	}
	else if (peptide_method == "all") {
		setkeyv(data.dt,c("run_id","protein_id","concentration"))
	
		if (peptide_summary == "mean") {
			data.dt<-data.dt[, list("response"=mean(peptide_intensity)), by=key(data.dt)]
		}
		else if (peptide_summary == "median") {
			data.dt<-data.dt[, list("response"=median(peptide_intensity)), by=key(data.dt)]
		}
		else if (peptide_summary == "sum") {
			data.dt<-data.dt[, list("response"=sum(peptide_intensity)), by=key(data.dt)]
		}
	}
	else if (peptide_method == "iBAQ") {
    	if (is.na(fasta)) stop("This peptide_method requries a FASTA file.")
		peptides<-ldply(peptide_sequences.fasta, function(x) ldply(x))
		names(peptides)<-c("protein_id","peptide_sequence")

		data.dt<-subset(data.dt,peptide_sequence %in% peptides$peptide_sequence)

		setkeyv(data.dt,c("run_id","protein_id","concentration"))

		data.dt<-data.dt[, list("response"=sum(peptide_intensity)/length(unlist((lapply(as.list(peptide_sequences.fasta[[protein_id[1]]]),function(X){if(nchar(X)>=6 && nchar(X)<=30){return(X)}else{return(NA)}})))[!is.na(unlist((lapply(as.list(peptide_sequences.fasta[[protein_id[1]]]),function(X){if(nchar(X)>=6 && nchar(X)<=30){return(X)}else{return(NA)}}))))])), by=key(data.dt)]
	}
	else if (peptide_method == "APEX") {
    	if (is.na(fasta)) stop("This peptide_method requries a FASTA file.")
		peptides<-ldply(peptide_sequences.fasta, function(x) ldply(x))
		names(peptides)<-c("protein_id","peptide_sequence")

		data.dt<-subset(data.dt,peptide_sequence %in% peptides$peptide_sequence)

		peptide_sequences.af <- apexFeatures(data.frame("peptide_sequence" = unique(as.vector(unlist(peptide_sequences.fasta[data.dt$protein_id]))), "apex"=NA, stringsAsFactors=FALSE))

        peptide_sequences.apex <- predict(apex_model,peptide_sequences.af)$prediction[,c("peptide_sequence","apex")]

        setkeyv(data.dt,c("run_id","protein_id","concentration"))
        data.dt<-data.dt[, list("response"=sum(peptide_intensity)/sum(peptide_sequences.apex$apex[which(peptide_sequences.apex$peptide_sequence %in% unlist(unname(peptide_sequences.fasta[protein_id[1]])))])), by=key(data.dt)]
	}
	else if (peptide_method == "NSAF") {
    	if (is.na(fasta)) stop("This peptide_method requries a FASTA file.")
		peptides<-ldply(peptide_sequences.fasta, function(x) ldply(x))
		names(peptides)<-c("protein_id","peptide_sequence")

		data.dt<-subset(data.dt,peptide_sequence %in% peptides$peptide_sequence)

		setkeyv(data.dt,c("run_id","protein_id","concentration"))

		proteins <- ldply(read.fasta(file = fasta, seqtype = "AA", as.string = TRUE, seqonly = FALSE, strip.desc = TRUE), function(x) ldply(x))
		names(proteins) <- c("protein_id","protein_sequence")
		proteins<-data.table(proteins)
		setkey(proteins,protein_id)
		proteins<-proteins[,list("protein_sequence_length"=nchar(protein_sequence)), by=key(proteins)]

		data.dt<-merge(data.dt,proteins, by=c("protein_id"))
		setkeyv(data.dt,c("run_id","protein_id","concentration"))

		data.dt<-data.dt[, list("peptide_intensity"=sum(peptide_intensity)/protein_sequence_length), by=key(data.dt)]

		setkey(data.dt,"run_id")
		data.dt<-data.dt[, list("protein_id"=protein_id,"concentration"=concentration,"response"=peptide_intensity/sum(peptide_intensity)), by=key(data.dt)]

		setkeyv(data.dt,c("run_id","protein_id","concentration"))

		data.dt<-unique(data.dt)
	}
	else if (peptide_method == "SCAMPI") {
    	if (is.na(fasta)) stop("This peptide_method requries a FASTA file.")
		peptides<-ldply(peptide_sequences.fasta, function(x) ldply(x))
		names(peptides)<-c("protein_id","peptide_sequence")

		data.dt<-subset(data.dt,peptide_sequence %in% peptides$peptide_sequence)

     	if (!combine_precursors) stop("This peptide_method requries to combine precursors. Set combine_precursors=TRUE.")
    	if (!(scampi_method %in% c("LSE","MLE"))) stop("Select a valid SCAMPI method: LSE or MLE.")

    	data.dt<-data.table(ddply(as.data.frame(data.dt),.(run_id),function(X){scampi.ProteinInference(X, peptide_sequences.fasta, scampi_method, scampi_iterations, scampi_outliers, scampi_outliers_iterations, scampi_outliers_threshold)}))
		
		setkeyv(data.dt,c("run_id","protein_id","concentration"))
	}
	
	return(data.dt)
}

# Peptide inference for aLFQ import data frame
PeptideInference <- function(data, ...)  UseMethod("PeptideInference")

PeptideInference.default <- function(data, transition_topx = 3, transition_strictness = "strict", transition_summary = "sum", consensus_proteins = TRUE, consensus_transitions = TRUE, ...) {
	run_id <- protein_id <- peptide_id <- precursor_charge <- transition_id <- omni <- omni_reference <- transition_intensity <- mean_transition_intensity <- min_mean_transition_intensity <- peptide_intensity <- NULL

	if ("data.table" %in% class(data)) {
		data.dt<-data
	}
	else {
		data.dt<-data.table(data)
	}

	# consensus filter
	if (consensus_transitions){
		# only use transition_ids that occur in all run_ids
		setkeyv(data.dt,c("protein_id","peptide_id","precursor_charge","transition_id"))

		omni.dt <- data.dt[, list("omni" = length(run_id)), by = key(data.dt)]

		data.dt<-data.dt[omni.dt]

		if (consensus_proteins) {
			data.dt<-subset(data.dt, omni == length(unique(data.dt$run_id)))
		}
		else {
			omni_reference.dt <- data.dt[, list("omni_reference" = length(unique(run_id))), by = protein_id]
			data.dt<-data.dt[omni_reference.dt][omni == omni_reference,]
		}

		# calculate mean of per transition_id
		mean_transitions.dt <- data.dt[, list("mean_transition_intensity"=mean(transition_intensity)), by=key(data.dt)]

		data.dt<-data.dt[mean_transitions.dt]
		setkeyv(data.dt,c("run_id","protein_id","peptide_id","precursor_charge"))

		# select top consensus transitions
		data.dt<-data.dt[data.dt[, list("min_mean_transition_intensity"=min(strictnessfilter.ProteinInference(sort(mean_transition_intensity,decreasing=TRUE)[1:transition_topx],strictness=transition_strictness)),na.rm=T), by=key(data.dt)]]

		data.dt<-subset(data.dt, mean_transition_intensity >= min_mean_transition_intensity)
	}

	setkeyv(data.dt,c("run_id","protein_id","peptide_id","peptide_sequence","precursor_charge","concentration"))

	if (transition_summary == "mean") {
		data.dt<-data.dt[, list("peptide_intensity"=mean(strictnessfilter.ProteinInference(sort(transition_intensity,decreasing=TRUE)[1:transition_topx],strictness=transition_strictness))), by=key(data.dt)]

	} else if (transition_summary == "median") {
		data.dt<-data.dt[, list("peptide_intensity"=median(strictnessfilter.ProteinInference(sort(transition_intensity,decreasing=TRUE)[1:transition_topx],strictness=transition_strictness))), by=key(data.dt)]
	} else if (transition_summary == "sum") {
		data.dt<-data.dt[, list("peptide_intensity"=sum(strictnessfilter.ProteinInference(sort(transition_intensity,decreasing=TRUE)[1:transition_topx],strictness=transition_strictness))), by=key(data.dt)]
	}
	data.dt<-subset(data.dt,!is.na(peptide_intensity))

	if ("data.table" %in% class(data)) {
		return(data.dt)
	}
	else {
		return(data.frame(data.dt))
	}
}

scampi.ProteinInference <- function(data.df, peptide_sequences.fasta, scampi_method, scampi_iterations, scampi_outliers, scampi_outliers_iterations, scampi_outliers_threshold) {
	peptides<-ldply(peptide_sequences.fasta, function(x) ldply(x))
	names(peptides)<-c("protein_id","peptide_sequence")

	scampi_peptides<-unique(data.df[,c("peptide_id","peptide_sequence","peptide_intensity")])
	scampi_peptides$pepId<-1:nrow(scampi_peptides)
	names(scampi_peptides)<-c("peptide_id","pepSeq","pepQty","pepId")

	scampi_proteins<-unique(data.df[,c("protein_id","protein_id","concentration")])
	scampi_proteins$protId<-(nrow(scampi_peptides)+1):(nrow(scampi_peptides)+nrow(scampi_proteins))
	names(scampi_proteins)<-c("protein_id","protName","concentration","protId")

	scampi_edgespp<-merge(merge(unique(merge(data.df,peptides)[,c("peptide_id","protein_id")]),scampi_peptides),scampi_proteins)[,c("pepId","protId")]

	if (scampi_outliers) {
		scampires<-iterateScampi(peptides=scampi_peptides,proteins=scampi_proteins,edgespp=scampi_edgespp,rescaling=TRUE,method="all",numIter=scampi_outliers_iterations,numMLEIter=scampi_iterations,thresh=scampi_outliers_threshold,verbose=FALSE)
	}
	else {
		scampires<-runScampi(peptides=scampi_peptides,proteins=scampi_proteins,edgespp=scampi_edgespp,rescaling=TRUE,method="all",quantifyPeptides=FALSE,numIter=10,verbose=FALSE)
	}

	data.res<-scampires@proteins[,c("protein_id","concentration","LSEScore","MLEScore")]

	if (scampi_method=="LSE") {
		data.res$response<-10^data.res$LSEScore
	}
	else if (scampi_method=="MLE") {
		data.res$response<-10^data.res$MLEScore
	}

	return(data.res[,c("protein_id","concentration","response")])
}

strictnessfilter.ProteinInference <- function(data, strictness="loose") {
	if (NA %in% data && strictness == "loose") {
		return(as.vector(na.omit(data)))
	}
	else if (NA %in% data && strictness == "strict") {
		return(as.numeric(NA))
	}
	else {
		return(as.vector(data))
	}
}

trypsin <- function(fasta, ...)  UseMethod("trypsin")

trypsin.default <- function(fasta, ...) {
	proteins <- read.fasta(file = fasta, seqtype = "AA", as.string = TRUE, seqonly = FALSE, strip.desc = TRUE)

	sequences <- sapply(proteins, strsplit, "(?!P)(?<=[RK])", perl = TRUE)
	return(sequences)
}