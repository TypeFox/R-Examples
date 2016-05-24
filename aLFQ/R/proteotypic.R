proteotypic <- function(fasta, ...) UseMethod("proteotypic")

proteotypic.default <- function(fasta, apex_model, min_aa=4 , max_aa=20, ...) {
	peptide_sequence_length <- peptide_sequence <- NULL

	peptides.fasta <- trypsin(fasta)

	peptides <- ldply(peptides.fasta, function(x) ldply(x))
	names(peptides) <- c("protein_id","peptide_sequence")

	nonproteotypic <- peptides[duplicated(peptides$peptide_sequence),]$peptide_sequence
	peptides$peptide_sequence_length<-nchar(peptides$peptide_sequence)
	peptides.filt <- subset(peptides,peptide_sequence_length >= min_aa & peptide_sequence_length <= max_aa & !(peptide_sequence %in% nonproteotypic))

	peptides.natfilt <- peptides.filt[!grepl("B", peptides.filt$peptide_sequence),]
	peptides.natfilt <- peptides.natfilt[!grepl("J", peptides.natfilt$peptide_sequence),]
	peptides.natfilt <- peptides.natfilt[!grepl("O", peptides.natfilt$peptide_sequence),]
	peptides.natfilt <- peptides.natfilt[!grepl("U", peptides.natfilt$peptide_sequence),]
	peptides.natfilt <- peptides.natfilt[!grepl("X", peptides.natfilt$peptide_sequence),]
	peptides.natfilt <- peptides.natfilt[!grepl("Z", peptides.natfilt$peptide_sequence),]

	peptide_sequences.af <- apexFeatures(data.frame("peptide_sequence" = peptides.natfilt$peptide_sequence, "apex"=NA, stringsAsFactors=FALSE))
    peptide_sequences.apex <- predict(apex_model,peptide_sequences.af)

    peptide_apex <- peptide_sequences.apex$prediction[,c("peptide_sequence","apex")]

    proteotypic <- merge(peptide_apex, peptides, all.x=T, all.y =F)[,c("protein_id","peptide_sequence","apex")]

    return(proteotypic)
}
