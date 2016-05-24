vlfFun <-
function(x, p=0.001, seqlength=648, own = NULL){
	species.names <- x[,2]

	specimen.Number <- nrow(x)

	rownames(x) <- species.names

	Nuc.count <- count.function(x, specimen.Number,seqlength)

	frequency.matrix <- ffrequency.matrix.function(Nuc.count,seqlength)

	spec.freq <- specimen.frequencies(frequency.matrix, x, specimen.Number, species.names,seqlength)

	nucleotide.modalSequence <- MODE(frequency.matrix,seqlength)

	first.modal.frequencies <- MODE.freq(frequency.matrix,seqlength)

	second.modal.frequencies <- MODE.second.freq(frequency.matrix,seqlength)

	First_conserved_100 <- conservation_first(first.modal.frequencies, 1,seqlength)

	First_conserved_99.9 <- conservation_first(first.modal.frequencies, (1-p),seqlength)

	FirstAndSecond_conserved_99.9 <- conservation_two(first.modal.frequencies, second.modal.frequencies, (1-p),seqlength)

	specimen_VLFcount <- VLF.count.spec(spec.freq, p,seqlength)

	position_VLFcount <- VLF.count.pos(spec.freq, p,seqlength)

	VLFconvert <- VLF.convert.matrix(x, spec.freq, p,seqlength)

	VLFnuc <- VLF.nucleotides(VLFconvert, x,seqlength)

	VLFreduced <- VLF.reduced(VLFnuc, specimen_VLFcount, seqlength)

	species <- separate(VLFreduced)

	singleAndShared <- find.singles(species,seqlength)
	
	if(is.null(own)){
		foo<-list(modal=nucleotide.modalSequence, con100=First_conserved_100, conp=First_conserved_99.9,combine=FirstAndSecond_conserved_99.9,specimen=specimen_VLFcount, position=position_VLFcount, sas=singleAndShared, VLFmatrix = VLFreduced)
		class(foo)<-"vlf"
		foo
	}
	else{
		ownspec.freq <- specimen.frequencies(frequency.matrix, own, nrow(own), own[,2], seqlength)
		ownspec.VLFcount <- VLF.count.spec(ownspec.freq, p, seqlength)
		ownpos.VLFcount <- VLF.count.pos(ownspec.freq, p, seqlength)
		own.VLFconvert <- VLF.convert.matrix(own, ownspec.freq, p, seqlength)
		own.VLFnuc <- VLF.nucleotides(own.VLFconvert, own, seqlength)
		own.VLFreduced <- VLF.reduced(own.VLFnuc, ownspec.VLFcount, seqlength)
		
		foo<-list(modal=nucleotide.modalSequence, con100=First_conserved_100, conp=First_conserved_99.9,combine=FirstAndSecond_conserved_99.9,specimen=specimen_VLFcount, position=position_VLFcount, sas=singleAndShared, VLFmatrix = VLFreduced, ownSpecCount = ownspec.VLFcount, ownPosCount = ownpos.VLFcount, ownVLFMatrix = own.VLFnuc, ownVLFreduced = own.VLFreduced)
		class(foo)<-"vlf"
		foo
	}	
}
