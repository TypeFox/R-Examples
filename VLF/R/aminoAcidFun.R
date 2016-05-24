aminoAcidFun <-
function(x, p=0.001, seqlength = 216, own = NULL){
	species.names <- x[,2]
	
	specimen.Number <- nrow(x)
	
	rownames(x) <- species.names
	
	aminoAcid_count <- aa.count.function(x, seqlength)

	aminoAcid_frequency.Matrix <- aa.frequency.matrix.function(aminoAcid_count, seqlength)

	aminoAcid_specfrequencies <- aa.specimen.frequencies(aminoAcid_frequency.Matrix, x, species.names, seqlength)

	aminoAcid_Modal <- aa.MODE(aminoAcid_frequency.Matrix, seqlength)

	aminoAcid_firstModalFreq <- aa.MODE.freq(aminoAcid_frequency.Matrix, seqlength)

	aminoAcid_secondModalFreq <- aa.MODE.second.freq(aminoAcid_frequency.Matrix, seqlength)

	aminoAcid_firstConservation_100 <- aa.conservation_first(aminoAcid_firstModalFreq, 1, seqlength)

	aminoAcid_first_Conservation_99.9 <- aa.conservation_first(aminoAcid_firstModalFreq, 1-p, seqlength)

	aminoAcid_secondConservation_99.9 <- aa.conservation_two(aminoAcid_firstModalFreq, aminoAcid_secondModalFreq, 1-p, seqlength)

	aminoAcid_specimenVLFcount <- aa.VLF.count.spec(aminoAcid_specfrequencies, p, seqlength)

	aminoAcid_positionVLFcount <- aa.VLF.count.pos(aminoAcid_specfrequencies, p, seqlength)

	aaVLFconvert <- aa.VLF.convert.matrix(x, aminoAcid_specfrequencies, p, seqlength)

	aminoAcidVLFs <- VLF.aminoAcids(aaVLFconvert, x, seqlength)

	aaVLFreduced <- aa.VLF.reduced(aminoAcidVLFs, aminoAcid_specimenVLFcount, seqlength)

	aaSpecies <- separate(aaVLFreduced)

	aminoAcid_singleAndShared <- aa.find.singles(aaSpecies, seqlength)
	
	if(is.null(own)){
		foo <- list(modal = aminoAcid_Modal, con100 = aminoAcid_firstConservation_100, conp = aminoAcid_first_Conservation_99.9, combine = aminoAcid_secondConservation_99.9, specimen = aminoAcid_specimenVLFcount, position = aminoAcid_positionVLFcount, sas = aminoAcid_singleAndShared, VLFmatrix = aaVLFreduced) #ADD OWN VARIABLE INFO HERE
		class(foo) <- "aaVLF"
		foo
	}
	else{
		ownaa.spec.freq <- aa.specimen.frequencies(aminoAcid_frequency.Matrix, own, own[,2], seqlength)
		ownaa.spec.VLFcount <- aa.VLF.count.spec(ownaa.spec.freq, p, seqlength)
		ownpos.aaVLFcount <- aa.VLF.count.pos(ownaa.spec.freq, p, seqlength)
		own.aaVLFconvert <- aa.VLF.convert.matrix(own, ownaa.spec.freq, p, seqlength)
		own.aaVLFnuc <- VLF.aminoAcids(own.aaVLFconvert, own, seqlength)
		own.aareduced <- aa.VLF.reduced(own.aaVLFnuc,ownaa.spec.VLFcount, seqlength)
		
		foo <- list(modal = aminoAcid_Modal, con100 = aminoAcid_firstConservation_100, conp = aminoAcid_first_Conservation_99.9, combine = aminoAcid_secondConservation_99.9, specimen = aminoAcid_specimenVLFcount, position = aminoAcid_positionVLFcount, sas = aminoAcid_singleAndShared, VLFmatrix = aaVLFreduced, ownSpecCount = ownaa.spec.VLFcount, ownPosCount = ownpos.aaVLFcount, ownVLFMatrix = own.aaVLFnuc, ownVLFreduced = own.aareduced)
		class(foo) <- "aaVLF"
		foo
	}	
}
