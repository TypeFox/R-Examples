#R Script written by Bettina Budeus, rokkaku-fight@gmx.de
#from 2012-06-22 to 2012-?-?
#checks if the found comutation between two AA is in the possible epitopes found for this Allel

#search for "change" for points where changes are possible
#please check before executing

#libs
library(Biostrings)

#filenames, change as required...
### consensus sequence
#consensus <- 0
#sequences <- c()
#pos_epi <- c()
#mut_core <- c()

#function to create the same output as the old readFASTA from new read.AAStringSet
create_correct_FASTA_input <- function(sequences){
	liste <- as.character(sequences, use.names=TRUE)
	result_list <- c()
	for (i in 1:length(liste)){
		seq=liste[i][[1]]

		name = names(liste[i])
		entry <- list (name, seq)
		names(entry) <- c("desc", "seq")
		result_list[[i]] = entry
	}
	return (result_list)
}

pos_epi_get_input_file_sequences <- structure(function(
	### set the input file (sequence).
	path_to_file
	### file with sequence data in FASTA format.  For reference please look in example file.
	){
	### the sequences from the FASTA file.
	sequences <- readAAStringSet(path_to_file)
	#c_sequences <- create_correct_FASTA_input(sequences)
	.GlobalEnv[["sequences"]] <- sequences
},ex=function(){
	pos_epi_get_input_file_sequences("SeqFeatR/extdata/Example.fasta")
})


pos_epi_get_input_file_pos_epi <- structure(function(
	### set the input file (csv result from epitope analysis).
	path_to_file
	### file with possible epitope result data. For reference please look in example file.
	){
	### the possible epitope data.
	pos_epi_p <- read.csv2(path_to_file, na.strings=" ", colClasses = c("character"))
	.GlobalEnv[["pos_epi"]] <- remove_double_positions(pos_epi_p)
},ex=function(){
	pos_epi_get_input_file_pos_epi("SeqFeatR/extdata/epitope_results.csv")
})


pos_epi_get_input_file_mut_core <- structure(function(
	### set the input file (csv result from co mutation analysis).
	path_to_file
	### file with results from co mutation analysis. For reference please look in example file.
	, p.value
	### p-value to be defined as significant.
	){
	### the co mutation analysis data.
	mut_core_p <- read.csv2(path_to_file, na.strings=" ", colClasses = c("character"))
	.GlobalEnv[["mut_core"]] <- mut_core_p[which(abs(as.numeric(mut_core_p[,6])) < p.value),]
},ex=function(){
	pos_epi_get_input_file_mut_core("SeqFeatR/extdata/co_mutation_results.csv", 0.05)
})


#functions

is_in_pos_epi <- structure(function(row){
	pos_epi <- .GlobalEnv[["pos_epi"]]
	allel <- as.character(row[2][[1]]) ##<<
	first <- as.numeric(row[4][[1]]) ##<<
	second <- as.numeric(row[5][[1]]) ##<<
	allels_epi <- as.character(pos_epi[,6]) ##<<
	part_pos_epi_allel <- pos_epi[which(allel==allels_epi),] ##<<
	if (length(part_pos_epi_allel[[1]])==0){
		return (FALSE)
	}
	else if(is_in_epitope(first,part_pos_epi_allel)){
		return (TRUE)
	}
	else if (is_in_epitope(second,part_pos_epi_allel)){
		return (TRUE)
	}
	else {
		return (FALSE)
	}
},ex=function(){
	Require(tcltk)
})


is_in_epitope <- structure(function(position, part_pos_epi){
	result <- c()
	for (row in 1:nrow(part_pos_epi)){
		#print (part_pos_epi)
		range_of_epi <- c(part_pos_epi[row, 2]:part_pos_epi[row, 4])
		#print (range_of_epi)
		if (position %in% range_of_epi){
			return (TRUE)
		}else{
			result <- FALSE
		}
	}
	return (result)
},ex=function(){
	Require(tcltk)
})


remove_double_positions <- structure(function(list){
	list_2 <- subset(list,!duplicated(list$start))
	return (list_2)
},ex=function(){
	Require(tcltk)
})


#tests if given character is consensus on that position
is_consensus <- structure(function(letter, position){
	consensus <- .GlobalEnv[["consensus"]]
	if (letter == consensus[position]){
		return ("_")
	}
	else{
		return("Y")
	}
},ex=function(){
	Require(tcltk)
})

#-----------------------------------------------------------------------------
assocpointpair <- structure(function(# Compare results of co-mutation and epitopes
	### Checks which positions in the results of epitope finder are also in the results of given co-mutation. Therefore just compares certain columns of the results, so they have to be in the output format from epitope finder and co-mutation.
	##seealso<< \code{\link{find_possible_epitopes}}
	##seealso<< \code{\link{test_for_comutation}}
	##seealso<< \code{\link{test_for_comutation_without_allel}}
	##details<< This function takes the results from epitope finder and co-mutation analysis and tries to combine them into one result. It checks if a possible epitope, as it is an output of epitope finder, is also one of the positions which are found in the co mutation. If this is the case, they may be a point of compensatory mutation.
	path_to_file_sequence_alignment = NULL,
	### file with sequence data in FASTA format.  For reference please look in example file.
	path_to_file_assocpoint_csv_result = NULL,
	### file with possible epitope result data. For reference please look in example file.
	path_to_file_assocpairfeat_csv_result = NULL,
	### file with results from co mutation analysis. For reference please look in example file.
	significance_level = 0.05,
	### p-value to be defined as significant.
	save_name_csv,
	### the file name of the result file.
	save_name_pos
	### the file name of the possible compensatory file. Only written if there are such possible compensatory mutations.
	){
	result <- check_for_pos_epi_mut_core_inner(path_to_file_sequence_alignment, path_to_file_assocpoint_csv_result, path_to_file_assocpairfeat_csv_result, significance_level, save_name_csv, save_name_pos)
	return (result)
},ex=function(){
	ex <- system.file("extdata", "Example.fasta", package="SeqFeatR")
	ep <- system.file("extdata", "epitope_results.csv", package="SeqFeatR")
	co <- system.file("extdata", "co_mutation_results.csv", package="SeqFeatR")
	assocpointtuple(ex, ep, co, 0.05,
	 "co_mut_in_epitopes_result.csv",
	 "possible_compensatory_mutation.csv")
})

check_for_pos_epi_mut_core_inner <- function(path_to_file_s = NULL,path_to_file_e = NULL,path_to_file_m = NULL,p.value = 0.05,save_name_result,save_name_pos){
	
	if (is.null(path_to_file_s)==FALSE){
		pos_epi_get_input_file_sequences(path_to_file_s)
	}
	if (is.null(path_to_file_e)==FALSE){
		pos_epi_get_input_file_pos_epi(path_to_file_e)
	}
	if (is.null(path_to_file_m)==FALSE){
		pos_epi_get_input_file_mut_core(path_to_file_m, p.value)
	}

	sequences <- .GlobalEnv[["sequences"]]
	mut_core <- .GlobalEnv[["mut_core"]] 

	.GlobalEnv[["consensus"]] <- consensusString(sequences, ambiguityMap="?", threshold=0.6)

	lines <- nrow(mut_core)

	result <- array(rep(0, lines*9), dim=c(lines, 9))

	i <- 1

	for (entry in 1:nrow(mut_core)){
		cat (entry, "/", lines, "\n")
		if (is_in_pos_epi(mut_core[entry, ])){
			line <- list(mut_core[entry,])
			first_letter <- substr(mut_core[entry,3],1,1)
			first_position <- as.numeric(mut_core[entry,4])
			second_letter <- substr(mut_core[entry,3],2,2)
			second_position <- as.numeric(mut_core[entry,5])
			res <- c(as.matrix(mut_core[entry,]),TRUE,is_consensus(first_letter, first_position),is_consensus(second_letter, second_position))
			result[i,] <- res
			i <- i+1
		}
	}
	### the result of the calculation if a given position in epitope finder results are also in co-mutation results.
	result <- result[1:(i-1),2:ncol(result)]
	if (is.null(nrow(result))){
		result <- "Nothing found"
	}else {
		colnames(result) <- c("allel","AA_pair","f.pos", "s.pos","p-value", "is_true", "first_letter_consensus", "second_letter_consensus")

		possible_compensatory_mutation <- array(dim=c(1,8))

		i <- 0

		for (row in 1:nrow(result)){
			if(result[row,7]=="Y" && result[row,8]=="Y"){
				i <- i + 1
				possible_compensatory_mutation <- rbind(possible_compensatory_mutation, result[row,])
			}
		}
		write.csv2(possible_compensatory_mutation, paste(save_name_pos))
	}
	write.csv2(result, paste(save_name_result))
	return (result)
	### the corresponding table with an added column with q-values
}

#---------------------------------------------------------------------
#assocpointpair("../inst/extdata/Example_aa.fasta", "../inst/extdata/assocpoint_results.csv", "../inst/extdata/assocpairfeat_results.csv", significance_level=0.05, save_name_csv="assocpointpair_result.csv", save_name_pos="possible_compensatory_mutation.csv")
