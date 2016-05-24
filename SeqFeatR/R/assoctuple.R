require(Biostrings)
require(parallel)

sm_wo_set_input_file_ori_sequences <- structure(function(
	### set the input file (known epitopes).
	path_to_file
	### a csv file with known epitopes. For reference please look in example file.
	){
	### known epitopes
	if (path_to_file != ""){
		.GlobalEnv[["sequences"]] <- readAAStringSet(path_to_file)
	}
	else {
		.GlobalEnv[["sequences"]] <- c()
	}
},ex=function(){
	
})

sm_wo_set_input_file_epi_results <- structure(function(
	###read aa sequences from fasta file - they have to be aligned in order to work!
	path_to_file
	### a FASTA file with sequence data. For reference please look in example file.
	){
	### the sequences from the FASTA file.
	data <- read.csv2(path_to_file, stringsAsFactors = FALSE, sep=";")
	.GlobalEnv[["result"]] <- data
},ex=function(){
		
})

#removes "0"s from a string
remove.zeroes <- function(s){
	if(substr(s, 1, 1) == "0")
  		a <- remove.zeroes(substr(s, 2, nchar(s)))  
	else
  		a <- s
}

get_seqs_for_certain_position <- function(position, StringSequences, AA){

	singleAA <- as.character(subseq(StringSequences, position, position))

	with <- which(singleAA==AA)
	seq_with_AA <- StringSequences[with]

	return (seq_with_AA)
	
}

get_identifier <- function(xstringset, idet){
	names_of <- names(xstringset)
	trop_1 <- length(grep(idet, names_of, fixed = TRUE))##all seqs which have certain identifier
	return (trop_1)
}

get_intres <- function(xstringset, intres, A11, A12, A21, A22, B11, B12, B21, B22){
	names_of <- names(xstringset)
	number_of_intres <- 0
	intres <- "A2"
	for (i in 1:length(xstringset)){
		allel.A1 <- paste("A", remove.zeroes(substr(names_of[i], A11, A12)), sep="")
    	  	allel.A2 <- paste("A", remove.zeroes(substr(names_of[i], A21, A22)), sep="")
   	   	allel.B1 <- paste("B", remove.zeroes(substr(names_of[i], B11, B12)), sep="")
  	    	allel.B2 <- paste("B", remove.zeroes(substr(names_of[i], B21, B22)), sep="")
		if(allel.A1 == intres || allel.A2 == intres || allel.B1 == intres || allel.B2 == intres){
			number_of_intres <- number_of_intres+1
		}
	}	
	return (number_of_intres)
}

get_p_fish <- function(a,b,c,d){
	small.table = matrix(c(a,b,c,d),nrow=2)
	#only calculate p-values from contingency tables with row sums above threshold or zero
	if(sum(small.table[1,])*sum(small.table[2,])!=0 & (sum(small.table[1,])<=0 | sum(small.table[2,])<=0)){
		return ("F")
	}
        #calculate fisher's exact test for each pair (allel, acid)    
        test.result <- fisher.test(small.table)
	return (test.result$p.value)
}

assoctuple <- structure(function(
	### Calculates the...
	### Details: Please be aware that it just useses the position in your csv file. If this is NOT the correct position, bechause you removed some positions, then please correct them first in the csv file or keep in mind, that the result file will also be just the relative position from your fasta file. 
	### Please use the same fasta file you had for epitope analysis!
	### The values for the analysis to estimate which positions should be included are normaly the corrected p-values but can be anything else. Even an own added column.
	path_to_file_sequence_alignment,
	###
	path_to_file_assocpoint_csv_result,
	###
	threshold,
	###
	min_number_of_elements_in_tuple,
	###
	max_number_of_elements_in_tuple,
	###
	save_name_csv,
	###
	column_of_feature,
	### the column in which the type is located for which the analysis should be done
	column_of_position,
	#### NEEEEWWW
	column_of_p_values,
	### the column from which the (p) values should be taken. See details.
	column_of_aa,
	A11,
	### the position of the start of the first HLA A Allel in the description block of the FASTA file.
	A12, 
	### the position of the end of the first HLA A Allel in the description block of the FASTA file.
	A21, 
	### the position of the start of the second HLA A Allel in the description block of the FASTA file.
	A22, 
	### the position of the end of the second HLA A Allel in the description block of the FASTA file.
	B11, 
	### the position of the start of the first HLA B Allel in the description block of the FASTA file.
	B12, 
	### the position of the end of the first HLA B Allel in the description block of the FASTA file.
	B21, 
	### the position of the start of the second HLA B Allel in the description block of the FASTA file.
	B22,
	### the position of the end of the second HLA B Allel in the description block of the FASTA file.
	one_feature,
	### if there is only one identifier
	feature
	### the identifier which should be analzed. See details. 	
	### idet is the one identifier inside the comment line of the fasta file which should be compared with all other identifiers. E.g X4 for HIV tropism.
	){

	result <- get_shared_mutations_inside(path_to_file_sequence_alignment, path_to_file_assocpoint_csv_result, threshold, min_number_of_elements_in_tuple, max_number_of_elements_in_tuple, save_name_csv, column_of_feature, column_of_position, column_of_p_values, column_of_aa, A11, A12, A21, A22, B11, B12, B21, B22, one_feature, feature)
	
	return (result)

},ex=function(){
	ex <- system.file("extdata", "Example.fasta", package="SeqFeatR")
	ep <- system.file("extdata", "epitope_results.csv", package="SeqFeatR")
	assoctuple(ex,
	ep,
	threshold,
	min_number_of_ele_in_tupel,
	max_number_of_ele_in_tuple,
	save_name_csv,
	column,
	column_of_position,
	column_of_values,
	column_of_aas,
	A11,
	A12, 
	A21, 
	A22, 
	B11, 
	B12, 
	B21, 
	B22,
	one_identifier,
	identifier
	)
})

get_shared_mutations_inside <- function(original_fasta_seq, epi_results, threshold, min_number_of_ele_in_tupel, max_number_of_ele_in_tuple, result_filename, column, column_of_position, column_of_values, column_of_aas, A11, A12, A21, A22, B11, B12, B21, B22, one_ident, idet){

	if (is.null(original_fasta_seq)==FALSE){
		sm_wo_set_input_file_ori_sequences(original_fasta_seq)
	}
	if (is.null(epi_results)==FALSE){
		sm_wo_set_input_file_epi_results(epi_results)
	}

	StringSequences <- .GlobalEnv[["sequences"]]
	data <- .GlobalEnv[["result"]]
	#print (data)
	#print (threshold)
	#print (column_of_values)
	
	#1. Get all input from csv

	#data2 <- read.csv2("Tropismus_results_test_clus.csv", stringsAsFactors = FALSE)
	#2. Get all AAs below a certain limit
	#print (data[,column_of_values])
	below_threshold_f <- data[which(as.numeric(data[,column_of_values]) < threshold),]
	#print (below_threshold_f)
	#mir <- occ
	AA <- data[which(as.numeric(data[,column_of_values]) < threshold), column_of_aas]
	below_threshold_number_f <- data[which(as.numeric(data[,column_of_values]) < threshold), column_of_position]
	below_threshold_number <- c()

	if(length(below_threshold_number_f) < 2){
		stop("There is no or only one position with a value below your threshold", "\n")
	}

	#3. Get List of doubles, triples, quadrupels, etc of positions
	all_result <- list()
	for (i in 1:length(below_threshold_number_f)){
		below_threshold_number <- c(below_threshold_number, below_threshold_number_f[i])
	}
	cat ("The following tuples will be examinated: ", below_threshold_number, "\n")

	#for (m in 2:length(below_threshold_number)){
	for (m in min_number_of_ele_in_tupel:max_number_of_ele_in_tuple){
		result <- combn(below_threshold_number, m, FUN = NULL, simplify = TRUE)
		all_result[[length(all_result) + 1L]] <- result
	}

	#3.5 make set of data for the x positions and mutations,
	below_threshold <- data[below_threshold_number,]
	all_aa_seqs <- list()
	for (i in 1:length(below_threshold_number)){
		### first: the position of the AA, seqs, the AA with the lowest p-value
		result <- get_seqs_for_certain_position(below_threshold_number[i], StringSequences, AA[i])
		all_aa_seqs[[length(all_aa_seqs) + 1L]] <- result
	}

	#3.6 get Number of interessting one and others
	if (one_ident){
		identifier <- get_identifier(StringSequences, idet)
		rest <- length(StringSequences)-identifier
	}else{
		intres <- colnames(data)[column]
		cat("You choose to analyse ", intres, "\n")
		n_o_intres <- get_intres(StringSequences, intres, A11, A12, A21, A22, B11, B12, B21, B22)
		other <- length(StringSequences)-n_o_intres
	}

	list_of_26 <- c()

	#4. For each of entrys in above List:  create set for the x single entrys and test if some are the same
	for (i in 1:length(all_result)){
		cat (i, "/", length(all_result), "\n")
		list_of_tuple <- all_result[[i]]
		for (j in 1:ncol(list_of_tuple)){
			cat (j, "/", ncol(list_of_tuple), "\n")
			number_in_list <- which(list_of_tuple[1,j]==below_threshold_number)
			final_set <- all_aa_seqs[[number_in_list]]
			spec_re <- 0
			cat ("first set ", length(final_set), "\n")
			for (k in 2:nrow(list_of_tuple)){
				#print (final_set)
				next_number <- which(list_of_tuple[k,j]==below_threshold_number)
				next_set <- all_aa_seqs[[next_number]]
				cat ("next", length(next_set), "\n")
				both_together <- append(final_set, next_set)
				final_set <- intersect(final_set, next_set)
				if (length(final_set) > 0){
					names_final_set <- c()
					for (seq in 1:length(final_set)){
						names_final_set <- c(names_final_set, (names(both_together)[which(final_set[seq] == both_together)][1]))
					}
					names(final_set) <- names_final_set
					cat ("end final", length(final_set), "\n")
				}
				else{
					break
				}
			}
			if (length(final_set) > 0){
				result <- length(final_set)
				#5. Make Fish for tupel and intres
				if (one_ident){
					identifier_res <- get_identifier(final_set, idet)
					rest_res <- length(final_set)-identifier_res
					a <- identifier_res ##identifier in result
					b <- rest_res ##!identifier in result
					c <- identifier-identifier_res##X4 in rest
					d <- rest-rest_res##!identifier in rest
				}else{
					intresres <- get_intres(final_set, n_o_intres, A11, A12, A21, A22, B11, B12, B21, B22)
					other_res <- length(final_set)-intresres
					a <- intresres ##intres in result
					b <- other_res ##!intres in result
					c <- n_o_intres-intresres##X4 in rest
					d <- other-other_res##!intres in rest
				}
				#cat (a,b,c,d,"\n")
				p_value <- get_p_fish (a,b,c,d)
			}
			else {
				result <- 0
				p_value <- "F"
				a <- "F"
				b <- "F"
				c <- "F"
				d <- "F"
			}
			pos <- strsplit(toString(list_of_tuple[,j]), ", ")
			positions_list <- c()
			for (posi in 1:length(pos[[1]])){
				positions_list <- paste(positions_list, pos[[1]][posi], "\t", sep="")
			}
			#6. write answer in file: position combo, AA combo, in x seqs together
			cat(j, "\tPosition combination: \t", positions_list, "were in\t", result, "\tsequences  together\t", "\twith p-value of\t", p_value, "\t",a, "\t",b, "\t",c, "\t",d,"\n",file = result_filename, append = TRUE)
			gc()
		}
	}
return (result_filename)


}
