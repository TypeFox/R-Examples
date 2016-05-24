#R Script written by Bettina Budeus, rokkaku-fight@gmx.de
#from 2012-08-02 to 2012-?-?
#removes the significant co mutation pairs, which are within a branch of the tree for the given sequences

#search for "change" for points where changes are possible
#please check before executing

#libs
library(Biostrings)
library(phangorn)

sequences <- c()
co_matrix_pre <- c()
tree <- c()
same <- c()
co_matrix <- c()
max_numbers <- c()
max_rows <- c()
result_list <- c()
distances <- c()
list_of_subtrees <- c()
l_list_of_subtrees <- c()

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

#read aa sequences from fasta file - they have to be aligned in order to work!
fe_wo_get_input_file_known_sequences <- structure(function(
	### set the input file (sequences).
	path_to_file
	### a FASTA file with sequence data. For reference please look in example file.
	){
	### the sequences from the FASTA file.
	sequences <- readAAStringSet(path_to_file)
	c_sequences <- create_correct_FASTA_input(sequences)
	.GlobalEnv[["sequences"]] <- c_sequences
	#print (length(sequenceB))
},ex=function(){
	fe_wo_get_input_file_known_sequences("SeqFeatR/extdata/Example.fasta")
})

#read co mutation results from csv file
fe_wo_get_input_file_co_mut_results <- structure(function(
	### set the input file (results from co-mutation analysis).
	path_to_file
	### a csv file with the results of the co-mutation analysis. For reference please look in example file.
	){
	### the results from co-mutation.
	.GlobalEnv[["co_matrix_pre"]] <- read.csv2(path_to_file, na.strings = "", colClasses = "character")
},ex=function(){
	fe_wo_get_input_file_co_mut_results("SeqFeatR/extdata/co_mutation_results_wo_allels.csv")

})

#read tree file
fe_wo_get_input_file_tree <- structure(function(
	### set the input file (tree file of the given sequences).
	path_to_file
	### a nexus file with tree data. For reference please look in example file.
	){
	### the tree.
	.GlobalEnv[["tree"]] <- read.nexus(path_to_file)
},ex=function(){
	fe_wo_get_input_file_tree("SeqFeatR/extdata/Example_tree.nh")
})

#functions
test_if_same_branch <- function(name_1, name_2){
	disctances <- .GlobalEnv[["distances"]]
	pos_1 <- which(dimnames(distances)[[2]]==name_1)
	pos_2 <- which(dimnames(distances)[[1]]==name_2)
	distance <- distances[pos_2, pos_1]
	return (distance)
}

# removes the rows in which our pair only occurs once
remove_only_ones <- function(list){
	result_list <- c()
	for (row in 1:nrow(list)){
		if (list[row,7] != 1){
			result_list <- rbind(result_list, list[row,])
		}
	}
	return (result_list)
}

#test a list of tips with the same pair mutation for occurency in a sub tree
test_if_same_subtree <- function(list){
	list_of_subtrees <- .GlobalEnv[["list_of_subtrees"]]
	sequences <- .GlobalEnv[["sequences"]]
	l_list_of_subtrees <- .GlobalEnv[["l_list_of_subtrees"]] 
	same <- .GlobalEnv[["same"]]
	length_of_list <- length(list)
	which_tree_numbers <- c()
	#generates a list of subtrees with enough tips as in the list
	for (tree_number in 1:l_list_of_subtrees){
		number_of_t <- length(list_of_subtrees[[tree_number]]$tip.label)
		if (number_of_t >= length_of_list){
			which_tree_numbers <- c(which_tree_numbers, tree_number)
		}
	}
	#calculates the list of subtrees with all tips from list
	subtrees_which_are_possible <- c()
	for (subtree in 1:length(which_tree_numbers)){
		is_in_subtree_list <- c()
		list_of_labels <- list_of_subtrees[[which_tree_numbers[subtree]]]$tip.label
		for (actual_tip in 1:length_of_list){
			is_in_subtree_list <- c(is_in_subtree_list, which(list_of_labels == list[actual_tip]))
			if (length(is_in_subtree_list) == length_of_list){
				subtrees_which_are_possible <- c(subtrees_which_are_possible, which_tree_numbers[subtree])
			}
		}
	}
	is_sub <- FALSE
	lowest_number_of_tips <- length(sequences)
	if(!is.null(subtrees_which_are_possible)){
		for (pos_tree in 1:length(subtrees_which_are_possible)){
			number_of_tips <- list_of_subtrees[[subtrees_which_are_possible[pos_tree]]]$Ntip
			if (number_of_tips >= same*length_of_list && !is_sub){
				is_sub <- FALSE
				lowest_number_of_tips <- number_of_tips
			}
			else{
				is_sub <- TRUE
				if (number_of_tips < lowest_number_of_tips){
					lowest_number_of_tips <- number_of_tips
				}
			}
		}
	}
	return (c(lowest_number_of_tips, is_sub))
}


#Main______________________________________________________________________________

foundereffectfinder <- structure(function(# Founder Eliminator
	### takes the result of co-mutation analysis and corrects them with an added column with a note if the row is in the same branch.
	path_to_file_sequence_alignment = NULL,
	### a FASTA file with sequence data. For reference please look in example file.
	path_to_file_assocpair_csv_result = NULL,
	### a csv file with the results of the co-mutation analysis. For reference please look in example file.
	path_to_file_nexus_tree = NULL,
	### a nexus file with tree data. For reference please look in example file.
	save_name_csv,
	### the file name of the result file.
	threshold = 7
	### number of tips so that the whole bulk is considered as in one branch.
	##details<< Takes the positions of the co-mutation results and a tree of the sequences with which the co-mutation was created and counts how often the pair of AAs is in one branch. If there are more than 'value' times the number of found co-mutations in one generated subtree, they are considered to be because of the founder effect.
	##seealso<< \code{\link{test_for_comutation_without_allel}}
	##note<< Only use files generated without allel usage!
	##references<< Mayr, Ernst (1954). "Change of genetic environment and evolution". In Julian Huxley. 
	##Evolution as a Process. London: George Allen & Unwin. OCLC 974739
	){
	result <- founder_eliminator_main_inner(path_to_file_sequence_alignment, path_to_file_assocpair_csv_result, path_to_file_nexus_tree, save_name_csv, threshold)
	return (result)
	
},ex=function(){
	ex <- system.file("extdata", "Example.fasta", package="SeqFeatR")
	ep <- system.file("extdata", "co_mutation_results_wo_allels.csv", package="SeqFeatR")
	co <- system.file("extdata", "Example_tree.nh", package="SeqFeatR")
	founderelim(ex,
 ep,
 co,
 "d.csv",
 7)
})

founder_eliminator_main_inner <- function(path_to_file_s = NULL, path_to_file_p = NULL, path_to_file_m = NULL, save_name, value = 7){
	fe_wo_get_input_file_known_sequences(path_to_file_s)
	fe_wo_get_input_file_co_mut_results(path_to_file_p)
	fe_wo_get_input_file_tree(path_to_file_m)

	sequences <- .GlobalEnv[["sequences"]]
	co_matrix_pre <- .GlobalEnv[["co_matrix_pre"]]
	tree <- .GlobalEnv[["tree"]]

	.GlobalEnv[["same"]] <- value
	co_matrix <- remove_only_ones(co_matrix_pre)
	max_numbers <- max(as.numeric(co_matrix[,7]))
	max_rows  <- nrow(co_matrix)
	result_list <- array(rep(0, max_numbers, max_rows), dim=(c(max_rows,max_numbers)))
	.GlobalEnv[["distances"]] <- cophenetic.phylo(tree)
	list_of_subtrees <- subtrees(tree)
	.GlobalEnv[["list_of_subtrees"]] <- list_of_subtrees
	.GlobalEnv[["l_list_of_subtrees"]] <- length(list_of_subtrees)

	for (row in 1:nrow(co_matrix)){
		position_1 <- co_matrix[row,2]
		position_1_letter <- co_matrix[row,4]
		position_2 <- co_matrix[row,3]
		position_2_letter <- co_matrix[row,5]
		i <- 1
		for (seq in sequences){
			letter_one_in_sequence <- substr(seq$seq, position_1, position_1)
			letter_two_in_sequence <- substr(seq$seq, position_2, position_2)
			if (letter_one_in_sequence == position_1_letter && letter_two_in_sequence == position_2_letter){
				result_list[row,i] <- seq$desc
				i <- i + 1
			}
		}
	}


	result_col <- c()
	nr_of_tips <- c()

	for (row in 1:nrow(result_list)){
		liney <- result_list[row,which(result_list[row,]!="0")]
		little_result_list <- c()
		if (length(liney) > 0){
			is_in_same_subtree <- test_if_same_subtree(liney)
			nr_of_tips <- c(nr_of_tips, paste(is_in_same_subtree[1], "/", length(liney)))
			if (is_in_same_subtree[2]){
				result_col <- c(result_col, "in subtree")
			}else{
				result_col <- c(result_col, "not in subtree")
			}
		}else{
			result_col <- c(result_col, "error")
		}
	}

	co_matrix_new <- cbind(co_matrix, nr_of_tips, result_col)

	write.csv2(co_matrix_new, save_name)

	return (co_matrix_new)
	### the same file as the co muation input file with an added column in which the user can read if the row is found possibly due to the founder effect.

}

#ex <- "../inst/extdata/Example_aa.fasta"
#ep <- "../inst/extdata/co_mutation_results_wo_allels.csv"
#co <- "../inst/extdata/Example_tree.nh"
#founderelim(ex,
# ep,
# co,
# "d.csv",
# 7)
