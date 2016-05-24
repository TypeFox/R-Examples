#R Script written by Bettina Budeus, rokkaku-fight@gmx.de
#from 2012-05-16 to 2012-?-?
#test for pair mutations without allel

#search for "change" for points where changes are possible
#please check before executing

#libs
library(Biostrings)

rowsum.threshold <- -1

#some global variables
consensus <- 0

min_FASTA_length <- 0

par.aacids <- 0

length_of_pairs <- 0

char_array <- 0

HLA_array <- 0

pair_array <- 0

length_allels <- 0

counter <- 0

result_matrix <- 0

result_matrix <- 0

results <- 0

sequences <- c()

sequences.count <- c()

multicores <- c()

thr.sig.fi <- c()

#some functions needed in the following

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

cm_wo_get_input_file_sequences <- structure(function(
	###read aa sequences from fasta file - they have to be aligned in order to work!
	path_to_file
	### a FASTA file with sequence data. For reference please look in example file.
	){
	sequences <- readAAStringSet(path_to_file)
	c_sequences <- create_correct_FASTA_input(sequences)
	.GlobalEnv[["sequences"]] <- c_sequences
	.GlobalEnv[["sequences.count"]] <- length(c_sequences)
},ex=function(){
	cm_get_input_file("SeqFeatR/extdata/Example.fasta")
})

cm_wo_get_input_file_consensus <- structure(function(
	### set the input file (consensus sequences).
	path_to_file
	### a FASTA file with consensus data. For reference please look in example file.
	){
	consensus_sequences <- readAAStringSet(path_to_file)
	c_consensus_sequences <- create_correct_FASTA_input(consensus_sequences)
	.GlobalEnv[["consensus"]] <- c_consensus_sequences[[1]]$seq
},ex=function(){
	cm_get_input_file("SeqFeatR/extdata/Example_Consensus.fasta")
})

#removes "0"s from a string
remove.zeroes <- function(s){
	if(substr(s, 1, 1) == "0")
  		a <- remove.zeroes(substr(s, 2, nchar(s)))  
	else
  		a <- s
	}

#get indices of non-zero-values
find.neqzero <- function(array){
	r <- c()
	for(i in 1:length(array)){
  		if(array[i]!=0)
    			r <- c(r,i)
	}
	r
}

get_AA_from_consensus <- function(i){
	consensus <- .GlobalEnv[["consensus"]]
	return (substr(consensus, i, i))

}

#checks if both pairs are totally different
both_are_different <- function(pair_1, pair_2){
	letter_1a <- substr(pair_1, 1, 1)
	letter_1b <- substr(pair_1, 2, 2)
	letter_2a <- substr(pair_2, 1, 1)
	letter_2b <- substr(pair_2, 2, 2)
	if ((letter_1a != letter_2a) && (letter_1b != letter_2b)){
		return (TRUE)
	}
	else{
		return (FALSE)
	}
}

#calculates the length of result list
calculate_length <- function(list){
	counter <- 0
	print (list)
	for (i in 1:length(list)){
		counter <- counter + nrow(list[[i]])
	}
	return (counter)
}

#function do make char array from the sequences
cm_wo_make_char_array <- function(){
	sequences <- .GlobalEnv[["sequences"]]
	min_FASTA_length <- .GlobalEnv[["min_FASTA_length"]]
	char_array <- array(rep(0,min_FASTA_length*length(sequences)),c(length(sequences),min_FASTA_length))
	for (sequence in 1:length(sequences)){
		for (i in 1:min_FASTA_length){
			letter <- (substr(sequences[[sequence]]$seq, i, i))
			char_array[sequence,i] <- letter
		}
		
	}
	return (char_array)
}

#function to make list of HLA-types for each sequence
#so that all are the same
cm_wo_make_HLA_types_array <- function(){
	sequences <- .GlobalEnv[["sequences"]]
	HLA_array <- array(rep(0,length(sequences)*4),c(length(sequences),4))
	for (sequence in 1:length(sequences)){
		a1 <- "A02"
      		a2 <- "A02"
      		b1 <- "A02"
      		b2 <- "A02"
		HLA_array[sequence,1] <- a1
		HLA_array[sequence,2] <- a2
		HLA_array[sequence,3] <- b1
		HLA_array[sequence,4] <- b2
	}
	return (HLA_array)
}

#function to create table with all possible pairs
cm_wo_make_pair_array <- function(){
	length_of_pairs <- .GlobalEnv[["length_of_pairs"]]
	min_FASTA_length <- .GlobalEnv[["min_FASTA_length"]]
	pair_array <- array(rep(0,length_of_pairs*2),c(2,length_of_pairs))
	counter <- 1
	for (i in 1:(min_FASTA_length-1)){
		for (j in (i+1):min_FASTA_length){
			col <- c(i,j)
			pair_array[,counter] <- col
			counter <- counter+1
		}
	}
	pair_array <- split(pair_array, col(pair_array))
	return (pair_array)
}

make_special_pair_array <- function(){
	pair_array <- array(rep(0,1),c(2,1))
	counter <- 1
	for (i in 4:4){
		for (j in 66:66){
			col <- c(i,j)
			pair_array[,counter] <- col
			counter <- counter+1
		}
	}
	pair_array <- split(pair_array, col(pair_array))
	return (pair_array)
}

#for controll and test purpose only!!!!
make_fisher_test <- function(a,b,c,d){
	small.table = matrix(c(a,b-a,c-a,d-b-c+a),nrow=2)
	if(sum(small.table[1,])*sum(small.table[2,])!=0 & (sum(small.table[1,])<=rowsum.threshold | sum(small.table[2,])<=rowsum.threshold))
	next
	#calculate fisher's exact test for each pair (allel, acid)
	test.result <- fisher.test(small.table)
	return (test.result$p.value)
}

#create mini_table for calculating all relevant numbers for one position pair i,j
#@result: a table with all relevant data aka the sum for all occurings
cm_wo_create_mini_table <- function(i,j){
	last_col <- c("sum")
	last_row <- c("sum")
	sequences <- .GlobalEnv[["sequences"]]
	char_array <- .GlobalEnv[["char_array"]]
	length_allels <- .GlobalEnv[["length_allels"]]
	mini_table <- array(rep(0, (length_allels+2)),c(length_allels+2))
	mini_table <- rbind(mini_table, c(rep(0, (length_allels+2))))
	for (sequence in 1:length(sequences)){		
		l1 <- char_array[sequence,i]
		l2 <- char_array[sequence,j]
		pos.l1 = which(mini_table[,1] == l1)
		pos.l2 = which(mini_table[1,] == l2)
		if (!length(pos.l1)){
			mini_table <- rbind(mini_table, c(0))
			rows = nrow(mini_table)
			mini_table[rows,1] = l1
		}
		if (!length(pos.l2)){
			mini_table <- cbind(mini_table, c(0))
			cols = ncol(mini_table)
			mini_table[1, cols] = l2
		}
		pos.l1 = which(mini_table[,1] == l1)
		pos.l2 = which(mini_table[1,] == l2)
		mini_table[pos.l1,pos.l2] = as.numeric(mini_table[pos.l1,pos.l2]) + 1
	}
	for (col in 2:ncol(mini_table)){
		last_col <- c(last_col, sum(as.numeric(mini_table[(2:nrow(mini_table)),col])))
	}
	mini_table <- rbind(mini_table, last_col)
	for (row in 2:nrow(mini_table)){
		last_row <- c(last_row, sum(as.numeric(mini_table[row,(2:ncol(mini_table))])))
	}
	mini_table <- cbind(mini_table, last_row)
	mini_table_r <- mini_table[-2,-(2:3)]
	cons_i <- get_AA_from_consensus(i)
	cons_j <- get_AA_from_consensus(j)
	cons_row <- rep(cons_i,(ncol(mini_table_r)))
	cons_col <- rep(cons_j,(nrow(mini_table_r)+1))

	mini_table_r <- rbind(mini_table_r, cons_row)	
	mini_table_r <- cbind(mini_table_r, cons_col)
	mini_table_r[nrow(mini_table_r), ncol(mini_table_r)] <- paste(cons_i, cons_j, sep="")
	return (mini_table_r)
}

#main function to call:
# for each position pair:
# generate mini_table
# check if entrys in mini_table are inside threshold for significant F Test
cm_wo_main <- function (){
	pair_array <- .GlobalEnv[["pair_array"]]
	counter <- .GlobalEnv[["counter"]]
	result <- lapply(pair_array, FUN=cm_wo_inside_main)
	result_wo_zero <- result[!sapply(result, is.null)]
	length_of_all <- calculate_length(result_wo_zero)
	print (length_of_all)
	result_matrix <- matrix(rep(0,length_of_all*10),ncol = 10,dimnames = list(NULL, c("First Position", "Second Position", "First A", "Second A", "Consensus", "#of_first_and_second","#of_first_w_o_test_pair","#of_second_w_o_test_pair","#of_all","p-value")))	
	counter <- 1
	for (i in 1:length(result_wo_zero)){
		for (j in 1:nrow(result_wo_zero[[i]])){
			result_matrix[counter,] <- result_wo_zero[[i]][j,]
			counter <- counter+1
		 }
	}
	.GlobalEnv[["result_matrix"]] <- result_matrix
}

#inside main function with use of apply
cm_wo_inside_main <- function(col){
	thr.sig.fi <- .GlobalEnv[["thr.sig.fi"]]
	i = col[1]
	j = col[2]
	result <- c()
	cat (i, j, "\n")
	mini_table <- cm_wo_create_mini_table(i,j)
	for (col in 2:(ncol(mini_table)-1)){
		if (mini_table[1,col]!=mini_table[1,ncol(mini_table)]){
			for (row in 2:(nrow(mini_table)-1)){
				if (mini_table[row,1]!=mini_table[nrow(mini_table),1]){
					a <- as.numeric(mini_table[row,col])
					if (a != 0){
						b <- as.numeric(mini_table[row,(ncol(mini_table)-1)])
						c <- as.numeric(mini_table[(nrow(mini_table)-1),col])
						d <- as.numeric(mini_table[(nrow(mini_table)-1),(ncol(mini_table)-1)])
						p.value <- make_fisher_test(a,b,c,d)
						if (p.value <= thr.sig.fi){
							f.pos <- i
							s.pos <- j
							f.a <- mini_table[row,1]
							s.a <- mini_table[1,col]
							cons <- mini_table[nrow(mini_table),ncol(mini_table)]
							result <- rbind(result, c(f.pos, s.pos, f.a, s.a, cons, a,b-a,c-a,d, p.value))
						}
					}
				}
			}
		}
	}
	return (result)
}

#adds a column with corrected p-values
add_correction <- function(matrix){
	statistical_correction <- .GlobalEnv[["correction"]]
	p_values <- abs(as.numeric(matrix[,10]))
	all.pvals.corrected <- p.adjust(p_values, method = statistical_correction)
	matrix_c <- cbind(matrix, all.pvals.corrected)
	return (matrix_c)
}

#___________________________________________________________

assocpair <- structure(function(#Test for co-mutation without HLA types
	### Calculates if there are two positions in the given sequence which may have a co-mutation.
	##details<< For every position in the sequence given within the FASTA file a fishers exact test is made with every other position in the sequence and the Consensus sequence. The result p-values are collected in one big table if they are below a given value (thr.sig.f). If there are any HLA types in the FASTA description they are simply ignored. Also uses p.adjust from stats package to calculate some p-value correction additionally as an extra column in the csv output.
	##note<< If you want an graphical output, please use  \code{\link{test_for_comutation_only_graphics_wo_allels}}.
	##seealso<< \code{\link{test_for_comutation_only_graphics_wo_allels}}
	path_to_file_sequence_alignment = NULL,
	### a FASTA file with sequence data. For reference please look in example file.
	path_to_file_consensus = NULL,
	### a FASTA file with consensus data. For reference please look in example file.
	save_name_csv, 
	### the file name of the result file in csv format
	dna = FALSE,
	### if the data is in DNA or amino Acid Code
	significance_level = 0.05, 
	### p-value threshold below which the results from fishers exact test should be added to output.
	multiple_testing_correction = "bonferroni"
	### the statistical correction applied to the p-values. Input can be: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
	){
	result <- test_for_comutation_without_allel_inner(path_to_file_sequence_alignment, path_to_file_consensus, save_name_csv, dna, significance_level, multiple_testing_correction)
	return (result)
	
},ex=function(){
	ex <- system.file("extdata", "Example.fasta", package="SeqFeatR")
	ep <- system.file("extdata", "Example_Consensus.fasta", package="SeqFeatR")
	assocpair(ex,
 ep,
 "co_mutation_results_wo_allels.csv",
 FALSE,
 0.05,
 "bonferroni")	
})

test_for_comutation_without_allel_inner <- function(path_to_file_s = NULL, path_to_file_c = NULL, save_name, dna = FALSE, thr.sig.f = 0.05, statistical_correction = "bonferroni"){
	if (is.null(path_to_file_s)==FALSE){
		cm_wo_get_input_file_sequences(path_to_file_s)
	}
	if (is.null(path_to_file_c)==FALSE){
		cm_wo_get_input_file_consensus(path_to_file_c)
	}
	sequences <- .GlobalEnv[["sequences"]]

	.GlobalEnv[["thr.sig.fi"]] <- thr.sig.f

	#find the shortest length of FASTA String, if the sequences are of different length
	#min_FASTA_length holds this value
	min_FASTA_length <- nchar(sequences[[1]]$seq)
	for (i in 1:length(sequences)){
      		if(nchar(sequences[[i]]$seq)<=min_FASTA_length) min_FASTA_length <- nchar(sequences[[i]]$seq)
	}
	.GlobalEnv[["min_FASTA_length"]] <- min_FASTA_length

	#find out all allels that occur in given sequences
	allels.all <- c()
	allels.count <- c()

	for(i in 1:length(sequences))
  	#the used values depend on fasta file; they possibly have to be changed - the same holds for the block some lines below
  	#allels.all holds all allels from the FASTA files, even duplicates
  	#allels holds no duplicates and none which are only a letter, without number (if the type is 00)
  	#allels.count counts how often the certain alles is there
  	allels.all <- c(allels.all, "A02", "A02", "A02", "A02")

	allels <- sort(unique(allels.all[allels.all!="A" & allels.all!="B"]))

	for(i in 1:length(allels)){
  	allels.count <- c(allels.count, sum(allels.all==allels[i]))
	}

	#possible one letter codes for amino acids
	if (!dna){
		aacids <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-", "B", "X")
	}else {
		aacids <- c("A","C","G","T","-","R", "Y", "M", "K", "S", "W", "N")
	}

	#possible parings of amino acids
	par.aacids <- c()
	names.par.aacids <- c()
	for (i in 1:length(aacids)){
		for (j in 1:length(aacids)){
			both <- data.frame(c(aacids[i], aacids[j]))
			names.par.aacids <- c(names.par.aacids, paste(aacids[i],aacids[j], sep = ""))
			#print (both)
			par.aacids <- c(par.aacids, both)
		}
	}

	par.aacids <- data.frame(lapply(par.aacids, as.character), stringsAsFactors=FALSE)
	.GlobalEnv[["length_of_pairs"]] <- (min_FASTA_length*(min_FASTA_length-1))/2
	.GlobalEnv[["char_array"]] <- cm_wo_make_char_array()
	.GlobalEnv[["HLA_array"]] <- cm_wo_make_HLA_types_array()
	.GlobalEnv[["pair_array"]] <- cm_wo_make_pair_array()
	.GlobalEnv[["length_allels"]] <- length(allels.count)
	length_allels <- .GlobalEnv[["length_allels"]]
	counter_a <- array(rep(1, length(allels)), c(length_allels), list(allels))

	sequences.count <- .GlobalEnv[["sequences.count"]]
	cat (sequences.count, length(par.aacids), sum(allels.count), "\n")

	allels.counter <- 0

	allels.list <- c()
	.GlobalEnv[["counter"]] <- 1
	result_matrix_p <- cm_wo_main()
	result_matrix <- subset(result_matrix_p, !duplicated(result_matrix_p))

	results <- add_correction(result_matrix)

	dimnames(results) <- list(NULL, c("First Position", "Second Position", "First A", "Second A", "Consensus", "#of_first_and_second","#of_first_w_o_test_pair","#of_second_w_o_test_pair","#of_all","p_value", "corrected p_values"))

	write.csv2(results, paste(save_name, sep=""))

	return (results)
	### a table with every possible co-mutation below the given p-value.
}

#ex <- "../inst/extdata/Example_aa.fasta"
#ep <- "../inst/extdata/Example_Consensus_aa.fasta"
#	assocpair(ex,
# ep,
# "co_mutation_results_wo_allels.csv",
# FALSE,
# 0.05,
# "bonferroni")
