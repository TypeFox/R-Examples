#R Script written by Bettina Budeus, rokkaku-fight@gmx.de
#from 2012-05-16 to 2012-?-?
#test for co-mutations in csv with p-values

#search for "change" for points where changes are possible
#please check before executing

#libs
library(Biostrings)

rowsum.threshold <- -1

min_FASTA_length <- 0

par.aacids <- 0

length_of_pairs <- 0

char_array <- 0

HLA_array <- 0

pair_array <- 0

length_allels <- 0

counter_a <- 0

allels.list <- 0

counter <- 0

result_matrix <- 0

result_matrix <- 0

results <- 0

A11 <- 0

A12 <- 0

A21 <- 0

A22 <- 0

B11 <- 0

B12 <- 0

B21 <- 0

B22 <- 0

thr.sig.fi <- c()

patnum.threshold <- c()

allels <- c()

threshold <- c()

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

cm_get_input_file <- structure(function(
	###read aa sequences from fasta file - they have to be aligned in order to work!
	path_to_file
	### a FASTA file with sequence data. For reference please look in example file.
	){
	sequences <- readAAStringSet(path_to_file)
	.GlobalEnv[["StringSetsequences"]] <- sequences
	c_sequences <- create_correct_FASTA_input(sequences)
	.GlobalEnv[["sequences"]] <- c_sequences
	.GlobalEnv[["sequences.count"]] <- length(c_sequences)
},ex=function(){
	cm_get_input_file("SeqFeatR/extdata/Example.fasta")
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

#calculates the length of result list
calculate_length <- function(list){
	counter <- 0
	for (i in 1:length(list)){
		counter <- counter + nrow(list[[i]])
	}
	return (counter)
}

#function do make char array from the sequences
make_char_array <- function(){
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

#function to create table with all possible pairs
make_pair_array <- function(){
	min_FASTA_length <- .GlobalEnv[["min_FASTA_length"]]
	length_of_pairs <- .GlobalEnv[["length_of_pairs"]]
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
	for (i in 3:3){
		for (j in 23:23){
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
	#print (small.table)
	if(sum(small.table[1,])*sum(small.table[2,])!=0 & (sum(small.table[1,])<=rowsum.threshold | sum(small.table[2,])<=rowsum.threshold))
	next
	#calculate fisher's exact test for each pair (allel, acid)    
	test.result <- fisher.test(small.table)
	return (test.result$p.value)
}


#main function to call:
# for each position pair:
# generate mini_table
# check if entrys in mini_table are inside threshold for significant F Test
main <- function (){
	pair_array <- .GlobalEnv[["pair_array"]]
	length_allels <- .GlobalEnv[["length_allels"]]
	result <- lapply(pair_array, FUN=inside_main)
	result_wo_zero <- result[!sapply(result, is.null)]
	length_of_all <- calculate_length(result_wo_zero)
	print (length_of_all)
	result_matrix <- matrix(rep(0,length_of_all*5),ncol = 5,dimnames = list(NULL, c("allel", "AAPair", "First Position", "Second Position", "p-value")))	
	counter <- 1
	for (i in 1:length(result_wo_zero)){
		for (j in 1:nrow(result_wo_zero[[i]])){
			result_matrix[counter,] <- result_wo_zero[[i]][j,]
			counter <- counter+1
		}
	}
	return (result_matrix[order(result_matrix[,1]),])
}

#inside main function with use of apply
inside_main <- function(col){
	threshold <- .GlobalEnv[["threshold"]]
	allels <- .GlobalEnv[["allels"]]
	allels.count <- .GlobalEnv[["allels.count"]]
	patnum.threshold <- .GlobalEnv[["patnum.threshold"]]
	sequencesSet <- .GlobalEnv[["StringSetsequences"]]
	sequences <- .GlobalEnv[["sequences"]]
	patnum.threshold <- .GlobalEnv[["patnum.threshold"]]
	length_allels <- .GlobalEnv[["length_allels"]]
	char_array <- .GlobalEnv[["char_array"]]
	allels <- .GlobalEnv[["allels"]]
	new_sequence_names <- .GlobalEnv[["new_sequence_names"]]
	thr.sig.fi <- .GlobalEnv[["thr.sig.fi"]]
	i = col[1]
	j = col[2]
	result <- c()
	cat (i, j, "\n")
	for(allel.num in 1:length(allels)){
  		if(allels.count[allel.num]>patnum.threshold){
			current_allel <- allels[allel.num]
			has_allel <- grep(current_allel, new_sequence_names)
			has_no_allel <- grep(current_allel, new_sequence_names, invert = T)

			just_this_position1 <- subseq(sequencesSet, i, i)
			mat1 <- as.matrix(just_this_position1)
			just_this_position2 <- subseq(sequencesSet, j, j)
			mat2 <- as.matrix(just_this_position2)
			mat <- mat1[,1] <- as.matrix(paste(mat1[,1], mat2[,1], sep=""))

			subset_has_allel <- mat[has_allel]
			subset_not_has_allel <- mat[has_no_allel]

			for (pair in 1:nrow(mat)){
				ap <- mat[pair, 1]
				ap_and_h <- length(which(subset_has_allel == ap))
				ap_and_n_h <- length(which(subset_not_has_allel == ap))
				n_ap_and_h <- length(which(subset_has_allel != ap))
				n_ap_and_n_h <- length(which(subset_not_has_allel != ap))
				small.table = matrix(c(ap_and_h, n_ap_and_h, ap_and_n_h, n_ap_and_n_h),nrow=2)
	  			#only calculate p-values from contingency tables with row sums above threshold or zero
	  			if(sum(small.table[1,])*sum(small.table[2,])!=0 & (sum(small.table[1,])<=rowsum.threshold | sum(small.table[2,])<=rowsum.threshold))
				next
        			#print(small.table)
				#calculate fisher's exact test for each pair (tropis, acid)    
        			test.result <- fisher.test(small.table)
				p.value <- test.result$p.value
				if (p.value < thr.sig.fi){
					result <- rbind(result,c(allels[allel.num],ap,i,j,p.value))
				}
			}
		}
	}
	return (result)
	### a list with co-mutation entrys
}

#adds a column with corrected p-values
add_correction_cm <- function(matrix){
	statistical_correction <- .GlobalEnv[["correction"]]
	p_values <- abs(as.numeric(matrix[,5]))
	all.pvals.corrected <- p.adjust(p_values, method = statistical_correction)
	matrix_c <- cbind(matrix, all.pvals.corrected)
	return (matrix_c)
}

#___________________________________________________________
assocpairfeat <- structure(function( #Test for co-mutation
	### Calculates if there are two positions in the given sequence which may have a co-mutation.
	##details<< For every position in the sequence given within the FASTA file a fishers exact test is made with every other position in the sequence and every HLA-type. The result p-values are collected in one big table if they are below a given value (thr.sig.f). Also uses p.adjust from stats package to calculate some p-value correction additionally as an extra column in the csv output.
	##seealso<< \code{\link{test_for_comutation_only_graphics}}
	##note<< If you want an graphical output, please use \code{\link{test_for_comutation_only_graphics}}.
	## If a patient has a homozygot HLA Allel, then please change the second one to "00" (without ") instead!
	path_to_file_sequence_alignment = NULL,
	### a FASTA file with sequence data. For reference please look in example file.
	save_name_csv, 
	### the file name of the result file in csv format
	dna = FALSE,
	### if the data is in DNA or amino Acid Code
	patnum_threshold = 1, 
	### the minimum number of patients of one HLA type to consider in the calculation.
	significance_level = 0.05, 
	### p-value threshold below which the results from fishers exact test should be added to output.
	A11a,
	### the position of the start of the first HLA A Allel in the description block of the FASTA file.
	A12a, 
	### the position of the end of the first HLA A Allel in the description block of the FASTA file.
	A21a, 
	### the position of the start of the second HLA A Allel in the description block of the FASTA file.
	A22a, 
	### the position of the end of the second HLA A Allel in the description block of the FASTA file.
	B11a, 
	### the position of the start of the first HLA B Allel in the description block of the FASTA file.
	B12a, 
	### the position of the end of the first HLA B Allel in the description block of the FASTA file.
	B21a, 
	### the position of the start of the second HLA B Allel in the description block of the FASTA file.
	B22a,
	### the position of the end of the second HLA B Allel in the description block of the FASTA file.
	multiple_testing_correction = "bonferroni"
	### the statistical correction applied to the p-values. Input can be: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
	){
	result <- test_for_comutation_inner(path_to_file_sequence_alignment, save_name_csv, dna, patnum_threshold, significance_level, A11a, A12a, A21a, A22a, B11a, B12a, B21a, B22a, multiple_testing_correction)
	return (result)
	
},ex=function(){
	ex <- system.file("extdata", "Example.fasta", package="SeqFeatR")
	assocpairfeat(ex,
 "co_mutation_results.csv",
 FALSE,
 1,
 0.05,
 22,
 23,
 25,
 26,
 29,
 30,
 32,
 33,
 "bonferroni")	
})	
#___________________________________________________________

test_for_comutation_inner <- function(path_to_file = NULL, save_name, dna = FALSE, patnum.threshol = 1, thr.sig.f = 0.05, A11a, A12a, A21a, A22a, B11a, B12a, B21a, B22a, statistical_correction = "bonferroni"){
	if (is.null(path_to_file)==FALSE){
		cm_get_input_file(path_to_file)
	}

	sequences <- .GlobalEnv[["sequences"]]
	sequences.count <- .GlobalEnv[["sequences.count"]]

	.GlobalEnv[["thr.sig.fi"]] <- thr.sig.f
	.GlobalEnv[["patnum.threshold"]] <- patnum.threshol
	.GlobalEnv[["A11"]] <- A11a
	.GlobalEnv[["A12"]] <- A12a
	.GlobalEnv[["A21"]] <- A21a
	.GlobalEnv[["A22"]] <- A22a
	.GlobalEnv[["B11"]] <- B11a
	.GlobalEnv[["B12"]] <- B12a
	.GlobalEnv[["B21"]] <- B21a
	.GlobalEnv[["B22"]] <- B22a
	.GlobalEnv[["correction"]] <- statistical_correction
	min_FASTA_length <- nchar(sequences[[1]]$seq)

	A11 <- .GlobalEnv[["A11"]]
	A12 <- .GlobalEnv[["A12"]]
	A21 <- .GlobalEnv[["A21"]]
	A22 <- .GlobalEnv[["A22"]]
	B11 <- .GlobalEnv[["B11"]]
	B12 <- .GlobalEnv[["B12"]]
	B21 <- .GlobalEnv[["B21"]]
	B22 <- .GlobalEnv[["B22"]]

	for (i in 1:length(sequences)){
      		if(nchar(sequences[[i]]$seq)<=min_FASTA_length) min_FASTA_length <- nchar(sequences[[i]]$seq)
	}
	.GlobalEnv[["min_FASTA_length"]] <- min_FASTA_length

	#find out all allels that occur in given sequences
	allels.all <<- c()
	allels.count <<- c()
	new_sequence_names <- c()

	for(i in 1:length(sequences)){
 	#the used values depend on fasta file; they possibly have to be changed - the same holds for the block some lines below
 	#allels.all holds all allels from the FASTA files, even duplicates
 	#allels holds no duplicates and none which are only a letter, without number (if the type is 00)
	#allels.count counts how often the certain alles is there
  	        f <- paste("A", remove.zeroes(substr(sequences[[i]]$desc, A11, A12)), sep="")
                s <- paste("A", remove.zeroes(substr(sequences[[i]]$desc, A21, A22)), sep="")
                t <- paste("B", remove.zeroes(substr(sequences[[i]]$desc, B11, B12)), sep="")
                o <- paste("B", remove.zeroes(substr(sequences[[i]]$desc, B21, B22)), sep="")
		allels.all <- c(allels.all, f, s, t, o)
		new_sequence_names <- c(new_sequence_names, paste(f, s, t, o))
	}
	allels <- sort(unique(allels.all[allels.all!="A" & allels.all!="B"]))
	.GlobalEnv[["allels"]] <- allels
	.GlobalEnv[["new_sequence_names"]] <- new_sequence_names

	for(i in 1:length(allels)){
  	allels.count <- c(allels.count, sum(allels.all==allels[i]))
	}

	#possible one letter codes for amino acids
	if (!dna){
		aacids <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-","B","X")
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
	.GlobalEnv[["char_array"]] <- make_char_array()
	.GlobalEnv[["pair_array"]] <- make_pair_array()
	.GlobalEnv[["length_allels"]] <- length(allels.count)
	.GlobalEnv[["allels.count"]] <- allels.count
	length_allels <- .GlobalEnv[["length_allels"]]
	counter_a <- array(rep(1, length(allels)), c(length_allels), list(allels))

	cat (sequences.count, length(par.aacids), sum(allels.count), "\n")
	
	allels.counter <- 0

	allels.list <- c()

	counter <- 1

	result_matrix <- main()

	result_matrix <- subset(result_matrix, !duplicated(result_matrix))

	results <- add_correction_cm(result_matrix)

	dimnames(results) <- list(NULL, c("allel", "AAPair", "First Position", "Second Position", "p_value", "corrected p_values"))

	write.csv2(results, paste(save_name, sep=""))

	return (results)
	### a table with every possible co-mutation below the given p-value.
}

#ex <- "../inst/extdata/Example_aa.fasta"
#	assocpairfeat(ex,
# "co_mutation_results.csv",
# FALSE,
# 1,
# 0.2,
# 10,
# 11,
# 13,
# 14,
# 17,
# 18,
# 20,
# 21,
# "bonferroni")
