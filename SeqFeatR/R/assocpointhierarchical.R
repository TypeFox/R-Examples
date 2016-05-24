#http://www.r-bloggers.com/developing-a-user-friendly-regular-expression-function-easygregexpr/
easyGregexpr <- function(pattern, charvector, ...) {
  #require(plyr)

  if (storage.mode(charvector) != "character") stop("easyGregexpr expects charvector to be a character vector.")

  #identify all matches
  regexpMatches <- gregexpr(pattern, charvector, ...)
  
  convertMatches <- c()
  for (i in 1:length(regexpMatches)) {
    thisLine <- regexpMatches[[i]]
    #only append if there is at least one match on this line
    if (thisLine[1] != -1) {
      convertMatches <- rbind(convertMatches, data.frame(element=i, start=thisLine, end=thisLine + attr(thisLine, "match.length") - 1))
    }
  }
  
  #if no matches exist, return null (otherwise, will break adply)
  if (is.null(convertMatches)) return(NULL)
  
  #We now have a data frame with the line, starting position, and ending position of every match
  #Add the matched string to the data.frame
  #Use adply to iterate over rows and apply substr func
  convertMatches <- adply(convertMatches, 1, function(row) {
        row$match <- substr(charvector[row$element], row$start, row$end)
        return(as.data.frame(row))
      })

  #need to convert from factor to character because adply uses stringsAsFactors=TRUE even when overridden
  convertMatches$match <- as.character(convertMatches$match)
  return(convertMatches)
}

#search for "change" for points where changes are possible
#please check before executing

#libs
library(Biostrings)
library(plyr)
library(R2jags)
library(coda)
set.seed(12345)

rowsum.threshold <- -1
#min_FASTA_length <- 0
#epitopes <- c()
#sequences <- c()
#pos.epitopes <- c()
#pos.epi.matrix <- c()

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

#read epitopes from csv file
ep_wo_set_input_file_known_epitopes <- structure(function(
	### set the input file (known epitopes).
	path_to_file
	### a csv file with known epitopes. For reference please look in example file.
	){
	### known epitopes
	if (path_to_file != ""){
		.GlobalEnv[["epitopes"]] <- read.csv2(path_to_file)
	}
	else {
		.GlobalEnv[["epitopes"]] <- c()
	}
},ex=function(){
	ep_wo_set_input_file_known_epitopes("SeqFeatR/extdata/Example_epitopes.csv")
})

ep_wo_set_input_file_known_sequences <- structure(function(
	###read aa sequences from fasta file - they have to be aligned in order to work!
	path_to_file
	### a FASTA file with sequence data. For reference please look in example file.
	){
	### the sequences from the FASTA file.
	sequences <- readAAStringSet(path_to_file)
	consensus <- consensusString(sequences)
	.GlobalEnv[["StringSetsequences"]] <- sequences
	.GlobalEnv[["consensus"]] <- consensus
	c_sequences <- create_correct_FASTA_input(sequences)
	.GlobalEnv[["sequences"]] <- sequences
},ex=function(){
	ep_wo_set_input_file_known_sequences("SeqFeatR/extdata/Example.fasta")	
})

#read possible epitopes from csv file
ep_wo_set_input_file_known_pos_epi <- structure(function(
	### set the input file (HLA binding motifs).
	path_to_file
	### a csv file with HLA binding motifs if available. For reference please look in example file.
	){
	### the HLA binding motifs.
	if (path_to_file != ""){
		.GlobalEnv[["pos.epitopes"]] <- read.csv2(path_to_file,stringsAsFactors = FALSE)
	}
	else {
		.GlobalEnv[["pos.epitopes"]] <- c()
	}
},ex=function(){
	ep_wo_set_input_file_known_pos_epi("SeqFeatR/extdata/Example_HLA_binding_motifs.csv")
})

#some functions needed in the following

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

#return a vector from 0.1 till lowest 10 based potence of input
get_10_based <- function(value){
	vector <- c()
	potence <- 0
	while (round(value, digits = potence) == 0){
		potence <- potence + 1
	}
	lowest_value <- round(value, digits = potence)
	for (i in 1:(potence+1)){
		vector <- c(vector, 1*10^-(i))
	}
	#print (vector)
	#miracle <- occure
	return (vector)	
}

#check if given possible epitope is anywhere in range of possible epitope position
find_epitope_in_given_sequence <- function(pos.epitope, matrix.only.pvals.above.threshold, sequences, allel){
	#print (matrix.only.pvals.above.threshold)
	#print (pos.epitope)	
	min_FASTA_length <- .GlobalEnv[["min_FASTA_length"]]
	count.letter.inside <- -2
	count.praenthesis <- 0
	true.false.letter <- FALSE	
	for (letter in 1:nchar(pos.epitope)){
		if (substr(pos.epitope, letter, letter)=="["){
			count.praenthesis <- count.praenthesis + 1
			true.false.letter <- TRUE
		}
		else if (substr(pos.epitope, letter, letter)=="]"){
			true.false.letter <- FALSE	
		}
		if (true.false.letter && substr(pos.epitope, letter, letter)!="("){
			if (true.false.letter && substr(pos.epitope, letter, letter)!=")"){
				count.letter.inside <- count.letter.inside+1
			}
		}
	}
	sub <- gsub("[[|]|(|A-Z|)]", "", pos.epitope)
	length.of.epitope <- (nchar(gsub("[(|)]", "", sub)) - count.letter.inside + count.praenthesis)-1
	is_in_seq <- c()
	if(class(matrix.only.pvals.above.threshold)==class(2)){
		entry <- matrix.only.pvals.above.threshold[1]
		min <- max(1, entry-length.of.epitope)
		max <- min(entry+length.of.epitope, min_FASTA_length)
		for (single_seq in 1:length(sequences)){
			is_in_seq <- rbind(is_in_seq, test_for_occurence(pos.epitope, min, max, sequences[single_seq], allel, entry))
		}
	}else{
		for (entry in matrix.only.pvals.above.threshold[,1]){
			min <- max(1, entry-length.of.epitope)
			max <- min(entry+length.of.epitope, min_FASTA_length)
			for (single_seq in 1:length(sequences)){
				is_in_seq <- rbind(is_in_seq, test_for_occurence(pos.epitope, min, max, sequences[single_seq], allel, entry))
			}
		}
	}
	return (is_in_seq)
}

#tests if given sequence is in certain sequence
#min = min position, max = max position, original = where p-value is
test_for_occurence <- function(pos.epitope, min, max, seq, allel, original_position){
	eppi <- gsub("[(|)]", "", pos.epitope)
	result_array <- c()
	test.str <- as.character(seq)#as.character(subseq(seq, min, max))
	corrected.pos.epitope <- gsub("[x]","[A-Z]",eppi)
	if (!is.null(easyGregexpr(corrected.pos.epitope,test.str))){
		list.of.matches <- easyGregexpr(corrected.pos.epitope,test.str)
	#	print (list.of.matches)
		for (i in 1:nrow(list.of.matches)){
			start <- list.of.matches[i, 2]
			end <- list.of.matches[i, 3]
			result_array <- rbind(result_array, start, end)
		}
		#cat(min, max, allel, original_position,"\n")
		create_pos_epi_csv(result_array, seq, pos.epitope, allel, original_position, test.str)
		
	}
	return (result_array)
}

#creates plot data for possible epitopes from list with min, max values
create_pos_epi_plot <- function(pos.epitope_list, min_FASTA_length){
	plot <- c()
	starting <- FALSE
	for (i in 1:min_FASTA_length){
		if (any(pos.epitope_list == i) & !starting){
			starting <- TRUE
			plot <- c(plot, 1)
		}
		else if (any(pos.epitope_list == i) & starting){
			plot <- c(plot, 1)
			starting <- FALSE
		}
		else if (starting){
			plot <- c(plot, 1)
		}
		else{
			plot <- c(plot, 0)
		}
	}
	return (plot)
}


create_pos_epi_csv <- function(pos.epitope_list, seq, epitope, allel, original_position, test.str){
	pos.epi.matrix <- .GlobalEnv[["pos.epi.matrix"]]
	for (entry in 1:(length(pos.epitope_list)/2)){
		x <- 1
		if (!is.null(pos.epitope_list[x])){
			s <- pos.epitope_list[x]-1
			e <- pos.epitope_list[x+1]-1
			sequence <- test.str
			#cat (s, e, sequence, "\n")
			#cat (pos.epitope_list[x]-1, original_position, pos.epitope_list[x+1]-1, epitope, allel, sequence, "\n")
			pos.epi.matrix_p <- rbind(pos.epi.matrix, c(pos.epitope_list[x]-1, original_position, pos.epitope_list[x+1]-1, epitope, allel, sequence))
			.GlobalEnv[["pos.epi.matrix"]] <- pos.epi.matrix_p
		}
		x <- x+2
	}
}

find_smallest_cluster_number <- function(cluster, seqs, seqs_to_comp){
	if (length(seqs_to_comp) > 1){
		for (i in seq(from = length(seqs), to = 1)){
			groups <- cutree(cluster, k=i)
			seqs_in_same_cluster <- test_if_in_one_cluster(groups, seqs_to_comp)		
			if(seqs_in_same_cluster){
				return (i/20)
			}
		}
	}
	else{
		return (1/20)
	}
}

test_if_in_one_cluster <- function(group, seqs){
	special_groups <- group[1:length(seqs)]
	identical <- all(special_groups[1] == special_groups)
	return (identical)
}


#---------------------------------------------------------------------------------------------
assocpointhierarchical <- structure(function(#Find possible epitopes
	### Searches in a given sequence for possible epitopes.
	##details<< For every position in the sequences a fisher's exact test of (Certain HLA type, Other HLA types, AA and not this AA) 	 is calculated and the p-value is pasted in a table.
	## Every HLA type has a column, every position a row. There is also a graphical output generated in which every HLA type has it's own page. In this output, a line is inserted, which shows the user where the p-value limit is considered above to be significant. 
	## Integrated in this graphic are the known epitopes if there are any and the calculated possible epitopes, which are calculated with the known HLA binding motifs. As an extras information a second graphical output is added in which the user can see the level of significant muations in a 9er sliding window.
	## Each point in this graphic is the number of mutations in the next nine sequence positions.
	## If pylogenetic comparison was choosen, the programm compares for every lowest p.value below a certain user given value the sequence similarity of the sequences with this trai with the similarity of all sequences with a wilcoxon test and adds the p.values in another column.
	## Uses EasyGregExpr from http://www.r-bloggers.com/developing-a-user-friendly-regular-expression-function-easygregexpr/
	## Also uses p.adjust from stats package to calculate some p-value correction additionally for the graphical output and as an extra column in the csv output.
	## Please consider that you have to use the dna/aa switch depending on your data!
	##seealso<< \code{\link{check_for_pos_epi_mut_core}}
	##note<< If a patient has a homozygot HLA Allel, then please change the second one to "00" (without ") instead!
	path_to_file_sequence_alignment = NULL,
	### a FASTA file with sequence data. Homzygot patiens have to have a 00 instead of the HLA typeFor reference please look in
	### example file.
	path_to_file_known_epitopes = NULL,
	### a csv file with known epitopes. For reference please look in example file.
	path_to_file_binding_motifs = NULL,
	### a csv file with HLA binding motifs if available. For reference please look in example file.
	save_name_pdf,
	### the file name of the result file in pdf format
	save_name_csv, 
	### the file name of the result file in csv format
	dna = FALSE,
	### if the data is in DNA or amino Acid Code
	patnum_threshold = 1,
	### the minimum number of patients of one HLA type to consider in the calculation.
	optical_significance_level = 0.05, 
	### the height of the optical horizontal line in the graphical output. It should be yout chosen p-value limit, like 0.05.
	star_significance_level = 0.001, 
	### the height of the invisible horizontal line above which alle points are marked as stars. Should be a high p-value limit,
	### like 0,001.
	### Should be the same as optical.significance_leve. 
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
	path_to_file_reference_sequence = NULL,
	### Reference sequence for better comparison of results
	one_feature=FALSE,
	### if there is only one identifier
	window_size=9,
	### size for the sliding window
	epi_plot=FALSE,
	### if epi plot
	own_model,
	###
	number_of_simulations, 
	###
	number_of_burnin, 
	### 
	number_of_chains, 
	### 
	response_inits, 
	### 
	further_response_data, 
	### 
	response_parameters 
	###
	){
	#print (reference_sequence)
	#print (path_to_file_p)
	#print (path_to_file_m)
	result <- find_possible_epitopes_inner_hierach(path_to_file_sequence_alignment, path_to_file_known_epitopes, path_to_file_binding_motifs, save_name_pdf, save_name_csv, dna, patnum_threshold, optical_significance_level, star_significance_level, optical_significance_level, optical_significance_level, A11, A12, A21, A22, B11, B12, B21, B22, path_to_file_reference_sequence, one_feature, window_size, epi_plot, own_model, number_of_simulations, number_of_burnin, number_of_chains, response_inits, further_response_data, response_parameters)
	return(result)

},ex=function(){
	ex <- system.file("extdata", "Example.fasta", package="SeqFeatR")
	ep <- system.file("extdata", "Example_epitopes.csv", package="SeqFeatR")
	co <- system.file("extdata", "Example_HLA_binding_motifs.csv", package="SeqFeatR")
	assocpointhierarchical(ex,
 ep,
 co,
 "epitope_results.pdf",
 "epitope_results.csv",
 FALSE,
 1,
 0.05,
 0.001,
 0.05,
 0.05,
 22,
 23,
 25,
 26,
 29,
 30,
 32,
 33,
 "bonferroni",
 FALSE,
 matrix_for_phylo = NULL,
 reference_sequence = NULL,
 one_feature=FALSE,
 window_size=9,
 pos_epi_plot=FALSE)
})

find_possible_epitopes_inner_hierach <- function(path_to_file_s = NULL, path_to_file_p = NULL, path_to_file_m = NULL, save_name, save_name_csv, dna = FALSE, patnum.threshold = 1, optical.significance_level = 0.05, star.significance_level = 0.001, pot.epitope_level = 0.05, window_threshold = 0.05, A11, A12, A21, A22, B11, B12, B21, B22, reference_sequence = NULL, one_ident, window_size, pos_epi_plot=FALSE, own_model, number_of_simulations, number_of_burnin, number_of_chains, response_inits, further_response_data, response_parameters){

	if (!is.null(reference_sequence) && reference_sequence == 0){
		reference_sequence <- NULL
	}
	if (!is.null(reference_sequence) && !(reference_sequence == "") ){
		ref_seq <- readAAStringSet(reference_sequence)
	}

	if (!is.null(path_to_file_s) && path_to_file_s == 0){
		path_to_file_s <- NULL
	}
	if (!is.null(path_to_file_s) && !(path_to_file_s == "") ){
		ep_wo_set_input_file_known_sequences(path_to_file_s)
	}

	if (!is.null(path_to_file_p) && path_to_file_p == 0){
		path_to_file_p <- NULL
	}
	if (!is.null(path_to_file_p) && !(path_to_file_p == "") ){
		ep_wo_set_input_file_known_epitopes(path_to_file_p)
	}

	if (!is.null(path_to_file_m) && path_to_file_m == 0){
		path_to_file_m <- NULL
	}
	if (!is.null(path_to_file_m) && !(path_to_file_m == "") ){
		ep_wo_set_input_file_known_pos_epi(path_to_file_m)
	}
	if(window_size <=2){
		print ("Beware! Window size smaller than 2! Will be set to 2")
		window_size <- 2
	}

	sequences <- .GlobalEnv[["StringSetsequences"]]
	epitopes <- .GlobalEnv[["epitopes"]]
	pos.epitopes <- .GlobalEnv[["pos.epitopes"]]
	consensus <- .GlobalEnv[["consensus"]]

	number_of_sequences <- length(sequences)

#find the shortest length of FASTA String, if the sequences are of different length
#min_FASTA_length holds this value
.GlobalEnv[["min_FASTA_length"]] <- min(width(sequences))

	min_FASTA_length <- .GlobalEnv[["min_FASTA_length"]]

#find out all allels that occur in given sequences
allels.all <- c()
allels.count <- c()

new_sequence_names <- c()

if (one_ident){
	for(i in 1:length(sequences)){
  	#the used values depend on fasta file; they possibly have to be changed - the same holds for the block some lines below
  	#tropismus.all holds all tropismus from the FASTA files, even duplicates
  	#tropismus holds no duplicates and none which are only a letter, without number (if the type is 00)
  	#tropismus.count counts how often the certain alles is there
  	allels.all <- c(allels.all, strsplit(names(sequences)[i], ";")[[1]][2])
	new_sequence_names <- c(new_sequence_names, strsplit(names(sequences)[i], ";")[[1]][2])
	}
}else{
	for(i in 1:length(sequences)){
	 	#the used values depend on fasta file; they possibly have to be changed - the same holds for the block some lines below
	  	#allels.all holds all allels from the FASTA files, even duplicates
	  	#allels holds no duplicates and none which are only a letter, without number (if the type is 00)
	  	#allels.count counts how often the certain alles is there
		f <- paste("A", remove.zeroes(substr(names(sequences)[i], A11, A12)), sep="")
		s <- paste("A", remove.zeroes(substr(names(sequences)[i], A21, A22)), sep="")
		t <- paste("B", remove.zeroes(substr(names(sequences)[i], B11, B12)), sep="")
		o <- paste("B", remove.zeroes(substr(names(sequences)[i], B21, B22)), sep="")
	  	allels.all <- c(allels.all, f, s, t, o)
		new_sequence_names <- c(new_sequence_names, paste(f, s, t, o))
	}
}

allels <- sort(unique(allels.all[allels.all!="A" & allels.all!="B"]))

#print (allels)

for(i in 1:length(allels)){
  allels.count <- c(allels.count, sum(allels.all==allels[i]))
}

if (dna){
	all_of_letters <- c("A","C","G","T", "-", "X")
}else{
	all_of_letters <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-", "X")
}

#save all p-values here
all.pvals <- c()

#save all aas here
all.aas <- c()

#save the number of aas here
aa.all_level <- c()

##save fishertable for best AA here
all.aa_1 <- c()
all.aa_2 <- c()
all.aa_3 <- c()
all.aa_4 <- c()
all.odds <- c()

#save stattest
all.stattest <- c()

#save all above star level here
above.star.all <- c()

#save all pos.epi here
#pos.epi.matrix <- data.frame(start=numeric(), original_position=numeric(), stop=numeric(), epitope=character(), allel=character(), protein_sequence=character())
pos.epi.matrix <- matrix(nrow = 1, ncol = 6)
colnames(pos.epi.matrix) <- (c("start", "original_position", "stop", "epitope", "allel", "protein_sequence") )
.GlobalEnv[["pos.epi.matrix"]] <- pos.epi.matrix
all.pos.epi <- c()

#do the following for each allel with a number of patients over a certain threshold
pdf(paste(save_name, sep=""), width = 11.69, height = 18.27)

# for each allel (where allel.num is the index of the allel in allels):
# if the count of this allel is above patnum.threshold.
# allels.counter counts the allels which are abouve the patnum.threshold
allels.counter <- 0

sequences_names <- names(sequences)

p.values.all <- c()

for(allel.num in 1:length(allels)){
  	if(allels.count[allel.num]>patnum.threshold){
		probs.max <- c()
		probs.norm <- c()
		aa_values_m <- c()

		##### CREATE HUGE MATRIX (# HLA vs not THIS HLA, positions as columns, amino acids/nucleotides as sheets) WITH ALL INFOS FOR THE MODEL!
		#### Reminder: All HLA in one huge table is a big problem with dopuble data, because an amino acid is sorted into more than one bin, which violates the multinominal distribution!

		huge_matrix <- array(NA, dim=c(2, length(all_of_letters), min_FASTA_length), dimnames=list(c(allels[allel.num], paste("!", allels[allel.num])), all_of_letters, c(1:min_FASTA_length)))

		print (allel.num)
		current_allel <- allels[allel.num]
		has_allel <- grep(current_allel, new_sequence_names)
		has_no_allel <- grep(current_allel, new_sequence_names, invert = T)
		#table for illustration
		allels.counter <- allels.counter + 1
		# for each letter in FASTA:
	    	for (j in 1:min_FASTA_length){
			just_this_position <- subseq(sequences, j, j)
			mat <- as.matrix(just_this_position)
			subset_has_allel <- mat[has_allel]
			subset_not_has_allel <- mat[has_no_allel]
			for(k in 1:length(all_of_letters)){
				aa <- all_of_letters[k]
				allel_has_aa <- length(which(subset_has_allel == aa))
				not_allel_has_aa <- length(which(subset_not_has_allel == aa))
				huge_matrix[1,k,j] <- allel_has_aa
				huge_matrix[2,k,j] <- not_allel_has_aa
			}
		}
		
		### Save result from following run:
		feat_results <- array(0, dim=c(length(all_of_letters), min_FASTA_length), dimnames=list(all_of_letters, c(1:min_FASTA_length)))
		### CREATE MODEL
		if(length(own_model) > 0){

			model <- readLines(own_model)
			model <- gsub("min_FASTA_length", min_FASTA_length, model, fixed=T)
			model <- gsub("length(all_of_letters)", length(all_of_letters), model, fixed=T)

			matrix <- huge_matrix
			nS <- t(apply(matrix, 3, rowSums))

			command <- parse(text=paste0("list(", response_inits, ")"))

			response.inits <- function(){
				eval(command)
			}
			
			response.data <- c("matrix", "nS", further_response_data)

			response.fit <- jags(response.data, inits=response.inits, response_parameters, model.file=textConnection(model), n.iter=number_of_simulations, n.burnin=number_of_burnin, n.chains=number_of_chains)

		}else{
			model_string <- paste("model{\n", "for (j in 1:",min_FASTA_length,"){\n", "for (i in 1:2){\n", "matrix[i,1:",length(all_of_letters),", j] ~ dmulti(theta[i,1:",length(all_of_letters),", j], nS[j, i])\n", "theta[i,1:",length(all_of_letters),",j] ~ ddirch(vcounts[1:",length(all_of_letters),"] + 0.03)\n","}\n","}\n", "for (i in 1:",length(all_of_letters),"){ \n", "vcounts[i] ~ dgamma(0.1,0.1)", "}\n", "}\n",sep="")
			matrix <- huge_matrix
			nS=t(apply(matrix, 3, rowSums))

			response.data=c("matrix", "nS")
			response.inits=function(){
				list("vcounts"=rep(1, length(all_of_letters)), "theta"=array(1/length(all_of_letters),dim=c(2, length(all_of_letters), min_FASTA_length)))
			}
			response.params=c("theta", "vcounts")

			response.fit <- jags(response.data, inits=response.inits, response.params, model.file=textConnection(model_string), n.iter=number_of_simulations, n.burnin=number_of_burnin, n.chains=number_of_chains)#, progress.bar = "none")
		}
		response.mcmc <- as.mcmc(response.fit)

		samples <- as.matrix(response.mcmc)

		for (amino in 1:length(all_of_letters)){
			for (sam in 1:min_FASTA_length){
				feat_results[amino, sam] <- mean(samples[,paste("theta[1,",amino,",",sam,"]", sep="")] > samples[,paste("theta[2,",amino,",",sam,"]", sep="")])
			}
		}
		
		data.df = data.frame(t(feat_results))
		probs.max <- c(probs.max, do.call(pmax, data.df))
		position <- apply(feat_results,2,which.max)
		probs.norm <- c(probs.norm, 100*probs.max)
		aa_values_m <- c(aa_values_m, all_of_letters[position])

		#evaluate
    		#get epitopes of this allel
		#1:2 is column, first number in allel.epitopes is row
    		allel.epitopes <- epitopes[epitopes[,3]==allels[allel.num],1:2]

    		#only evaluate if any epitopes have already been discovered for this allel
		#fill allel.epitopes2 with zeros for FASTA-length, then put max of p.values.min for all values between epitope border	
		allel.epitopes2 <- rep(0,min_FASTA_length)
		
    		if(length(allel.epitopes[,1])>0){
			for (i in 1:length(allel.epitopes[,1])){
        			allel.epitopes2[allel.epitopes[i,1]:allel.epitopes[i,2]] <- max(probs.max)
			}
		}
		all.pvals_one_allel<-c(probs.max)
		matrix.all.pvals_pre <- matrix(all.pvals_one_allel, ncol=1, dimnames=list(1:length(probs.max),allels[allel.num]))

		#post_production for every allel above threshold:
		#first: create a matrix with a column for position,
		#second: remove all rows with values beneath a certain threshold
		#third: calculate with the rest possible epitopes

		matrix.all.pvals <- cbind(1:min_FASTA_length, matrix.all.pvals_pre)

		matrix.one.allel <- c()
		#columns <- c(1, i+1)
		matrix.one.allel <- matrix.all.pvals
		
		rowsums.all.pvals <- rowSums(matrix.one.allel)

		removal_counter <- c() #for pval selection
		remove_counter <- c() # for sequence clean up

		for (i in 1:min_FASTA_length){
			remove_counter <- c(remove_counter,0)
			if (rowsums.all.pvals[i]<=(matrix.one.allel[i,1]+(1-pot.epitope_level))){ # this is becaus of sum of number + value
				removal_counter <- c(removal_counter, 1)
			}
			else{
				removal_counter <- c(removal_counter, 0)
			}
		}

		matrix.only.pvals.above.threshold <- matrix.one.allel[!removal_counter, ]
	
		if (is.null(names(matrix.only.pvals.above.threshold))){
			name <- colnames(matrix.only.pvals.above.threshold)
		}else{
			name <- names(matrix.only.pvals.above.threshold)
		}
		allel <- name[2]

		if(one_ident || is.null(pos.epitopes)){
			pos.epitope.plot <- c(0)
		}

		if (!is.null(name) && !is.null(pos.epitopes)){		
			pos.epitopes.f.allel <- c()
			for (i in 1:length(pos.epitopes[,1])){
				pos.epi <- gsub("A[*]0*", "A", substr(pos.epitopes[i,1],1, 4)) # replaces the *and leading 0s in HLA type for possible epitope
				pos.epi <- gsub("B[*]0*", "B", pos.epi)		
				if (pos.epi == allel){
					pos.epitopes.f.allel <- c(pos.epitopes.f.allel, pos.epitopes[i,2])
				}
			}

			pos.epitope_list <- c()
			if ((!is.null(pos.epitopes.f.allel)) && (length(matrix.only.pvals.above.threshold)>0)){
				for (pos.epitope in pos.epitopes.f.allel){
					pos.epitope_list <- c(pos.epitope_list, find_epitope_in_given_sequence(pos.epitope, matrix.only.pvals.above.threshold, sequences[has_allel], allel))
				}
			}
			pos.epitope.plot <- create_pos_epi_plot(pos.epitope_list, min_FASTA_length)
		}

		#window-slinding-thingy
		number_of_mutations_in_sequence <- c()
		thresh <- window_threshold/100
		for (i in 1:(length(all.pvals_one_allel)-(window_size-1))){
			number_of_mutations_in_window <- 0
			for (j in 0:(window_size-1)){
				value <- all.pvals_one_allel[i+j]
				if (value > thresh){
					number_of_mutations_in_window <- (number_of_mutations_in_window + 1)
				}
			}
			number_of_mutations_in_sequence <- c(number_of_mutations_in_sequence, ((1.0/window_size)*number_of_mutations_in_window))
		}

		#### GRAPHICS

		if (pos_epi_plot){
			par(fig=c(0,0.75,0.55,1))#, new=TRUE)
			if (is.null(path_to_file_m) & !is.null(path_to_file_p)){
				plot (panel.first=
				c(lines(allel.epitopes2, type="l", col="red")), #known epitopes in red
				number_of_mutations_in_sequence, type="l", col = "black", main=paste(allels[allel.num],"(",allels.count[allel.num]," patient(s))",sep=""),ylab="fraction of significant Probability", xlab="", xlim=c(0, min_FASTA_length), ylim=c(0,1)) # plot for percental number of mutations in 9er
				par(fig=c(0,1,0.6,1))
				legend("right", c("known epitopes", paste("# mutations in ", window_size,"er \nsliding window", sep="")),lty=c(1,1,1),lwd=c(1.5,1.5,1.5),col=c("red","black"),bty="n", cex=0.8, xjust=0)
			}else if (!is.null(path_to_file_m) & is.null(path_to_file_p)){
				plot (panel.first=
				c(lines(pos.epitope.plot, type="l", col = "darkgoldenrod1")), #known epitopes in red
				number_of_mutations_in_sequence, type="l", col = "black", main=paste(allels[allel.num],"(",allels.count[allel.num]," patient(s))",sep=""),ylab="fraction of significant Probability", xlab="", xlim=c(0, min_FASTA_length), ylim=c(0,1)) # plot for percental number of mutations in 9er
				par(fig=c(0,1,0.6,1))
				legend("right", c("possible new epitopes", paste("# mutations in ", window_size,"er \nsliding window", sep="")),lty=c(1,1,1),lwd=c(1.5,1.5,1.5),col=c("darkgoldenrod1","black"),bty="n", cex=0.8, xjust=0)
			}else if (is.null(path_to_file_p) & is.null(path_to_file_m)){
				plot (number_of_mutations_in_sequence, type="l", col = "black", main=paste(allels[allel.num],"(",allels.count[allel.num]," patient(s))",sep=""),ylab="fraction of significant Probability", xlab="", xlim=c(0, min_FASTA_length), ylim=c(0,1)) # plot for percental number of mutations in 9er
				par(fig=c(0,1,0.6,1))
				legend("right", c(paste("# mutations in ", window_size,"er \nsliding window", sep="")),lty=c(1,1,1),lwd=c(1.5,1.5,1.5),col=c("black"),bty="n", cex=0.8, xjust=0)

			}else{
				plot (panel.first=
				c(lines(allel.epitopes2, type="l", col="red"), #known epitopes in red
				lines(pos.epitope.plot, type="l", col = "darkgoldenrod1")), # possible epitopes in darkgoldenrod1
				number_of_mutations_in_sequence, type="l", col = "black", main=paste(allels[allel.num],"(",allels.count[allel.num]," patient(s))",sep=""),ylab="fraction of significant Probability", xlab="", xlim=c(0, min_FASTA_length), ylim=c(0,1)) # plot for percental number of mutations in 9er
				par(fig=c(0,1,0.6,1))
				legend("right", c("known epitopes","possible new epitopes", paste("# mutations in ", window_size,"er \nsliding window", sep="")),lty=c(1,1,1),lwd=c(1.5,1.5,1.5),col=c("red","darkgoldenrod1","black"),bty="n", cex=0.8, xjust=0)
			}
			par(fig=c(0,0.75,0.0,0.6), new=TRUE)
		}	
		plot(probs.norm, type="p", col="black", pch=20, ylab="Probability[%]", xlab="Alignment position", xlim=c(0, min_FASTA_length))
		if (!pos_epi_plot){
			mtext(paste(allels[allel.num],"(",allels.count[allel.num]," patient(s))",sep=""))
		}
		corrected_opti <- optical.significance_level
		corrected_star <- star.significance_level
		abline(h=optical.significance_level, col = "gray60", lty=2)
		probs.above_sig_level_norm <- probs.norm
		probs.above_star_level_norm <- probs.norm
		for (i in 1:length(probs.norm)){
			if (probs.norm[i] >= (corrected_opti)){
				if (probs.norm[i] >= (corrected_star)){
					probs.above_star_level_norm[i] <- probs.norm[i]
				}
				else{
					probs.above_sig_level_norm[i] <- probs.norm[i]
					probs.above_star_level_norm[i] <- NA
				}
			}
			else {
				probs.above_sig_level_norm[i] <- NA
				probs.above_star_level_norm[i] <- NA
			}
		}
		vert_lines <- which(probs.above_star_level_norm > 0)
		for (sp in 1:length(vert_lines)){
			abline(v=vert_lines[sp], col = "gray60", lty=2)
		}
		lines(probs.above_sig_level_norm, type="p", col = "red", pch=20) #points above a certain significance level
		lines(probs.above_star_level_norm, type="p", col = "red", pch=8) #stars above a certain significance level

		all.pvals<-c(all.pvals,probs.norm)
		all.aas <- cbind(all.aas, aa_values_m)
		above.star.all <- cbind(above.star.all, probs.above_star_level_norm)
	}
}

all.pvals.corrected <- all.pvals

#forth: write this epitope in a csv
all.pos.epi <- subset(pos.epi.matrix, !duplicated(pos.epi.matrix))
#all.pos.epi <- pos.epi.matrix[!duplicated(pos.epi.matrix)]
#print (all.pos.epi)
if (nrow(all.pos.epi)>1){
	all.pos.epi <- all.pos.epi[2:nrow(all.pos.epi),]
	#print(all.pos.epi)
}

result_output <- (matrix(all.pvals, ncol=sum(allels.count>patnum.threshold), dimnames=list(1:length(probs.max),allels[allels.count>patnum.threshold])))

all.pvals.c_output <- (matrix(all.pvals.corrected, ncol=sum(allels.count>patnum.threshold), dimnames=list(1:length(probs.max),allels[allels.count>patnum.threshold])))


result <- c()

if (!is.null(reference_sequence) && !(reference_sequence == "") ){
	ref_row <- as.character(subseq(ref_seq, 1, 1))

	for (j in 2:width(ref_seq)){
		ref_row <- rbind(ref_row, as.character(subseq(ref_seq, j, j)))
	}
}



if (!is.null(reference_sequence) && !(reference_sequence == "") ){
	if (is.null(dim(all.aas))){
		result <- cbind(result_output, all.aas, as.character(subseq(ref_seq, j, j)))
	}
	else {	
		for (i in 1:ncol(result_output)){
			result <- cbind(result, result_output[,i], scale(result_output[,i],center=TRUE,scale=TRUE)[,1], all.aas[,i], ref_row)
		}
	}
}else{
	if (is.null(dim(all.aas))){
		result <- cbind(result_output, all.aas)
	}
	else {	
		for (i in 1:ncol(result_output)){
			result <- cbind(result, result_output[,i], scale(result_output[,i],center=TRUE,scale=TRUE)[,1], all.aas[,i])
		}
	}


}

colname <- c()
if (!is.null(reference_sequence) && !(reference_sequence == "") ){
	for (i in 1:(ncol(result)/4)){
		colname <- c(colname, allels[allels.count>patnum.threshold][i], "z-scores", "AA with lowest p value", "ref_seq")
	}
}else{
	for (i in 1:(ncol(result)/3)){
		colname <- c(colname, allels[allels.count>patnum.threshold][i], "z-scores", "AA with lowest p value")
	}
}

colnames(result) <- colname

write.csv2(result, paste(save_name_csv, sep=""))

#write all pvalues in csv file

#we have finished our pdf writing
dev.off()

#print (result)
### returns the result as a table of HLA type vs sequence position
return (result)
}


#assocpointhierarchical(
#"HBV core 12-8-2014 complete final.fasta",#../inst/extdata/Example_aa.fasta",
#"",#../inst/extdata/Example_epitopes_aa.csv",
#"",#../inst/extdata/Example_HLA_binding_motifs_aa.csv",
#save_name_pdf = "HBVhi_results.pdf",
#save_name_csv = "HBVhi_results.csv",
#dna = FALSE,
#patnum_threshold = 1,
#optical_significance_level = 97, 
#star_significance_level = 99,
#A11 = 12,
#A12 = 13,
#A21 = 15,
#A22 = 16,
#B11 = 19,
#B12 = 20,
#B21 = 22,
#B22 = 23,
#path_to_file_reference_sequence = "",#../inst/extdata/Example_reference_aa.fasta",
#one_feature=FALSE,
#window_size=9,
#epi_plot=TRUE,
#own_model="../inst/extdata/Example_model_for_bayes.txt",
#number_of_simulations=20,
#number_of_burnin=floor(20/2),
#number_of_chains=2,
#response_inits="'vcounts' = rep(1, length(all_of_letters)), 'theta' = array(1/length(all_of_letters),dim=c(2, length(all_of_letters), min_FASTA_length))",
#further_response_data=c(),
#response_parameters=c("theta", "vcounts")
#)
