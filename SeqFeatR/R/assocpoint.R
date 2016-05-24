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

bayes_fac <- function(y, K, m){

    rdirichlet = function(n, alpha) {
        l <- length(alpha)
        x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
        sm <- x %*% rep(1, l)
        return(x/as.vector(sm))
    }
    ldirichlet = function(alpha) {
        return(rowSums(lgamma(alpha)) - lgamma(rowSums(alpha)))
    }
    yc = colSums(y)
    yr = rowSums(y)
    n = sum(yc)
    d = dim(y)
    I = d[1]
    J = d[2]
    etaA = rdirichlet(m, yr + 1)
    etaB = rdirichlet(m, yc + 1)
    Keta = c()
    KetaY = c()
    for (i in 1:I) {
        for (j in 1:J) {
            Keta = cbind(Keta, K * etaA[, i] * etaB[, j])
            KetaY = cbind(KetaY, K * etaA[, i] * etaB[, j] + 
                y[i, j])
        }
    }
    logint = ldirichlet(KetaY) - ldirichlet(Keta)
    for (i in 1:I) logint = logint - yr[i] * log(etaA[, i])
    for (j in 1:J) logint = logint - yc[j] * log(etaB[, j])
    int = exp(logint)
    return(list(bf = mean(int), nse = sd(int)/sqrt(m)))
}



#search for "change" for points where changes are possible
#please check before executing

#libs
library(Biostrings)
library(plyr)

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

KD <- function(mat) {
        rs = rowSums(mat)
        cs = colSums(mat)
        eta = rs %*% t(cs)/sum(mat)
        Kd = sum(mat) - sum(abs(eta - mat))
        Kd
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
		groups <- cutree(cluster, k=seq(from = (length(seqs)-1), to = 1))
		result_vectors <- c()
		for (i in 1:ncol(groups)){
			seqs_in_same_cluster <- test_if_in_one_cluster(groups[,i], seqs_to_comp)
			result_vectors <- rbind(result_vectors, seqs_in_same_cluster)
			#if(seqs_in_same_cluster){
			#	return ((ncol(groups)-i)/length(seqs))
			#}
		}
		rS <- rowSums(result_vectors)
		positions <- seq(1:ncol(groups))
		return (lm(rS ~ positions)$coefficients[1]*-1)
	}
	else{
		return (NA)
	}
}

test_if_in_one_cluster <- function(group, seqs_to_comp){
	group_number <- 0
	number_of_member <- 0
	result_vector <- c(0,0,0,0,0,0)
	for(i in 1:max(group)){
		group_with_max <- group[group == i]
		names_of <- names(group_with_max)
		is_inside <- as.numeric(names(seqs_to_comp) %in% names_of)
		#print (is_inside)
		if(100/length(is_inside)*sum(is_inside) >= 40){
			result_vector[1] <- 1
		}
		if(100/length(is_inside)*sum(is_inside) >= 50){
			result_vector[2] <- 1
		}
		if(100/length(is_inside)*sum(is_inside) >= 60){
			result_vector[3] <- 1
		}
		if(100/length(is_inside)*sum(is_inside) >= 70){
			result_vector[4] <- 1
		}
		if(100/length(is_inside)*sum(is_inside) >= 80){
			result_vector[5] <- 1
		}
		if(100/length(is_inside)*sum(is_inside) >= 90){
			result_vector[6] <- 1
		}
	}
	return (result_vector)
}

addtolist <- function(items, counter){

	p.values.all.all <- .GlobalEnv[["p.values.all.all"]]

	p.values.all.all[[counter]] <- items

	.GlobalEnv[["p.values.all.all"]] <- p.values.all.all
}


#---------------------------------------------------------------------------------------------
assocpoint <- structure(function(#Find possible epitopes
	### Searches in a given sequence for possible epitopes.
	##details<< For every position in the sequences a fisher's exact test of (Certain HLA type, Other HLA types, AA and not this AA) 	 is calculated and the p-value is pasted in a table.
	## Every HLA type has a column, every position a row. There is also a graphical output generated in which every HLA type has it's own page. In this output, a line is inserted, which shows the user where the p-value limit is considered above to be significant. 
	## Integrated in this graphic are the known epitopes if there are any and the calculated possible epitopes, which are calculated with the known HLA binding motifs. As an extras information a second graphical output is added in which the user can see the level of significant muations in a 9er sliding window.
	## Each point in this graphic is the number of mutations in the next nine sequence positions.
	## If pylogenetic comparison was choosen, the programm compares for every lowest p.value below a certain user given value the sequence similarity of the sequences with this trai with the similarity of all sequences with a Phylogenetic biason test and adds the p.values in another column.
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
	multiple_testing_correction = "bonferroni",
	### the statistical correction applied to the p-values. Input can be: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
	bayes_factor=FALSE,
	### if bayes factor should be applied instead of Fisher's exact test for association. See details.
	constant_dirichlet_precision_parameter=FALSE, 
	### if K is constant
	dirichlet_precision_parameter=20,
	### The Dirichlet precision parameter for evaluation with Bayes factor. See details.
	phylo_bias_check = FALSE,
	### if the sequences with the certain trait and Amino Acid with the lowest p.value should be compared in sequence similarity with all sequences. See details.
	path_to_file_reference_sequence = NULL,
	### Reference sequence for better comparison of results
	one_feature=FALSE,
	### if there is only one identifier
	window_size=9,
	### size for the sliding window
	epi_plot=FALSE
	){
	#print (path_to_file_reference_sequence)
	#print (multiple_testing_correction)
	#print (path_to_file_m)
	result <- find_possible_epitopes_inner(path_to_file_sequence_alignment, path_to_file_known_epitopes, path_to_file_binding_motifs, save_name_pdf, save_name_csv, dna, patnum_threshold, optical_significance_level, star_significance_level, optical_significance_level, optical_significance_level, A11, A12, A21, A22, B11, B12, B21, B22, multiple_testing_correction, bayes_factor, constant_dirichlet_precision_parameter, dirichlet_precision_parameter, phylo_bias_check, path_to_file_reference_sequence, one_feature, window_size, epi_plot)
	return(result)

},ex=function(){
	ex <- system.file("extdata", "Example.fasta", package="SeqFeatR")
	ep <- system.file("extdata", "Example_epitopes.csv", package="SeqFeatR")
	co <- system.file("extdata", "Example_HLA_binding_motifs.csv", package="SeqFeatR")
	assocpoint(ex,
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
 reference_sequence = NULL,
 one_feature=FALSE,
 window_size=9,
 pos_epi_plot=FALSE)
})

find_possible_epitopes_inner <- function(path_to_file_s = NULL, path_to_file_p = NULL, path_to_file_m = NULL, save_name, save_name_csv, dna = FALSE, patnum.threshold = 1, optical.significance_level = 0.05, star.significance_level = 0.001, pot.epitope_level = 0.05, window_threshold = 0.05, A11, A12, A21, A22, B11, B12, B21, B22, statistical_correction = "bonferroni", bayes_factor, constant_dirichlet_precision_parameter, dirichlet_precision_parameter=20, with_phylogenetic_comparison = FALSE, reference_sequence = NULL, one_ident, window_size, pos_epi_plot=FALSE){

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

	if (with_phylogenetic_comparison){
		StringSetsequences_p <- .GlobalEnv[["StringSetsequences"]]
		StringSetsequences <- StringSetsequences_p
		counter <- 0
		dist.vec <- dist.vec <- stringDist(StringSetsequences, method = "levenshtein")
		mean_dist_all <- mean(dist.vec)
	}

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

#possible one letter codes for amino acids B is for amBigious
#if (!dna){
#	aacids <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-", "B", "X")
#}else {
#	aacids <- c("A","C","G","T","-","R", "Y", "M", "K", "S", "W", "N")
#}

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

#all standard errors for bayes factor
all_sses <- c()

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

#all svings
p.values.all.all <- list(NULL)
length(p.values.all.all) <- length(allels)

.GlobalEnv[["p.values.all.all"]] <- p.values.all.all

for(allel.num in 1:length(allels)){
	#print (allels[allel.num])
  	if(allels.count[allel.num]>patnum.threshold){
		#if (allels[allel.num] == "A24"){
		print (allel.num)
		current_allel <- allels[allel.num]
		has_allel <- grep(current_allel, new_sequence_names)
		has_no_allel <- grep(current_allel, new_sequence_names, invert = T)
		#table for illustration
		allels.counter <- allels.counter + 1
		probs.max <- c()
		p.values.all <- c()
    		p.values.min <- c()
		odds.ratio.all <- c()
		aa_values_m <- c()
		aa_level <- c()
		aa_1 <- c()
		aa_2 <- c()
		aa_3 <- c()
		aa_4 <- c()
		p.values.norm <- c()
    		p.values2 <- c()
		probs.norm <- c()
    		probs2 <- c()
		sse_all <- c()
		stattest <- c()
		# for each letter in FASTA:
   		for (j in 1:min_FASTA_length){
			#print (j)
			#mir <- occ
      			#for each point in sequence do the following:
      			#for certain allel calculate fisher's exact test's minimum p-value among the tests for all amino acids (aa/!aa ; allel/!allel)
      			p.values <- c()
			probs <- c()
			sse <- c()
			odds.ratio <- c()
			aa_values <- c()
			aa_values2 <- c()
			aa_values3 <- c()
			aa_values4 <- c()
			# for each aa in length of possible aas:
			just_this_position <- subseq(sequences, j, j)
			mat <- as.matrix(just_this_position)
			#print (mat)
			aacids <- unique(mat[,1])
			subset_has_allel <- mat[has_allel]
			subset_not_has_allel <- mat[has_no_allel]
			#print (subset_has_allel)
			#print (subset_not_has_allel)
			#count from tabel for each FASTA letter and then each aa for the allel, the occurence
			# 1,1 is: occurency of aa in sequence in allel
			# 1,2 is: occurency of all aa in sequence in allel - 1,1
			# 2,1 is: occurency of aa in sequence in all allels - 1,1
			# 2,2 is: occurency of all aas in sequence - 1,2 - 2,1 + 1,1
			#@result: list of p values, sorted by aa from aacids
      			for(k in 1:length(aacids)){
				aa <- aacids[k]
				aa_and_h <- length(which(subset_has_allel == aa))
				aa_and_n_h <- length(which(subset_not_has_allel == aa))
				n_aa_and_h <- length(which(subset_has_allel != aa))
				n_aa_and_n_h <- length(which(subset_not_has_allel != aa))
        			small.table = matrix(c(aa_and_h, n_aa_and_h, aa_and_n_h, n_aa_and_n_h),nrow=2)
        			#print (small.table)
				####BAYES - FACTOR
				if(!constant_dirichlet_precision_parameter){
					dirichlet_precision_parameter <- KD(small.table)
					if (dirichlet_precision_parameter == 0){
						dirichlet_precision_parameter <- 0.01
					}
				}
				if (bayes_factor){					
					m=1000
					res <- bayes_fac(small.table, as.numeric(dirichlet_precision_parameter), m)
					probs <- c(probs, res$bf)
					odds.ratio <- c(odds.ratio, ((aa_and_h * n_aa_and_n_h)/(n_aa_and_h * aa_and_n_h)))
					sse <- c(sse, res$nse)
				}else{
		  			#only calculate p-values from contingency tables with row sums above threshold or zero
		  			if(sum(small.table[1,])*sum(small.table[2,])!=0 & (sum(small.table[1,])<=rowsum.threshold | sum(small.table[2,])<=rowsum.threshold))
					next
					#cat (j, k, "\n")
					#print(small.table)
					#calculate fisher's exact test for each pair (tropis, acid)    
					test.result <- fisher.test(small.table)
					#print (as.character(subseq(ref_seq, j, j))) 
	       				p.values <- c(p.values, test.result$p.value)
	       				#print(test.result$p.value)
	       				#print(test.result$estimate)
					if (!is.null(reference_sequence) && !(reference_sequence == "") ){
						cat((allel.num), "\t", j, "\t", as.character(subseq(ref_seq, j, j)), "\t", (small.table[1,1]), "\t", (small.table[2,1]), "\t", (small.table[1,2]), "\t", (small.table[2,2]), "\t", aacids[k], "\t", test.result$p.value, "\n", file= "logging.csv", append=TRUE)
					}else{
						cat((allel.num), "\t", j, "\t", (small.table[1,1]), "\t", (small.table[2,1]), "\t", (small.table[1,2]), "\t", (small.table[2,2]), "\t", aacids[k], "\t", test.result$p.value, "\n", file= "logging.csv", append=TRUE)
					}
					odds.ratio <- c(odds.ratio, test.result$estimate)
				}
       				#print(test.result)
				aa_values <- c(aa_values, small.table[1,1])
				aa_values2 <- c(aa_values2, small.table[2,1])
				aa_values3 <- c(aa_values3, small.table[1,2])
				aa_values4 <- c(aa_values4, small.table[2,2])
     			}
			####BAYES - FACTOR
			if (bayes_factor){	
				probs.max <- c(probs.max, max(probs))
				position <- which(probs == max(probs))
				probs.norm <- c(probs.norm, max(probs))
				sse_all <- c(sse_all, sse[position])
				p.values.all.little <- cbind(probs, NA)
			}else{
      				p.values.min <- c(p.values.min, 1-min(p.values))
				position <- which(p.values == min(p.values))
				p.values.norm <- c(p.values.norm, min(p.values))
				p.values.all.little <- cbind(p.values, NA)
			}
			if (length(position)>1){
				random_pos <- sample(position, 1)
				aa_1 <- c(aa_1, paste("/", aa_values[random_pos]))
				aa_2 <- c(aa_2, paste("/", aa_values2[random_pos]))
				aa_3 <- c(aa_3, paste("/", aa_values3[random_pos]))
				aa_4 <- c(aa_4, paste("/", aa_values4[random_pos]))
				aa_values_m <- c(aa_values_m, aacids[random_pos])
				aa_level <- c(aa_level, paste("/", aa_values[random_pos]))
				odds.ratio.all <- c(odds.ratio.all, odds.ratio[random_pos])
				p.values.all.little[random_pos, 2] <- "min"
			}
			else{
				aa_1 <- c(aa_1, aa_values[position])
				aa_2 <- c(aa_2, aa_values2[position])
				aa_3 <- c(aa_3, aa_values3[position])
				aa_4 <- c(aa_4, aa_values4[position])
				aa_values_m <- c(aa_values_m, aacids[position])
				aa_level <- c(aa_level, aa_values[position])
				odds.ratio.all <- c(odds.ratio.all, odds.ratio[position])
				p.values.all.little[position, 2] <- "min"
			}
			#print (p.values.all)
			#print (odds.ratio.all)
			p.values.all <- rbind(p.values.all, p.values.all.little)
			####BAYES - FACTOR
			if (bayes_factor & with_phylogenetic_comparison){	
				print ("Nothing here yet (phylogenetic bias)")
			}else{
				if (p.values.norm[j] <= optical.significance_level && with_phylogenetic_comparison){
					singleAA <- as.character(subseq(StringSetsequences_p, j, j))
					seq_number_with_AA <- which(singleAA==aa_values_m[j])
					seqs_to_comp <- StringSetsequences_p[seq_number_with_AA]
		  			labls <- which(labels(dist.vec) %in% names(seqs_to_comp))
		  			dist_mut <- as.matrix(dist.vec)[labls,labls]
					number_of_cluster <- 1-(mean(dist_mut)/mean_dist_all)
					stattest <- c(stattest, as.character(number_of_cluster))		
				}
				else{
					stattest <- c(stattest, "/")
				}
				#print (position)
			}
    		}
		addtolist(p.values.all, allel.num)
		###### BAYES FACTOR
		if (bayes_factor){
			all.pvals<-c(all.pvals,probs.norm)
			aa.all_level <- cbind(aa.all_level, aa_level)
			all.aas <- cbind(all.aas, aa_values_m)
			all.aa_1 <- cbind(all.aa_1, aa_1)
			all.aa_2 <- cbind(all.aa_2, aa_2)
			all.aa_3 <- cbind(all.aa_3, aa_3)
			all.aa_4 <- cbind(all.aa_4, aa_4)
			all.odds <- cbind(all.odds, odds.ratio.all)
			all_sses <- cbind(all_sses, sse_all)

		}else{
			all.pvals<-c(all.pvals,p.values.norm)
			aa.all_level <- cbind(aa.all_level, aa_level)
			all.aas <- cbind(all.aas, aa_values_m)
			all.aa_1 <- cbind(all.aa_1, aa_1)
			all.aa_2 <- cbind(all.aa_2, aa_2)
			all.aa_3 <- cbind(all.aa_3, aa_3)
			all.aa_4 <- cbind(all.aa_4, aa_4)
			all.odds <- cbind(all.odds, odds.ratio.all)
			all.stattest <- cbind(all.stattest, stattest)
		}
  	}
}

p.values.all.all <- .GlobalEnv[["p.values.all.all"]]

	####BAYES - FACTOR
if (bayes_factor){
	all.pvals.corrected <- c()
}else{
	all.pvals.corrected <- c()
	all.pvals.corrected.all.all <- lapply(p.values.all.all,function(r){p.adjust(r[,1], method = statistical_correction)})
}

for(allel.num in 1:length(allels)){
  	if(allels.count[allel.num]>patnum.threshold){
		#if (allels[allel.num] == "A24"){

		p.values.all <- p.values.all.all[[allel.num]]

		if(bayes_factor){
			probs.max <- as.numeric(p.values.all[which(p.values.all[,2]=="min")])
			probs.norm <- probs.max
			all.pvals.corrected <- c(all.pvals.corrected, probs.norm)
		}else{
			all.pvals.corrected.all <- all.pvals.corrected.all.all[[allel.num]]
			p.values.min <- 1-as.numeric(all.pvals.corrected.all[which(p.values.all[,2]=="min")])
			p.values.norm <- as.numeric(all.pvals.corrected.all[which(p.values.all[,2]=="min")])

			all.pvals.corrected <- c(all.pvals.corrected, p.values.norm)
		}
	#	print (p.values.min)
	#	print (p.values.norm)
	#	mir <- occ
    		#evaluate
    		#get epitopes of this allel
		#1:2 is column, first number in allel.epitopes is row
    		allel.epitopes <- epitopes[epitopes[,3]==allels[allel.num],1:2]

    		#only evaluate if any epitopes have already been discovered for this allel
		#fill allel.epitopes2 with zeros for FASTA-length, then put max of p.values.min for all values between epitope border	
		allel.epitopes2 <- rep(0,min_FASTA_length)
		
		####BAYES - FACTOR
		if (bayes_factor){
			if(length(allel.epitopes[,1])>0)
			for (i in 1:length(allel.epitopes[,1])){
        			allel.epitopes2[allel.epitopes[i,1]:allel.epitopes[i,2]] <- 1
			}
			all.pvals_one_allel<-c(probs.max)
			matrix.all.pvals_pre <- matrix(all.pvals_one_allel, ncol=1, dimnames=list(1:length(probs.max),allels[allel.num]))
		}else{
	    		if(length(allel.epitopes[,1])>0)
	      		for (i in 1:length(allel.epitopes[,1])){
				allel.epitopes2[allel.epitopes[i,1]:allel.epitopes[i,2]] <- 1
			}
			all.pvals_one_allel<-c(p.values.min)
			matrix.all.pvals_pre <- matrix(all.pvals_one_allel, ncol=1, dimnames=list(1:length(p.values.min),allels[allel.num]))
		}
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
		#print (name)
		#Get column names
		allel <- name[2]

		if(one_ident || is.null(pos.epitopes)){
			pos.epitope.plot <- c(0)
		}

		if (!is.null(name) && !is.null(pos.epitopes)){		
			#cat ("allel",allel,"\n")
			pos.epitopes.f.allel <- c()
			for (i in 1:length(pos.epitopes[,1])){
				pos.epi <- gsub("A[*]0*", "A", substr(pos.epitopes[i,1],1, 4)) # replaces the *and leading 0s in HLA type for possible epitope
				pos.epi <- gsub("B[*]0*", "B", pos.epi)
				if (pos.epi == allel){
					pos.epitopes.f.allel <- c(pos.epitopes.f.allel, pos.epitopes[i,2])
				}
			}

			pos.epitope_list <- c()
			#print (seq)
			if ((!is.null(pos.epitopes.f.allel)) && (length(matrix.only.pvals.above.threshold)>0)){
				for (pos.epitope in pos.epitopes.f.allel){
					pos.epitope_list <- c(pos.epitope_list, find_epitope_in_given_sequence(pos.epitope, matrix.only.pvals.above.threshold, sequences[has_allel], allel))
				}
			}
			pos.epitope.plot <- create_pos_epi_plot(pos.epitope_list, min_FASTA_length)
			#print (pos.epitope.plot)
		}

		#window-slinding-thingy
		if (bayes_factor){
			number_of_mutations_in_sequence <- c()
			for (i in 1:(length(all.pvals_one_allel)-(window_size-1))){
				number_of_mutations_in_window <- 0
				for (j in 0:(window_size-1)){
					value <- all.pvals_one_allel[i+j]
					thresh <- 10^window_threshold
					if (value > thresh){
						number_of_mutations_in_window <- (number_of_mutations_in_window + 1)
					}
				}
				number_of_mutations_in_sequence <- c(number_of_mutations_in_sequence, ((1.0/window_size)*number_of_mutations_in_window))
			}
		}else{
			number_of_mutations_in_sequence <- c()
			for (i in 1:(length(all.pvals_one_allel)-(window_size-1))){
				number_of_mutations_in_window <- 0
				for (j in 0:(window_size-1)){
					value <- all.pvals_one_allel[i+j]
					thresh <- 1-window_threshold
					if (value > thresh){
						number_of_mutations_in_window <- (number_of_mutations_in_window + 1)
					}
				}
				number_of_mutations_in_sequence <- c(number_of_mutations_in_sequence, ((1.0/window_size)*number_of_mutations_in_window))
			}
		}

		####BAYES - FACTOR
		if (bayes_factor){	
			highest_prob_value <- max(probs.norm)
			probs.norm <- log10(probs.norm)
		}else{
			lowest_p_value <- min(p.values.norm)
			vector_for_yaxis <- get_10_based(min(p.values.norm))
		}
		####BAYES - FACTOR
		if (bayes_factor){
			if (pos_epi_plot){
				par(fig=c(0,0.75,0.55,1))#, new=TRUE)
				if (is.null(path_to_file_m) & !is.null(path_to_file_p)){
					plot (panel.first=
					c(lines(allel.epitopes2, type="l", col="red")), #known epitopes in red
					number_of_mutations_in_sequence, type="l", col = "black", main=paste(allels[allel.num]," (",allels.count[allel.num]," sequences)",sep=""),ylab="fraction of significant Bayes Factors", xlab="", xlim=c(0, min_FASTA_length), ylim=c(0,1)) # plot for percental number of mutations in 9er
					par(fig=c(0,1,0.6,1))
					legend("right", c("known epitopes", paste("# mutations in ", window_size,"er \nsliding window", sep="")),lty=c(1,1,1),lwd=c(1.5,1.5,1.5),col=c("red","black"),bty="n", cex=0.8, xjust=0)
				}else if (!is.null(path_to_file_m) & is.null(path_to_file_p)){
					plot (panel.first=
					c(lines(pos.epitope.plot, type="l", col = "darkgoldenrod1")), #known epitopes in red
					number_of_mutations_in_sequence, type="l", col = "black", main=paste(allels[allel.num]," (",allels.count[allel.num]," sequences)",sep=""),ylab="fraction of significant Bayes Factors", xlab="", xlim=c(0, min_FASTA_length), ylim=c(0,1)) # plot for percental number of mutations in 9er
					par(fig=c(0,1,0.6,1))
					legend("right", c("possible new epitopes", paste("# mutations in ", window_size,"er \nsliding window", sep="")),lty=c(1,1,1),lwd=c(1.5,1.5,1.5),col=c("darkgoldenrod1","black"),bty="n", cex=0.8, xjust=0)
				}else if (is.null(path_to_file_p) & is.null(path_to_file_m)){
					plot (number_of_mutations_in_sequence, type="l", col = "black", main=paste(allels[allel.num]," (",allels.count[allel.num]," sequences)",sep=""),ylab="fraction of significant Bayes Factors", xlab="", xlim=c(0, min_FASTA_length), ylim=c(0,1)) # plot for percental number of mutations in 9er
					par(fig=c(0,1,0.6,1))
					legend("right", c(paste("# mutations in ", window_size,"er \nsliding window", sep="")),lty=c(1,1,1),lwd=c(1.5,1.5,1.5),col=c("black"),bty="n", cex=0.8, xjust=0)

				}else{
					plot (panel.first=
					c(lines(allel.epitopes2, type="l", col="red"), #known epitopes in red
					lines(pos.epitope.plot, type="l", col = "darkgoldenrod1")), # possible epitopes in darkgoldenrod1
					number_of_mutations_in_sequence, type="l", col = "black", main=paste(allels[allel.num]," (",allels.count[allel.num]," sequences)",sep=""),ylab="fraction of significant Bayes Factors", xlab="", xlim=c(0, min_FASTA_length), ylim=c(0,1)) # plot for percental number of mutations in 9er
					par(fig=c(0,1,0.6,1))
					legend("right", c("known epitopes","possible new epitopes", paste("# mutations in ", window_size,"er \nsliding window", sep="")),lty=c(1,1,1),lwd=c(1.5,1.5,1.5),col=c("red","darkgoldenrod1","black"),bty="n", cex=0.8, xjust=0)
				}
				par(fig=c(0,0.75,0.0,0.6), new=TRUE)
			}	
			plot(probs.norm, type="p", col="black", pch=20, ylab="log10(Bayes Factor)", xlab="Alignment position", xlim=c(0, min_FASTA_length))
			if (!pos_epi_plot){
				mtext(paste(allels[allel.num]," (",allels.count[allel.num]," sequences)",sep=""))
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

		}else{	
			if (pos_epi_plot){
				par(fig=c(0,0.75,0.55,1))#, new=TRUE)
				if (is.null(path_to_file_m) & !is.null(path_to_file_p)){
					plot (panel.first=
					c(lines(allel.epitopes2, type="l", col="red")), #known epitopes in red
					number_of_mutations_in_sequence, type="l", col = "black", main=paste(allels[allel.num]," (",allels.count[allel.num]," sequences)",sep=""),ylab="fraction of significant p-values", xlab="", xlim=c(0, min_FASTA_length), ylim=c(0,1)) # plot for percental number of mutations in 9er
					par(fig=c(0,1,0.6,1))
					legend("right", c("known epitopes", paste("# mutations in ", window_size,"er \nsliding window", sep="")),lty=c(1,1,1),lwd=c(1.5,1.5,1.5),col=c("red","black"),bty="n", cex=0.8, xjust=0)
				}else if (!is.null(path_to_file_m) & is.null(path_to_file_p)){
					plot (panel.first=
					c(lines(pos.epitope.plot, type="l", col = "darkgoldenrod1")), #known epitopes in red
					number_of_mutations_in_sequence, type="l", col = "black", main=paste(allels[allel.num]," (",allels.count[allel.num]," sequences)",sep=""),ylab="fraction of significant p-values", xlab="", xlim=c(0, min_FASTA_length), ylim=c(0,1)) # plot for percental number of mutations in 9er
					par(fig=c(0,1,0.6,1))
					legend("right", c("possible new epitopes", paste("# mutations in ", window_size,"er \nsliding window", sep="")),lty=c(1,1,1),lwd=c(1.5,1.5,1.5),col=c("darkgoldenrod1","black"),bty="n", cex=0.8, xjust=0)
				}else if (is.null(path_to_file_p) & is.null(path_to_file_m)){
					plot (number_of_mutations_in_sequence, type="l", col = "black", main=paste(allels[allel.num]," (",allels.count[allel.num]," sequences)",sep=""),ylab="fraction of significant p-values", xlab="", xlim=c(0, min_FASTA_length)) # plot for percental number of mutations in 9er
					par(fig=c(0,1,0.6,1))
					legend("right", c(paste("# mutations in ", window_size,"er \nsliding window", sep="")),lty=c(1,1,1),lwd=c(1.5,1.5,1.5),col=c("black"),bty="n", cex=0.8, xjust=0)

				}else{
					plot (panel.first=
					c(lines(allel.epitopes2, type="l", col="red"), #known epitopes in red
					lines(pos.epitope.plot, type="l", col = "darkgoldenrod1")), # possible epitopes in darkgoldenrod1
					number_of_mutations_in_sequence, type="l", col = "black", main=paste(allels[allel.num]," (",allels.count[allel.num]," sequences)",sep=""),ylab="fraction of significant p-values", xlab="", xlim=c(0, min_FASTA_length), ylim=c(0,1)) # plot for percental number of mutations in 9er
					par(fig=c(0,1,0.6,1))
					legend("right", c("known epitopes","possible new epitopes", paste("# mutations in ", window_size,"er \nsliding window", sep="")),lty=c(1,1,1),lwd=c(1.5,1.5,1.5),col=c("red","darkgoldenrod1","black"),bty="n", cex=0.8, xjust=0)
				}
				par(fig=c(0,0.75,0.0,0.6), new=TRUE)
			}
			plot(p.values.norm, type="p", ylim=c(1.0,lowest_p_value), log="y", col="black", yaxt="n", pch=20, ylab="p-values", xlab="Alignment position", xlim=c(0, min_FASTA_length))
			if (!pos_epi_plot){
				mtext(paste(allels[allel.num]," (",allels.count[allel.num]," sequences)",sep=""))
			}
			axis(side=2, at=vector_for_yaxis, labels=vector_for_yaxis)
			corrected_opti <- optical.significance_level
			corrected_star <- star.significance_level
			abline(h=optical.significance_level, col = "gray60", lty=2)
			p.values.above_sig_level_norm <- p.values.norm
			p.values.above_star_level_norm <- p.values.norm
			for (i in 1:length(p.values.norm)){
				if (p.values.norm[i] <= (corrected_opti)){
					if (p.values.norm[i] <= (corrected_star)){
						p.values.above_star_level_norm[i] <- p.values.norm[i]
					}
					else{
						p.values.above_sig_level_norm[i] <- p.values.norm[i]
						p.values.above_star_level_norm[i] <- 0
					}
				}
				else {
					p.values.above_sig_level_norm[i] <- 0
					p.values.above_star_level_norm[i] <- 0
				}
			}
			vert_lines <- which(p.values.above_star_level_norm > 0)
			for (sp in 1:length(vert_lines)){
				abline(v=vert_lines[sp], col = "gray60", lty=2)
			}
			lines(p.values.above_sig_level_norm, type="p", col = "red", pch=20) #points above a certain significance level
			lines(p.values.above_star_level_norm, type="p", col = "red", pch=8) #stars above a certain significance level
		}
  	}
}

#forth: write this epitope in a csv
all.pos.epi <- subset(pos.epi.matrix, !duplicated(pos.epi.matrix))
#all.pos.epi <- pos.epi.matrix[!duplicated(pos.epi.matrix)]
#print (all.pos.epi)
if (nrow(all.pos.epi)>1){
	all.pos.epi <- all.pos.epi[2:nrow(all.pos.epi),]
	#print(all.pos.epi)
}

		####BAYES - FACTOR
if (bayes_factor){	
	result_output <- (matrix(all.pvals, ncol=sum(allels.count>patnum.threshold), dimnames=list(1:length(probs.max),allels[allels.count>patnum.threshold])))

	all.pvals.c_output <- (matrix(all.pvals.corrected, ncol=sum(allels.count>patnum.threshold), dimnames=list(1:length(probs.max),allels[allels.count>patnum.threshold])))
}else{
	result_output <- (matrix(all.pvals, ncol=sum(allels.count>patnum.threshold), dimnames=list(1:length(p.values.min),allels[allels.count>patnum.threshold])))
	all.pvals.c_output <- (matrix(all.pvals.corrected, ncol=sum(allels.count>patnum.threshold), dimnames=list(1:length(p.values.min),allels[allels.count>patnum.threshold])))
}

result <- c()

if (!is.null(reference_sequence) && !(reference_sequence == "") ){
	ref_row <- as.character(subseq(ref_seq, 1, 1))

	for (j in 2:width(ref_seq)){
		ref_row <- rbind(ref_row, as.character(subseq(ref_seq, j, j)))
	}
}

		####BAYES - FACTOR
if (bayes_factor){	
	if (!is.null(reference_sequence) && !(reference_sequence == "") ){
		if (is.null(dim(all.aas))){
			result <- cbind(result_output, all_sses, all.aas, all.odds, as.character(subseq(ref_seq, j, j)))
		}
		else {	
			for (i in 1:ncol(result_output)){
				result <- cbind(result, result_output[,i], scale(result_output[,i],center=TRUE,scale=TRUE)[,1], all_sses[,i], all.aas[,i], all.odds[,i], ref_row)
			}
		}
	}else{
		if (is.null(dim(all.aas))){
			result <- cbind(result_output, all_sses, all.aas, all.odds)
		}
		else {	
			for (i in 1:ncol(result_output)){
				result <- cbind(result, result_output[,i], scale(result_output[,i],center=TRUE,scale=TRUE)[,1], all_sses[,i], all.aas[,i], all.odds[,i])
			}
		}
	}

	colname <- c()
	if (!is.null(reference_sequence) && !(reference_sequence == "") ){
		for (i in 1:(ncol(result)/6)){
			colname <- c(colname, allels[allels.count>patnum.threshold][i], "z-scores", "standard error BF", "AA with lowest p value", "Odds ratio", "ref_seq")
		}
	}else{
		for (i in 1:(ncol(result)/5)){
			colname <- c(colname, allels[allels.count>patnum.threshold][i], "z-scores", "standard error BF", "AA with lowest p value", "Odds ratio")
		}
	}
}else{

	if (!is.null(reference_sequence) && !(reference_sequence == "") ){
		if (is.null(dim(all.aas))){
			result <- cbind(result_output, all.aas, all.odds, all.stattest, as.character(subseq(ref_seq, j, j)))
		}
		else {	
			for (i in 1:ncol(result_output)){
				result <- cbind(result, result_output[,i], all.pvals.c_output[,i], scale(result_output[,i],center=TRUE,scale=TRUE)[,1], all.aas[,i], all.odds[,i], all.stattest[,i], 
ref_row)
			}
		}
	}else{
		if (is.null(dim(all.aas))){
			result <- cbind(result_output, all.aas, all.odds, all.stattest)
		}
		else {	
			for (i in 1:ncol(result_output)){
				result <- cbind(result, result_output[,i], all.pvals.c_output[,i], scale(result_output[,i],center=TRUE,scale=TRUE)[,1], all.aas[,i], all.odds[,i], all.stattest[,i])
			}
		}


	}
	colname <- c()

	if (!is.null(reference_sequence) && !(reference_sequence == "") ){
		for (i in 1:(ncol(result)/7)){
			colname <- c(colname, allels[allels.count>patnum.threshold][i],"Corrected p-value", "z-scores", "AA with lowest p value", "Odds ratio", "Phylogenetic bias", "ref_seq")
		}
	}else{
		for (i in 1:(ncol(result)/6)){
			colname <- c(colname, allels[allels.count>patnum.threshold][i],"Corrected p-value", "z-scores", "AA with lowest p value", "Odds ratio", "Phylogenetic bias")
		}
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


#assocpoint(
#"Input_Assocpoint.fasta",
#"",
#"",
#save_name_pdf = "results.pdf",
#save_name_csv = "results.csv",
#dna = FALSE,
#patnum_threshold = 1,
#optical_significance_level = 0.05, 
#star_significance_level = 0.001,
#A11 = 10,
#A12 = 11,
#A21 = 13,
#A22 = 14,
#B11 = 17,
#B12 = 18,
#B21 = 20,
#B22 = 21,
#multiple_testing_correction = "none",
#bayes_factor=FALSE,
#constant_dirichlet_precision_parameter=FALSE,
#dirichlet_precision_parameter=200,
#phylo_bias_check = FALSE,
#path_to_file_reference_sequence = "",
#one_feature=FALSE,
#window_size=9,
#epi_plot=TRUE
#)
