#R Script written by Bettina Budeus, rokkaku-fight@gmx.de
#from 2012-06-19 to 2012--
#analysis of correlation of HLA and ancestral sequence to look for escape mutations or reversions
#!!!! BEWARE!!!! You have to name the fasta sequences so that phylo can read and write them correctly. It will replace any space charakter with _ and ignore colon!!!!!!

#search for "change" for points where changes are possible
#please check before executing

#libs
library(Biostrings)
library(plyr)
library(plotrix)
library(ape)
library(phangorn)
library(qvalue)

#filenames, change as required...
files.dir <- "./"
#all patients (all clades w/o consensus!)
#files.sequences <- "Anti-D_m_woX_same_length.fasta"
files.sequences <- "Anti-D_test_hlong.fasta"
#change: p-values will only be calculated for allels with a number of patients above threshold
patnum.threshold <- 0
#change: p-values will only be calculated from contingency tables with row sums above threshold or zero
rowsum.threshold <- -1
#change: threshold level for significance in fisher test search
thr.sig.fi <- 0.05


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

ht_get_input_file_known_sequences <- structure(function(
	###read aa sequences from fasta file - they have to be aligned in order to work!
	path_to_file
	### a FASTA file with sequence data. For reference please look in example file.
	){
	### the sequences from the FASTA file.
	sequences <- readAAStringSet(path_to_file)
	#c_sequences <- create_correct_FASTA_input(sequences)
	.GlobalEnv[["sequences"]] <- sequences
},ex=function(){
	ht_get_input_file_known_sequences("SeqFeatR/extdata/Example_epitopes.csv")	
})

#read tree file
ht_get_input_file_tree <- structure(function(
	### set the input file (tree file of the given sequences).
	path_to_file
	### a nexus file with tree data. For reference please look in example file.
	){
	### the tree.
	.GlobalEnv[["tree"]] <- read.nexus(path_to_file)
},ex=function(){
	ht_wo_get_input_file_tree("SeqFeatR/extdata/Example_tree.nh")
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

#calculates from stupid phydat the original sequences without this patterns thingy
#@input a phydat file
#@return a list of sequences, each character one column
get_original_sequence <- function(phydat){
	weight <- attributes(phydat)$weight
	index <- attributes(phydat)$index
	levels <- attributes(phydat)$levels
	names <- attributes(phydat)$names
	number_of_sequences <- length(names)
	number_of_letters <- length(index)
	new_seqs <- array(rep(0,number_of_sequences*number_of_letters), c(number_of_sequences,number_of_letters))
	rownames(new_seqs) <- names
	for (entry in 1:length(index)){
		for (sequence in 1:number_of_sequences){
			sequence_pattern <- phydat[sequence][[1]]
			new_seqs[sequence,entry] <- sequence_pattern[index[entry]]
		}
	}
	for (row in 1:nrow(new_seqs)){
		for (col in 1:ncol(new_seqs)){
			AAnr <- as.numeric(new_seqs[row,col][[1]])
			new_seqs[row,col] <- levels[AAnr]
		}
	}
	return (new_seqs)
}

#calculates from stupid phydat-phangorn the original sequences without this patterns thingy (differnd because of matrix not single array)
#@input a phydat-phangorn file, such as from ancestral sequence reconstruction
#@return a list of sequences, each character one column
get_original_sequence_ph <- function(phydat){
	weight <- attributes(phydat)$weight
	index <- attributes(phydat)$index
	levels <- attributes(phydat)$levels
	names <- attributes(phydat)$names
	number_of_sequences <- length(names)
	number_of_letters <- length(index)
	new_seqs <- array(rep(0,number_of_sequences*number_of_letters), c(number_of_sequences,number_of_letters))
	rownames(new_seqs) <- names
	for (sequence in 1:number_of_sequences){
		for (entry in 1:length(index)){
			sequence_pattern <- subset(phydat, sequence)[[1]]
			sequence_pattern_names <- rownames(sequence_pattern)
			if (is.null(sequence_pattern_names)){
				subset <- sequence_pattern[index[entry],]
				max <- max(subset)
				which <- (which(max == subset))
				letter <- levels[which]
				if (length(letter) > 1){
					new_seqs[sequence,entry] <- letter[1]
				}else{
					new_seqs[sequence,entry] <- letter
				}
			}
			else {
				new_seqs[sequence,entry] <- sequence_pattern_names[index[entry]]
			}
		}
	}
	return (new_seqs)
}


#function do make char array from the sequences
make_char_array <- function(){
	sequences <- .GlobalEnv[["sequences"]]
	min_FASTA_length <- .GlobalEnv[["min_FASTA_length"]]
	char_array <- array(rep(0,min_FASTA_length*length(sequences)),c(length(sequences),min_FASTA_length))
	for (sequence in 1:length(sequences)){
		for (i in 1:min_FASTA_length){
			#letter <- (substr(sequences[[sequence]]$seq, i, i))
			letter <- subseq(sequences[sequence], i, i)
			char_array[sequence,i] <- as.character(letter)
		}
		
	}
	return (char_array)
}

#function to make list of HLA-types for each sequence
make_HLA_types_array <- function(){
	A11 <- .GlobalEnv[["A11"]]
	A12 <- .GlobalEnv[["A12"]]
	A21 <- .GlobalEnv[["A21"]]
	A22 <- .GlobalEnv[["A22"]]
	B11 <- .GlobalEnv[["B11"]]
	B12 <- .GlobalEnv[["B12"]]
	B21 <- .GlobalEnv[["B21"]]
	B22 <- .GlobalEnv[["B22"]]
	sequences <- .GlobalEnv[["sequences"]]
	HLA_array <- array(rep(0,length(sequences)*4),c(length(sequences),4))
	for (sequence in 1:length(sequences)){
		a1 <- paste("A", remove.zeroes(substr(names(sequences)[sequence], A11, A12)), sep="")
		a2 <- paste("A", remove.zeroes(substr(names(sequences)[sequence], A21, A22)), sep="")
		b1 <- paste("B", remove.zeroes(substr(names(sequences)[sequence], B11, B12)), sep="")
		b2 <- paste("B", remove.zeroes(substr(names(sequences)[sequence], B21, B22)), sep="")
		HLA_array[sequence,1] <- a1
		HLA_array[sequence,2] <- a2
		HLA_array[sequence,3] <- b1
		HLA_array[sequence,4] <- b2
	}
	return (HLA_array)
}

#function to return the potential ancestral AA for given position. It searches the tip.labels for the given sequence and calucalts
get_ancestral_for_this_position <- function(sequence, position){
	sequences <- .GlobalEnv[["sequences"]]
	tree <- .GlobalEnv[["tree"]]
	tree.tip.label <- .GlobalEnv[["tree.tip.label"]]
	a_sequences <- .GlobalEnv[["a_sequences"]]
	tip.N <- which(tree.tip.label==names(sequence))
	ancestor <- Ancestors(tree, tip.N, type=c("parent"))
	ancestral_letter <- a_sequences[ancestor, position]
	return (ancestral_letter)

}

#function to create table with all positions for better lapply use
make_array <- function(){
	min_FASTA_length <- .GlobalEnv[["min_FASTA_length"]]
	anc_new_array <- c()
	for (i in 1:(min_FASTA_length)){
		anc_new_array <- c(anc_new_array, i)
		
	}
	return (anc_new_array)
}

#function to find the minmal number of AA pairs at a certain allel to garanty significance
get.threshold.f.fisher.pair <- function(number.sequences, number.pairs, allels.count){
	thr.sig.fi <- .GlobalEnv[["thr.sig.fi"]]
	patnum.threshold <- .GlobalEnv[["patnum.threshold"]]
	length_allels <- .GlobalEnv[["length_allels"]]
	allels <- .GlobalEnv[["allels"]]
	threshold <- c()
	x <- 0
	for (allel in 1:length(allels)){
		if (allels.count[allel] > patnum.threshold){
			min_threshold <- 0
			for (x in 0:allels.count[allel]){
				y <- allels.count[allel]
				l <- sum(allels.count)
				for (z in seq(x,(l-y))){
					#cat (x,y-x,z-x,l-y-z+x,"\n")
					small.table = matrix(c(x,y-x,z-x,l-y-z+x),nrow=2)
					if(sum(small.table[1,])*sum(small.table[2,])!=0 & (sum(small.table[1,])<=rowsum.threshold | sum(small.table[2,])<=rowsum.threshold))
					next
	       				#calculate fisher's exact test for each pair (allel, acid)    
	       				test.result <- fisher.test(small.table)
					if (test.result$p.value < thr.sig.fi){
						threshold <- rbind(threshold, c(x,y-x,z-x,l-y-z+x,test.result$p.value))
					}
				}
			}
		}
	}
	threshold <- subset(threshold, !duplicated(threshold))
	#print (threshold)
	#miracle <- ocure
	return (threshold)
}

#calculates the length of result list
calculate_length <- function(list){
	counter <- 0
	for (i in 1:length(list)){
		counter <- counter + nrow(list[[i]])
	}
	return (counter)
}

#create mini_table for calculating all relevant numbers for one position pair i,j
#@result: a table with all relevant data aka the sum for all occurings
create_mini_table_tr <- function(i){
	sequences <- .GlobalEnv[["sequences"]]
	patnum.threshold <- .GlobalEnv[["patnum.threshold"]]
	char_array <- .GlobalEnv[["char_array"]]
	HLA_array <- .GlobalEnv[["HLA_array"]]
	length_allels <- .GlobalEnv[["length_allels"]]
	allels <- .GlobalEnv[["allels"]]
	last_line <- c("sum")
	mini_table <- array(rep(0, (length_allels+2)),c(length_allels+2), list(c("Anc_New_P",allels,"sum")))
	mini_table <- rbind(mini_table, c(rep(0, (length_allels+2))))
	for (sequence in 1:length(sequences)){
      		a1 <- HLA_array [sequence,1]
      		a2 <- HLA_array [sequence,2]
      		b1 <- HLA_array [sequence,3]
      		b2 <- HLA_array [sequence,4]
		l1 <- char_array[sequence,i]
		l2 <- get_ancestral_for_this_position(sequence, i)
		line <- paste(l1,l2,sep = "")
		pos.a1 = which(allels == a1)+1
		pos.a2 = which(allels == a2)+1
		pos.b1 = which(allels == b1)+1
		pos.b2 = which(allels == b2)+1
		pos.aa = which(mini_table[,1] == line)
		if (!length(pos.aa)){
			mini_table <- rbind(mini_table, c(0))
			rows = nrow(mini_table)
			mini_table[rows,1] = line
			mini_table[rows,pos.a1] = 1
			mini_table[rows,pos.a2] = 1
			mini_table[rows,pos.b1] = 1
			mini_table[rows,pos.b2] = 1
			mini_table[rows,ncol(mini_table)] = sum(as.numeric(mini_table[rows,2:(ncol(mini_table)-1)]))
		}
		else{
			mini_table[pos.aa,pos.a1] = as.numeric(mini_table[pos.aa,pos.a1]) + 1
			mini_table[pos.aa,pos.a2] = as.numeric(mini_table[pos.aa,pos.a2]) + 1
			mini_table[pos.aa,pos.b1] = as.numeric(mini_table[pos.aa,pos.b1]) + 1
			mini_table[pos.aa,pos.b2] = as.numeric(mini_table[pos.aa,pos.b2]) + 1
			mini_table[pos.aa,ncol(mini_table)] = sum(as.numeric(mini_table[pos.aa,2:(ncol(mini_table)-1)]))
		}
	}
	#print (mini_table)
	for (col in 2:ncol(mini_table)){
		last_line <- c(last_line, sum(as.numeric(mini_table[,col])))
	}
	mini_table <- rbind(mini_table, last_line)
	mini_table_r <- mini_table[,(which(as.numeric(mini_table[nrow(mini_table),2:ncol(mini_table)])>patnum.threshold)+1)]
	mini_table_r <- cbind (mini_table[,1],mini_table_r)
	mini_table_r <- mini_table_r[3:nrow(mini_table_r),]
	print (mini_table)	
	return (mini_table_r)
}

#for controll and test purpose only!!!!
make_fisher_test <- function(a,b,c,d){
	small.table = matrix(c(a,b-a,c-a,d-b-c+a),nrow=2)
	#print (small.table)
	if(sum(small.table[1,])*sum(small.table[2,])!=0 & (sum(small.table[1,])<=rowsum.threshold | sum(small.table[2,])<=rowsum.threshold))
	next
	if(small.table[1,1]*small.table[2,2] < small.table[1,2]*small.table[2,1]){
	        small.table[2,1] <- small.table[2,1]+1
	        small.table[1,2] <- small.table[1,2]+1
	}
	else if(small.table[1,1]*small.table[2,2] > small.table[1,2]*small.table[2,1]){
		small.table[1,1] <- small.table[1,1]+1
	        small.table[2,2] <- small.table[2,2]+1
	}
	#calculate fisher's exact test for each pair (allel, acid)    
	test.result <- fisher.test(small.table)
	return (test.result$p.value)
}

#main function to call:
# for each position pair:
# generate mini_table
# check if entrys in mini_table are inside threshold for significant F Test
mai <- function(){
	anc_new_array <- .GlobalEnv[["anc_new_array"]]
	length_allels <- .GlobalEnv[["length_allels"]]
	result <- lapply(anc_new_array, FUN=inside_main_tr)
	result_wo_zero <- result[!sapply(result, is.null)]
	length_of_all <- calculate_length(result_wo_zero)
	result_matrix <- matrix(rep(0,length_of_all*4),ncol = 4,dimnames = list(NULL, c("allel", "Anc_New", "First Position", "p-value")))	
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
inside_main_tr <- function(col){
	threshold <- .GlobalEnv[["threshold"]]
	i = col[1]
	result <- c()
	#cat ("col", i, "\n")
	mini_table <- create_mini_table_tr(i)
	#print (mini_table)
	for (col in 2:(ncol(mini_table)-1)){
		for (row in 1:(nrow(mini_table)-1)){
			a <- as.numeric(mini_table[row,col])
			b <- as.numeric(mini_table[nrow(mini_table),col])
			c <- as.numeric(mini_table[row,ncol(mini_table)])
			d <- as.numeric(mini_table[nrow(mini_table),ncol(mini_table)])
			#cat (a,b-a,c-a,d-b-c+a,"\n")
			#print (threshold)
			is_in <- threshold[which(threshold[,1]==a & threshold[,2]==(b-a) & threshold[,3]==(c-a) & threshold[,4]==(d-c-b+a)),]
			#print (is_in)
			if (length(is_in)>0){
				p_value <- is_in[5]
				#print (mini_table)
				#print (col)
				if (a == 0){
					result <- rbind(result,c(colnames(mini_table)[col],mini_table[row,1],i,-(p_value)))
				}else{
					line <- c(colnames(mini_table)[col],mini_table[row,1],i,p_value)
					#print (line)
					result <- rbind(result,line)
					#print (result)
				}
			}
		}
	#print (result)
	}
	return (result)
}

create_ancestral_sequences <- function(path_to_file_s){
	sequences_for_tree <- read.phyDat(path_to_file_s, type="AA", format="fasta")
	tree <- .GlobalEnv[["tree"]]
	#print (tree)
	#print (names(sequences_for_tree))
	fit = pml(tree, sequences_for_tree)
	#anc.ml = ancestral.pml(fit, type = "ml")
	anc.acctran = ancestral.pars(tree, sequences_for_tree, "ACCTRAN")
	#anc.mpr = ancestral.pars(tree, sequences_for_tree, "MPR")
	return (anc.acctran)
}


#get_length_of_nested_list <- function(nested_list){
#	leng <- 0
#	for (i in 1: length(nested_list)){
#		leng <- (leng + nrow(nested_list[[i]]))
#	}
#	return (leng)
#}

#___________________________________________________________

comparewithancestral <- structure(function(#Find possible escape mutations
	### Takes a tree and a set of sequences and searches for possible escape mutations
	##details<< 
	
	##seealso<< 
	##note<< 
	path_to_file_sequence_alignment = NULL,
	### a FASTA file with sequence data. Homzygot patiens have to have a 00 instead of the HLA typeFor reference please look in
	### example file.
	path_to_file_nexus_tree = NULL,
	### a tree file of the sequences. For reference please look in example file.
	save_name_csv, 
	### the file name of the result file in csv format
	patnum_threshold = 1, 
	### the minimum number of patients of one HLA type to consider in the calculation.
	significance_level = 0.05, 
	### p-value threshold below which the results from fishers exact test should be added to output.
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
	B22
	### the position of the end of the second HLA B Allel in the description block of the FASTA file.
	){
	HLA_and_tree_inner(path_to_file_sequence_alignment, path_to_file_nexus_tree, save_name_csv, patnum_threshold, significance_level, A11, A12, A21, A22, B11, B12, B21, B22)
	
},ex=function(){
	ex <- system.file("extdata", "Example.fasta", package="SeqFeatR")
	ep <- system.file("extdata", "Example_tree.nh", package="SeqFeatR")
	hlaTree(ex,
 ep,
 "ancestral_analysis_results.csv",
 1,
 0.05,
 22,
 23,
 25,
 26,
 29,
 30,
 32,
 33)	
})

HLA_and_tree_inner <- function(path_to_file_s = NULL, path_to_file_t = NULL, save_name_csv, patnum.threshol = 1, thr.sig.f = 0.05, A11, A12, A21, A22, B11, B12, B21, B22){
	#read aa sequences from fasta file - they have to be aligned in order to work!
	ht_get_input_file_known_sequences(path_to_file_s)
	ht_get_input_file_tree(path_to_file_t)

	sequences <- .GlobalEnv[["sequences"]]
	tree <- .GlobalEnv[["tree"]]
	.GlobalEnv[["thr.sig.fi"]] <- thr.sig.f
	.GlobalEnv[["patnum.threshold"]] <- patnum.threshol

	.GlobalEnv[["A11"]] <- A11
	.GlobalEnv[["A12"]] <- A12
	.GlobalEnv[["A21"]] <- A21
	.GlobalEnv[["A22"]] <- A22
	.GlobalEnv[["B11"]] <- B11
	.GlobalEnv[["B12"]] <- B12
	.GlobalEnv[["B21"]] <- B21
	.GlobalEnv[["B22"]] <- B22
	
	min_FASTA_length <- min(width(sequences))

	ancestral_sequences <- create_ancestral_sequences(path_to_file_s)
	#ancestral_sequences <- load(paste(files.dir, "tree_data.RData",sep=""))

	a_sequences <- get_original_sequence_ph(ancestral_sequences)
	.GlobalEnv[["a_sequences"]] <- a_sequences
	#a_sequences <- read.csv2(paste(files.dir, "tree.csv",sep=""))

	#list of the tipnames in correct order (hopefully) to get the correct tip number
	tree.tip.label <- tree$tip.label

	.GlobalEnv[["tree.tip.label"]] <- tree.tip.label

	sequences.count <- length(sequences)


	.GlobalEnv[["min_FASTA_length"]] <- min_FASTA_length

	#find out all allels that occur in given sequences
	allels.all <- c()
	allels.count <- c()

	new_sequence_names <- c()

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

	allels <- sort(unique(allels.all[allels.all!="A" & allels.all!="B"]))
	.GlobalEnv[["allels"]] <- allels
	for(i in 1:length(allels)){
	  allels.count <- c(allels.count, sum(allels.all==allels[i]))
	}
	#possible one letter codes for amino acids
	aacids <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-", "B", "X")

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
	.GlobalEnv[["HLA_array"]] <- make_HLA_types_array()
	.GlobalEnv[["anc_new_array"]] <- make_array()
	.GlobalEnv[["length_allels"]] <- length(allels.count)
	length_allels <- .GlobalEnv[["length_allels"]]

	counter_a <- array(rep(1, length(allels)), c(length_allels), list(allels))

	cat (sequences.count, length(par.aacids), sum(allels.count), "\n")

	threshold <- get.threshold.f.fisher.pair(sequences.count, length(par.aacids), allels.count)

	.GlobalEnv[["threshold"]] <- threshold

	allels.counter <- 0

	allels.list <- c()

	counter <- 1

	result_matrix <- mai()

	p_values <- result_matrix[,4]

	p_v <- c()

	for (entry in 1:length(p_values)){
		value <- abs(as.numeric(p_values[entry]))
		p_v <- c(p_v,value)
	}

	#write.csv2(p_values, paste(files.dir, "p_values.csv", sep=""))

	#q_value <- qvalue(p=p_v, lambda=0.04, pi0.method="smoother", fdr.level=NULL, robust=TRUE, gui=FALSE, smooth.df=3, smooth.log.pi0=FALSE)

	q_value <- qvalue(p=p_v, lambda=0, fdr.level=0.05)


	#print (q_value$qvalues)

	result_matrix <- cbind(result_matrix, q_value$qvalues)

	colnames(result_matrix)[5] <- "q-value"

	#print (result_matrix)

	write.csv2(result_matrix, save_name_csv)
}

#hlaTree("../inst/extdata/Example_aa.fasta",
# "../inst/extdata/Example_tree.nh",
# "test.csv",
# 1,
# 0.05,
# 10,
# 11,
# 13,
# 14,
# 17,
# 18,
# 20,
# 21)

