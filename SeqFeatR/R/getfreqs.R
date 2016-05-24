#libs
require(Biostrings)
require(ggplot2)

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


gF_wo_set_input_file_known_sequences <- structure(function(
	###read aa sequences from fasta file - they have to be aligned in order to work!
	path_to_file
	### a FASTA file with sequence data. For reference please look in example file.
	){
	### the sequences from the FASTA file.
	sequences <- readAAStringSet(path_to_file)
	.GlobalEnv[["StringSetsequences"]] <- sequences
	.GlobalEnv[["sequences"]] <- sequences
},ex=function(){
	ep_wo_set_input_file_known_sequences("SeqFeatR/extdata/Example.fasta")	
})

gF_wo_set_consensus_file_known_sequences <- structure(function(
	###read aa sequences from fasta file - they have to be aligned in order to work!
	path_to_file
	### a FASTA file with sequence data. For reference please look in example file.
	){
	### the sequences from the FASTA file.
	sequences <- readAAStringSet(path_to_file)
	.GlobalEnv[["epi_consensus"]] <- as.character(sequences)
	.GlobalEnv[["sequences"]] <- sequences
},ex=function(){
	ep_wo_set_input_file_known_sequences("SeqFeatR/extdata/Example.fasta")	
})

#removes "0"s from a string
remove.zeroes <- function(s){
	if(substr(s, 1, 1) == "0")
  		a <- remove.zeroes(substr(s, 2, nchar(s)))  
	else
  		a <- s
}


getfreqs <- structure(function(
	path_to_file_sequence_alignment=NULL,
	save_csv,
	save_png,
	epitope_position_start,
	epitope_position_end,
	path_to_file_consensus,
	patnum_threshold = 1,
	### the minimum number of patients of one HLA type to consider in the calculation.
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
	one_feature=FALSE
	### if there is only one identifier
	){
	get_freqs_inner(path_to_file_sequence_alignment, save_csv, save_png, epitope_position_start, epitope_position_end, path_to_file_consensus, patnum_threshold, A11, A12, A21, A22, B11, B12, B21, B22, one_feature)
},ex=function(){
	getfreqs("../inst/extdata/Example_aa.fasta",
	"../inst/extdata/Freqs.csv",
	"../inst/extdata/Freqs.png",
	1,
	8,
	"../inst/extdata/Example_Consensus_aa.fasta",
	patnum.threshold=1,
	A11 = 10,
	A12 = 11,
	A21 = 13,
	A22 = 14,
	B11 = 17,
	B12 = 18,
	B21 = 20,
	B22 = 21,
	one_ident = FALSE)
})


get_freqs_inner <- function(path_to_file_s, save_csv, save_png, epi_pos_start, epi_pos_end, epi_consensus, patnum.threshold, A11, A12, A21, A22, B11, B12, B21, B22, one_ident){
	
	ggplot <- NULL
	aes <- NULL
	kind <- NULL
	geom_bar <- NULL
	facet_grid <- NULL
	xlab <- NULL
	ylab <- NULL
	theme <- NULL
	element_text <- NULL
	ggsave <- NULL

	if (!is.null(path_to_file_s) && path_to_file_s == 0){
		path_to_file_s <- NULL
	}
	if (!is.null(path_to_file_s) && !(path_to_file_s == "") ){
		gF_wo_set_input_file_known_sequences(path_to_file_s)
	}

	if (!is.null(epi_consensus) && epi_consensus == 0){
		epi_consensus <- NULL
	}
	if (!is.null(epi_consensus) && !(epi_consensus == "") ){
		gF_wo_set_consensus_file_known_sequences(epi_consensus)
	}

	sequences <- .GlobalEnv[["StringSetsequences"]]
	epi_consensus <- .GlobalEnv[["epi_consensus"]]

	#find the shortest length of FASTA String, if the sequences are of different length
	#min_FASTA_length holds this value
	.GlobalEnv[["min_FASTA_length"]] <- min(width(sequences))

	min_FASTA_length <- .GlobalEnv[["min_FASTA_length"]]
	if (min_FASTA_length < epi_pos_start | min_FASTA_length < epi_pos_end){
		print ("Warning! Your alignment is shorter than your epitope end position!")
	}

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

	for(i in 1:length(allels)){
	  allels.count <- c(allels.count, sum(allels.all==allels[i]))
	}

	# allels.counter counts the allels which are abouve the patnum.threshold
	allels.counter <- 1

	sequences_names <- names(sequences)

	result_matrix_abs <- matrix(rep(0, length(allels.count)*(nchar(epi_consensus)+2)), ncol=nchar(epi_consensus)+2)
	colnames(result_matrix_abs) <- c("Allel", "Freq", strsplit(epi_consensus, "")[[1]])

	result_matrix_rel <- matrix(rep(0, length(allels.count)*(nchar(epi_consensus)+2)), ncol=nchar(epi_consensus)+2)
	colnames(result_matrix_rel) <- c("Allel", "Freq", strsplit(epi_consensus, "")[[1]])

	result_matrix_abs_not <- matrix(rep(0, length(allels.count)*(nchar(epi_consensus)+2)), ncol=nchar(epi_consensus)+2)
	colnames(result_matrix_abs) <- c("Allel", "Freq", strsplit(epi_consensus, "")[[1]])

	result_matrix_rel_not <- matrix(rep(0, length(allels.count)*(nchar(epi_consensus)+2)), ncol=nchar(epi_consensus)+2)
	colnames(result_matrix_rel) <- c("Allel", "Freq", strsplit(epi_consensus, "")[[1]])

	for(allel.num in 1:length(allels)){
  		if(allels.count[allel.num]>patnum.threshold){
			cat ("allel: ", allel.num, "\n")
			current_allel <- allels[allel.num]
			has_allel <- grep(current_allel, new_sequence_names)
			has_allel_not <- grep(current_allel, new_sequence_names, invert = T)

			result_matrix_abs[allels.counter,2] <- length(has_allel)
			result_matrix_abs_not[allels.counter,2] <- length(has_allel_not)

			result_matrix_abs[allels.counter,1] <- allels[allel.num]
			result_matrix_rel[allels.counter,1] <- allels[allel.num]
			result_matrix_abs_not[allels.counter,1] <- paste("!", allels[allel.num])
			result_matrix_rel_not[allels.counter,1] <- paste("!", allels[allel.num])
			# for each letter in FASTA:
			col <- 3
	    		for (j in epi_pos_start:epi_pos_end){
	      			#for each point in sequence do the following:
	      			# for each aa in length of possible aas:
				just_this_position <- subseq(sequences, j, j)
				mat <- as.matrix(just_this_position)
				aa_in_consensus <- substr(epi_consensus, j, j)
				subset_has_allel <- mat[has_allel]
				not_consensus <- length(which(subset_has_allel != aa_in_consensus))
				result_matrix_abs[allels.counter,col] <- not_consensus
				result_matrix_rel[allels.counter,col] <- not_consensus/length(subset_has_allel)

				subset_has_not_allel <- mat[has_allel_not]
				not_consensus_not <- length(which(subset_has_not_allel != aa_in_consensus))
				result_matrix_abs_not[allels.counter,col] <- not_consensus_not
				result_matrix_rel_not[allels.counter,col] <- not_consensus_not/length(subset_has_not_allel)

				col <- col + 1
				
			}
			allels.counter <- allels.counter + 1
		}
	}
	if (length(which(result_matrix_abs[,1]=="0")) > 0){
		result_matrix_abs <- result_matrix_abs[-(which(result_matrix_abs[,1]=="0")),]
		result_matrix_rel <- result_matrix_rel[-(which(result_matrix_rel[,1]=="0")),]
		result_matrix_abs_not <- result_matrix_abs_not[-(which(result_matrix_abs_not[,1]=="0")),]
		result_matrix_rel_not <- result_matrix_rel_not[-(which(result_matrix_rel_not[,1]=="0")),]
	}

	result_matrix <- result_matrix_abs[,c(1,2)]
	result_matrix_not <- result_matrix_abs_not[,c(1,2)]
	coln <- c("Allels", "Freqs")

	for (i in 3:ncol(result_matrix_rel)){
		result_matrix <- cbind(result_matrix, cbind(result_matrix_abs[,i], result_matrix_rel[,i]))
		coln <- c(coln, paste(colnames(result_matrix_abs)[i], "_abs", sep=""), paste(colnames(result_matrix_abs)[i], "_rel", sep=""))
	}
	colnames(result_matrix) <- coln

	coln <- c("Allels", "Freqs")
	for (i in 3:ncol(result_matrix_rel_not)){
		result_matrix_not <- cbind(result_matrix_not, cbind(result_matrix_abs_not[,i], result_matrix_rel_not[,i]))
		coln <- c(coln, paste(colnames(result_matrix_abs_not)[i], "_abs", sep=""), paste(colnames(result_matrix_abs_not)[i], "_rel", sep=""))
	}
	colnames(result_matrix_not) <- coln

	result_matrix <- rbind(result_matrix, result_matrix_not)
  
	grap_result <- c()
	for (i in 1:nrow(result_matrix)){
		pos <- 3
		for (j in 1:nchar(epi_consensus) ){
      			if (length(grep("!", result_matrix[i,1], fixed=TRUE)) > 0){
        			grap_result <- rbind(grap_result, c(unlist(strsplit(result_matrix[i,1], split=' ', fixed=TRUE))[2], letter=paste(j, "_", substr(epi_consensus, j, j), sep=""), result_matrix[i,pos], kind="not allel"))
     			}else{
				grap_result <- rbind(grap_result, c(result_matrix[i,1], letter=paste(j, "_", substr(epi_consensus, j, j), sep=""), result_matrix[i,pos], kind="allel"))
      			}
			pos <- pos+2
		}
	}
	grap_result <- as.data.frame(grap_result)
  print (grap_result)
	colnames(grap_result) <- c("allels", "letter", "value", "kind")

	grap_result$value <- as.numeric(as.character(grap_result$value))
	grap_result <- within(grap_result, letter <- factor(letter, levels = unique(grap_result$letter),ordered = TRUE))

	plot <- ggplot(data=grap_result, aes(x=letter, y=value, fill=kind)) + geom_bar(stat="identity", position="dodge") + facet_grid(allels ~.) + xlab("") +  ylab("not epitope amino acid[%]") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	ggsave(plot, file=save_png)

	write.csv2(result_matrix, file=save_csv)
}

#getfreqs("Example_aa.fasta", #<- Deine Fasta
# "Freqs.csv", 
#  "Freqs.svg",
# 1, #<- Kann so bleiben
# 40, #<- Laenge des Consensus (was ja gleichlang wie dein alignment sein sollte)
# "Example_Consensus_aa.fasta", # <- Dein Consensus
# patnum.threshold=1,
# A11 = 10, # <- HLA in Header
# A12 = 11, # <- HLA in Header
# A21 = 13, # <- HLA in Header
# A22 = 14, # <- HLA in Header
# B11 = 17, # <- HLA in Header
# B12 = 18, # <- HLA in Header
# B21 = 20, # <- HLA in Header
# B22 = 21, # <- HLA in Header
# one_ident = FALSE)
